/*
 * ParallelRunner.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "ParallelRunner.hpp"

namespace castle {
function<void(size_t, size_t)> ParallelRunner::empty_ranged_func = [&](size_t, size_t){};
ParallelRunner::ParallelRunner() {

}

ParallelRunner::~ParallelRunner() {
}

void ParallelRunner::create_bfi_index(vector<int64_t>& block_boundary, const string& input_path, const string& output_index_path, const int64_t block_size) {
	using namespace BamTools;
	BamReader serial_reader;
	if (!serial_reader.Open(input_path)) {
		cerr << (boost::format("[ParallelRunner.create_bfi_index] cannot open the BAM file: %s\n") % input_path).str();
		return;
	}
	auto& data = serial_reader.GetBGZF();
	BamAlignment al;
	int64_t n_entries = 0;
	string previous_name;
	block_boundary.push_back(0);
	int64_t global_n_entries = 0;
	cout << "[ParallelRunner.create_bfi_index] start indexing the BAM file\n";
	int64_t previous_pos = -1;
	vector<char> tmp_buf;
	while (serial_reader.ReadNextAlignment(tmp_buf)) {
		++global_n_entries;
		++n_entries;
		if (n_entries > block_size) {
			while (serial_reader.GetNextAlignmentWithName(al)) {
				if (!previous_name.empty() && previous_name != al.Name) {
					break;
				}
				previous_name = al.Name;
				previous_pos = data.Tell();
			}
			n_entries = 0;
			previous_name = "";
			block_boundary.push_back(previous_pos);
		}
	}
	serial_reader.Close();
	block_boundary.push_back(numeric_limits<int64_t>::max());
	ofstream out(output_index_path, ios::binary);
	for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
		out << block_boundary[block_id] << "\n";
	}
	cout << (boost::format("[ParallelRunner.create_bfi_index] %d entries were indexed\n") % global_n_entries).str();
}

void ParallelRunner::load_bfi_index(vector<int64_t>& block_boundary, const string& input_path) {
	string line;
	ifstream in(input_path);
	while (getline(in, line, '\n')) {
		block_boundary.push_back(boost::lexical_cast<int64_t>(line));
	}
}
void ParallelRunner::create_bni_index(const string& input_path, const string& output_index_path){
	castle::TimeChecker checker;
	checker.setTarget("ParallelRunner.create_bni_index");
	checker.start();
	BamTools::BamReader serial_reader;

	string bai_index_path = input_path + ".bai";
	if (!serial_reader.Open(input_path, bai_index_path)) {
		cerr << (boost::format("[ParallelRunner.create_bni_index] cannot open the BAM file: %s\n") % input_path).str();
		return;
	}

	vector<BamTools::BlockOffset> block_boundary;
	serial_reader.GetBlockOffsets(block_boundary);
	ofstream out(output_index_path, ios::binary);
	for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
		out << block_boundary[block_id].offset << "\t" << block_boundary[block_id].ref_id << "\t" << block_boundary[block_id].position << "\n";
	}
	cout << (boost::format("[ParallelRunner.create_bni_index] %d entries were indexed\n") % block_boundary.size()).str();
	cout << checker;
}


void ParallelRunner::load_bni_index(vector<BamTools::BlockOffset>& offset_blocks, const string& input_BAM_name, const string& output_bni_index_path) {
	string a_path(input_BAM_name);
	if (!boost::filesystem::exists(output_bni_index_path)) {
		create_bni_index(a_path, output_bni_index_path);
	}
	offset_blocks.clear();
	cout << (boost::format("[ParallelRunner.load_bni_index] bni: %s\n") % output_bni_index_path).str();
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in(output_bni_index_path, ios::binary);
	while (getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		BamTools::BlockOffset an_offset;
		an_offset.offset = boost::lexical_cast<int64_t>(data[0]);
		an_offset.ref_id = boost::lexical_cast<int32_t>(data[1]);
		an_offset.position = boost::lexical_cast<int64_t>(data[2]);
		offset_blocks.push_back(an_offset);
	}
	cout << (boost::format("[ParallelRunner.load_bni_index] # blocks: %d\n") % offset_blocks.size()).str();
}

void ParallelRunner::load_bai_bni_index(vector<int64_t>& block_boundary, const string& input_path, const string& output_bni_index_path, const int64_t block_size, const int64_t n_cores) {
	using namespace BamTools;
	string a_bai_index_path(input_path + ".bai");
	BamReader reader;
	if (!reader.Open(input_path, a_bai_index_path)) {
		cerr << (boost::format("[ParallelRunner] could not open BAM file: %s\n") % input_path).str();
		cerr << (boost::format("[ParallelRunner] relevant BAI file: %s\n") % a_bai_index_path).str();
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelRunner.load_bai_bni_index");
	checker.start();
	BamAlignment al;
	const RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	cout << (boost::format("[ParallelRunner.load_bai_bni_index] # refs: %d\n") % n_refs).str();
	int64_t size_total_ref = 0;
	for (int32_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		auto& a_ref_data = a_ref_vector[ref_id];
		size_total_ref += a_ref_data.RefLength;
	}

	int64_t estimated_n_blocks = size_total_ref / (double) block_size;
	++estimated_n_blocks;

	// calculate the positions of boundary entries.
	vector<pair<int32_t, int32_t>> boundary_positions;
	boundary_positions.push_back(make_pair(0, 0));
	int64_t boundary_id = 0;
	int32_t ref_id = 0;
	int64_t n_remaining_bases = a_ref_vector[ref_id].RefLength;
	while (n_remaining_bases >= 0) {
		auto& a_ref_data = a_ref_vector[ref_id];
		int64_t last_base_pos = boundary_positions[boundary_id].second;
		n_remaining_bases = a_ref_data.RefLength - last_base_pos;
		if (n_remaining_bases >= block_size) {
			boundary_positions.push_back(make_pair(ref_id, last_base_pos + block_size));
			++boundary_id;
		} else {
			if (a_ref_data.RefLength > block_size) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + block_size), a_ref_data.RefLength)));
				++boundary_id;
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			++boundary_id;
		}
	}

	int64_t calculated_n_blocks = boundary_positions.size();
	vector<meerkat::BlockBoundary> fixed_size_blocks(calculated_n_blocks);
	vector<function<void()> > tasks;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(input_path, a_bai_index_path)) {
				return;
			}
			auto& m_bgzf = local_reader.GetBGZF();
			BamTools::BamAlignment local_alignment_entry;
			bool success = local_reader.Jump(boundary_positions[block_id].first, boundary_positions[block_id].second);
			if(success) {
				string the_pos;
				string old_pos;
				do {
					int64_t cur_offset = m_bgzf.Tell();
					success = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
					if(success) {
						the_pos = (boost::format("%s:%s") % local_alignment_entry.RefID % local_alignment_entry.Position).str();
						fixed_size_blocks[block_id].read_name = local_alignment_entry.Name;
						fixed_size_blocks[block_id].ref_id = local_alignment_entry.RefID;
						fixed_size_blocks[block_id].offset = cur_offset;
						fixed_size_blocks[block_id].pos = local_alignment_entry.Position;
						fixed_size_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
						fixed_size_blocks[block_id].jump_pos = boundary_positions[block_id].second;
						if(!old_pos.empty() && old_pos != the_pos) {
							break;
						}
					}
					else {
						cout << (boost::format("[ParallelRunner.load_bai_bni_index] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					old_pos = the_pos;
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						break;
//						cout << (boost::format("[ParallelRunner.load_bai_bni_index] unaligned: %s %d:%d\n")
//								% local_alignment_entry.Name
//								% local_alignment_entry.RefID
//								% local_alignment_entry.Position).str();
					}
				}while(true);
			}
			local_reader.Close();
		});
	}
	vector<BlockOffset> unmapped_offsets;
	tasks.push_back([&] {
		load_bni_index(unmapped_offsets, input_path, output_bni_index_path);
	});
	run_unbalanced_load(n_cores, tasks);
	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
	if (fixed_size_blocks.size() > 0 && 0 == fixed_size_blocks[0].offset) {
		fixed_size_blocks.erase(fixed_size_blocks.begin());
	}
	fixed_size_blocks.erase(unique(fixed_size_blocks.begin(), fixed_size_blocks.end()), fixed_size_blocks.end());
	for(auto& an_offset : fixed_size_blocks) {
		block_boundary.push_back(an_offset.offset);
	}

	for (auto& an_offset : unmapped_offsets) {
		if (-1 == an_offset.ref_id) {
			block_boundary.push_back(an_offset.offset);
		}
	}
	cout << (boost::format("[ParallelRunner.load_bai_bni_index] actual # blocks: %d\n") % block_boundary.size()).str();
	cout << checker;
}

} /* namespace castle */
