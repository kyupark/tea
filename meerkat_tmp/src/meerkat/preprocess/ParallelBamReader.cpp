/*
 * ParallelBamReader.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#include "ParallelBamReader.hpp"

namespace meerkat {
ParallelBamReader::ParallelBamReader() :
		num_total(0), max_read_len(151), num_unmapped(0), num_clipped(0), unmapped_rejected(0), n_read_groups(0), qenc(0) {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
//	softclips.set_empty_key(" ");
}

ParallelBamReader::~ParallelBamReader() {
}
void ParallelBamReader::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}
void ParallelBamReader::collect_boundaries(const int64_t size_block) {
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.collect_boundaries");
	checker.start();
	BamTools::BamAlignment al;
	bool has_index = reader.HasIndex();
	if (!has_index) {
		cout << "[ParallelBamReader.collect_boundaries] no index is found so it is generating.";
		reader.CreateIndex();
	}
	fixed_size_blocks.clear();
	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	cout << (boost::format("[ParallelBamReader.collect_boundaries] # refs: %d\n") % n_refs).str();
	int64_t size_total_ref = 0;
	for (int32_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		auto& a_ref_data = a_ref_vector[ref_id];
//cout << a_ref_data.RefName << ": " << a_ref_data.RefLength << "\n";
		size_total_ref += a_ref_data.RefLength;
	}

	int64_t estimated_n_blocks = size_total_ref / (double) size_block;
	++estimated_n_blocks;

// calculate the positions of boundary entries.
	vector<pair<int32_t, int32_t>> boundary_positions;
	boundary_positions.push_back(make_pair(0, 0));
	if (verbose) {
		cout << "[ParallelBamReader.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
	}
	int64_t boundary_id = 0;
	int32_t ref_id = 0;
	int64_t n_remaining_bases = a_ref_vector[ref_id].RefLength;
	while (n_remaining_bases >= 0) {
		auto& a_ref_data = a_ref_vector[ref_id];
		int64_t last_base_pos = boundary_positions[boundary_id].second;
		n_remaining_bases = a_ref_data.RefLength - last_base_pos;
		if (n_remaining_bases >= size_block) {
			boundary_positions.push_back(make_pair(ref_id, last_base_pos + size_block));
			++boundary_id;
			if (verbose) {
				cout << "[ParallelBamReader.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
		} else {
			if (a_ref_data.RefLength > size_block) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
				++boundary_id;
				if (verbose) {
					cout << "[ParallelBamReader.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
				}
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			if (verbose) {
				cout << "[ParallelBamReader.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
			++boundary_id;
		}
	}

//for(uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
//auto& a_ref_data = a_ref_vector[ref_id];
//boundary_positions.push_back(make_pair(ref_id, a_ref_data.RefLength));
//}

	int64_t calculated_n_blocks = boundary_positions.size();

	cout << "[ParallelBamReader.collect_boundaries] estimated # blocks: " << estimated_n_blocks << "\n";
	cout << "[ParallelBamReader.collect_boundaries] calculated # blocks: " << calculated_n_blocks << "\n";
	cout << "[ParallelBamReader.collect_boundaries] total ref. size: " << size_total_ref << "\n";

	fixed_size_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	cout << "[ParallelBamReader.collect_boundaries] collect boundary positions\n";
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			auto& m_bgzf = local_reader.GetBGZF();
			BamTools::BamAlignment local_alignment_entry;
			bool success = local_reader.Jump(boundary_positions[block_id].first, boundary_positions[block_id].second);
			if(success) {
				do {
					int64_t cur_offset = m_bgzf.Tell();
					success = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
					if(success) {
						fixed_size_blocks[block_id].read_name = local_alignment_entry.Name;
						fixed_size_blocks[block_id].ref_id = local_alignment_entry.RefID;
						fixed_size_blocks[block_id].offset = cur_offset;
						fixed_size_blocks[block_id].pos = local_alignment_entry.Position;
						fixed_size_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
						fixed_size_blocks[block_id].jump_pos = boundary_positions[block_id].second;
					}
					else {
						cout << (boost::format("[ParallelBamReader.collect_boundaries] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						cout << (boost::format("[ParallelBamReader.collect_boundaries] unaligned: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
					}
				}while(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position);
			}
			local_reader.Close();
		});
	}
	vector<BlockOffset> unmapped_offsets;
	tasks.push_back([&]{
		collect_boundaries_alt(unmapped_offsets);
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
	if (fixed_size_blocks.size() > 0 && 0 == fixed_size_blocks[0].offset) {
		fixed_size_blocks.erase(fixed_size_blocks.begin());
	}
	fixed_size_blocks.erase(unique(fixed_size_blocks.begin(), fixed_size_blocks.end()), fixed_size_blocks.end());
	unmapped_included_blocks = fixed_size_blocks;
	BlockBoundary a_block_boundary;
	a_block_boundary.read_name = "last";
	a_block_boundary.ref_id = n_refs;
	a_block_boundary.offset = numeric_limits<int64_t>::max();
	a_block_boundary.pos = 0;
	a_block_boundary.jump_pos = -1;
	fixed_size_blocks.push_back(a_block_boundary);
	for(auto& an_offset : unmapped_offsets) {
		if(-1 == an_offset.ref_id) {
			BlockBoundary a_boundary;
			a_boundary.offset = an_offset.offset;
			a_boundary.ref_id = an_offset.ref_id;
			a_boundary.pos = an_offset.position;
			unmapped_included_blocks.push_back(a_boundary);
		}
	}

	cout << (boost::format("[ParallelBamReader.collect_boundaries] actual # blocks: %d\n") % fixed_size_blocks.size()).str();
	if (verbose) {
		for (uint64_t block_id = 0; block_id < fixed_size_blocks.size(); ++block_id) {
			auto& a_block = fixed_size_blocks[block_id];
			cout << (boost::format("[ParallelBamReader.collect_boundaries] BlockBoundary(%d): %d %s %d-%d(%d)\n") % block_id % a_block.offset % a_block.read_name % a_block.ref_id % a_block.pos % a_block.jump_pos).str();
		}
	}
	actual_blocks = fixed_size_blocks;
	cout << checker;
}

void ParallelBamReader::collect_chromosome_wide_boundaries() {

	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.collect_chromosome_wide_boundaries");

	if (silent) {
		checker.start_without_output();
	} else {
		checker.start();
	}
	chromosome_wide_blocks.clear();
	BamTools::BamAlignment al;
	bool has_index = reader.HasIndex();
	if (!has_index) {
		cout << "[ParallelBamReader.collect_chromosome_wide_boundaries] no index is found and generating.";
		reader.CreateIndex();
	}

	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	if (!silent) {
		cout << (boost::format("[ParallelBamReader.collect_chromosome_wide_boundaries] # refs: %d\n") % n_refs).str();
	}
	int64_t size_total_ref = 0;
	for (int32_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		auto& a_ref_data = a_ref_vector[ref_id];
		//		cout << a_ref_data.RefName << ": " << a_ref_data.RefLength << "\n";
		size_total_ref += a_ref_data.RefLength;
	}
//			const int64_t size_block = 4096000;
//			int64_t estimated_n_blocks = size_total_ref / (double) size_block;
//			++estimated_n_blocks;

// calculate the positions of boundary entries.
	vector<pair<int32_t, int32_t>> boundary_positions;
//			boundary_positions.push_back(make_pair(0, 0));
//			if (!silent) {
//				if (verbose) {
//					cout
//					<< "[ParallelBamReader.collect_chromosome_wide_boundaries] Inserted-"
//					<< (boundary_positions.size() - 1) << ":"
//					<< boundary_positions.back().first << "/"
//					<< boundary_positions.back().second << "\n";
//				}
//			}
	for (int64_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		boundary_positions.push_back(make_pair(ref_id, 0));
	}

	int64_t calculated_n_blocks = boundary_positions.size();

	if (!silent) {
//				cout
//				<< "[ParallelBamReader.collect_boundaries] estimated # blocks: "
//				<< estimated_n_blocks << "\n";
//				cout
//				<< "[ParallelBamReader.collect_chromosome_wide_boundaries] calculated # blocks: "
//				<< calculated_n_blocks << "\n";
		cout << "[ParallelBamReader.collect_chromosome_wide_boundaries] total ref. size: " << size_total_ref << "\n";
	}

	chromosome_wide_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	if (!silent) {
		cout << "[ParallelBamReader.collect_chromosome_wide_boundaries] collect boundary positions\n";
	}
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			BamTools::BamAlignment local_alignment_entry;
			bool success = local_reader.Jump(boundary_positions[block_id].first, boundary_positions[block_id].second);
			if(success) {
				do {
					success = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
					if(success) {
						chromosome_wide_blocks[block_id].read_name = local_alignment_entry.Name;
						chromosome_wide_blocks[block_id].ref_id = local_alignment_entry.RefID;
						chromosome_wide_blocks[block_id].pos = local_alignment_entry.Position;
						chromosome_wide_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
						chromosome_wide_blocks[block_id].jump_pos = boundary_positions[block_id].second;
					} else {
						cout << (boost::format("[ParallelBamReader.collect_chromosome_wide_boundaries] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						cout << (boost::format("[ParallelBamReader.collect_chromosome_wide_boundaries] unaligned: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
					}
				}while(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position);
			}
			local_reader.Close();
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//			sort(chromosome_wide_blocks.begin(), chromosome_wide_blocks.end());
//			independent_blocks.erase(unique(chromosome_wide_blocks.begin(), chromosome_wide_blocks.end()),
//					chromosome_wide_blocks.end());
	BlockBoundary a_block_boundary;
	a_block_boundary.read_name = "last";
	a_block_boundary.ref_id = n_refs;
	a_block_boundary.pos = 0;
	a_block_boundary.jump_pos = -1;
	chromosome_wide_blocks.push_back(a_block_boundary);

	if (!silent) {
		cout << (boost::format("[ParallelBamReader.collect_chromosome_wide_boundaries] actual # blocks: %d\n") % chromosome_wide_blocks.size()).str();
	}
	if (verbose) {
		for (uint64_t block_id = 0; block_id < chromosome_wide_blocks.size(); ++block_id) {
			auto a_block = chromosome_wide_blocks[block_id];
			cout << (boost::format("[ParallelBamReader.collect_chromosome_wide_boundaries] BlockBoundary(%d): %s %d-%d\n") % block_id % a_block.read_name % a_block.ref_id % a_block.pos).str();
		}
	}
	actual_blocks = chromosome_wide_blocks;
	if (!silent) {
		cout << checker;
	}

}

void ParallelBamReader::collect_boundaries_alt(vector<BlockOffset>& offset_blocks) {
	string a_path = options.prefix + ".bam";
	string bni_index_path;
	get_bni_index_path(a_path, bni_index_path);
	offset_blocks.clear();
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in(bni_index_path, ios::binary);
	while (getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		BamTools::BlockOffset an_offset;
		an_offset.offset = boost::lexical_cast<int64_t>(data[0]);
		an_offset.ref_id = boost::lexical_cast<int32_t>(data[1]);
		an_offset.position = boost::lexical_cast<int64_t>(data[2]);
		offset_blocks.push_back(an_offset);
	}
	cout << (boost::format("[ParallelBamReader.collect_boundaries_alt] # blocks: %d\n") % offset_blocks.size()).str();
}
void ParallelBamReader::preprocess() {
	string a_path = options.prefix + ".bam";

	string bni_path = a_path + ".bni";
	if(!options.working_dir.empty()) {
		bni_path = options.working_prefix + ".bam.bni";
	}

	create_bni_even_index(a_path, bni_path);

	int64_t size_block = 8192000;
	collect_boundaries(size_block);
	string rfile = options.prefix + ".insert.r";
	string pdffile = options.prefix + ".pdf";
	string resultfile = options.prefix + ".isinfo";
	string rdist_unmapfile = options.prefix + ".unmapped.rdist";
	string rdist_scfile = options.prefix + ".softclips.rdist";
	string sr1_file = options.prefix + ".sr.1.fq.gz";
	string sr2_file = options.prefix + ".sr.2.fq.gz";

	if (!options.working_dir.empty()) {
		rfile = options.working_prefix + ".insert.r";
		pdffile = options.working_prefix + ".pdf";
		resultfile = options.working_prefix + ".isinfo";
		rdist_unmapfile = options.working_prefix + ".unmapped.rdist";
		rdist_scfile = options.working_prefix + ".softclips.rdist";
		sr1_file = options.working_prefix + ".sr.1.fq.gz";
		sr2_file = options.working_prefix + ".sr.2.fq.gz";
	}

//	collect_boundaries_alt();

	if (!boost::filesystem::exists(rfile) || !boost::filesystem::exists(pdffile) || !boost::filesystem::exists(resultfile) || !boost::filesystem::exists(rdist_unmapfile) || !boost::filesystem::exists(rdist_scfile)
	|| !boost::filesystem::exists(sr1_file) || !boost::filesystem::exists(sr2_file)) {
		/*** First run through the BAM: record some basic statistics ***
		 *** and discover all soft clipped reads.  All of the soft   ***
		 *** clipped reads are stored in a map for the second run.   ***
		 *** GetNextAlignmentWithName() is a custom extension to the ***
		 *** BamTools library--it works nearly as fast as *Core, but ***
		 *** also fills in the read name.  We need to store the read ***
		 *** name so that we can identify mated pairs in the second  ***
		 *** pass.                                                   ***/

		collect_basic_statistics();
		/*** Second run through the file.  We need to output soft     ***
		 *** clips and their mates simultaneously to different files. ***
		 *** We use a map to buffer read mates which we know (using   ***
		 *** the previously constructed map) are either soft clipped  ***
		 *** or have soft clipped mates.  Once both mates are read,   ***
		 *** we output both reads and remove them from the buffer.    ***
		 *** All streams are gzipped streams to reduce I/O load.      ***/
		//	size_block = 8192000;
		//	collect_boundaries(size_block);
		//	collect_chromosome_wide_boundaries();
		collect_second_statistics();
		create_reports();
	}
	align_clipped_reads();
}

void ParallelBamReader::collect_basic_statistics() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.collect_basic_statistics");
	checker.start();
// only for debugging

	string first_stat_filename = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		first_stat_filename = options.working_prefix + ".firststat";
	}
	if (boost::filesystem::exists(first_stat_filename)) {
		load_first_stat();
		return;
	}
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);

//	vector<BlockBoundary> unmapped_included_blocks = unmapped_included_blocks;
	int64_t calculated_n_blocks = unmapped_included_blocks.size();
	vector<function<void()> > tasks;
	vector<BAMStatistics> local_BAMStatistics(calculated_n_blocks - 1);
	vector<vector<uint64_t>> local_quality_sample_spaces(calculated_n_blocks - 1);
//	vector<unordered_map<string, int>> local_softclips(calculated_n_blocks - 1);
//	vector<boost::unordered_set<string>> local_softclips(calculated_n_blocks - 1);
	vector<vector<string>> local_softclips(calculated_n_blocks - 1);
//	vector<StringDenseSet> local_softclips(calculated_n_blocks - 1);
//	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
//		local_softclips[block_id].set_empty_key(" ");
//	}
// only for debugging
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
	if (verbose) {
		for (auto itr = unmapped_included_blocks.begin(); unmapped_included_blocks.end() != itr; ++itr) {
			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
		}
	}
	vector<string> read_groups;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;
			int64_t num_unmapped = 0;
			int64_t num_clipped = 0;
			int64_t max_read_len = 0;
			int bps = min(options.big_s_bps, options.frag_size);
//			vector<string> the_local_softclips;
				auto& the_local_softclips = local_softclips[block_id];

				BamTools::BamAlignment local_alignment_entry;
				vector<uint64_t> local_quality_sample_space(255);
				int32_t the_current_ref_id = unmapped_included_blocks[block_id].ref_id;
				int32_t the_current_ref_pos = unmapped_included_blocks[block_id].pos;

				int32_t the_next_ref_id = unmapped_included_blocks[block_id + 1].ref_id;
				int32_t the_next_ref_pos = unmapped_included_blocks[block_id + 1].pos;
//				uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
				string str_block_id = boost::lexical_cast<string>(block_id);

				string the_next_block_read_name = unmapped_included_blocks[block_id + 1].read_name;

				int64_t the_current_ref_offset = unmapped_included_blocks[block_id].offset;
				int64_t the_next_ref_offset = unmapped_included_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;

//				auto& ref_vec = local_reader.GetReferenceData();

//string the_current_block_read_name = actual_blocks[block_id].read_name;
//				bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//bool region_success = local_reader.SetRegion(the_current_ref_id, the_current_ref_pos,
//the_next_ref_id, the_next_ref_pos);
//				if(!jump_success) {
//					cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (Jump fail): (%d/%d)-(%d/%d)\n")
//							% block_id % the_current_ref_id % the_current_ref_pos
//							% the_next_ref_id % the_next_ref_pos).str();
//					local_reader.Close();
//					return;
//				}

				if(verbose) {
					cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (start) (%d/%d)-(%d/%d)\n")
							% block_id % the_current_ref_id % the_current_ref_pos
							% the_next_ref_id % the_next_ref_pos).str();
				}
				while (local_reader.LoadNextAlignmentWithName(local_alignment_entry)) {
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;
//					if(cur_offset >= the_next_ref_offset) {
//						if(local_alignment_entry.RefID == the_next_ref_id
//								&& local_alignment_entry.Position == the_next_ref_pos
//								&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//								&& local_alignment_entry.Name == the_next_block_read_name
//						) {
//							cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (correct)\n")
//														% block_id).str();
//						} else {
//							cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (incorrect)\n")
//																					% block_id).str();
//						}
//						break;
//					}
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%d %s %d-%d %d")
//								% cur_offset % local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (first) %s\n")
//								% block_id % a_block_boundary_str).str();
//					}

//					if(local_alignment_entry.RefID == the_next_ref_id
//							&& local_alignment_entry.Position == the_next_ref_pos
//							&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//							&& local_alignment_entry.Name == the_next_block_read_name
//					) {
//						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (last) offset: jump: %d, calc: %d\n")
//																	% block_id % prev_offset % the_next_ref_offset).str();
//						break;
//					}

					++num_total;
					if (local_alignment_entry.Length > max_read_len) {
						max_read_len = local_alignment_entry.Length;
					}

					witness_qualities(local_alignment_entry, local_quality_sample_space);

					/* UM+UM pairs and singleton UMs aren't interesting */
					if(!local_alignment_entry.IsMapped() && (!local_alignment_entry.IsMateMapped()
									|| !local_alignment_entry.IsPaired())) {
						++num_unmapped;
						continue;
					}

					/* Save the names of interesting read pairs.  "Interesting"
					 * is defined as: (1) the mate is mapped and (2) either this
					 * read contains a "big" S operation or is unmapped. */
					/* S can only occur in the first or last cigar op */
					bool extract = ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.front(), bps)
					|| ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.back(), bps) || !local_alignment_entry.IsMapped();
					if (extract && local_alignment_entry.IsMateMapped()) {
//						bool debug = string::npos != local_alignment_entry.Name.find("ST-E00104:502:HFJN5CCXX:6:2101:30888:66004");
//						the_local_softclips[local_alignment_entry.Name] = ReadGroup::getMateNumber(local_alignment_entry);
						string a_key = local_alignment_entry.Name + boost::lexical_cast<string>(ReadGroup::getMateNumber(local_alignment_entry));
//						if(debug) {
//							cout << "[ParallelBamReader.collect_basic_statistics] key: " << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
//							cout << "[ParallelBamReader.collect_basic_statistics] key: " << a_key << "\n";
//						}
						the_local_softclips.push_back(a_key);
//						the_local_softclips.insert(a_key);
						++num_clipped;
					}
				}
//				cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (last) offset: jump: %d, cur: %d, next: %d\n")
//					% block_id % prev_offset % the_current_ref_offset % the_next_ref_offset).str();
				local_reader.Close();
				local_BAMStatistics[block_id].num_total = num_total;
				local_BAMStatistics[block_id].num_unmapped = num_unmapped;
				local_BAMStatistics[block_id].num_clipped = num_clipped;
				local_BAMStatistics[block_id].max_read_len = max_read_len;
				local_quality_sample_spaces[block_id] = local_quality_sample_space;
//				auto& the_local_softclips = local_softclips[block_id];
//				local_softclips[block_id].swap(the_local_softclips);
//				local_softclips[block_id] = the_local_softclips;
				done_vector[block_id] = 'D';

				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (last-wrong) %s\n")
								% block_id % a_block_boundary_str).str();
					} else {
						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (last) %s\n")
								% block_id % a_block_boundary_str).str();
					}
				} else {
					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
					double processed = n/(double)done_vector.size() * 100.0;
					cout << (boost::format("%.2f %%\n") % processed).str();
				}
			});
	}
// sometimes the unaligned reads are in the last portion of BAM file, hence
// changing the order.
//	swap(tasks[0], tasks.back());
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	cout << "[ParallelBamReader.collect_basic_statistics] gather scattered statistics\n";
	vector<uint64_t> quality_sample_space(255);
	int64_t n_soft_clips = 0;
	for (int64_t block_id = calculated_n_blocks - 2; block_id >= 0; --block_id) {
		num_total += local_BAMStatistics[block_id].num_total;
		num_unmapped += local_BAMStatistics[block_id].num_unmapped;
		num_clipped += local_BAMStatistics[block_id].num_clipped;
		max_read_len = max(local_BAMStatistics[block_id].max_read_len, max_read_len);
		for (int64_t quality_id = 0; quality_id < 255; ++quality_id) {
			auto& a_quality_sample_space = local_quality_sample_spaces[block_id];
			quality_sample_space[quality_id] += a_quality_sample_space[quality_id];
		}
		auto& current_soft_clips = local_softclips[block_id];
		n_soft_clips += current_soft_clips.size();
	}

	softclips.rehash(n_soft_clips);

// this runs in reverse order in order to correctly assign the softclip information
	for (int64_t block_id = calculated_n_blocks - 2; block_id >= 0; --block_id) {
		auto& current_soft_clips = local_softclips[block_id];
//		softclips.insert(current_soft_clips.begin(), current_soft_clips.end());
		for (auto itr = current_soft_clips.begin(); current_soft_clips.end() != itr; ++itr) {
			auto& the_key = *itr;
			char the_mate_num = the_key[the_key.size() - 1];
			if ('1' == the_mate_num) {
				the_key[the_key.size() - 1] = '2';
				if (softclips.end() == softclips.find(the_key)) {
					the_key[the_key.size() - 1] = '1';
//					softclips[itr->first] = itr->second;
					softclips.insert(the_key);
				}
			} else {
				the_key[the_key.size() - 1] = '1';
				if (softclips.end() == softclips.find(the_key)) {
					the_key[the_key.size() - 1] = '2';
//					softclips[itr->first] = itr->second;
					softclips.insert(the_key);
				}
			}
		}
	}

	if (verbose) {
//		for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
//			cout << itr->first << ":" << itr->second << "\n";
//		}
		for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
			cout << *itr << "\n";
		}
	}

	cout << options.fname << ": \n   " << num_total << " reads\n   " << num_unmapped << " unaligned\n   " << num_clipped << " clipped (" << softclips.size() << " hashed)\n   " << max_read_len << " bps in the longest read\n";
	report_qualities(quality_sample_space);
	store_first_stat();
	cout << checker;
}

void ParallelBamReader::load_first_stat() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.load_first_stat");
	checker.start();
	vector<function<void()> > tasks;
	string first_stat_filename = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		first_stat_filename = options.working_prefix + ".firststat";
	}
	{
		string line;
		ifstream in(first_stat_filename, ios::binary);
		getline(in, line, '\n');
		num_total = boost::lexical_cast<int64_t>(line);
		getline(in, line, '\n');
		num_unmapped = boost::lexical_cast<int64_t>(line);
		getline(in, line, '\n');
		num_clipped = boost::lexical_cast<int64_t>(line);
		getline(in, line, '\n');
		max_read_len = boost::lexical_cast<int64_t>(line);
		getline(in, line, '\n');
		qenc = boost::lexical_cast<int>(line);
	}

	int64_t file_size = castle::IOUtils::get_file_size(first_stat_filename);
	int64_t block_size = file_size / (double) n_cores;
	int64_t n_blocks = file_size / (double) block_size;
	vector<uint64_t> skip_points;
	castle::IOUtils::find_skip_points(skip_points, first_stat_filename, block_size, file_size, n_blocks, n_cores);
	cout << (boost::format("[ParallelBamReader.load_first_stat] # blocks: %d\n") % n_blocks).str();
	vector<vector<string>> the_lines(n_blocks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id, block_size] {
			string line;
			const int64_t beg_pos = skip_points[block_id];
			int64_t cur_pos = beg_pos;
			const int64_t next_pos = skip_points[block_id + 1];
			auto& local_lines = the_lines[block_id];
			local_lines.reserve(block_size);
			ifstream in(first_stat_filename, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			if(0 == block_id) {
				for(int64_t l_id = 0; l_id < 5; ++l_id) {
					getline(in, line, '\n');
				}
			}
			while (getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				local_lines.push_back(line);
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[ParallelBamReader.load_first_stat] start inserting soft-clipped read names\n";
	int64_t n_lines = 0;
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		n_lines += the_lines[block_id].size();
	}
	softclips.rehash(n_lines);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		auto& local_lines = the_lines[block_id];
		for (auto a_line : local_lines) {
			softclips.insert(a_line);
		}
	}
//	while (getline(in, line, '\n')) {
////			castle::StringUtils::c_string_multi_split(line, delims, a_cols);
////			softclips[a_cols[0]] = boost::lexical_cast<int>(a_cols[1]);
//		softclips.insert(line);
//	}
	cout << options.fname << ": \n   " << num_total << " reads\n   " << num_unmapped << " unaligned\n   " << num_clipped << " clipped (" << softclips.size() << " hashed)\n   " << max_read_len << " bps in the longest read\n";
	cout << checker;
}

void ParallelBamReader::store_first_stat() {
	string first_stat_filename = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		first_stat_filename = options.working_prefix + ".firststat";
	}
	ofstream scfile(first_stat_filename, ios::binary);
	scfile << num_total << "\n" //
			<< num_unmapped << "\n" //
			<< num_clipped << "\n" //
			<< max_read_len << "\n" //
			<< qenc << "\n";
	for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
//		scfile << itr->first << "\t" << static_cast<int64_t>(itr->second) << "\n";
		scfile << *itr << "\n";
	}
}

void ParallelBamReader::collect_basic_statistics_alt() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.collect_basic_statistics");
	checker.start();
	// only for debugging
//		string first_stat_filename = options.prefix + ".firststat";
//		if (boost::filesystem::exists(first_stat_filename)) {
//			string line;
//			ifstream in(first_stat_filename, ios::binary);
//			getline(in, line, '\n');
//			num_total = boost::lexical_cast<int64_t>(line);
//			getline(in, line, '\n');
//			num_unmapped = boost::lexical_cast<int64_t>(line);
//			getline(in, line, '\n');
//			num_clipped = boost::lexical_cast<int64_t>(line);
//			getline(in, line, '\n');
//			max_read_len = boost::lexical_cast<int64_t>(line);
//			getline(in, line, '\n');
//			qenc = boost::lexical_cast<int>(line);
//			vector<string> a_cols;
//			const char* delims = "\t";
//			while (getline(in, line, '\n')) {
//				castle::StringUtils::c_string_multi_split(line, delims, a_cols);
//				softclips[a_cols[0]] = boost::lexical_cast<int>(a_cols[1]);
//			}
//			cout << options.fname << ": \n   " << num_total << " reads\n   " << num_unmapped << " unaligned\n   " << num_clipped << " clipped ("
//					<< softclips.size() << " hashed)\n   " << max_read_len << " bps in the longest read\n";
//			cout << checker;
//			return;
//		}
	vector<uint64_t> quality_sample_space(255);
	find_quality_standard_and_max_read_length(quality_sample_space);
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);

	int64_t calculated_n_blocks = actual_blocks.size();

	vector<BAMStatistics> local_BAMStatistics(calculated_n_blocks - 1);
//	vector<StringInt8Map> local_softclips(calculated_n_blocks - 1);
	vector<StringDenseSet> local_softclips(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		local_softclips[block_id].set_empty_key(" ");
	}
	vector<int64_t> unmapped_rejected_lists(calculated_n_blocks - 1);
	vector<Histogram> uhist_lists(calculated_n_blocks - 1);
	vector<boost::unordered_map<string, BamAlignment>> umbufs_first(calculated_n_blocks - 1);
	vector<boost::unordered_map<string, BamAlignment>> umbufs_second(calculated_n_blocks - 1);
//	vector<map<string, ReadGroupAlt>> read_group_maps(calculated_n_blocks - 1);

// only for debugging
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
	if (verbose) {
		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
		}
	}
//		vector<string> read_groups;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;
			int64_t num_unmapped = 0;
			int64_t num_clipped = 0;
			int bps = min(options.big_s_bps, options.frag_size);
			auto& the_local_softclips = local_softclips[block_id];
			BamTools::BamAlignment local_alignment_entry;

			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
			int32_t the_current_ref_pos = actual_blocks[block_id].pos;

			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
			if(!jump_success) {
				cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (Jump fail): (%d/%d)-(%d/%d)\n")
						% block_id % the_current_ref_id % the_current_ref_pos
						% the_next_ref_id % the_next_ref_pos).str();
				local_reader.Close();
				return;
			}
			if(verbose) {
				cout << (boost::format("[ParallelBamReader.collect_basic_statistics] block-%d (start) (%d/%d)-(%d/%d)\n")
						% block_id % the_current_ref_id % the_current_ref_pos
						% the_next_ref_id % the_next_ref_pos).str();
			}
			int64_t *cov = 0;
			int64_t covsize = sizeof(int64_t) * max_read_len;
			if (options.coverage_cutoff > 0) {
				cov = (int64_t *) malloc(covsize);
				bzero((char *) cov, covsize);
			}

			int64_t cur_pos = the_current_ref_pos;
			int64_t cur_chrom = the_current_ref_id;
			auto& refnames = local_reader.GetReferenceData();
//			auto& local_read_group_map = read_group_maps[block_id];

				string rg_name;
				int64_t local_unmapped_rejected = 0;
				const int64_t min_frag_size = options.frag_size << 1;
				auto& local_uhist = uhist_lists[block_id];
				auto& local_umbuf_first = umbufs_first[block_id];
				auto& local_umbuf_second = umbufs_second[block_id];
				ofstream umfile(options.umfname + "." + str_block_id, ios::binary);
				ofstream blist(options.blistname + "." + str_block_id, ios::binary);

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					if(verbose && 0 == num_total) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag).str();
						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
							cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (first-wrong) %s\n")
									% block_id % a_block_boundary_str).str();
						} else {
							cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (first) %s\n")
									% block_id % a_block_boundary_str).str();
						}
					}
					//if((local_alignment_entry.Name == the_next_block_read_name
					//&& the_next_ref_id <= local_alignment_entry.RefID)) {
					if(local_alignment_entry.RefID == the_next_ref_id
							&& local_alignment_entry.Position == the_next_ref_pos
							&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
							&& local_alignment_entry.Name == the_next_block_read_name
					) {
						break;
					}
					++num_total;

					/* UM+UM pairs and singleton UMs aren't interesting */
					if(!local_alignment_entry.IsMapped() && (!local_alignment_entry.IsMateMapped()
									|| !local_alignment_entry.IsPaired())) {
						++num_unmapped;
					}

					/* Save the names of interesting read pairs.  "Interesting"
					 * is defined as: (1) the mate is mapped and (2) either this
					 * read contains a "big" S operation or is unmapped. */
					/* S can only occur in the first or last cigar op */
					bool extract = ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.front(), bps)
					|| ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.back(), bps) || !local_alignment_entry.IsMapped();
					if (extract && local_alignment_entry.IsMateMapped()) {
						string key = local_alignment_entry.Name + boost::lexical_cast<string>(ReadGroup::getMateNumber(local_alignment_entry));
//						the_local_softclips[local_alignment_entry.Name] = ReadGroup::getMateNumber(local_alignment_entry);
						the_local_softclips.insert(key);
						++num_clipped;
					}
					if(!local_alignment_entry.GetReadGroup(rg_name)) {
						rg_name = "none";
					}
					if(black_listed.end() != black_listed.find(rg_name)) {
						continue;
					}

					if(!local_alignment_entry.IsMapped()) {
						//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped trimming) before: %d\n")
						//% block_id % local_alignment_entry.Length).str();
						ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
						//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped trimming) after: %d\n")
						//% block_id % local_alignment_entry.Length).str();
						local_uhist.add(local_alignment_entry);
						//string a_block_boundary_str = (boost::format("%s %d-%d %d")
						//% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
						//% local_alignment_entry.AlignmentFlag).str();
						//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped) %s\n")
						//% block_id % a_block_boundary_str).str();
						if (local_alignment_entry.Length >= min_frag_size) {
							BamAlignment copy(local_alignment_entry);
							if (local_alignment_entry.IsPaired()) {
								copy.Name += copy.IsFirstMate() ? "_1" : "_2";
							}
							ReadGroup::writeFQ(umfile, copy);
							//ReadGroup::split_read(copy, split1, split2, options.frag_size, options.n_cutoff);
						} else {
							++local_unmapped_rejected;
						}
					}
					/* 2. Handle UU pairs */
					if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
						if(local_alignment_entry.IsFirstMate()) {
							auto x_first = local_umbuf_first.find(local_alignment_entry.Name);
							if (x_first == local_umbuf_first.end()) {
								local_umbuf_first[local_alignment_entry.Name] = local_alignment_entry;
							}
						} else if(local_alignment_entry.IsSecondMate()) {
							auto x = local_umbuf_first.find(local_alignment_entry.Name);
							if(x == local_umbuf_second.end()) {
								local_umbuf_second[local_alignment_entry.Name] = local_alignment_entry;
							}
						}
					}
					/* XXX: does not account for CIGAR string data like calDepth */
					/* these 4 Is* conditions match calDepth's default mask */
					if (options.coverage_cutoff > 0 && local_alignment_entry.IsMapped()
							&& local_alignment_entry.IsPrimaryAlignment()
							&& !local_alignment_entry.IsFailedQC() && !local_alignment_entry.IsDuplicate()) {
						if (local_alignment_entry.RefID != cur_chrom && cur_chrom != -1) {
							//if("MT" == refnames[cur_chrom].RefName && 0 == cur_pos) {
							//string a_block_boundary_str = (boost::format("%s (%d-%d)->(%d-%d) %d")
							//% local_alignment_entry.Name % cur_chrom % cur_pos % local_alignment_entry.RefID % local_alignment_entry.Position
							//% local_alignment_entry.AlignmentFlag).str();
							//cout << (boost::format("[ParallelBamReader.output_cov_and_readgroup_alt] Block-%d (writing case-1) %s\n")
							//% block_id % a_block_boundary_str).str();
							//}
							for (int i = 0; i < max_read_len; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
								cov[i] = 0;
							}
						} else if (local_alignment_entry.Position != cur_pos && cur_pos != -1) {
							/* compute the number of positions in our cov
							 * buffer that need to be kicked out */
							//if("MT" == refnames[cur_chrom].RefName && 0 == cur_pos) {
							//string a_block_boundary_str = (boost::format("%s (%d-%d)->(%d-%d) %d")
							//% local_alignment_entry.Name % cur_chrom % cur_pos % local_alignment_entry.RefID % local_alignment_entry.Position
							//% local_alignment_entry.AlignmentFlag).str();
							//cout << (boost::format("[ParallelBamReader.output_cov_and_readgroup_alt] Block-%d (writing case-2) %s\n")
							//% block_id % a_block_boundary_str).str();
							//}
							int x = min(local_alignment_entry.Position - cur_pos, max_read_len);
							for (int i = 0; i < x; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
							}
							memmove(cov, cov + x, covsize - sizeof(int64_t) * x);
							bzero((char *) (cov + max_read_len - x), sizeof(int64_t) * x);
						}
						cur_pos = local_alignment_entry.Position;
						cur_chrom = local_alignment_entry.RefID;
						for (int i = 0; i < local_alignment_entry.Length; ++i) {
							++cov[i];
						}
					}
				}
				// since it is not possible to recognize overlaps at the block boundary with only above loop,
				// the additional check will be performed.
				int64_t the_next_max_pos = cur_pos + (max_read_len << 1);
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					if(local_alignment_entry.RefID != cur_chrom) {
						cur_pos = local_alignment_entry.Position;
						the_next_max_pos = cur_pos + (max_read_len << 1);
						cur_chrom = local_alignment_entry.RefID;
					}
					if(cur_pos > the_next_max_pos) {
						break;
					}
					/* XXX: does not account for CIGAR string data like calDepth */
					/* these 4 Is* conditions match calDepth's default mask */
					if (options.coverage_cutoff > 0 && local_alignment_entry.IsMapped()
							&& local_alignment_entry.IsPrimaryAlignment()
							&& !local_alignment_entry.IsFailedQC() && !local_alignment_entry.IsDuplicate()) {
						if (local_alignment_entry.RefID != cur_chrom && cur_chrom != -1) {
							//if("MT" == refnames[cur_chrom].RefName && 0 == cur_pos) {
							//string a_block_boundary_str = (boost::format("%s (%d-%d)->(%d-%d) %d")
							//% local_alignment_entry.Name % cur_chrom % cur_pos % local_alignment_entry.RefID % local_alignment_entry.Position
							//% local_alignment_entry.AlignmentFlag).str();
							//cout << (boost::format("[ParallelBamReader.output_cov_and_readgroup_alt] Block-%d (writing case-1) %s\n")
							//% block_id % a_block_boundary_str).str();
							//}
							for (int i = 0; i < max_read_len; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
								cov[i] = 0;
							}
						} else if (local_alignment_entry.Position != cur_pos && cur_pos != -1) {
							/* compute the number of positions in our cov
							 * buffer that need to be kicked out */
							//if("MT" == refnames[cur_chrom].RefName && 0 == cur_pos) {
							//string a_block_boundary_str = (boost::format("%s (%d-%d)->(%d-%d) %d")
							//% local_alignment_entry.Name % cur_chrom % cur_pos % local_alignment_entry.RefID % local_alignment_entry.Position
							//% local_alignment_entry.AlignmentFlag).str();
							//cout << (boost::format("[ParallelBamReader.output_cov_and_readgroup_alt] Block-%d (writing case-2) %s\n")
							//% block_id % a_block_boundary_str).str();
							//}
							int x = min(local_alignment_entry.Position - cur_pos, max_read_len);
							for (int i = 0; i < x; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
							}
							memmove(cov, cov + x, covsize - sizeof(int64_t) * x);
							bzero((char *) (cov + max_read_len - x), sizeof(int64_t) * x);
						}
						cur_pos = local_alignment_entry.Position;
						cur_chrom = local_alignment_entry.RefID;
						for (int i = 0; i < local_alignment_entry.Length; ++i) {
							++cov[i];
						}
					}
				}
				local_reader.Close();
				free(cov);
				local_BAMStatistics[block_id].num_total = num_total;
				local_BAMStatistics[block_id].num_unmapped = num_unmapped;
				local_BAMStatistics[block_id].num_clipped = num_clipped;
				local_BAMStatistics[block_id].max_read_len = max_read_len;
				unmapped_rejected_lists[block_id] = local_unmapped_rejected;

				done_vector[block_id] = 'D';

				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (last-wrong) %s\n")
								% block_id % a_block_boundary_str).str();
					} else {
						cout << (boost::format("[ParallelBamReader.collect_basic_statistics] Block-%d (last) %s\n")
								% block_id % a_block_boundary_str).str();
					}
				}
//			else {
//				size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//				double processed = n/(double)done_vector.size() * 100.0;
//				cout << (boost::format("%.2f %%\n") % processed).str();
//			}
			});
	}
	// sometimes the unaligned reads are in the last portion of BAM file, hence
	// changing the order.
	swap(tasks[0], tasks.back());
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	cout << "[ParallelBamReader.collect_basic_statistics] gather scattered statistics\n";
	int64_t n_soft_clips = 0;
	for (int64_t block_id = calculated_n_blocks - 2; block_id >= 0; --block_id) {
		auto& current_soft_clips = local_softclips[block_id];
		n_soft_clips += current_soft_clips.size();
	}
//	softclips.resize(n_soft_clips);

	// this runs in reverse order in order to correctly assign the softclip information
	for (int64_t block_id = calculated_n_blocks - 2; block_id >= 0; --block_id) {
		num_total += local_BAMStatistics[block_id].num_total;
		num_unmapped += local_BAMStatistics[block_id].num_unmapped;
		num_clipped += local_BAMStatistics[block_id].num_clipped;
		max_read_len = max(local_BAMStatistics[block_id].max_read_len, max_read_len);
		auto& current_soft_clips = local_softclips[block_id];
		softclips.insert(current_soft_clips.begin(), current_soft_clips.end());
//		for (auto itr = current_soft_clips.begin(); current_soft_clips.end() != itr; ++itr) {
//			if (softclips.end() == softclips.find(itr->first)) {
//				softclips[itr->first] = itr->second;
//			}
//		}
	}

	if (verbose) {
//		for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
//			cout << itr->first << ":" << itr->second << "\n";
//		}
		for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
			cout << *itr << "\n";
		}
	}

	cout << options.fname << ": \n   " << num_total << " reads\n   " << num_unmapped << " unaligned\n   " << num_clipped << " clipped (" << softclips.size() << " hashed)\n   " << max_read_len << " bps in the longest read\n";

	report_qualities(quality_sample_space);
//		ofstream scfile(first_stat_filename, ios::binary);
//		scfile << num_total << "\n" //
//				<< num_unmapped << "\n" //
//				<< num_clipped << "\n" //
//				<< max_read_len << "\n" //
//				<< qenc << "\n";
//		for (auto itr = softclips.begin(); softclips.end() != itr; ++itr) {
//			scfile << itr->first << "\t" << static_cast<int64_t>(itr->second) << "\n";
//		}

	boost::unordered_map<string, boost::unordered_map<string, PairedAlignment>> paired_alns;
	tasks.push_back([&, calculated_n_blocks]
	{
		if(!options.processUU) {
			return;
		}
		string rg_name;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1;
				++block_id) {
			auto& local_umbuf_first = umbufs_first[block_id];
			for(auto itr_first = local_umbuf_first.begin(); local_umbuf_first.end() != itr_first; ++itr_first) {
				if(!itr_first->second.GetReadGroup(rg_name)) {
					rg_name = "none";
				}
				paired_alns[rg_name][itr_first->first].pair_1 = itr_first->second;
			}
			local_umbuf_first.clear();
			auto& local_umbuf_second = umbufs_second[block_id];
			for(auto itr_second = local_umbuf_second.begin(); local_umbuf_second.end() != itr_second; ++itr_second) {
				if(!itr_second->second.GetReadGroup(rg_name)) {
					rg_name = "none";
				}
				paired_alns[rg_name][itr_second->first].pair_2 = itr_second->second;
			}
			local_umbuf_second.clear();
		}
		umbufs_first.clear();
		umbufs_second.clear();
		for(auto itr_rg_groups = paired_alns.begin(); paired_alns.end() != itr_rg_groups; ++itr_rg_groups) {
			auto an_itr_map = itr_rg_groups->second;
			for(auto itr_aln = an_itr_map.begin(); an_itr_map.end() != itr_aln; ++itr_aln) {
				itr_aln->second.should_process = ReadGroup::isMatePair(itr_aln->second.pair_1, itr_aln->second.pair_2);
			}
		}
	});

	Histogram uhist;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		unmapped_rejected = unmapped_rejected_lists[block_id];
		uhist += uhist_lists[block_id];
	}
	ofstream umrdist(options.umrdistname, ios::binary);
	uhist.print(umrdist);
	vector<string> umf_file_names;
	vector<string> blist_file_names;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string a_umf_file_name = options.umfname + "." + str_block_id;
		string a_blist_file_name = options.blistname + "." + str_block_id;
		umf_file_names.push_back(a_umf_file_name);
		blist_file_names.push_back(a_blist_file_name);
	}
	// since in the next steps of meerkat, entries in the black list file are stored in a map, hence redundant keys due to the boundary loop
	// would not be a problem. The code below reduces such redundant key, but not so efficient.
//		tasks.push_back([&] {
//			ogzstream blist_final(options.blistname.c_str());
//			string line;
//			string prev_line;
//			for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				istringstream in(castle::IOUtils::read_fully(options.blistname + "." + str_block_id));
//	// at the boundary of each parallel block, some redundant entry would be generated
//	// hence a duplicated entry will be ignored.
//				while(getline(in, line, '\n')) {
//					if(line != prev_line) {
//						blist_final << line << "\n";
//						prev_line = line;
//					}
//				}
//			}
//			blist_final.close();
//		});

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_compress_and_merge(options.umfname, umf_file_names, n_cores, false);
	castle::IOUtils::plain_file_compress_and_merge(options.blistname, blist_file_names, n_cores, true);
	cout << checker;
}

void ParallelBamReader::output_blacklist_file() {
	if (boost::filesystem::exists(options.blistname)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.output_blacklist_file");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	vector<BlockBoundary> actual_blocks = unmapped_included_blocks;
	int64_t calculated_n_blocks = actual_blocks.size();
// only for debugging
	string done_vector(calculated_n_blocks - 1, 'U');
// find read group
//vector<string> read_groups;
//{
//BamTools::BamReader local_reader;
//if (!local_reader.Open(a_path, an_index_path)) {
//return;
//}
//istringstream in(local_reader.GetHeaderText());
//string line;
//const char* delim = "\t:";
//vector<string> a_cols;
//while(getline(in, line, '\n')) {
//if(string::npos == line.find("@RG")) {
//continue;
//}
//castle::StringUtils::c_string_multi_split(line, delim, a_cols);
//if(a_cols.size() < 2) {
//continue;
//}
//read_groups.push_back(a_cols[2]);
//}
//local_reader.Close();
//}
	vector<map<string, ReadGroupAlt>> read_group_maps(calculated_n_blocks - 1);
//max_read_len = 151;
	vector<BlockBoundary> black_list_block_boundaries(calculated_n_blocks);
	for (int64_t block_id = 1; block_id < calculated_n_blocks; ++block_id) {
		black_list_block_boundaries[block_id - 1] = actual_blocks[block_id];
	}
// calculate boundary strings
	cout << "[ParallelBamReader.output_blacklist_file] start finding boundary\n";
	for (int64_t block_id = 1; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}

			BamTools::BamAlignment local_alignment_entry;
			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
			if(-1 == the_current_ref_id && -1 == the_next_ref_id) {
				local_reader.Close();
				return;
			}
//uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
				string str_block_id = boost::lexical_cast<string>(block_id);

//ofstream blist(options.blistname + "." + str_block_id, ios::binary);
				string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//				bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//				if(!jump_success) {
//					cout << (boost::format("[ParallelBamReader.output_blacklist_file] block-%d (Jump fail): (%d/%d)-(%d/%d)\n")
//							% block_id % the_current_ref_id % the_current_ref_pos
//							% the_next_ref_id % the_next_ref_pos).str();
//					local_reader.Close();
//					return;
//				}

				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
//				int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}

				long *cov = 0;
				long covsize = sizeof(long) * max_read_len;
				if (options.coverage_cutoff > 0) {
					cov = (long *) malloc(covsize);
					bzero((char *) cov, covsize);
				}

				long cur_pos = the_current_ref_pos, cur_chrom = the_current_ref_id;
//auto& refnames = local_reader.GetReferenceData();
//auto& local_read_group_map = read_group_maps[block_id];
				bool terminate_search = false;
				int64_t num_total = 0;
				int32_t the_current_max_pos = -1;
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;
				if(cur_offset != the_current_ref_offset) {
				cout << (boost::format("[ParallelBamReader.output_blacklist_file] block-%d (start) offset: jump: %d, calc: %d\n")
				% block_id % cur_offset % the_current_ref_offset).str();
				}

				while (local_reader.LoadNextAlignmentWithName(local_alignment_entry)) {

//if(local_alignment_entry.Name == the_next_block_read_name
//&& local_alignment_entry.RefID == the_next_ref_id
//&& local_alignment_entry.Position == the_next_ref_pos
//&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//) {
//break;
//}
					if( -1 == local_alignment_entry.Position) {
						continue;
					}
					cur_offset = m_bgzf.Tell();
					if(-1 == the_current_max_pos) {
						the_current_max_pos = local_alignment_entry.Position + (max_read_len << 3);
					}

					if(verbose && 0 == num_total) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag).str();
						cout << (boost::format("[ParallelBamReader.output_blacklist_file] Block-%d (boundary first) %s\n")
								% block_id % a_block_boundary_str).str();
					}

					++num_total;
					/* XXX: does not account for CIGAR string data like calDepth */
					/* these 4 Is* conditions match calDepth's default mask */
					if (options.coverage_cutoff > 0 && local_alignment_entry.IsMapped()
							&& local_alignment_entry.IsPrimaryAlignment()
							&& !local_alignment_entry.IsFailedQC() && !local_alignment_entry.IsDuplicate()) {
						if (local_alignment_entry.RefID != cur_chrom && cur_chrom != -1) {
							for (int i = 0; i < max_read_len; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
//string a_boundary_str =
//(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
//% (cur_pos + i + 1) % cov[i]).str();
									black_list_block_boundaries[block_id - 1].read_name = local_alignment_entry.Name;
									black_list_block_boundaries[block_id - 1].ref_id = local_alignment_entry.RefID;
									black_list_block_boundaries[block_id - 1].offset = prev_offset;
									black_list_block_boundaries[block_id - 1].pos = local_alignment_entry.Position;
									black_list_block_boundaries[block_id - 1].aln_flag = local_alignment_entry.AlignmentFlag;
//block_boundaries[block_id - 1] = a_block_boundary_str;
									terminate_search = true;
									break;
								}
							}
						} else if (local_alignment_entry.Position != cur_pos && cur_pos != -1) {
							/* compute the number of positions in our cov
							 * buffer that need to be kicked out */
							int x = min(local_alignment_entry.Position - cur_pos, max_read_len);
							for (int i = 0; i < x; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
//string a_boundary_str =
//(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
//% (cur_pos + i + 1) % cov[i]).str();
									black_list_block_boundaries[block_id - 1].read_name = local_alignment_entry.Name;
									black_list_block_boundaries[block_id - 1].ref_id = local_alignment_entry.RefID;
									black_list_block_boundaries[block_id - 1].offset = prev_offset;
									black_list_block_boundaries[block_id - 1].pos = local_alignment_entry.Position;
									black_list_block_boundaries[block_id - 1].aln_flag = local_alignment_entry.AlignmentFlag;
									terminate_search = true;
									break;
								}
							}

							if(!terminate_search) {
								memmove(cov, cov + x, covsize - sizeof(long) * x);
								bzero((char *) (cov + max_read_len - x), sizeof(long) * x);
							}
						}
						if(terminate_search) {
							break;
						}
						cur_pos = local_alignment_entry.Position;
						cur_chrom = local_alignment_entry.RefID;
						for (int i = 0; i < local_alignment_entry.Length; ++i) {
							++cov[i];
						}
					}
					prev_offset = cur_offset;

					if((-1 != local_alignment_entry.RefID && the_current_ref_id != local_alignment_entry.RefID) ||
							(-1 != local_alignment_entry.Position && the_current_max_pos < local_alignment_entry.Position)) {
						break;
					}
				}

				local_reader.Close();
				free(cov);
				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d<%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% the_current_max_pos
							% local_alignment_entry.AlignmentFlag).str();
					cout << (boost::format("[ParallelBamReader.output_blacklist_file] Block-%d (boundary last) %s\n")
							% block_id % a_block_boundary_str).str();
				}
//done_vector[block_id] = 'D';

//else {
//size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//double processed = n/(double)done_vector.size() * 100.0;
//cout << (boost::format("%.2f %%\n") % processed).str();
//}
			});
	}
//	if(tasks.size() > 0) {
//		swap(tasks[0], tasks.back());
//	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//block_boundaries[calculated_n_blocks - 2].read_name = "last";
	if (verbose) {
		cout << "[ParallelBamReader.output_blacklist_file] print boundaries\n";
		for (uint64_t boundary_id = 0; boundary_id < black_list_block_boundaries.size(); ++boundary_id) {
//if(block_boundaries[boundary_id].empty()) {
//continue;
//}
			cout << (boost::format("[ParallelBamReader.output_blacklist_file] Block-%d (boundary str) %s\n") % boundary_id % black_list_block_boundaries[boundary_id].str()).str();
		}
	}

	cout << "[ParallelBamReader.output_blacklist_file] start writing\n";
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			BamTools::BamAlignment local_alignment_entry;
			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
			int32_t the_current_ref_pos = actual_blocks[block_id].pos;

			int32_t the_next_ref_id = black_list_block_boundaries[block_id].ref_id;
//			int32_t the_next_ref_pos = black_list_block_boundaries[block_id].pos;
//			uint32_t the_next_aln_flag = black_list_block_boundaries[block_id].aln_flag;
			if(-1 == the_current_ref_id && -1 == the_next_ref_id) {
				local_reader.Close();
				return;
			}
				string the_next_block_read_name = black_list_block_boundaries[block_id].read_name;
				string str_block_id = boost::lexical_cast<string>(block_id);

				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
				int64_t the_next_ref_offset = black_list_block_boundaries[block_id].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}

				ofstream blist(options.blistname + "." + str_block_id, ios::binary);

//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[ParallelBamReader.output_blacklist_file] block-%d (Jump fail): (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
				long *cov = 0;
				long covsize = sizeof(long) * max_read_len;
				if (options.coverage_cutoff > 0) {
					cov = (long *) malloc(covsize);
					bzero((char *) cov, covsize);
				}

				long cur_pos = the_current_ref_pos, cur_chrom = the_current_ref_id;
				auto& refnames = local_reader.GetReferenceData();
				auto& local_read_group_map = read_group_maps[block_id];
//int32_t state = 0;
				int64_t num_total = 0;
//int64_t the_current_max_pos = -1;
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;
				if(cur_offset != the_current_ref_offset) {
					cout << (boost::format("[ParallelBamReader.output_blacklist_file] block-%d (start) offset: jump: %d, calc: %d\n")
							% block_id % cur_offset % the_current_ref_offset).str();
				}

				while (local_reader.LoadNextAlignmentWithName(local_alignment_entry)) {
					if(verbose && 0 == num_total) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag).str();
						cout << (boost::format("[ParallelBamReader.output_blacklist_file] Block-%d (writing first) %s\n")
								% block_id % a_block_boundary_str).str();
					}

					++num_total;
					/* XXX: does not account for CIGAR string data like calDepth */
					/* these 4 Is* conditions match calDepth's default mask */
					if (options.coverage_cutoff > 0 && local_alignment_entry.IsMapped()
							&& local_alignment_entry.IsPrimaryAlignment()
							&& !local_alignment_entry.IsFailedQC() && !local_alignment_entry.IsDuplicate()) {
						if (local_alignment_entry.RefID != cur_chrom && cur_chrom != -1) {
							for (int i = 0; i < max_read_len; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
								cov[i] = 0;
							}
						} else if (local_alignment_entry.Position != cur_pos && cur_pos != -1) {
							/* compute the number of positions in our cov
							 * buffer that need to be kicked out */
//if("MT" == refnames[cur_chrom].RefName && 0 == cur_pos) {
//string a_block_boundary_str = (boost::format("%s (%d-%d)->(%d-%d) %d")
//% local_alignment_entry.Name % cur_chrom % cur_pos % local_alignment_entry.RefID % local_alignment_entry.Position
//% local_alignment_entry.AlignmentFlag).str();
//cout << (boost::format("[ParallelBamReader.output_cov_and_readgroup_alt] Block-%d (writing case-2) %s\n")
//% block_id % a_block_boundary_str).str();
//}
							int x = min(local_alignment_entry.Position - cur_pos, max_read_len);
							for (int i = 0; i < x; ++i) {
								if (cov[i] >= options.coverage_cutoff) {
									string a_boundary_str =
									(boost::format("%s\t%d\t%d\n") % refnames[cur_chrom].RefName
											% (cur_pos + i + 1) % cov[i]).str();
									blist << a_boundary_str;
								}
							}
							memmove(cov, cov + x, covsize - sizeof(long) * x);
							bzero((char *) (cov + max_read_len - x), sizeof(long) * x);
						}
						cur_pos = local_alignment_entry.Position;
						cur_chrom = local_alignment_entry.RefID;
						for (int i = 0; i < local_alignment_entry.Length; ++i) {
							++cov[i];
						}
					}
					cur_offset = m_bgzf.Tell();
//					if(local_alignment_entry.RefID == the_next_ref_id
//							&& local_alignment_entry.Position == the_next_ref_pos
//							&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//							&& local_alignment_entry.Name == the_next_block_read_name
//					) {
//						break;
//					}

					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;
///* Determine this alignment's readgroup */
					string rg_name;
					if(!local_alignment_entry.GetReadGroup(rg_name)) {
						rg_name = "none";
					}
					if(black_listed.end() != black_listed.find(rg_name)) {
						continue;
					}
					local_read_group_map[rg_name].witness(local_alignment_entry, options.max_isize, options.isize_samples);
				}
				free(cov);
				blist.close();
				local_reader.Close();
				done_vector[block_id] = 'D';
				if(prev_offset != the_next_ref_offset) {
					cout << (boost::format("[ParallelBamReader.output_blacklist_file] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n") % block_id % prev_offset % cur_offset % the_next_ref_offset).str();
				}
				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					cout << (boost::format("[ParallelBamReader.output_blacklist_file] Block-%d (writing last) %s\n")
							% block_id % a_block_boundary_str).str();
				}
//				else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
			});
	}
//	if(tasks.size() > 0) {
//		swap(tasks[0], tasks.back());
//	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[ParallelBamReader.output_blacklist_file] gathers scattered data\n";
	vector<string> blist_file_names;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string a_blist_file_name = options.blistname + "." + str_block_id;
		blist_file_names.push_back(a_blist_file_name);
	}
	// since in the next steps of meerkat, entries in the black list file are stored in a map, hence redundant keys due to the boundary loop
	// would not be a problem. The code below reduces such redundant key, but not so efficient.
	//		tasks.push_back([&] {
	//			ogzstream blist_final(options.blistname.c_str());
	//			string line;
	//			string prev_line;
	//			for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
	//				string str_block_id = boost::lexical_cast<string>(block_id);
	//				istringstream in(castle::IOUtils::read_fully(options.blistname + "." + str_block_id));
	//	// at the boundary of each parallel block, some redundant entry would be generated
	//	// hence a duplicated entry will be ignored.
	//				while(getline(in, line, '\n')) {
	//					if(line != prev_line) {
	//						blist_final << line << "\n";
	//						prev_line = line;
	//					}
	//				}
	//			}
	//			blist_final.close();
	//		});

	castle::IOUtils::plain_file_compress_and_merge(options.blistname, blist_file_names, n_cores, true);
//	tasks.push_back([&] {
//		ogzstream blist_final(options.blistname.c_str());
//		string line;
//		string prev_line;
//		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			istringstream in(castle::IOUtils::read_fully(options.blistname + "." + str_block_id));
//// at the boundary of each parallel block, some redundant entry would be generated
//// hence a duplicated entry will be ignored.
//			while(getline(in, line, '\n')) {
//				if(line != prev_line) {
//					blist_final << line << "\n";
//					prev_line = line;
//				}
//			}
//		}
//		blist_final.close();
//	});
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	{
//		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
//			tasks.push_back([&, block_id] {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				boost::filesystem::remove(options.blistname + "." + str_block_id);
//			});
//		}
//		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	}
	cout << checker;
}

void ParallelBamReader::find_quality_standard_and_max_read_length(vector<uint64_t>& quality_sample_space) {
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.find_quality_standard_and_max_read_length");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = actual_blocks.size();
	vector<vector<uint64_t>> quality_sample_space_lists(calculated_n_blocks - 1);
	vector<uint64_t> max_read_length_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}

//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				return;
//			}
				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}
				int64_t local_max_read_len = 0;
				string str_block_id = boost::lexical_cast<string>(block_id);
				int64_t n_samples = 0;
				const int64_t N_MAX_SAMPLES = 600000 /(double) actual_blocks.size();
				BamTools::BamAlignment local_alignment_entry;
				auto& local_quality_sample_space = quality_sample_space_lists[block_id];
				local_quality_sample_space.resize(255);
				while (local_reader.LoadNextAlignmentWithName(local_alignment_entry)) {
					witness_qualities(local_alignment_entry, local_quality_sample_space);
					if (local_alignment_entry.Length > local_max_read_len) {
						local_max_read_len = local_alignment_entry.Length;
					}
					++n_samples;
					if(n_samples > N_MAX_SAMPLES) {
						break;
					}
				}
				max_read_length_lists[block_id] = local_max_read_len;
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	uint64_t temp_max_read_len = 0;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		auto& local_quality_sample_spaces = quality_sample_space_lists[block_id];
		for (int64_t quality_id = 0; quality_id < 255; ++quality_id) {
			quality_sample_space[quality_id] += local_quality_sample_spaces[quality_id];
		}
		if (max_read_length_lists[block_id] > temp_max_read_len) {
			temp_max_read_len = max_read_length_lists[block_id];
		}
	}
	max_read_len = temp_max_read_len;
	int a = 0;
	int b = 0;
	for (int i = 0; i < 255; ++i) {
		if (quality_sample_space[i] != 0) {
			a = i;
			break;
		}
	}
	for (int i = 255 - 1; i >= 0; --i) {
		if (quality_sample_space[i] != 0) {
			b = i;
			break;
		}
	}

	/* check min max range */
	set<int> ignoring_standards;

	for (unsigned int i = 0; i < sizeof(ReadGroup::known_qencodings) / sizeof(QENCODING); ++i) {
		if (ReadGroup::known_qencodings[i].min > a) {
			ignoring_standards.insert(i);
			continue;
		}
		if(ReadGroup::known_qencodings[i].max < b) {
			ignoring_standards.insert(i);
		}
	}

	/* Simple guess at the encoding type */
	int best = 0, best_diff = INT_MAX;
	for (unsigned int i = 0; i < sizeof(ReadGroup::known_qencodings) / sizeof(QENCODING); ++i) {
		if(ignoring_standards.end() != ignoring_standards.find(i)) {
			continue;
		}
		int d = abs(ReadGroup::known_qencodings[i].min - a) + abs(ReadGroup::known_qencodings[i].max - b);
		if (d < best_diff) {
			best = i;
			best_diff = d;
		}
	}
	qenc = best;
	cout << checker;
}
void ParallelBamReader::output_unmapped() {
	if (boost::filesystem::exists(options.umfname)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.output_unmapped");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
//	vector<BlockBoundary> unmapped_included_blocks = actual_blocks;
	int64_t calculated_n_blocks = unmapped_included_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');

	vector<int64_t> unmapped_rejected_lists(calculated_n_blocks - 1);
	vector<Histogram> uhist_lists(calculated_n_blocks - 1);
//	vector<map<string, BamAlignment>> umbufs_first(calculated_n_blocks - 1);
//	vector<map<string, BamAlignment>> umbufs_second(calculated_n_blocks - 1);
	cout << (boost::format("[ParallelBamReader.output_unmapped] min frag. size: %d\n") % (options.frag_size << 1)).str();
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t n_alignments = 0;
			int64_t local_unmapped_rejected = 0;
			const int64_t min_frag_size = options.frag_size << 1;
			Histogram local_uhist;
			BamTools::BamAlignment local_alignment_entry;
			int32_t the_current_ref_id = unmapped_included_blocks[block_id].ref_id;
			int32_t the_current_ref_pos = unmapped_included_blocks[block_id].pos;

			int32_t the_next_ref_id = unmapped_included_blocks[block_id + 1].ref_id;
			int32_t the_next_ref_pos = unmapped_included_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
				string str_block_id = boost::lexical_cast<string>(block_id);

				int64_t the_current_ref_offset = unmapped_included_blocks[block_id].offset;
				int64_t the_next_ref_offset = unmapped_included_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}

//			auto& local_umbuf_first = umbufs_first[block_id];
//			auto& local_umbuf_second = umbufs_second[block_id];
				ofstream umfile(options.umfname + "." + str_block_id, ios::binary);
//ofstream split1(options.split1name + "." + str_block_id, ios::binary);
//ofstream split2(options.split2name + "." + str_block_id, ios::binary);
				string the_next_block_read_name = unmapped_included_blocks[block_id + 1].read_name;
//string the_current_block_read_name = actual_blocks[block_id].read_name;
//				bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//				if(!jump_success) {
//					cout << (boost::format("[ParallelBamReader.output_unmapped] block-%d (Jump fail): (%d/%d)-(%d/%d)\n")
//							% block_id % the_current_ref_id % the_current_ref_pos
//							% the_next_ref_id % the_next_ref_pos).str();
//					local_reader.Close();
//					return;
//				}
				if(verbose) {
					cout << (boost::format("[ParallelBamReader.output_unmapped] block-%d (start) (%d/%d)-(%d/%d)\n")
							% block_id % the_current_ref_id % the_current_ref_pos
							% the_next_ref_id % the_next_ref_pos).str();
				}
				string rg_name;
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;
//				if(cur_offset != the_current_ref_offset) {
//					cout << (boost::format("[ParallelBamReader.output_unmapped] block-%d (start) offset: jump: %d, calc: %d\n")
//							% block_id % cur_offset % the_current_ref_offset).str();
//				}
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}

					if(verbose && 0 == n_alignments) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag).str();
						cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (first) %s\n")
								% block_id % a_block_boundary_str).str();
					}

//					if(local_alignment_entry.RefID == the_next_ref_id
//							&& local_alignment_entry.Position == the_next_ref_pos
//							&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//							&& local_alignment_entry.Name == the_next_block_read_name
//					) {
//						break;
//					}
					prev_offset = cur_offset;
					++n_alignments;
					if(!local_alignment_entry.GetReadGroup(rg_name)) {
						rg_name = "none";
					}
					if(black_listed.end() != black_listed.find(rg_name)) {
						continue;
					}

					if(!local_alignment_entry.IsMapped()) {
//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped trimming) before: %d\n")
//% block_id % local_alignment_entry.Length).str();
						ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped trimming) after: %d\n")
//% block_id % local_alignment_entry.Length).str();
						local_uhist.add(local_alignment_entry);
//string a_block_boundary_str = (boost::format("%s %d-%d %d")
//% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//% local_alignment_entry.AlignmentFlag).str();
//cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (unmapped) %s\n")
//% block_id % a_block_boundary_str).str();
						if (local_alignment_entry.Length >= min_frag_size) {
							BamAlignment copy(local_alignment_entry);
							if (local_alignment_entry.IsPaired()) {
								copy.Name += copy.IsFirstMate() ? "_1" : "_2";
							}
							ReadGroup::writeFQ(umfile, copy);
//ReadGroup::split_read(copy, split1, split2, options.frag_size, options.n_cutoff);
						} else {
							++local_unmapped_rejected;
						}
					}
					/* 2. Handle UU pairs */
				}
//				if(prev_offset != the_next_ref_offset) {
//					cout << (boost::format("[ParallelBamReader.output_unmapped] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n")
//							% block_id % prev_offset % cur_offset % the_next_ref_offset).str();
//				}
				local_reader.Close();
				unmapped_rejected_lists[block_id] = local_unmapped_rejected;
				uhist_lists[block_id] = local_uhist;
				done_vector[block_id] = 'D';

				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					cout << (boost::format("[ParallelBamReader.output_unmapped] Block-%d (last) %s\n")
							% block_id % a_block_boundary_str).str();
				}
				else {
					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
					double processed = n/(double)done_vector.size() * 100.0;
					cout << (boost::format("%.2f %%\n") % processed).str();
				}
			});
	}
// sometimes the unaligned reads are in the last portion of BAM file, hence
// changing the order.
//	swap(tasks[0], tasks.back());
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	int64_t local_unmapped_rejected = 0;
	Histogram uhist;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		local_unmapped_rejected += unmapped_rejected_lists[block_id];
		uhist += uhist_lists[block_id];
	}
	unmapped_rejected = local_unmapped_rejected;
	cout << "Final read counts after rejection due to quality (-q " << options.q << " -m " << options.min_read_len << "):\n   unmapped reads: "
	<< (num_unmapped - unmapped_rejected) << " (" << unmapped_rejected << " rejected)\n";
	ofstream umrdist(options.umrdistname, ios::binary);
	uhist.print(umrdist);
	vector<string> file_names;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string file_name = options.umfname + "." + str_block_id;
		file_names.push_back(file_name);
	}
	castle::IOUtils::plain_file_compress_and_merge(options.umfname, file_names, n_cores, true);
	cout << checker;
}

void ParallelBamReader::collect_second_statistics() {
	output_unmapped();
	output_blacklist_file();
// dependency on softclips
	output_softclips();
	output_read_groups_alt();
}

void ParallelBamReader::report_qualities(vector<uint64_t>& sample_space) {
	int a = 0;
	int b = 0;
	for (int i = 0; i < 255; ++i) {
		if (sample_space[i] != 0) {
			a = i;
			break;
		}
	}
	for (int i = 255 - 1; i >= 0; --i) {
		if (sample_space[i] != 0) {
			b = i;
			break;
		}
	}
	if (options.verbose) {
		cout << "Quality score distribution:\n";
		for (int i = 0; i < 255; ++i)
			cout << i << "\t" << sample_space[i] << "\n";
	}

	cout << "   Quality string info:\n";
	cout << "      quality score range: [" << a << ", " << b << "]\n";

	/* check min max range */
	set<int> ignoring_standards;

	for (unsigned int i = 0; i < sizeof(ReadGroup::known_qencodings) / sizeof(QENCODING); ++i) {
		if (ReadGroup::known_qencodings[i].min > a) {
			ignoring_standards.insert(i);
			continue;
		}
		if(ReadGroup::known_qencodings[i].max < b) {
			ignoring_standards.insert(i);
		}
	}
	/* Simple guess at the encoding type */
	int best = 0, best_diff = INT_MAX;
	for (unsigned int i = 0; i < sizeof(ReadGroup::known_qencodings) / sizeof(QENCODING); ++i) {
		if(ignoring_standards.end() != ignoring_standards.find(i)) {
			continue;
		}
		int d = abs(ReadGroup::known_qencodings[i].min - a) + abs(ReadGroup::known_qencodings[i].max - b);
		if (d < best_diff) {
			best = i;
			best_diff = d;
		}
	}

	qenc = best;
	for (int i = 0; i < static_cast<int>(sizeof(ReadGroup::known_qencodings) / sizeof(QENCODING)); ++i) {
		if (qenc == i)
			cout << "    * ";
		else
			cout << "      ";

		cout.width(25);
		cout << left << ReadGroup::known_qencodings[i].name << "[" << ReadGroup::known_qencodings[i].min << ", " << ReadGroup::known_qencodings[i].max << "]\n";
	}
}

/* TEMPLATE
 castle::TimeChecker checker;
 checker.setTarget("ParallelBamReader.output_unmapped");
 checker.start();
 const bool verbose = false;
 vector<function<void()> > tasks;
 string a_path(options.fname);
 string an_index_path(a_path);
 an_index_path += ".bai";
 int64_t calculated_n_blocks = actual_blocks.size();
 string done_vector(calculated_n_blocks - 1, 'U');

 vector<int64_t> unmapped_rejected_lists(calculated_n_blocks - 1);
 for (int64_t block_id = 0; block_id < calculated_n_blocks - 1;
 ++block_id) {
 tasks.push_back(
 [&, block_id] {
 BamTools::BamReader local_reader;
 if (!local_reader.Open(a_path, an_index_path)) {
 return;
 }
 int64_t num_total = 0;
 int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
 int32_t the_current_ref_pos = actual_blocks[block_id].pos;

 int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
 int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
 uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
 string the_next_block_read_name = actual_blocks[block_id + 1].read_name;

 string str_block_id = boost::lexical_cast<string>(block_id);
 auto& local_umbuf_first = unmapped_rejected_lists[block_id];

 bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
 if(!jump_success) {
 cout << (boost::format("[ParallelBamReader.output_unmapped] (Jump fail) Block-%d (%d/%d)-(%d/%d)\n")
 % block_id % the_current_ref_id % the_current_ref_pos
 % the_next_ref_id % the_next_ref_pos).str();
 exit(1);
 }
 if(verbose) {
 cout << (boost::format("[ParallelBamReader.output_unmapped] (start) Block-%d (%d/%d)-(%d/%d)\n")
 % block_id % the_current_ref_id % the_current_ref_pos
 % the_next_ref_id % the_next_ref_pos).str();
 }
 string rg_name;
 BamTools::BamAlignment local_alignment_entry;

 while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
 if(verbose && 0 == num_total) {
 string a_block_boundary_str = (boost::format("%s %d-%d %d")
 % local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
 % local_alignment_entry.AlignmentFlag).str();
 cout << (boost::format("[ParallelBamReader.output_unmapped] (first) Block-%d %s\n")
 % block_id % a_block_boundary_str).str();
 }
 if(local_alignment_entry.RefID == the_next_ref_id
 && local_alignment_entry.Position == the_next_ref_pos
 && local_alignment_entry.AlignmentFlag == the_next_aln_flag
 && local_alignment_entry.Name == the_next_block_read_name
 ) {
 break;
 }
 ++num_total;
 if(!local_alignment_entry.GetReadGroup(rg_name)) {
 rg_name = "none";
 }
 if(black_listed.end() != black_listed.find(rg_name)) {
 continue;
 }

 if(!local_alignment_entry.IsMapped()) {
 ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
 }
 }
 local_reader.Close();
 done_vector[block_id] = 'D';
 if(verbose) {
 string a_block_boundary_str = (boost::format("%s %d-%d %d")
 % local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
 % local_alignment_entry.AlignmentFlag).str();
 cout << (boost::format("[ParallelBamReader.output_raw_softclips] (last) Block-%d %s\n")
 % block_id % a_block_boundary_str).str();
 } else {
 size_t n = count(done_vector.begin(), done_vector.end(), 'D');
 double processed = n/(double)done_vector.size() * 100.0;
 cout << (boost::format("%.2f %%\n") % processed).str();
 }
 });
 }
 castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
 cout << checker;
 */

void ParallelBamReader::output_softclips() {
	if (boost::filesystem::exists(options.clipname) && boost::filesystem::exists(options.split1name) && boost::filesystem::exists(options.split2name)) {
		if (options.generate_mapped_sc_um_file) {
			string mapped_um_name = options.prefix + ".mapped_um.bam";
			string mapped_sc_name = options.prefix + ".mapped_sc.bam";
			if (!options.working_dir.empty()) {
				mapped_um_name = options.working_prefix + ".mapped_um.bam";
				mapped_sc_name = options.working_prefix + ".mapped_sc.bam";
			}
//			if (boost::filesystem::exists(mapped_um_name) && boost::filesystem::exists(mapped_sc_name)) {
//				return;
//			}
		}
//		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.output_softclips");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
//	vector<BlockBoundary> unmapped_included_blocks = actual_blocks;
	int64_t calculated_n_blocks = unmapped_included_blocks.size();

	string done_vector(calculated_n_blocks - 1, 'U');
	boost::unordered_map<string, PairedAlignmentWithTwoIds> paired_alns;
	vector<Histogram> chist_lists(calculated_n_blocks - 1);
	vector<int64_t> clipped_rejected_list(calculated_n_blocks - 1);
	vector<int64_t> unmapped_rejected_list(calculated_n_blocks - 1);
	Histogram chist;
	int64_t unmapped_rejected = 0;
	int64_t clipped_rejected = 0;

//	string target_key = "2_1308_1023_1_0_0_0_0:0:0_0:0:0_cda8f9_ABSL";
	{
		vector<boost::unordered_map<string, vector<BamAlignment>>> albuf_lists(calculated_n_blocks - 1);
//		mutex a_merge_mutex;

		vector<string> read_groups;
		map<string, int64_t> read_groups_reverse_index;
		map<string, int64_t> ref_reverse_index;
		BamTools::BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
		const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
		for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
			auto& a_ref = a_ref_vector[ref_id];
			ref_reverse_index[a_ref.RefName] = ref_id;
		}
		{
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}

			istringstream in(local_reader.GetHeaderText());
			string line;
			const char* delim = "\t:";
			vector<string> a_cols;
			while (getline(in, line, '\n')) {
				if (string::npos == line.find("@RG")) {
					continue;
				}
				castle::StringUtils::c_string_multi_split(line, delim, a_cols);
				if (a_cols.size() < 2) {
					continue;
				}
				read_groups_reverse_index[a_cols[2]] = read_groups.size();
				read_groups.push_back(a_cols[2]);
			}
			string last_key = "none";
			read_groups_reverse_index[last_key] = read_groups.size();
			read_groups.push_back(last_key);
			local_reader.Close();
		}
		n_read_groups = read_groups.size();

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				BamTools::BamReader local_reader;
				if (!local_reader.Open(a_path, an_index_path)) {
					return;
				}
//				if(1 != block_id){
//					return;
//				}
//				if(225 != block_id && 20 != block_id) {
//					return;
//				}
					int64_t num_total = 0;
					int64_t local_clipped_rejected = 0;
					int64_t local_unmapped_rejected = 0;
//					int32_t the_current_ref_id = unmapped_included_blocks[block_id].ref_id;
//					int32_t the_current_ref_pos = unmapped_included_blocks[block_id].pos;
//
//					int32_t the_next_ref_id = unmapped_included_blocks[block_id + 1].ref_id;
//					int32_t the_next_ref_pos = unmapped_included_blocks[block_id + 1].pos;
//					uint32_t the_next_aln_flag = unmapped_included_blocks[block_id + 1].aln_flag;
					string the_next_block_read_name = unmapped_included_blocks[block_id + 1].read_name;

					//int bps = min(options.big_s_bps, options.frag_size);

					string str_block_id = boost::lexical_cast<string>(block_id);

					int64_t the_current_ref_offset = unmapped_included_blocks[block_id].offset;
					int64_t the_next_ref_offset = unmapped_included_blocks[block_id + 1].offset;
					auto& m_bgzf = local_reader.GetBGZF();
					if(0 != block_id) {
						if(!m_bgzf.Seek(the_current_ref_offset)) {
							local_reader.Close();
							return;
						}
					}

//					bool jump_success = local_reader.Jump(unmapped_included_blocks[block_id].ref_id, unmapped_included_blocks[block_id].jump_pos);
//					if(!jump_success) {
//						cout << (boost::format("[ParallelBamReader.output_softclips] (Jump fail) Block-%d (%d/%d)-(%d/%d)\n")
//								% block_id % the_current_ref_id % the_current_ref_pos
//								% the_next_ref_id % the_next_ref_pos).str();
//						local_reader.Close();
//						return;
//					}
//					if(verbose) {
//						cout << (boost::format("[ParallelBamReader.output_softclips] (start) Block-%d (%d/%d)-(%d/%d)\n")
//								% block_id % the_current_ref_id % the_current_ref_pos
//								% the_next_ref_id % the_next_ref_pos).str();
//					}
					string rg_name;
					BamTools::BamAlignment local_alignment_entry;

					ReadGroup rg;
					auto& local_albuf = albuf_lists[block_id];
					auto& local_chist = chist_lists[block_id];
					const int frag_size_doubled = options.frag_size << 1;
//					int bps = min(options.big_s_bps, options.frag_size);
					boost::unordered_set<string> resolved;
//					int32_t decoy_ref_id = -8;
//					auto decoy_itr = ref_reverse_index.find("hs37d5");
//					if(ref_reverse_index.end() != decoy_itr) {
//						decoy_ref_id = decoy_itr->second;
//					}
//					cout << (boost::format("[ParallelBamReader.output_softclips] Decoy Ref Id: %d \n") % decoy_ref_id).str();

					const RefVector& refnames = local_reader.GetReferenceData();
					string mapped_um_name = options.prefix + ".mapped_um.sam." + str_block_id;
					string mapped_sc_name = options.prefix + ".mapped_sc.sam." + str_block_id;
					if(!options.working_dir.empty()) {
						mapped_um_name = options.working_prefix + ".mapped_um.sam." + str_block_id;
						mapped_sc_name = options.working_prefix + ".mapped_sc.sam." + str_block_id;
					}
					ofstream clipfile(options.clipname + "." + str_block_id, ios::binary);
					ofstream split1(options.split1name + "." + str_block_id, ios::binary);
					ofstream split2(options.split2name + "." + str_block_id, ios::binary);
					BamWriter mapped_sc_file;
					BamWriter mapped_um_file;

					if (options.generate_mapped_sc_um_file) {
						if(0 == block_id) {
							if (!mapped_sc_file.SAMOpen(mapped_sc_name,
											local_reader.GetHeaderText(), refnames)) {
								cout << "ERROR: could not open output SAM file '"
								<< mapped_sc_name << "' for writing" << endl;
								exit(1);
							}
							if (!mapped_um_file.SAMOpen(mapped_um_name,
											local_reader.GetHeaderText(), refnames)) {
								cout << "ERROR: could not open output SAM file '"
								<< mapped_um_name << "' for writing" << endl;
								exit(1);
							}
						} else {
							if (!mapped_sc_file.SAMOpenNoHeader(mapped_sc_name, refnames)) {
								cout << "ERROR: could not open output SAM file '"
								<< mapped_sc_name << "' for writing" << endl;
								exit(1);
							}
							if (!mapped_um_file.SAMOpenNoHeader(mapped_um_name, refnames)) {
								cout << "ERROR: could not open output SAM file '"
								<< mapped_um_name << "' for writing" << endl;
								exit(1);
							}
						}
					}
					int64_t cur_offset = m_bgzf.Tell();
					int64_t prev_offset = cur_offset;
//					if(cur_offset != the_current_ref_offset) {
//						cout << (boost::format("[ParallelBamReader.output_softclips] block-%d (start) offset: jump: %d, calc: %d\n")
//								% block_id % cur_offset % the_current_ref_offset).str();
//					}

					while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {

//						if(verbose && 0 == num_total) {
//							string a_block_boundary_str = (boost::format("%s %d-%d %d")
//									% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//									% local_alignment_entry.AlignmentFlag).str();
//							cout << (boost::format("[ParallelBamReader.output_softclips] (first) Block-%d %s\n")
//									% block_id % a_block_boundary_str).str();
//						}

//						if(local_alignment_entry.RefID == the_next_ref_id
//								&& local_alignment_entry.Position == the_next_ref_pos
//								&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//								&& local_alignment_entry.Name == the_next_block_read_name
//						) {
//							break;
//						}
//						auto& cur_name = local_alignment_entry.Name;
						const bool debug = false;
//						const bool debug = cur_name == target_key;
						cur_offset = m_bgzf.Tell();
						if(prev_offset >= the_next_ref_offset) {
							break;
						}
						prev_offset = cur_offset;
//						bool debug = string::npos != local_alignment_entry.Name.find("ST-E00104:502:HFJN5CCXX:7:2122:32400:13527");
//						if(the_next_ref_id == decoy_ref_id && local_alignment_entry.RefID == decoy_ref_id) {
//							break;
//						}
						++num_total;
						if(!local_alignment_entry.GetReadGroup(rg_name)) {
							rg_name = "none";
						}
//						if((decoy_ref_id == local_alignment_entry.RefID || decoy_ref_id == local_alignment_entry.MateRefID)) {
//							continue;
//						}

						if(black_listed.end() != black_listed.find(rg_name) || resolved.end() != resolved.find(local_alignment_entry.Name)) {
							continue;
						}
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] pass blacklisted\n";
						}
						if(!local_alignment_entry.IsMapped()) {
							if(debug) {
								cout << "[ParallelBamReader.output_softclips] unmapped\n";
							}
							ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
							if (local_alignment_entry.Length >= frag_size_doubled) {
								if(debug) {
									cout << "[ParallelBamReader.output_softclips] larger than frag_size_doubled\n";
								}
								BamAlignment copy(local_alignment_entry);
								if (local_alignment_entry.IsPaired()) {
									copy.Name += copy.IsFirstMate() ? "_1" : "_2";
								}
								ReadGroup::split_read(copy, split1, split2, options.frag_size, options.n_cutoff);
							} else {
								++local_unmapped_rejected;
							}
						}
						/* 2. Handle UU pairs */
						if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
							if(debug) {
								cout << "[ParallelBamReader.output_softclips] UU filtered\n";
							}
							continue;
						}
						/* 3a. Hardclipped reads aren't interesting.  In the new BWA-MEM,
						 * softclipped read records are often accompanied by a hard clipped
						 * record, resulting in 3 reads belonging to the same read pair to be
						 * reported in the BAM.  E.g.:
						 * ar3118743    pr1    chr13    18133804    0    48H27M    chr21    14051775    0    CCCAAGAGTGTGGAACTTTCCTTGGGC
						 * ar3118743    pR2    chr21    14051775    12    75M    =    14051906    179    ACTTCACAAAGCTAAGAAATACATATGCATATACCAAGCAAAACACACCATAAGGGTAAAATGATGACTTTTTTG
						 * ar3118743    pr1    chr21    14051906    27    48M27S    =    14051775    -179    TTAATAGGTCCAAATAACAGGTTTATGCTTTTGATTTTGCAGTGGAAGCCCAAGAGTGTGGAACTTTCCTTGGGC
						 */
						if (ReadGroup::isBigH(local_alignment_entry)) {
							if(debug) {
								cout << "[ParallelBamReader.output_softclips] has big H\n";
							}
							continue;
						}
						/* 3. Handle mapped+[softclipped|unmapped] read pairs */
//						int mate_num = 1;
						auto it = softclips.find(local_alignment_entry.Name+"1");
						if(it == softclips.end()) {
							it = softclips.find(local_alignment_entry.Name+"2");
//							mate_num = 2;
						}

//						if(debug) {
//							cout << "here-5\n";
//						}
						if (it != softclips.end()) {
							if(debug) {
								cout << "[ParallelBamReader.output_softclips] has softclip entry\n";
							}
							auto x = local_albuf.find(local_alignment_entry.Name);
							if (x == local_albuf.end()) {
								if(debug) {
									cout << "[ParallelBamReader.output_softclips] first hit for the local hash\n";
								}
								local_albuf[local_alignment_entry.Name].resize(2);
								if(local_alignment_entry.IsFirstMate()) {
									local_albuf[local_alignment_entry.Name][0] = local_alignment_entry;
								} else if(local_alignment_entry.IsSecondMate()) {
									local_albuf[local_alignment_entry.Name][1] = local_alignment_entry;
								}
							} else {
								if(local_alignment_entry.IsFirstMate()) {
									if(local_albuf[local_alignment_entry.Name][0].QueryBases.size() < local_alignment_entry.QueryBases.size()) {
										local_albuf[local_alignment_entry.Name][0] = local_alignment_entry;
									}
								} else if(local_alignment_entry.IsSecondMate()) {
									if(local_albuf[local_alignment_entry.Name][1].QueryBases.size() < local_alignment_entry.QueryBases.size()) {
										local_albuf[local_alignment_entry.Name][1] = local_alignment_entry;
									}
								}
								if(local_albuf[local_alignment_entry.Name][0].QueryBases.size() > 0 && local_albuf[local_alignment_entry.Name][1].QueryBases.size() > 0) {
									if(debug) {
										cout << "[ParallelBamReader.output_softclips] second hit for the local hash\n";
									}
									/* If we're here, then we've got both reads for
									 * some pair that was stored in softclips */
									BamAlignment untrimmed_alignment_entry = local_albuf[local_alignment_entry.Name][0];
									BamAlignment mate = local_albuf[local_alignment_entry.Name][1];
									if(debug) {
										cout << "here-8-cur: " << BamWriter::GetSAMAlignment(local_alignment_entry, refnames) << "\n";
										cout << "here-8-mate: " << BamWriter::GetSAMAlignment(mate, refnames) << "\n";
									}
									ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
									ReadGroup::trim_read(mate, options.q, qenc);

									boost::unordered_map<int8_t, CigarOp> local_map;

									if(local_alignment_entry.CigarData.size() > 0) {
										int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
										// current front
										if('S' == local_alignment_entry.CigarData.front().Type
											&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.front().Length)) >= options.frag_size) {
											if(debug) {
												cout << "[ParallelBamReader.output_softclips] local major front S - 0\n";
											}
											local_map[local_mate_num] = local_alignment_entry.CigarData.front();
										}
										// current back
										if('S' == local_alignment_entry.CigarData.back().Type
											&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.back().Length)) >= options.frag_size) {
											if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > local_alignment_entry.CigarData.back().Length) {
												if(debug) {
													cout << "[ParallelBamReader.output_softclips] local major back S - 0\n";
												}
												local_map[local_mate_num] = local_alignment_entry.CigarData.back();
											}
										}
									}
									if(mate.CigarData.size() > 0) {
										int8_t local_mate_num = mate.IsFirstMate() ? 1 : 2;
										// current front
										if('S' == mate.CigarData.front().Type
											&& min(mate.Length, static_cast<int32_t>(mate.CigarData.front().Length)) >= options.frag_size) {
											if(debug) {
												cout << "[ParallelBamReader.output_softclips] local mate front S - 1\n";
											}
											local_map[local_mate_num] = mate.CigarData.front();
										}
										// current back
										if('S' == mate.CigarData.back().Type
											&& min(mate.Length, static_cast<int32_t>(mate.CigarData.back().Length)) >= options.frag_size) {
											if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > mate.CigarData.back().Length) {
												if(debug) {
													cout << "[ParallelBamReader.output_softclips] local mate back S - 1\n";
												}
												local_map[local_mate_num] = mate.CigarData.back();
											}
										}
									}
									int8_t rep_mate = local_map[1].Length >= local_map[2].Length ? 1 : 2;

	//								ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.front(), bps);
	//								ReadGroup::isBigS(local_alignment_entry, local_alignment_entry.CigarData.back(), bps))

									if (!local_alignment_entry.IsMapped() && mate.IsMapped()) {
										local_chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);
									}
									if (local_alignment_entry.IsMapped() && mate.IsMapped()) {
										if(debug) {
											cout << "here-9-cur: " << BamWriter::GetSAMAlignment(local_alignment_entry, refnames) << "\n";
											cout << "here-9-mate: " << BamWriter::GetSAMAlignment(mate, refnames) << "\n";
											cout << "here-9-rep_mate: " << static_cast<int32_t>(rep_mate) << "\n";
											cout << "here-9-options.big_s_bps: " << options.big_s_bps << "\n";
											cout << "here-9-options.frag_size: " << options.frag_size << "\n";
											cout << "here-9-options.n_cutoff: " << options.n_cutoff << "\n";
										}

										int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
										int8_t mate_mate_num = mate.IsFirstMate() ? 1 : 2;

										BamAlignment local_left = ReadGroup::snip(local_alignment_entry, 0, options.frag_size);
										BamAlignment local_right = ReadGroup::snip(local_alignment_entry, local_alignment_entry.Length - options.frag_size, options.frag_size);
										size_t n_local_amb = count(local_left.QueryBases.begin(), local_left.QueryBases.end(), 'N');
										n_local_amb += count(local_right.QueryBases.begin(), local_right.QueryBases.end(), 'N');
										BamAlignment mate_left = ReadGroup::snip(mate, 0, options.frag_size);
										BamAlignment mate_right = ReadGroup::snip(mate, mate.Length - options.frag_size, options.frag_size);
										size_t n_mate_amb = count(mate_left.QueryBases.begin(), mate_left.QueryBases.end(), 'N');
										n_mate_amb += count(mate_right.QueryBases.begin(), mate_right.QueryBases.end(), 'N');
										if(rep_mate == local_mate_num && n_local_amb > static_cast<size_t>(options.n_cutoff)) {
											rep_mate = mate_mate_num;
										} else if(rep_mate == mate_mate_num && n_mate_amb > static_cast<size_t>(options.n_cutoff)) {
											rep_mate = local_mate_num;
										}

										if(debug) {
											cout << "here-0-local_left: " << local_left.QueryBases << "\n";
											cout << "here-0-local_right: " << local_right.QueryBases << "\n";
											cout << "here-0-mate_left: " << mate_left.QueryBases << "\n";
											cout << "here-0-mate_right: " << mate_right.QueryBases << "\n";
										}
										bool found_zero = false;
										for(auto& a_cigar : untrimmed_alignment_entry.CigarData) {
											if(0 == a_cigar.Length) {
												if(debug) {
													cout << "found zero 0\n";
												}
												found_zero = true;
												break;
											}
										}
										if(!found_zero) {
											for(auto& a_cigar : local_albuf[local_alignment_entry.Name][1].CigarData) {
												if(0 == a_cigar.Length) {
													if(debug) {
														cout << "found zero 1\n";
													}
													found_zero = true;
													break;
												}
											}
										}
										if(!found_zero) {
											for(auto& a_cigar : local_alignment_entry.CigarData) {
												if(0 == a_cigar.Length) {
													if(debug) {
														cout << "found zero 2\n";
													}
													found_zero = true;
													break;
												}
											}
										}
										if(!found_zero) {
											for(auto& a_cigar : mate.CigarData) {
												if(0 == a_cigar.Length) {
													if(debug) {
														cout << "found zero 3\n";
													}
													found_zero = true;
													break;
												}
											}
										}

										if(!found_zero) {
											local_chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);
											if (!rg.recordSCAltNoRG(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate,
													ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? mate : local_alignment_entry, rep_mate, clipfile, split1, split2, options.big_s_bps,
															options.frag_size, options.n_cutoff)) {
												if(debug) {
													cout << "[ParallelBamReader.output_softclips] clipped rejected\n";
												}
												++local_clipped_rejected;
											} else {
												if(options.generate_mapped_sc_um_file) {
													if(debug) {
														cout << "[ParallelBamReader.output_softclips] recorded softclipped\n";
													}
													mapped_sc_file.SaveSAMAlignment(untrimmed_alignment_entry);
													mapped_sc_file.SaveSAMAlignment(local_albuf[local_alignment_entry.Name][1]);
												}
											}
										}
									} else {
										if (!rg.recordUMAltNoRG(local_alignment_entry, mate, options.big_s_bps, options.n_cutoff)) {
											if(debug) {
												cout << "[ParallelBamReader.output_softclips] not recorded UM\n";
											}
											++local_clipped_rejected;
										}
										if (options.generate_mapped_sc_um_file &&
												local_alignment_entry.Length >= frag_size_doubled &&
												mate.Length >= frag_size_doubled) {
											if(debug) {
												cout << "[ParallelBamReader.output_softclips] recorded UM\n";
											}
											mapped_um_file.SaveSAMAlignment(untrimmed_alignment_entry);
											mapped_um_file.SaveSAMAlignment(local_albuf[local_alignment_entry.Name][1]);
										}
									}
									resolved.insert(local_alignment_entry.Name);
									local_albuf.erase(local_alignment_entry.Name);
								}
							}
						}
					}
					local_reader.Close();
					if (options.generate_mapped_sc_um_file) {
						mapped_sc_file.Close();
						mapped_um_file.Close();
					}
					clipped_rejected_list[block_id] = local_clipped_rejected;
					unmapped_rejected_list[block_id] = local_unmapped_rejected;
//					if(prev_offset != the_next_ref_offset) {
//						cout << (boost::format("[ParallelBamReader.output_softclips] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n") % block_id % prev_offset % cur_offset % the_next_ref_offset).str();
//					}

					done_vector[block_id] = 'D';
//					cout << (boost::format("[ParallelBamReader.output_softclips] (last) Block-%d: al_buf(%d)\n")
//							% block_id % local_albuf.size()).str();

					if(verbose) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d (%d)")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag % local_clipped_rejected).str();
						cout << (boost::format("[ParallelBamReader.output_softclips] (last) Block-%d %s\n")
								% block_id % a_block_boundary_str).str();
					}
					else {
						size_t n = count(done_vector.begin(), done_vector.end(), 'D');
						double processed = n/(double)done_vector.size() * 100.0;
						cout << (boost::format("%.2f %%\n") % processed).str();
					}
				});
		}

//		int64_t limit_swap = min(static_cast<int64_t>(10), static_cast<int64_t>(tasks.size()));
//		for (int64_t t_id = 0; t_id < limit_swap; ++t_id) {
//			swap(tasks[t_id], *(tasks.rbegin() + t_id));
//		}

		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			chist += chist_lists[block_id];
			unmapped_rejected += unmapped_rejected_list[block_id];
			clipped_rejected += clipped_rejected_list[block_id];
		}

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_albuf = albuf_lists[block_id];
			for (auto itr_aln = local_albuf.begin(); local_albuf.end() != itr_aln; ++itr_aln) {
				if(paired_alns[itr_aln->first].pair_1.QueryBases.size() < itr_aln->second[0].QueryBases.size()) {
					paired_alns[itr_aln->first].pair_1 = itr_aln->second[0];
					paired_alns[itr_aln->first].id_1 = block_id;
				}
				if(paired_alns[itr_aln->first].pair_2.QueryBases.size() < itr_aln->second[1].QueryBases.size()) {
					paired_alns[itr_aln->first].pair_2 = itr_aln->second[1];
					paired_alns[itr_aln->first].id_2 = block_id;
				}
			}
			local_albuf.clear();
		}
		cout << (boost::format("[ParallelBamReader.output_softclips] # uncommitted: %d\n") % paired_alns.size()).str();
	}
	{
		vector<shared_ptr<ofstream>> clipfiles;
		vector<shared_ptr<ofstream>> split1files;
		vector<shared_ptr<ofstream>> split2files;

		vector<shared_ptr<stringstream>> strstreams_clip;
		vector<shared_ptr<stringstream>> strstreams_split1;
		vector<shared_ptr<stringstream>> strstreams_split2;

		vector<uint64_t> clip_buf_sizes(calculated_n_blocks - 1);
		vector<uint64_t> split1_buf_sizes(calculated_n_blocks - 1);
		vector<uint64_t> split2_buf_sizes(calculated_n_blocks - 1);

		clipfiles.reserve(calculated_n_blocks - 1);
		split1files.reserve(calculated_n_blocks - 1);
		split2files.reserve(calculated_n_blocks - 1);

		strstreams_clip.reserve(calculated_n_blocks - 1);
		strstreams_split1.reserve(calculated_n_blocks - 1);
		strstreams_split2.reserve(calculated_n_blocks - 1);

//		const int frag_size_doubled = options.frag_size << 1;
		BamTools::BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
//		const RefVector& refnames = local_reader.GetReferenceData();

		// write additional entries to clipfile, split1 file and split2 files

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string clipfile_name(options.clipname + "." + str_block_id);
			string split1file_name(options.split1name + "." + str_block_id);
			string split2file_name(options.split2name + "." + str_block_id);
			clipfiles.push_back(make_shared<ofstream>(clipfile_name, ios::binary | ios::app));
			split1files.push_back(make_shared<ofstream>(split1file_name, ios::binary | ios::app));
			split2files.push_back(make_shared<ofstream>(split2file_name, ios::binary | ios::app));
			strstreams_clip.push_back(make_shared<stringstream>());
			strstreams_split1.push_back(make_shared<stringstream>());
			strstreams_split2.push_back(make_shared<stringstream>());
		}
		local_reader.Close();
		string rg_name;
		//int frag_size_doubled = options.frag_size << 1;
//		int bps = min(options.big_s_bps, options.frag_size);
		ReadGroup rg;


		auto pair_aln_itr = paired_alns.begin();

		while (paired_alns.end() != pair_aln_itr) {
			//for (auto pair_aln_itr = paired_alns.begin();
			//paired_alns.end() != pair_aln_itr; ++pair_aln_itr) {
			if (-1 == pair_aln_itr->second.id_1 || -1 == pair_aln_itr->second.id_2) {
				++pair_aln_itr;
				continue;
			}

			int64_t writing_block_id = max(pair_aln_itr->second.id_1, pair_aln_itr->second.id_2);
			auto& local_alignment_entry = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_1 : pair_aln_itr->second.pair_2;
//			const bool debug = local_alignment_entry.Name == target_key;
			const bool debug = false;
			if (!local_alignment_entry.GetReadGroup(rg_name)) {
				rg_name = "none";
			}

			if (black_listed.end() != black_listed.find(rg_name)) {
				if(debug) {
					cout << (boost::format("[ParallelBamReader.output_softclips] continue due to rg_name: %s\n") % rg_name).str();
				}
				++pair_aln_itr;
				continue;
			}

			/* 2. Handle UU pairs */
			if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] continue due to !al.IsMapped(), !al.IsMateMapped(), options.processUU\n";
				}
				++pair_aln_itr;
				continue;
			}

			auto& mate = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_2 : pair_aln_itr->second.pair_1;
//			int mate_num = 1;
			auto it = softclips.find(local_alignment_entry.Name + "1");
			if (softclips.end() == it) {
				it = softclips.find(local_alignment_entry.Name + "2");
//				mate_num = 2;
			}
			if (softclips.end() == it) {
				++pair_aln_itr;
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] continue since the iterator is missing in the softclips\n";
				}
				continue;
			}
			//bool debug = ("ST-E00104:502:HFJN5CCXX:7:2118:2615:30527"
			//== local_alignment_entry.Name);

			stringstream clipfile;
			stringstream split1;
			stringstream split2;
//			ofstream& clipfile = *clipfiles[writing_block_id].get();
//			ofstream& split1 = *split1files[writing_block_id].get();
//			ofstream& split2 = *split2files[writing_block_id].get();

			BamAlignment untrimmed_alignment_entry = local_alignment_entry;
			BamAlignment untrimmed_mate = mate;
			ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
			ReadGroup::trim_read(mate, options.q, qenc);
			boost::unordered_map<int8_t, CigarOp> local_map;

			if(local_alignment_entry.CigarData.size() > 0) {
				int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
				if(debug) {
					cout << (boost::format("[ParallelBamReader.output_softclips] local_mate mate num: %d\n") % local_mate_num).str();
				}
				// current front
				if('S' == local_alignment_entry.CigarData.front().Type
					&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.front().Length)) >= options.frag_size) {
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] front S-0\n";
					}
					local_map[local_mate_num] = local_alignment_entry.CigarData.front();
				}
				// current back
				if('S' == local_alignment_entry.CigarData.back().Type
					&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.back().Length)) >= options.frag_size) {
					if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > local_alignment_entry.CigarData.back().Length) {
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] back S-0\n";
						}
						local_map[local_mate_num] = local_alignment_entry.CigarData.back();
					}
				}
			}
			if(mate.CigarData.size() > 0) {
				int8_t local_mate_num = mate.IsFirstMate() ? 1 : 2;
				// current front
				if('S' == mate.CigarData.front().Type
					&& min(mate.Length, static_cast<int32_t>(mate.CigarData.front().Length)) >= options.frag_size) {
					local_map[local_mate_num] = mate.CigarData.front();
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] front S-1\n";
					}
				}
				// current back
				if('S' == mate.CigarData.back().Type
					&& min(mate.Length, static_cast<int32_t>(mate.CigarData.back().Length)) >= options.frag_size) {
					if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > mate.CigarData.back().Length) {
						local_map[local_mate_num] = mate.CigarData.back();
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] back S-1\n";
						}
					}
				}
			}
			int8_t rep_mate = local_map[1].Length >= local_map[2].Length ? 1 : 2;
//			uint32_t max_soft_clipped_len = 0;
//			for(auto& an_op : local_map) {
//				if(an_op.second.Length > max_soft_clipped_len) {
//					rep_mate = an_op.first;
//					max_soft_clipped_len = an_op.second.Length;
//				}
//			}

			if (!(local_alignment_entry.IsMapped() && mate.IsMapped())) {
				chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);
			}


			if (local_alignment_entry.IsMapped() && mate.IsMapped()) {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] al.IsMapped() && mate.IsMapped()\n";
				}
				int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
				int8_t mate_mate_num = mate.IsFirstMate() ? 1 : 2;

				BamAlignment local_left = ReadGroup::snip(local_alignment_entry, 0, options.frag_size);
				BamAlignment local_right = ReadGroup::snip(local_alignment_entry, local_alignment_entry.Length - options.frag_size, options.frag_size);
				size_t n_local_amb = count(local_left.QueryBases.begin(), local_left.QueryBases.end(), 'N');
				n_local_amb += count(local_right.QueryBases.begin(), local_right.QueryBases.end(), 'N');
				BamAlignment mate_left = ReadGroup::snip(mate, 0, options.frag_size);
				BamAlignment mate_right = ReadGroup::snip(mate, mate.Length - options.frag_size, options.frag_size);
				size_t n_mate_amb = count(mate_left.QueryBases.begin(), mate_left.QueryBases.end(), 'N');
				n_mate_amb += count(mate_right.QueryBases.begin(), mate_right.QueryBases.end(), 'N');
				if(rep_mate == local_mate_num && n_local_amb > static_cast<size_t>(options.n_cutoff)) {
					rep_mate = mate_mate_num;
				} else if(rep_mate == mate_mate_num && n_mate_amb > static_cast<size_t>(options.n_cutoff)) {
					rep_mate = local_mate_num;
				}
				chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);
				if (!rg.recordSCAltNoRG(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate,
						ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? mate : local_alignment_entry, rep_mate, clipfile, split1, split2, options.big_s_bps, options.frag_size, options.n_cutoff)) {
					++clipped_rejected;
				}
			} else {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] !al.IsMapped() || !mate.IsMapped()\n";
				}
				if (!rg.recordUMAltNoRG(local_alignment_entry, mate, options.big_s_bps, options.n_cutoff)) {
					++clipped_rejected;
				}
			}
			++pair_aln_itr;
			string clip_str = clipfile.str();
			string split1_str = split1.str();
			string split2_str = split2.str();
			stringstream& clipfile_mem = *strstreams_clip[writing_block_id].get();
			stringstream& split1_mem = *strstreams_split1[writing_block_id].get();
			stringstream& split2_mem = *strstreams_split2[writing_block_id].get();
			clip_buf_sizes[writing_block_id] += clip_str.size();
			split1_buf_sizes[writing_block_id] += split1_str.size();
			split2_buf_sizes[writing_block_id] += split2_str.size();
			clipfile_mem << clip_str;
			split1_mem << split1_str;
			split2_mem << split2_str;
			if (clip_buf_sizes[writing_block_id] > castle::IOUtils::SMALL_WRITE_BUFFER_SIZE) {
				ofstream& the_clipfile = *clipfiles[writing_block_id].get();
				the_clipfile << clipfile_mem.str();
				clipfile_mem.str(string());
				clip_buf_sizes[writing_block_id] = 0;
			}
			if (split1_buf_sizes[writing_block_id] > castle::IOUtils::SMALL_WRITE_BUFFER_SIZE) {
				ofstream& the_split1_file = *split1files[writing_block_id].get();
				the_split1_file << split1_mem.str();
				split1_mem.str(string());
				split1_buf_sizes[writing_block_id] = 0;
			}
			if (split1_buf_sizes[writing_block_id] > castle::IOUtils::SMALL_WRITE_BUFFER_SIZE) {
				ofstream& the_split2_file = *split2files[writing_block_id].get();
				the_split2_file << split2_mem.str();
				split2_mem.str(string());
				split2_buf_sizes[writing_block_id] = 0;
			}
		}
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			stringstream& clipfile_mem = *strstreams_clip[block_id].get();
			stringstream& split1_mem = *strstreams_split1[block_id].get();
			stringstream& split2_mem = *strstreams_split2[block_id].get();
			string clip_str = clipfile_mem.str();
			string split1_str = split1_mem.str();
			string split2_str = split2_mem.str();

			if (!clip_str.empty()) {
				ofstream& the_clipfile = *clipfiles[block_id].get();
				the_clipfile << clip_str;
			}
			if (!split1_str.empty()) {
				ofstream& the_split1_file = *split1files[block_id].get();
				the_split1_file << split1_str;
			}
			if (!split2_str.empty()) {
				ofstream& the_split2_file = *split2files[block_id].get();
				the_split2_file << split2_str;
			}
		}
	}

	{
		// write additional entries to mapped_sc and mapped_um
		vector<shared_ptr<BamWriter>> mapped_sc_files;
		vector<shared_ptr<BamWriter>> mapped_um_files;

		mapped_sc_files.reserve(calculated_n_blocks - 1);
		mapped_um_files.reserve(calculated_n_blocks - 1);


		const int frag_size_doubled = options.frag_size << 1;
		BamTools::BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
		const RefVector& refnames = local_reader.GetReferenceData();

		// write additional entries to clipfile, split1 file and split2 files

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			mapped_sc_files.push_back(make_shared<BamWriter>());
			mapped_um_files.push_back(make_shared<BamWriter>());
			auto& mapped_sc_file = *mapped_sc_files.back().get();
			auto& mapped_um_file = *mapped_um_files.back().get();
			string mapped_um_name = options.prefix + ".mapped_um.rest.sam." + str_block_id;
			string mapped_sc_name = options.prefix + ".mapped_sc.rest.sam." + str_block_id;
			if (!options.working_dir.empty()) {
				mapped_um_name = options.working_prefix + ".mapped_um.rest.sam." + str_block_id;
				mapped_sc_name = options.working_prefix + ".mapped_sc.rest.sam." + str_block_id;
			}
			if (options.generate_mapped_sc_um_file) {
				if (!mapped_sc_file.SAMOpenNoHeader(mapped_sc_name, refnames)) {
					cout << "ERROR: could not open output SAM file '" << mapped_sc_name << "' for writing" << endl;
					exit(1);
				}
				if (!mapped_um_file.SAMOpenNoHeader(mapped_um_name, refnames)) {
					cout << "ERROR: could not open output SAM file '" << mapped_um_name << "' for writing" << endl;
					exit(1);
				}
			}
		}
		local_reader.Close();
		string rg_name;
		//int frag_size_doubled = options.frag_size << 1;
//		int bps = min(options.big_s_bps, options.frag_size);
		ReadGroup rg;

		auto pair_aln_itr = paired_alns.begin();
		while (paired_alns.end() != pair_aln_itr) {
			//for (auto pair_aln_itr = paired_alns.begin();
			//paired_alns.end() != pair_aln_itr; ++pair_aln_itr) {
			if (-1 == pair_aln_itr->second.id_1 || -1 == pair_aln_itr->second.id_2) {
				++pair_aln_itr;
				continue;
			}

			int64_t writing_block_id = max(pair_aln_itr->second.id_1, pair_aln_itr->second.id_2);
			auto& local_alignment_entry = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_1 : pair_aln_itr->second.pair_2;

//			const bool debug = local_alignment_entry.Name == target_key;
			const bool debug = false;
			if (!local_alignment_entry.GetReadGroup(rg_name)) {
				rg_name = "none";
			}

			if (black_listed.end() != black_listed.find(rg_name)) {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] 2 - continue due to missing id 1 or 2\n";
				}
				++pair_aln_itr;
				continue;
			}

			/* 2. Handle UU pairs */
			if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] 2 - continue since al.IsMapped(), al.IsMateMapped, options.processUU\n";
				}
				++pair_aln_itr;
				continue;
			}

			auto& mate = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_2 : pair_aln_itr->second.pair_1;
//			int mate_num = 1;
			auto it = softclips.find(local_alignment_entry.Name + "1");
			if (softclips.end() == it) {
				it = softclips.find(local_alignment_entry.Name + "2");
//				mate_num = 2;
			}
			if (softclips.end() == it) {
				if(debug) {
					cout << "[ParallelBamReader.output_softclips] 2 - missing softclips\n";
				}
				++pair_aln_itr;
				continue;
			}
			//bool debug = ("ST-E00104:502:HFJN5CCXX:7:2118:2615:30527"
			//== local_alignment_entry.Name);

			stringstream clipfile;
			stringstream split1;
			stringstream split2;
			BamWriter& mapped_sc_file = *mapped_sc_files[writing_block_id].get();
			BamWriter& mapped_um_file = *mapped_um_files[writing_block_id].get();

			BamAlignment untrimmed_alignment_entry = local_alignment_entry;
			BamAlignment untrimmed_mate = mate;
			ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
			ReadGroup::trim_read(mate, options.q, qenc);
			boost::unordered_map<int8_t, CigarOp> local_map;

			if(local_alignment_entry.CigarData.size() > 0) {
				int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
				// current front
				if('S' == local_alignment_entry.CigarData.front().Type
					&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.front().Length)) >= options.frag_size) {
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] 2 - front S - 1\n";
					}
					local_map[local_mate_num] = local_alignment_entry.CigarData.front();
				}
				// current back
				if('S' == local_alignment_entry.CigarData.back().Type
					&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.back().Length)) >= options.frag_size) {
					if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > local_alignment_entry.CigarData.back().Length) {
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] 2 - back S - 1\n";
						}
						local_map[local_mate_num] = local_alignment_entry.CigarData.back();
					}
				}
			}
			if(mate.CigarData.size() > 0) {
				int8_t local_mate_num = mate.IsFirstMate() ? 1 : 2;
				// current front
				if('S' == mate.CigarData.front().Type
					&& min(mate.Length, static_cast<int32_t>(mate.CigarData.front().Length)) >= options.frag_size) {
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] 2 - front S - 2\n";
					}
					local_map[local_mate_num] = mate.CigarData.front();
				}
				// current back
				if('S' == mate.CigarData.back().Type
					&& min(mate.Length, static_cast<int32_t>(mate.CigarData.back().Length)) >= options.frag_size) {
					if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > mate.CigarData.back().Length) {
						local_map[local_mate_num] = mate.CigarData.back();
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] 2 - back S - 2\n";
						}
					}
				}
			}
			int8_t rep_mate = local_map[1].Length >= local_map[2].Length ? 1 : 2;
//			uint32_t max_soft_clipped_len = 0;
//			for(auto& an_op : local_map) {
//				if(an_op.second.Length > max_soft_clipped_len) {
//					rep_mate = an_op.first;
//					max_soft_clipped_len = an_op.second.Length;
//				}
//			}

//			if (!(local_alignment_entry.IsMapped() && mate.IsMapped())) {
//				chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);
//			}

			if (local_alignment_entry.IsMapped() && mate.IsMapped()) {
				int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
				int8_t mate_mate_num = mate.IsFirstMate() ? 1 : 2;

				BamAlignment local_left = ReadGroup::snip(local_alignment_entry, 0, options.frag_size);
				BamAlignment local_right = ReadGroup::snip(local_alignment_entry, local_alignment_entry.Length - options.frag_size, options.frag_size);
				size_t n_local_amb = count(local_left.QueryBases.begin(), local_left.QueryBases.end(), 'N');
				n_local_amb += count(local_right.QueryBases.begin(), local_right.QueryBases.end(), 'N');
				BamAlignment mate_left = ReadGroup::snip(mate, 0, options.frag_size);
				BamAlignment mate_right = ReadGroup::snip(mate, mate.Length - options.frag_size, options.frag_size);
				size_t n_mate_amb = count(mate_left.QueryBases.begin(), mate_left.QueryBases.end(), 'N');
				n_mate_amb += count(mate_right.QueryBases.begin(), mate_right.QueryBases.end(), 'N');
				if(rep_mate == local_mate_num && n_local_amb > static_cast<size_t>(options.n_cutoff)) {
					rep_mate = mate_mate_num;
				} else if(rep_mate == mate_mate_num && n_mate_amb > static_cast<size_t>(options.n_cutoff)) {
					rep_mate = local_mate_num;
				}
//				chist.add(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate);

				bool found_zero = false;
				for(auto& a_cigar : untrimmed_alignment_entry.CigarData) {
					if(0 == a_cigar.Length) {
						found_zero = true;
						break;
					}
				}
				if(!found_zero) {
					for(auto& a_cigar : untrimmed_mate.CigarData) {
						if(0 == a_cigar.Length) {
							found_zero = true;
							break;
						}
					}
				}
				if(!found_zero) {
					for(auto& a_cigar : local_alignment_entry.CigarData) {
						if(0 == a_cigar.Length) {
							found_zero = true;
							break;
						}
					}
				}
				if(!found_zero) {
					for(auto& a_cigar : mate.CigarData) {
						if(0 == a_cigar.Length) {
							found_zero = true;
							break;
						}
					}
				}
				if(debug) {
					if(found_zero) {
						cout << "[ParallelBamReader.output_softclips] found zero\n";
					} else {
						cout << "[ParallelBamReader.output_softclips] not found zero\n";
					}
				}
				if(!found_zero) {
					if (!rg.recordSCAltNoRG(ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate,
							ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? mate : local_alignment_entry, rep_mate, clipfile, split1, split2, options.big_s_bps, options.frag_size, options.n_cutoff)) {
						;
						if(debug) {
							cout << "[ParallelBamReader.output_softclips] not recorded softclipped reads\n";
						}
	//					++clipped_rejected;
					} else {
						if (options.generate_mapped_sc_um_file) {
							if(debug) {
								cout << "[ParallelBamReader.output_softclips] not recorded softclipped reads\n";
							}
							mapped_sc_file.SaveSAMAlignment(untrimmed_alignment_entry);
							mapped_sc_file.SaveSAMAlignment(untrimmed_mate);
						}
					}
				}
			} else {
				if (!rg.recordUMAltNoRG(local_alignment_entry, mate, options.big_s_bps, options.n_cutoff)) {
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] not recorded unmapped reads\n";
					}
					;
				}
				if (options.generate_mapped_sc_um_file && local_alignment_entry.Length >= frag_size_doubled && mate.Length >= frag_size_doubled) {
					if(debug) {
						cout << "[ParallelBamReader.output_softclips] recorded unmapped reads\n";
					}
					mapped_um_file.SaveSAMAlignment(untrimmed_alignment_entry);
					mapped_um_file.SaveSAMAlignment(untrimmed_mate);
				}
			}
			pair_aln_itr = paired_alns.erase(pair_aln_itr);
		}
	}
	{
		vector<string> file_names;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string a_file_name(options.clipname + "." + str_block_id);
			file_names.push_back(a_file_name);
		}
		castle::IOUtils::plain_file_compress_and_merge(options.clipname, file_names, n_cores, true);
	}
	{
		vector<string> file_names;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string a_file_name(options.split1name + "." + str_block_id);
			file_names.push_back(a_file_name);
		}
		// the partial files will be used in the alg stage.
		castle::IOUtils::plain_file_compress_and_merge(options.split1name, file_names, n_cores, true);
	}
	{
		vector<string> file_names;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string a_file_name(options.split2name + "." + str_block_id);
			file_names.push_back(a_file_name);
		}
		castle::IOUtils::plain_file_compress_and_merge(options.split2name, file_names, n_cores, true);
	}
	if (options.generate_mapped_sc_um_file) {
		string mapped_um_sam_name = options.prefix + ".mapped_um.sam";
		string mapped_sc_sam_name = options.prefix + ".mapped_sc.sam";

//		string mapped_um_tmp_name = options.prefix + ".mapped_um.tmp.bam";
//		string mapped_sc_tmp_name = options.prefix + ".mapped_sc.tmp.bam";

		string mapped_um_bam_name = options.prefix + ".mapped_um.bam";
		string mapped_sc_bam_name = options.prefix + ".mapped_sc.bam";

		string mapped_um_tmp_bam_name = options.prefix + ".mapped_um.tmp.bam";
		string mapped_sc_tmp_bam_name = options.prefix + ".mapped_sc.tmp.bam";
		if (!options.working_dir.empty()) {
			mapped_um_sam_name = options.working_prefix + ".mapped_um.sam";
			mapped_sc_sam_name = options.working_prefix + ".mapped_sc.sam";

			mapped_um_tmp_bam_name = options.working_prefix + ".mapped_um.tmp.bam";
			mapped_sc_tmp_bam_name = options.working_prefix + ".mapped_sc.tmp.bam";

			mapped_um_bam_name = options.working_prefix + ".mapped_um.bam";
			mapped_sc_bam_name = options.working_prefix + ".mapped_sc.bam";

		}
		vector<string> mapped_um_names(calculated_n_blocks - 1);
		vector<string> mapped_sc_names(calculated_n_blocks - 1);

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				string str_block_id = boost::lexical_cast<string>(block_id);
				string local_mapped_um_name = options.prefix + ".mapped_um.sam." + str_block_id;
				string local_mapped_um_rest_name = options.prefix + ".mapped_um.rest.sam." + str_block_id;
				string local_mapped_sc_name = options.prefix + ".mapped_sc.sam." + str_block_id;
				string local_mapped_sc_rest_name = options.prefix + ".mapped_sc.rest.sam." + str_block_id;

//				string local_mapped_um_bam_name = options.prefix + ".mapped_um.bam." + str_block_id;
//				string local_mapped_um_sorted_bam_name = options.prefix + ".mapped_um.sorted.bam." + str_block_id;
//				string local_mapped_sc_bam_name = options.prefix + ".mapped_sc.bam." + str_block_id;
//				string local_mapped_sc_sorted_bam_name = options.prefix + ".mapped_sc.sorted.bam." + str_block_id;
					if(!options.working_dir.empty()) {
						local_mapped_um_name = options.working_prefix + ".mapped_um.sam." + str_block_id;
						local_mapped_um_rest_name = options.working_prefix + ".mapped_um.rest.sam." + str_block_id;
						local_mapped_sc_name = options.working_prefix + ".mapped_sc.sam." + str_block_id;
						local_mapped_sc_rest_name = options.working_prefix + ".mapped_sc.rest.sam." + str_block_id;

//					local_mapped_um_bam_name = options.working_prefix + ".mapped_um.bam." + str_block_id;
//					local_mapped_um_sorted_bam_name = options.working_prefix + ".mapped_um.sorted.bam." + str_block_id;
//					local_mapped_sc_bam_name = options.working_prefix + ".mapped_sc.bam." + str_block_id;
//					local_mapped_sc_sorted_bam_name = options.working_prefix + ".mapped_sc.sorted.bam." + str_block_id;
					}
					string line;
					{
						ofstream out_um_sam(local_mapped_um_name, ios::binary | ios::app);
						ifstream in_um_sam(local_mapped_um_rest_name, ios::binary);
						while(getline(in_um_sam, line, '\n')) {
//						if('@' == line[0]) {
//							continue;
//						}
							out_um_sam << line << "\n";
						}
					}
					{
						ofstream out_sc_sam(local_mapped_sc_name, ios::binary | ios::app);
						ifstream in_sc_sam(local_mapped_sc_rest_name, ios::binary);
						while(getline(in_sc_sam, line, '\n')) {
//						if('@' == line[0]) {
//							continue;
//						}
							out_sc_sam << line << "\n";
						}
					}
					boost::filesystem::remove(local_mapped_um_rest_name);
					boost::filesystem::remove(local_mapped_sc_rest_name);
//				string sam_to_bam_um_cmd = (boost::format("samtools view -o %s -Sb %s") % local_mapped_um_bam_name % local_mapped_um_name).str();
//				string sam_to_bam_sc_cmd = (boost::format("samtools view -o %s -Sb %s") % local_mapped_sc_bam_name % local_mapped_sc_name).str();
//				string um_sort_cmd = (boost::format("sambamba sort -l 1 -t 1 -o %s %s") % local_mapped_um_sorted_bam_name % local_mapped_um_bam_name).str();
//				string sc_sort_cmd = (boost::format("sambamba sort -l 1 -t 1 -o %s %s") % local_mapped_sc_sorted_bam_name % local_mapped_sc_bam_name).str();
//				system(sam_to_bam_um_cmd.c_str());
//				system(sam_to_bam_sc_cmd.c_str());
//				system(um_sort_cmd.c_str());
//				system(sc_sort_cmd.c_str());
//				mapped_um_names[block_id] = local_mapped_um_sorted_bam_name;
//				mapped_sc_names[block_id] = local_mapped_sc_sorted_bam_name;
					mapped_um_names[block_id] = local_mapped_um_name;
					mapped_sc_names[block_id] = local_mapped_sc_name;
				});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

		castle::IOUtils::plain_file_merge(mapped_um_sam_name, mapped_um_names, n_cores, true);
		castle::IOUtils::plain_file_merge(mapped_sc_sam_name, mapped_sc_names, n_cores, true);
//		tasks.push_back([&] {
		string sam_to_bam_um_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % mapped_um_tmp_bam_name % mapped_um_sam_name).str();
		system(sam_to_bam_um_cmd.c_str());
//		});
//		tasks.push_back([&] {
		string sam_to_bam_sc_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % mapped_sc_tmp_bam_name % mapped_sc_sam_name).str();
		system(sam_to_bam_sc_cmd.c_str());
//		});
//		tasks.push_back([&]{
//			ofstream out_um_sam(mapped_um_sam_name, ios::binary);
//			string local_mapped_um_name = options.prefix + ".mapped_um.sam.0";
//			out_um_sam << castle::IOUtils::read_fully(local_mapped_um_name);
//			for (int64_t block_id = 1; block_id < calculated_n_blocks - 1; ++block_id) {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				string local_mapped_um_name = options.prefix + ".mapped_um.sam." + str_block_id;
//				if(!options.working_dir.empty()) {
//					local_mapped_um_name = options.working_prefix + ".mapped_um.sam." + str_block_id;
//				}
//				out_um_sam << castle::IOUtils::read_fully(local_mapped_um_name);
//			}
//		});
//		tasks.push_back([&]{
//			ofstream out_sc_sam(mapped_sc_sam_name, ios::binary);
//			string local_mapped_sc_name = options.prefix + ".mapped_sc.sam.0";
//			out_sc_sam << castle::IOUtils::read_fully(local_mapped_sc_name);
//			for (int64_t block_id = 1; block_id < calculated_n_blocks - 1; ++block_id) {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				string local_mapped_sc_name = options.prefix + ".mapped_sc.sam." + str_block_id;
//				if(!options.working_dir.empty()) {
//					local_mapped_sc_name = options.working_prefix + ".mapped_sc.sam." + str_block_id;
//				}
//				out_sc_sam << castle::IOUtils::read_fully(local_mapped_sc_name);
//			}
//		});

//		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//		string sambamba_merge_um_cmd = (boost::format("sambamba merge -l 1 -t %d %s %s") % n_cores % mapped_um_tmp_name % castle::StringUtils::join(mapped_um_names, " ")).str();
//		string sambamba_merge_sc_cmd = (boost::format("sambamba merge -l 1 -t %d %s %s") % n_cores % mapped_sc_tmp_name % castle::StringUtils::join(mapped_sc_names, " ")).str();
		string sambamba_sort_um_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % mapped_um_bam_name % mapped_um_tmp_bam_name).str();
		string sambamba_sort_sc_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % mapped_sc_bam_name % mapped_sc_tmp_bam_name).str();
//		system(sambamba_merge_um_cmd.c_str());
//		system(sambamba_merge_sc_cmd.c_str());
		system(sambamba_sort_um_cmd.c_str());
		system(sambamba_sort_sc_cmd.c_str());
	}

	{
		ofstream scrdist(options.scrdistname, ios::binary);
		chist.print(scrdist);
	}

	cout << "Final read counts after rejection due to quality (-q " << options.q << " -m " << options.min_read_len << "):\n   unmapped reads: "

	<< (num_unmapped - unmapped_rejected) << " (" << unmapped_rejected << " rejected)\n   clipped reads:  " << (num_clipped - clipped_rejected) << " (" << clipped_rejected << " rejected)\n";

	cout << checker;

}

void ParallelBamReader::output_read_groups_alt() {
	string pdffile = options.prefix + ".pdf";
	string resultfile = options.prefix + ".isinfo";
	string rdist_unmapfile = options.prefix + ".unmapped.rdist";
	string rdist_scfile = options.prefix + ".softclips.rdist";
	if (!options.working_dir.empty()) {
		pdffile = options.working_prefix + ".pdf";
		resultfile = options.working_prefix + ".isinfo";
		rdist_unmapfile = options.working_prefix + ".unmapped.rdist";
		rdist_scfile = options.working_prefix + ".softclips.rdist";
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.output_read_groups_alt");
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = unmapped_included_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');

	vector<string> read_groups;
	map<string, int64_t> read_groups_reverse_index;
	map<string, int64_t> ref_reverse_index;
	{
		BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
		const RefVector& a_ref_vector = local_reader.GetReferenceData();
		for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
			auto& a_ref = a_ref_vector[ref_id];
			ref_reverse_index[a_ref.RefName] = ref_id;
		}
		istringstream in(local_reader.GetHeaderText());
		string line;
		const char* delim = "\t:";
		vector<string> a_cols;
		while (getline(in, line, '\n')) {
			if (string::npos == line.find("@RG")) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim, a_cols);
			if (a_cols.size() < 2) {
				continue;
			}
			read_groups_reverse_index[a_cols[2]] = read_groups.size();
			read_groups.push_back(a_cols[2]);
		}
		string last_key = "none";
		read_groups_reverse_index[last_key] = read_groups.size();
		read_groups.push_back(last_key);
		bool all_processed = true;
		for (uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
			string rg_name = read_groups[rg_id];
			string filename_1 = options.prefix + "/" + rg_name + "_1.fq";
			string filename_2 = options.prefix + "/" + rg_name + "_2.fq";
			if (!options.working_dir.empty()) {
				filename_1 = options.working_prefix + "/" + rg_name + "_1.fq";
				filename_2 = options.working_prefix + "/" + rg_name + "_2.fq";
			}
			if (!boost::filesystem::exists(filename_1) || !boost::filesystem::exists(filename_2)) {
				all_processed = false;
				break;
			}
		}
		if (all_processed) {
			return;
		}
		checker.start();
		local_reader.Close();
	}
	n_read_groups = read_groups.size();
	vector<boost::unordered_map<string, BamAlignment>> albuf_lists(calculated_n_blocks);
	vector<boost::unordered_map<string, BamAlignment>> unmapped_lists(calculated_n_blocks);
	vector<boost::unordered_map<string, ReadGroupAlt>> rgs_lists(calculated_n_blocks - 1);

	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;
//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string the_next_block_read_name = unmapped_included_blocks[block_id + 1].read_name;
			string str_block_id = boost::lexical_cast<string>(block_id);

			vector<shared_ptr<ofstream>> read_groups_f1;
			vector<shared_ptr<ofstream>> read_groups_f2;

			for(uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
				string name = read_groups[rg_id];
				string filename_1 = options.prefix + "/" + name + "_1."
				+ str_block_id + ".fq";
				string filename_2 = options.prefix + "/" + name + "_2."
				+ str_block_id + ".fq";
				if(!options.working_dir.empty()) {
					filename_1 = options.working_prefix + "/" + name + "_1."
					+ str_block_id + ".fq";
					filename_2 = options.working_prefix + "/" + name + "_2."
					+ str_block_id + ".fq";
				}
				read_groups_f1.push_back(
						make_shared < ofstream > (filename_1, ios::binary));
				read_groups_f2.push_back(
						make_shared < ofstream > (filename_2, ios::binary));
			}

			int64_t the_current_ref_offset = unmapped_included_blocks[block_id].offset;
			int64_t the_next_ref_offset = unmapped_included_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[ParallelBamReader.output_read_groups_alt] (Jump fail) Block-%d (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//if(verbose) {
//cout << (boost::format("[ParallelBamReader.output_read_groups] (start) Block-%d (%d/%d)-(%d/%d)\n")
//% block_id % the_current_ref_id % the_current_ref_pos
//% the_next_ref_id % the_next_ref_pos).str();
//}

//				int32_t decoy_ref_id = -8;
//				auto decoy_itr = ref_reverse_index.find("hs37d5");
//				if(ref_reverse_index.end() != decoy_itr) {
//					decoy_ref_id = decoy_itr->second;
//				}
				string rg_name;
//				const RefVector& refnames = local_reader.GetReferenceData();
				BamAlignment local_alignment_entry;
				ReadGroup rg;
				auto& rgs = rgs_lists[block_id];
				auto& local_albuf = albuf_lists[block_id];
				auto& umbuf = unmapped_lists[block_id];

				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;
//				if(cur_offset != the_current_ref_offset) {
//					cout << (boost::format("[ParallelBamReader.output_read_groups_alt] block-%d (start) offset: jump: %d, calc: %d\n")
//				% block_id % cur_offset % the_current_ref_offset).str();
//				}

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						cout << (boost::format("[ParallelBamReader.output_read_groups_alt] (first) Block-%d %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//					if(local_alignment_entry.RefID == the_next_ref_id
//							&& local_alignment_entry.Position == the_next_ref_pos
//							&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//							&& local_alignment_entry.Name == the_next_block_read_name
//					) {
//						break;
//					}
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;
					//					if(the_next_ref_id == decoy_ref_id && local_alignment_entry.RefID == decoy_ref_id) {
//						break;
//					}
//					bool debug = string::npos != local_alignment_entry.Name.find("ST-E00104:502:HFJN5CCXX:1:1108:32826:10890");
					++num_total;
					if(!local_alignment_entry.GetReadGroup(rg_name)) {
						rg_name = "none";
					}
					rgs[rg_name].witness(local_alignment_entry, options.max_isize, options.isize_samples);
					if(black_listed.end() != black_listed.find(rg_name)) {
						continue;
					}

					if(!local_alignment_entry.IsMapped()) {
						ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
					}
//					if(debug) {
//						cout << "rg: here-0: " << block_id << "\n";
//					}
					int64_t local_rg_id = read_groups_reverse_index[rg_name];
					ofstream& f1 = *read_groups_f1[local_rg_id].get();
					ofstream& f2 = *read_groups_f2[local_rg_id].get();
					/* 2. Handle UU pairs */
					if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
//						if(debug) {
//							cout << "rg: here-1\n";
//						}
						auto x = umbuf.find(local_alignment_entry.Name);
						if (x == umbuf.end()) {
							umbuf[local_alignment_entry.Name] = local_alignment_entry;
							continue;
						}
//						if(debug) {
//							cout << "rg: here-2\n";
//						}
						BamAlignment mate = umbuf[local_alignment_entry.Name];
						if (!ReadGroup::isMatePair(local_alignment_entry, mate)) {
							continue;
						}

//						if(debug) {
//							cout << "rg: here-3\n";
//						}

						/* Both reads are trimmed by (1) before we get here */
						rg.recordUUAlt(f1, f2, local_alignment_entry, mate, options.big_s_bps, options.n_cutoff);
						umbuf.erase(x);
						continue;
					}
					/* 3a. Hardclipped reads aren't interesting.  In the new BWA-MEM,
					 * softclipped read records are often accompanied by a hard clipped
					 * record, resulting in 3 reads belonging to the same read pair to be
					 * reported in the BAM.  E.g.:
					 * ar3118743    pr1    chr13    18133804    0    48H27M    chr21    14051775    0    CCCAAGAGTGTGGAACTTTCCTTGGGC
					 * ar3118743    pR2    chr21    14051775    12    75M    =    14051906    179    ACTTCACAAAGCTAAGAAATACATATGCATATACCAAGCAAAACACACCATAAGGGTAAAATGATGACTTTTTTG
					 * ar3118743    pr1    chr21    14051906    27    48M27S    =    14051775    -179    TTAATAGGTCCAAATAACAGGTTTATGCTTTTGATTTTGCAGTGGAAGCCCAAGAGTGTGGAACTTTCCTTGGGC
					 */
//					if(debug) {
//						cout << "rg: here-4\n";
//					}
					if (ReadGroup::isBigH(local_alignment_entry)) {
						continue;
					}
//					if(debug) {
//						cout << "rg: here-5\n";
//					}
					/* 3. Handle mapped+[softclipped|unmapped] read pairs */
					int mate_num = 1;
					auto it = softclips.find(local_alignment_entry.Name+"1");

					if(softclips.end() == it) {
//						if(debug) {
//							cout << "rg: here-6\n";
//						}
						it = softclips.find(local_alignment_entry.Name+"2");
						mate_num = 2;
					}
					if (it != softclips.end()) {
//						if(debug) {
//							cout << "rg: here-7\n";
//						}
						auto x = local_albuf.find(local_alignment_entry.Name);
						if (x == local_albuf.end()) {
//							if(debug) {
//								cout << "rg: here-8\n";
//							}
							local_albuf[local_alignment_entry.Name] = local_alignment_entry;
						} else if (ReadGroup::isMatePair(local_alignment_entry, x->second)) {
//							if(debug) {
//								cout << "rg: here-9\n";
//							}
							/* If we're here, then we've got both reads for
							 * some pair that was stored in softclips */
							BamAlignment mate = x->second;
//							if(debug) {
//								cout << "rg: here-9-cur: " << BamWriter::GetSAMAlignment(local_alignment_entry, refnames) << "\n";
//								cout << "rg: here-9-mate: " << BamWriter::GetSAMAlignment(mate, refnames) << "\n";
//							}
							ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
							ReadGroup::trim_read(mate, options.q, qenc);

							boost::unordered_map<int8_t, CigarOp> local_map;

							if(local_alignment_entry.CigarData.size() > 0) {
								int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
								// current front
								if('S' == local_alignment_entry.CigarData.front().Type
									&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.front().Length)) >= options.frag_size) {
									local_map[local_mate_num] = local_alignment_entry.CigarData.front();
								}
								// current back
								if('S' == local_alignment_entry.CigarData.back().Type
									&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.back().Length)) >= options.frag_size) {
									if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > local_alignment_entry.CigarData.back().Length) {
										local_map[local_mate_num] = local_alignment_entry.CigarData.back();
									}
								}
							}
							if(mate.CigarData.size() > 0) {
								int8_t local_mate_num = mate.IsFirstMate() ? 1 : 2;
								// current front
								if('S' == mate.CigarData.front().Type
									&& min(mate.Length, static_cast<int32_t>(mate.CigarData.front().Length)) >= options.frag_size) {
									local_map[local_mate_num] = mate.CigarData.front();
								}
								// current back
								if('S' == mate.CigarData.back().Type
									&& min(mate.Length, static_cast<int32_t>(mate.CigarData.back().Length)) >= options.frag_size) {
									if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > mate.CigarData.back().Length) {
										local_map[local_mate_num] = mate.CigarData.back();
									}
								}
							}
							int8_t rep_mate = local_map[1].Length >= local_map[2].Length ? 1 : 2;
							int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
							int8_t mate_mate_num = mate.IsFirstMate() ? 1 : 2;

							BamAlignment local_left = ReadGroup::snip(local_alignment_entry, 0, options.frag_size);
							BamAlignment local_right = ReadGroup::snip(local_alignment_entry, local_alignment_entry.Length - options.frag_size, options.frag_size);
							size_t n_local_amb = count(local_left.QueryBases.begin(), local_left.QueryBases.end(), 'N');
							n_local_amb += count(local_right.QueryBases.begin(), local_right.QueryBases.end(), 'N');
							BamAlignment mate_left = ReadGroup::snip(mate, 0, options.frag_size);
							BamAlignment mate_right = ReadGroup::snip(mate, mate.Length - options.frag_size, options.frag_size);
							size_t n_mate_amb = count(mate_left.QueryBases.begin(), mate_left.QueryBases.end(), 'N');
							n_mate_amb += count(mate_right.QueryBases.begin(), mate_right.QueryBases.end(), 'N');
							if(rep_mate == local_mate_num && n_local_amb > static_cast<size_t>(options.n_cutoff)) {
								rep_mate = mate_mate_num;
							} else if(rep_mate == mate_mate_num && n_mate_amb > static_cast<size_t>(options.n_cutoff)) {
								rep_mate = local_mate_num;
							}

							if (local_alignment_entry.IsMapped() && mate.IsMapped()) {
//								if(debug) {
//									cout << "rg: here-10-mate_num: " << static_cast<int32_t>(mate_num) << "\n";
//									cout << "rg: here-10-cur: " << BamWriter::GetSAMAlignment(local_alignment_entry, refnames) << "\n";
//									cout << "rg: here-10-mate: " << BamWriter::GetSAMAlignment(mate, refnames) << "\n";
//									cout << "rg: here-10-rep_mate: " << static_cast<int32_t>(rep_mate) << "\n";
//									cout << "rg: here-10-options.big_s_bps: " << options.big_s_bps << "\n";
//									cout << "rg: here-10-options.frag_size: " << options.frag_size << "\n";
//									cout << "rg: here-10-options.n_cutoff: " << options.n_cutoff << "\n";
//								}
								rg.recordSCAltRG(f1, f2, ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate,
										ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? mate : local_alignment_entry, mate_num, options.big_s_bps,
										options.frag_size, options.n_cutoff);
							} else {
								rg.recordUMAltRG(f1, f2, local_alignment_entry, mate, options.big_s_bps, options.n_cutoff);
							}
							local_albuf.erase(x);
						}
					}
				}
//				if(prev_offset != the_next_ref_offset) {
//					cout << (boost::format("[ParallelBamReader.output_read_groups_alt] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n") % block_id % prev_offset % cur_offset % the_next_ref_offset).str();
//				}

				local_reader.Close();
				done_vector[block_id] = 'D';
				if(!silent) {
					if(verbose) {
						string a_block_boundary_str = (boost::format("%s %d-%d %d")
								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
								% local_alignment_entry.AlignmentFlag).str();
						cout << (boost::format("[ParallelBamReader.output_read_groups_alt] (last) Block-%d %s\n")
								% block_id % a_block_boundary_str).str();
					}
					else {
						size_t n = count(done_vector.begin(), done_vector.end(), 'D');
						double processed = n/(double)done_vector.size() * 100.0;
						cout << (boost::format("%.2f %%\n") % processed).str();
					}
				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	{
		boost::unordered_map<string, PairedAlignmentWithTwoIds> paired_alns;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_albuf = albuf_lists[block_id];
			for (auto itr_aln = local_albuf.begin(); local_albuf.end() != itr_aln; ++itr_aln) {
				if (itr_aln->second.IsFirstMate()) {
					paired_alns[itr_aln->first].pair_1 = itr_aln->second;
					paired_alns[itr_aln->first].id_1 = block_id;
				} else if (itr_aln->second.IsSecondMate()) {
					paired_alns[itr_aln->first].pair_2 = itr_aln->second;
					paired_alns[itr_aln->first].id_2 = block_id;
				}
			}
			local_albuf.clear();
		}
		cout << (boost::format("[ParallelBamReader.output_read_groups_alt] output unprocessed soft-clipped & UM reads: %d\n") % paired_alns.size()).str();
		output_unprocessed_read_groups(paired_alns, read_groups, read_groups_reverse_index);
	}

	tasks.push_back([&, calculated_n_blocks] {
		map<string, ReadGroupAlt> all_rgs;
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1;
				++block_id) {
			auto& local_rgs = rgs_lists[block_id];
			for (auto itr_a_map = local_rgs.begin();
					local_rgs.end() != itr_a_map; ++itr_a_map) {
				all_rgs[itr_a_map->first] += itr_a_map->second;
			}
		}
		ofstream isinfo(options.isinfoname, ios::binary);
		for (auto rgit = all_rgs.begin(); all_rgs.end() != rgit; ++rgit) {
			auto& x = rgit->second;
			cout << "Readgroup [" + rgit->first + "]: " << x.nreads
			<< " reads\n";
			isinfo << "Read length:\t" << rgit->first << "\n";
			isinfo << x.readlens[0] << "\n";
			/* Write the insert size distribution to file */
			string an_is_file_name = (options.prefix + "/" + rgit->first + ".is");
			if(!options.working_dir.empty()) {
				an_is_file_name = (options.working_prefix + "/" + rgit->first + ".is");
			}
			ofstream f(an_is_file_name.c_str());
			int64_t n_inserts = min(static_cast<int64_t>(options.isize_samples),
					static_cast<int64_t>(x.inserts.size()));
			for (int64_t i = 0; i < n_inserts; ++i) {
				f << x.inserts[i] << "\n";
			}
		}
	});

	for (uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
		tasks.push_back([&, rg_id, calculated_n_blocks] {
			string rg_name = read_groups[rg_id];
			string root_filename_1 = options.prefix + "/" + rg_name + "_1.fq";
			if(!options.working_dir.empty()) {
				root_filename_1 = options.working_prefix + "/" + rg_name + "_1.fq";
			}
			ofstream f1(root_filename_1, ios::binary);
			for(int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
				string str_block_id = boost::lexical_cast<string>(block_id);
				string filename_1 = options.prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
				if(!options.working_dir.empty()) {
					filename_1 = options.working_prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
				}
				f1 << castle::IOUtils::read_fully(filename_1);
			}
		});
		tasks.push_back([&, rg_id, calculated_n_blocks] {
			string rg_name = read_groups[rg_id];
			string root_filename_2 = options.prefix + "/" + rg_name + "_2.fq";
			if(!options.working_dir.empty()) {
				root_filename_2 = options.working_prefix + "/" + rg_name + "_2.fq";
			}
			ofstream f2(root_filename_2, ios::binary);
			for(int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
				string str_block_id = boost::lexical_cast<string>(block_id);
				string filename_2 = options.prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
				if(!options.working_dir.empty()) {
					filename_2 = options.working_prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
				}
				f2 << castle::IOUtils::read_fully(filename_2);
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	for (uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
//		tasks.push_back([&, rg_id, calculated_n_blocks] {
//			string rg_name = read_groups[rg_id];
//			for(int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				string filename_1 = options.prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
//				if(!options.working_dir.empty()) {
//					filename_1 = options.working_prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
//				}
//				if(boost::filesystem::exists(filename_1)) {
//					boost::filesystem::remove(filename_1);
//				}
//			}
//		});
//		tasks.push_back([&, rg_id, calculated_n_blocks] {
//			string rg_name = read_groups[rg_id];
//			for(int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
//				string str_block_id = boost::lexical_cast<string>(block_id);
//				string filename_2 = options.prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
//				if(!options.working_dir.empty()) {
//					filename_2 = options.working_prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
//				}
//				if(boost::filesystem::exists(filename_2)) {
//					boost::filesystem::remove(filename_2);
//				}
//			}
//		});
//	}
//	tasks.push_back([&] {
//		string filename_1 = options.prefix + "/none_1.fq";
//		if(!options.working_dir.empty()) {
//			filename_1 = options.working_prefix + "/none_1.fq";
//		}
//		if(boost::filesystem::exists(filename_1)) {
//			boost::filesystem::remove(filename_1);
//		}
//		string filename_2 = options.prefix + "/none_2.fq";
//		if(!options.working_dir.empty()) {
//			filename_2 = options.working_prefix + "/none_2.fq";
//		}
//		if(boost::filesystem::exists(filename_2)) {
//			boost::filesystem::remove(filename_2);
//		}
//	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void ParallelBamReader::output_unprocessed_read_groups(boost::unordered_map<string, PairedAlignmentWithTwoIds>& paired_alns, vector<string>& read_groups, map<string, int64_t>& read_groups_reverse_index) {
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.output_unprocessed_read_groups");
	checker.start();
	int64_t calculated_n_blocks = actual_blocks.size();
	vector<vector<shared_ptr<stringstream>>> read_groups_f1_lists(calculated_n_blocks - 1);
	vector<vector<shared_ptr<stringstream>>> read_groups_f2_lists(calculated_n_blocks - 1);
	vector<vector<uint64_t>> read_groups_f1_buf_size_lists(calculated_n_blocks - 1);
	vector<vector<uint64_t>> read_groups_f2_buf_size_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		//			string str_block_id = boost::lexical_cast<string>(block_id);
		auto& read_groups_f1 = read_groups_f1_lists[block_id];
		auto& read_groups_f2 = read_groups_f2_lists[block_id];
		auto& read_groups_buf_size_f1 = read_groups_f1_buf_size_lists[block_id];
		auto& read_groups_buf_size_f2 = read_groups_f2_buf_size_lists[block_id];
		//			vector<shared_ptr<stringstream>> read_groups_f2;
		for (uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
			//				string name = read_groups[rg_id];
			//				string filename_1 = options.prefix + "/" + name + "_1."
			//				+ str_block_id + ".fq";
			//				string filename_2 = options.prefix + "/" + name + "_2."
			//				+ str_block_id + ".fq";
			//				if(!options.working_dir.empty()) {
			//					filename_1 = options.working_prefix + "/" + name + "_1."
			//					+ str_block_id + ".fq";
			//					filename_2 = options.working_prefix + "/" + name + "_2."
			//					+ str_block_id + ".fq";
			//				}
			read_groups_f1.push_back(make_shared<stringstream>());
			read_groups_f2.push_back(make_shared<stringstream>());
			read_groups_buf_size_f1.push_back(0);
			read_groups_buf_size_f2.push_back(0);
		}
	}

	string rg_name;
	ReadGroup rg;
	auto pair_aln_itr = paired_alns.begin();
	while (paired_alns.end() != pair_aln_itr) {
		if (-1 == pair_aln_itr->second.id_1 || -1 == pair_aln_itr->second.id_2) {
			++pair_aln_itr;
			continue;
		}

		int64_t writing_block_id = max(pair_aln_itr->second.id_1, pair_aln_itr->second.id_2);
		auto& local_alignment_entry = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_1 : pair_aln_itr->second.pair_2;
		if (!local_alignment_entry.GetReadGroup(rg_name)) {
			rg_name = "none";
		}

		if (black_listed.end() != black_listed.find(rg_name)) {
			++pair_aln_itr;
			continue;
		}

		/* 2. Handle UU pairs */
		if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped() && options.processUU) {
			++pair_aln_itr;
			continue;
		}

		auto& mate = (writing_block_id == pair_aln_itr->second.id_1) ? pair_aln_itr->second.pair_2 : pair_aln_itr->second.pair_1;

		int mate_num = 1;
		auto it = softclips.find(local_alignment_entry.Name + "1");
		if (softclips.end() == it) {
			it = softclips.find(local_alignment_entry.Name + "2");
			mate_num = 2;
		}
		if (softclips.end() == it) {
			++pair_aln_itr;
			continue;
		}

		ReadGroup::trim_read(local_alignment_entry, options.q, qenc);
		ReadGroup::trim_read(mate, options.q, qenc);
		//			string str_block_id = boost::lexical_cast<string>(writing_block_id);
		//			string filename_1 = options.prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
		//			string filename_2 = options.prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
		//			if(!options.working_dir.empty()) {
		//				filename_1 = options.working_prefix + "/" + rg_name + "_1."
		//				+ str_block_id + ".fq";
		//				filename_2 = options.working_prefix + "/" + rg_name + "_2."
		//				+ str_block_id + ".fq";
		//			}
		auto& read_groups_f1 = read_groups_f1_lists[writing_block_id];
		auto& read_groups_f2 = read_groups_f2_lists[writing_block_id];
		auto& local_buf_size_f1 = read_groups_f1_buf_size_lists[writing_block_id];
		auto& local_buf_size_f2 = read_groups_f1_buf_size_lists[writing_block_id];

		stringstream f1;
		stringstream f2;

		boost::unordered_map<int8_t, CigarOp> local_map;

		if(local_alignment_entry.CigarData.size() > 0) {
			int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
			// current front
			if('S' == local_alignment_entry.CigarData.front().Type
				&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.front().Length)) >= options.frag_size) {
				local_map[local_mate_num] = local_alignment_entry.CigarData.front();
			}
			// current back
			if('S' == local_alignment_entry.CigarData.back().Type
				&& min(local_alignment_entry.Length, static_cast<int32_t>(local_alignment_entry.CigarData.back().Length)) >= options.frag_size) {
				if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > local_alignment_entry.CigarData.back().Length) {
					local_map[local_mate_num] = local_alignment_entry.CigarData.back();
				}
			}
		}
		if(mate.CigarData.size() > 0) {
			int8_t local_mate_num = mate.IsFirstMate() ? 1 : 2;
			// current front
			if('S' == mate.CigarData.front().Type
				&& min(mate.Length, static_cast<int32_t>(mate.CigarData.front().Length)) >= options.frag_size) {
				local_map[local_mate_num] = mate.CigarData.front();
			}
			// current back
			if('S' == mate.CigarData.back().Type
				&& min(mate.Length, static_cast<int32_t>(mate.CigarData.back().Length)) >= options.frag_size) {
				if(local_map.end() == local_map.find(local_mate_num) || local_map[local_mate_num].Length > mate.CigarData.back().Length) {
					local_map[local_mate_num] = mate.CigarData.back();
				}
			}
		}
		int8_t rep_mate = local_map[1].Length >= local_map[2].Length ? 1 : 2;
		int8_t local_mate_num = local_alignment_entry.IsFirstMate() ? 1 : 2;
		int8_t mate_mate_num = mate.IsFirstMate() ? 1 : 2;

		BamAlignment local_left = ReadGroup::snip(local_alignment_entry, 0, options.frag_size);
		BamAlignment local_right = ReadGroup::snip(local_alignment_entry, local_alignment_entry.Length - options.frag_size, options.frag_size);
		size_t n_local_amb = count(local_left.QueryBases.begin(), local_left.QueryBases.end(), 'N');
		n_local_amb += count(local_right.QueryBases.begin(), local_right.QueryBases.end(), 'N');
		BamAlignment mate_left = ReadGroup::snip(mate, 0, options.frag_size);
		BamAlignment mate_right = ReadGroup::snip(mate, mate.Length - options.frag_size, options.frag_size);
		size_t n_mate_amb = count(mate_left.QueryBases.begin(), mate_left.QueryBases.end(), 'N');
		n_mate_amb += count(mate_right.QueryBases.begin(), mate_right.QueryBases.end(), 'N');
		if(rep_mate == local_mate_num && n_local_amb > static_cast<size_t>(options.n_cutoff)) {
			rep_mate = mate_mate_num;
		} else if(rep_mate == mate_mate_num && n_mate_amb > static_cast<size_t>(options.n_cutoff)) {
			rep_mate = local_mate_num;
		}

		if (local_alignment_entry.IsMapped() && mate.IsMapped()) {
			rg.recordSCAltRG(f1, f2, ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? local_alignment_entry : mate,
					ReadGroup::getMateNumber(local_alignment_entry) == rep_mate ? mate : local_alignment_entry, mate_num, options.big_s_bps, options.frag_size, options.n_cutoff);
		} else {
			rg.recordUMAltRG(f1, f2, local_alignment_entry, mate, options.big_s_bps, options.n_cutoff);
		}

		int64_t local_rg_id = read_groups_reverse_index[rg_name];
		string f1_str = f1.str();
		string f2_str = f2.str();
		local_buf_size_f1[local_rg_id] += f1_str.size();
		local_buf_size_f2[local_rg_id] += f2_str.size();

		// f1_mem and f2_mem are in fact stringstream objects
		ostream& f1_mem = *read_groups_f1[local_rg_id].get();
		ostream& f2_mem = *read_groups_f2[local_rg_id].get();
		f1_mem << f1_str;
		f2_mem << f2_str;
		if (local_buf_size_f1[local_rg_id] > castle::IOUtils::WRITE_BUFFER_SIZE) {
			string str_block_id = boost::lexical_cast<string>(writing_block_id);
			string filename_1 = options.prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
			if (!options.working_dir.empty()) {
				filename_1 = options.working_prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
			}
			ofstream out_f1(filename_1, ios::binary | ios::app);
			stringstream& local_f1_mem = *read_groups_f1[local_rg_id].get();
			out_f1 << local_f1_mem.str();
			local_f1_mem.str(string());
			local_buf_size_f1[local_rg_id] = 0;
		}

		if (local_buf_size_f2[local_rg_id] > castle::IOUtils::WRITE_BUFFER_SIZE) {
			string str_block_id = boost::lexical_cast<string>(writing_block_id);
			string filename_2 = options.prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
			if (!options.working_dir.empty()) {
				filename_2 = options.working_prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
			}
			ofstream out_f2(filename_2, ios::binary | ios::app);
			stringstream& local_f2_mem = *read_groups_f2[local_rg_id].get();
			out_f2 << local_f2_mem.str();
			local_f2_mem.str(string());
			local_buf_size_f2[local_rg_id] = 0;
		}
		pair_aln_itr = paired_alns.erase(pair_aln_itr);
	}
	// write remaining FASTQ entries.
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		auto& read_groups_f1 = read_groups_f1_lists[block_id];
		auto& read_groups_f2 = read_groups_f2_lists[block_id];
		for (uint64_t rg_id = 0; rg_id < read_groups.size(); ++rg_id) {
			stringstream& f1_mem = *read_groups_f1[rg_id].get();
			stringstream& f2_mem = *read_groups_f2[rg_id].get();
			string f1_str = f1_mem.str();
			string f2_str = f2_mem.str();
			string rg_name = read_groups[rg_id];
			string filename_1 = options.prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
			string filename_2 = options.prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
			if (!options.working_dir.empty()) {
				filename_1 = options.working_prefix + "/" + rg_name + "_1." + str_block_id + ".fq";
				filename_2 = options.working_prefix + "/" + rg_name + "_2." + str_block_id + ".fq";
			}
			if (!f1_str.empty()) {
				ofstream out_f1(filename_1, ios::binary | ios::app);
				out_f1 << f1_str;
			}
			if (!f2_str.empty()) {
				ofstream out_f2(filename_2, ios::binary | ios::app);
				out_f2 << f2_str;
			}
		}
	}
	cout << checker;
}

void ParallelBamReader::create_reports() {
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = actual_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');

	vector<string> read_groups;
	map<string, int64_t> read_groups_reverse_index;
	map<string, int64_t> ref_reverse_index;
	{
		BamTools::BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
		const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
		for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
			auto& a_ref = a_ref_vector[ref_id];
			ref_reverse_index[a_ref.RefName] = ref_id;
		}
		istringstream in(local_reader.GetHeaderText());
		string line;
		const char* delim = "\t:";
		vector<string> a_cols;
		while (getline(in, line, '\n')) {
			if (string::npos == line.find("@RG")) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim, a_cols);
			if (a_cols.size() < 2) {
				continue;
			}
			read_groups_reverse_index[a_cols[2]] = read_groups.size();
			read_groups.push_back(a_cols[2]);
		}
		string last_key = "none";
		read_groups_reverse_index[last_key] = read_groups.size();
		read_groups.push_back(last_key);
		local_reader.Close();
	}
	n_read_groups = read_groups.size();
	string rfile = options.prefix + ".insert.r";
	string pdffile = options.prefix + ".pdf";
	string resultfile = options.prefix + ".isinfo";
	string rdist_unmapfile = options.prefix + ".unmapped.rdist";
	string rdist_scfile = options.prefix + ".softclips.rdist";
	if (!options.working_dir.empty()) {
		rfile = options.working_prefix + ".insert.r";
		pdffile = options.working_prefix + ".pdf";
		resultfile = options.working_prefix + ".isinfo";
		rdist_unmapfile = options.working_prefix + ".unmapped.rdist";
		rdist_scfile = options.working_prefix + ".softclips.rdist";
	}
	{
		int64_t rgi = 0;
		ofstream RFILE(rfile, ios::binary);
		RFILE << (boost::format("pdf(\"%s\")\n") % pdffile).str();
		for (auto& rg : read_groups) {
			string insertsizefile = (options.prefix + "/" + rg + ".is");
			if (!options.working_dir.empty()) {
				insertsizefile = (options.working_prefix + "/" + rg + ".is");
			}
			int64_t filesize = castle::IOUtils::get_file_size(insertsizefile);
			if (0 < filesize) {
				string rgname = "rg" + boost::lexical_cast<string>(rgi);
				RFILE << (boost::format("%s=read.table(\"%s\")\n") % rgname % insertsizefile).str();
				if (rgi) {
					RFILE << (boost::format("data=c(data,%s$V1)\n") % rgname).str();
				} else {
					RFILE << (boost::format("data=%s$V1\n") % rgname).str();
				}
				RFILE << (boost::format("m=mean(%s[,1])\n"
						"n=median(%s[,1])\n"
						"sd=sd(%s[,1])\n"
						"write(\"Mean insert size:\t%s\", append=T, file=\"%s\")\n"
						"write(m, append=T, file=\"%s\")\n"
						"write(\"Median insert size:\t%s\", append=T, file=\"%s\")\n"
						"write(n, append=T, file=\"%s\")\n"
						"write(\"Standard deviation of insert size:\t%s\", append=T, file=\"%s\")\n"
						"write(sd, append=T, file=\"%s\")\n") % rgname % rgname % rgname % rg % resultfile % resultfile % rg % resultfile % resultfile % rg % resultfile % resultfile).str();
			}
			++rgi;
		}
		RFILE << (boost::format("hist(data,xlim=c(0,%d),breaks=1000,xlab=\"Insert size\",main=\"\")\n"
				"unmap=read.table(\"%s\")\n"
				"sc=read.table(\"%s\")\n"
				"barplot(unmap[,3],col=1,space=0,names.arg=unmap[,1],main=\"Read length distribution of unmapped reads\",xlab=\"Read length (bp)\",ylab=\"Fraction\")\n"
				"barplot(sc[,3],col=1,space=0,names.arg=sc[,1],main=\"Read length distribution of soft clipped reads\",xlab=\"Read length (bp)\",ylab=\"Fraction\")\n") % options.plotrange % rdist_unmapfile % rdist_scfile).str();
	}
	string r_cmd = (boost::format("Rscript %s") % rfile).str();
	system(r_cmd.c_str());
}

void ParallelBamReader::align_clipped_reads() {
	if (!options.clip) {
		return;
	}
	string cl_sortbam = options.prefix + ".cl.sorted.bam";
	if (!options.working_dir.empty()) {
		cl_sortbam = options.working_prefix + ".cl.sorted.bam";
	}
//	if (boost::filesystem::exists(cl_sortbam)) {
//		return;
//	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.align_clipped_reads");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.fname);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);

	int64_t calculated_n_blocks = actual_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');

	cout << "[ParallelBamReader.align_clipped_reads] collect read group data\n";
	vector<string> read_groups;
	map<string, int64_t> read_groups_reverse_index;
	map<string, int64_t> ref_reverse_index;
	{
		BamTools::BamReader local_reader;
		if (!local_reader.Open(a_path, an_index_path)) {
			return;
		}
		const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
		for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
			auto& a_ref = a_ref_vector[ref_id];
			ref_reverse_index[a_ref.RefName] = ref_id;
		}
		istringstream in(local_reader.GetHeaderText());
		string line;
		const char* delim = "\t:";
		vector<string> a_cols;
		while (getline(in, line, '\n')) {
			if (string::npos == line.find("@RG")) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim, a_cols);
			if (a_cols.size() < 2) {
				continue;
			}
			read_groups_reverse_index[a_cols[2]] = read_groups.size();
			read_groups.push_back(a_cols[2]);
		}
		string cl_fq1 = options.prefix + "/none_1.fq";
		if (!options.working_dir.empty()) {
			cl_fq1 = options.working_prefix + "/none_1.fq";
		}
		if (castle::IOUtils::get_file_size(cl_fq1) > 0) {
			string last_key = "none";
			read_groups_reverse_index[last_key] = read_groups.size();
			read_groups.push_back(last_key);
		}
		local_reader.Close();
	}
	n_read_groups = read_groups.size();
//	if (boost::filesystem::exists("bwa.err")) {
//		boost::filesystem::remove("bwa.err");
//	}
	cout << "[ParallelBamReader.align_clipped_reads] align clipped reads\n";
	castle::TimeChecker sub_checker;
	sub_checker.setTarget("ParallelBamReader.align_clipped_reads - align");
	sub_checker.start();
	auto& rg_blacklist = options.rg_blacklist;
	for (auto& rgname : read_groups) {
		if (rg_blacklist.end() != rg_blacklist.find(rgname)) {
			continue;
		}

		string cl_fq1 = options.prefix + "/" + rgname + "_1.fq";
		string cl_fq2 = options.prefix + "/" + rgname + "_2.fq";
//		string cl_sai1 = options.prefix + rgname + "_1.sai";
//		string cl_sai2 = options.prefix + rgname + "_2.sai";
		string cl_sam = options.prefix + "/" + rgname + ".sam";
//		string cl_tmp_sam = options.prefix + rgname + ".tmp.sam";
		if (!options.working_dir.empty()) {
			cl_fq1 = options.working_prefix + "/" + rgname + "_1.fq";
			cl_fq2 = options.working_prefix + "/" + rgname + "_2.fq";
//			cl_sai1 = options.working_prefix + "/" + rgname + "_1.sai";
//			cl_sai2 = options.working_prefix + "/" + rgname + "_2.sai";
//			cl_tmp_sam = options.working_prefix + "/" + rgname + ".tmp.sam";
			cl_sam = options.working_prefix + "/" + rgname + ".sam";
		}
		string aln_param = "-l 40 -k 2";
		string sampe_param = (boost::format("-P -N %d") % options.alt_map_max).str();
		BWACaller bc;
		bc.set_n_cores(n_cores);
		int64_t n_blocks = bc.split_FASTQ_alt(cl_fq1, cl_fq2, 262144);
		bc.collect_align_tasks_alt(tasks, cl_sam, options.reference_path, aln_param, sampe_param, cl_fq1, cl_fq2, n_blocks);
//		bc.collect_align_tasks(tasks, cl_sam, options.reference_path, aln_param, sampe_param, cl_fq1, cl_fq2);
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << sub_checker;
	cout << "[ParallelBamReader.align_clipped_reads] add read groups\n";
	sub_checker.setTarget("ParallelBamReader.align_clipped_reads - add rg");
	sub_checker.start();

	vector<string> all_removal_names;
	vector<string> all_SAM_file_names;

	for (auto& rgname : read_groups) {
		if (rg_blacklist.end() != rg_blacklist.find(rgname)) {
			continue;
		}
		string fq1_file = options.prefix + "/" + rgname + "_1.fq";
		string fq2_file = options.prefix + "/" + rgname + "_2.fq";
		string is_file = options.prefix + "/" + rgname + ".is";
		all_removal_names.push_back(fq1_file);
		all_removal_names.push_back(fq2_file);
		all_removal_names.push_back(is_file);

		for (int64_t block_id = 0;; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string cl_fq1 = options.prefix + "/" + rgname + "_1.fq." + str_block_id;
			string cl_fq2 = options.prefix + "/" + rgname + "_2.fq." + str_block_id;
			string cl_sai1 = options.prefix + "/" + rgname + "_1.fq.sai." + str_block_id;
			string cl_sai2 = options.prefix + "/" + rgname + "_2.fq.sai." + str_block_id;

			string cl_tmp_sam = options.prefix + "/" + rgname + ".tmp.sam." + str_block_id;
			string cl_sam = options.prefix + "/" + rgname + ".sam." + str_block_id;
			string bwa_err = options.prefix + "/" + rgname + ".fq.bwa.err." + str_block_id;
			string cl_fq1_bwa_err = options.prefix + "/" + rgname + "_1.fq.bwa.err." + str_block_id;
			string cl_fq2_bwa_err = options.prefix + "/" + rgname + "_2.fq.bwa.err." + str_block_id;
			if (!options.working_dir.empty()) {
				cl_fq1 = options.working_prefix + "/" + rgname + "_1.fq." + str_block_id;
				cl_fq2 = options.working_prefix + "/" + rgname + "_2.fq." + str_block_id;
				cl_sai1 = options.working_prefix + "/" + rgname + "_1.fq.sai." + str_block_id;
				cl_sai2 = options.working_prefix + "/" + rgname + "_2.fq.sai." + str_block_id;
				cl_tmp_sam = options.working_prefix + "/" + rgname + ".tmp.sam." + str_block_id;
				cl_sam = options.working_prefix + "/" + rgname + ".sam." + str_block_id;
				bwa_err = options.working_prefix + "/" + rgname + ".fq.bwa.err." + str_block_id;
				cl_fq1_bwa_err = options.working_prefix + "/" + rgname + "_1.fq.bwa.err." + str_block_id;
				cl_fq2_bwa_err = options.working_prefix + "/" + rgname + "_2.fq.bwa.err." + str_block_id;
			}
			if (!boost::filesystem::exists(cl_sam)) {
				break;
			}
			all_SAM_file_names.push_back(cl_tmp_sam);
			all_removal_names.push_back(cl_fq1);
			all_removal_names.push_back(cl_fq2);
			all_removal_names.push_back(cl_sai1);
			all_removal_names.push_back(cl_sai2);
			all_removal_names.push_back(cl_tmp_sam);
			all_removal_names.push_back(bwa_err);
			all_removal_names.push_back(cl_fq1_bwa_err);
			all_removal_names.push_back(cl_fq2_bwa_err);
			all_removal_names.push_back(cl_sam);
		}
	}
	string the_first_rg;
	if (all_SAM_file_names.size() > 0) {
		the_first_rg = all_SAM_file_names[0];
	}
//	cout << (boost::format("[ParallelBamReader.align_clipped_reads] First: %s\n") % the_first_rg).str();
	for (auto& rgname : read_groups) {
		if (rg_blacklist.end() != rg_blacklist.find(rgname)) {
			continue;
		}
		for (int64_t block_id = 0;; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string cl_sam = options.prefix + "/" + rgname + ".sam." + str_block_id;
			if (!options.working_dir.empty()) {
				cl_sam = options.working_prefix + "/" + rgname + ".sam." + str_block_id;
			}
			if (!boost::filesystem::exists(cl_sam)) {
				break;
			}
			tasks.push_back([&, the_first_rg, block_id, str_block_id] {

				const char* delim_tab = "\t";
				vector<string> a;
				string cl_sam = options.prefix + "/" +rgname + ".sam." + str_block_id;
				string cl_tmp_sam = options.prefix + "/" + rgname + ".tmp.sam." + str_block_id;
				string cl_bam = options.prefix + "/" + rgname + ".bam." + str_block_id;
				string cl_sorted_bam = options.prefix + "/" + rgname + ".sorted.bam." + str_block_id;
				if (!options.working_dir.empty()) {
					cl_sam = options.working_prefix + "/" + rgname + ".sam." + str_block_id;
					cl_tmp_sam = options.working_prefix + "/" + rgname + ".tmp.sam." + str_block_id;
					cl_bam = options.working_prefix + "/" + rgname + ".bam." + str_block_id;
					cl_sorted_bam = options.working_prefix + "/" + rgname + ".sorted.bam." + str_block_id;
				}
				bool is_first_rg = the_first_rg == cl_tmp_sam;
				string line;
				ofstream out(cl_tmp_sam, ios::binary);
				ifstream in(cl_sam, ios::binary);
				// write the SAM header.
					if(is_first_rg) {
						while (getline(in, line, '\n')) {
							if(line.empty()) {
								continue;
							}
							if('@' == line[0]) {
								out << line << "\n";
								continue;
							}
							castle::StringUtils::tokenize(line, delim_tab, a);
							if (a.size() > 11) {
								a[11] = (boost::format("RG:Z:%s\t%s") % rgname % a[11]).str();
							} else {
								a.resize(12);
								a[11] = "RG:Z:" + rgname;
							}
							out << castle::StringUtils::join(a, "\t") << "\n";
						}
					} else {
						while (getline(in, line, '\n')) {
							if(line.empty()) {
								continue;
							}
							if('@' == line[0]) {
								continue;
							}
							castle::StringUtils::tokenize(line, delim_tab, a);
							if (a.size() > 11) {
								a[11] = (boost::format("RG:Z:%s\t%s") % rgname % a[11]).str();
							} else {
								a.resize(12);
								a[11] = "RG:Z:" + rgname;
							}
							out << castle::StringUtils::join(a, "\t") << "\n";
						}
					}
//				boost::filesystem::remove(cl_sam);
//				boost::filesystem::rename(cl_tmp_sam, cl_sam);
//			string samtools_cmd = (boost::format("samtools view -bt %s -o %s %s") % reference_fai % cl_bam % cl_tmp_sam).str();
//			string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t 1 -o %s %s") % cl_sorted_bam % cl_bam).str();
//			system(samtools_cmd.c_str());
//			system(sambamba_sort_cmd .c_str());
				});
		}
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << sub_checker;
	sub_checker.setTarget("ParallelBamReader.align_clipped_reads - merge SAM");
	sub_checker.start();
	cout << "[ParallelBamReader.align_clipped_reads] merge SAM files\n";
	string merged_sam = options.prefix + ".cl.sam";
	if (!options.working_dir.empty()) {
		merged_sam = options.working_prefix + ".cl.sam";
	}
	castle::IOUtils::plain_file_merge(merged_sam, all_SAM_file_names, n_cores, false);
	cout << sub_checker;
	cout << "[ParallelBamReader.align_clipped_reads] create BAM file\n";
	sub_checker.setTarget("ParallelBamReader.align_clipped_reads - create BAM");
	sub_checker.start();
	string cl_bam = options.prefix + ".cl.bam";
	if (!options.working_dir.empty()) {
		cl_bam = options.working_prefix + ".cl.bam";
	}
	string samtools_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % cl_bam % merged_sam).str();
	system(samtools_cmd.c_str());
	cout << sub_checker;

	sub_checker.setTarget("ParallelBamReader.align_clipped_reads - sort BAM");
	sub_checker.start();
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % cl_sortbam % cl_bam).str();
	system(sambamba_sort_cmd.c_str());
	cout << sub_checker;

//	for(auto& a_path : all_removal_names) {
//		cout << a_path << "\n";
//	}

//	castle::IOUtils::remove_files(all_removal_names, n_cores);

	cout << checker;
}


void ParallelBamReader::create_bni_even_index(const string& a_path, string bni_index_path) {
//	get_bni_index_path(a_path, bni_index_path);
	if (boost::filesystem::exists(bni_index_path)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.create_even_index");
	checker.start();

	string bai_index_path;
	get_bai_index_path(a_path, bai_index_path);
	BamTools::BamReader serial_reader;
	if (!serial_reader.Open(a_path, bai_index_path)) {
		cout << checker;
		serial_reader.Close();
		return;
	}

	vector<BamTools::BlockOffset> block_boundary;
	serial_reader.GetBlockOffsets(block_boundary);
	ofstream out(bni_index_path, ios::binary);
	for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
		out << block_boundary[block_id].offset << "\t" << block_boundary[block_id].ref_id << "\t" << block_boundary[block_id].position << "\n";
	}
	cout << (boost::format("[PreProcessor.create_even_index] %d entries were indexed\n") % block_boundary.size()).str();
	cout << checker;
}

} /* namespace meerkat */

