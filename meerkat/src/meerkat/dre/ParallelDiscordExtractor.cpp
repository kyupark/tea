/*
 * ParallelDiscordExtractor.cpp
 *
 *  Created on: Jun 14, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#include "ParallelDiscordExtractor.hpp"

namespace meerkat {

ParallelDiscordExtractor::ParallelDiscordExtractor() :
		not_paired(0), not_both_mapped(0), not_same_chr(0), same_pos(0), not_diff_strands(0), weird_insert_size(0), big_insert_size(0), num_memorized_clipped_discs(0) {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
	//		memo_clipped_disc.set_empty_key("");
}

ParallelDiscordExtractor::~ParallelDiscordExtractor() {
}

void ParallelDiscordExtractor::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	options.read_blacklist();
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}

void ParallelDiscordExtractor::extract_discordant_reads() {
	string out_sorted_bam_name = options.prefix + ".disc.sorted.bam";
	string out_sorted_num_bam_name = options.prefix + ".disc.num.bam";

	if(!options.working_dir.empty()) {
		out_sorted_bam_name = options.working_prefix + ".disc.sorted.bam";
		out_sorted_num_bam_name = options.working_prefix + ".disc.num.bam";
	}
	const bool write_dup_file = !options.dupfile_name.empty();
	if (boost::filesystem::exists(out_sorted_bam_name) && boost::filesystem::exists(out_sorted_num_bam_name)) {
		if (write_dup_file) {
			string out_sorted_dup_bam_name = options.prefix + ".dup.sorted.bam";
			string out_sorted_dup_sorted_num_bam_name = options.prefix + ".dup.num.bam";
			if(!options.working_dir.empty()) {
				out_sorted_dup_bam_name = options.working_prefix + ".dup.sorted.bam";
				out_sorted_dup_sorted_num_bam_name = options.working_prefix + ".dup.num.bam";
			}
			if (boost::filesystem::exists(out_sorted_dup_bam_name) && boost::filesystem::exists(out_sorted_dup_sorted_num_bam_name)) {
				return;
			}
		}
		else {
			return;
		}
	}

	if (options.verbose) {
		cout << "Memorizing clipped discordant reads...\n";
	}
	string a_path(options.infile_name);
	string bni_path;
	create_bni_even_index(a_path, bni_path);
	int64_t size_block = 8192000;
	collect_boundaries(size_block);
//		collect_independent_boundaries();
//		collect_chromosome_wide_boundaries();
	/* Build a set of read names that will be considered discordant.  This
	 * is specifically to handle clipped reads which, under BWA-MEM, sometimes
	 * generate 3 records per read pair rather than 2.  In this case, one of
	 * the 2 hard- and softclipped read records may be discordant, but that
	 * information may not be reflected in the chosen read mate.  See
	 * bamreader.cpp's source for an example.
	 */
	collect_clipped_discordants();
	write_discordants_alt();
	collect_discordants_sambamba();
	//		collect_discordants();
}

void ParallelDiscordExtractor::extract_discordant_reads_serial() {
	string out_sorted_bam_name = options.prefix + ".disc.sorted.bam";
	string out_sorted_num_bam_name = options.prefix + ".disc.num.bam";

	if(!options.working_dir.empty()) {
		out_sorted_bam_name = options.working_prefix + ".disc.sorted.bam";
		out_sorted_num_bam_name = options.working_prefix + ".disc.num.bam";
	}
	const bool write_dup_file = !options.dupfile_name.empty();
	if (boost::filesystem::exists(out_sorted_bam_name) && boost::filesystem::exists(out_sorted_num_bam_name)) {
		if (write_dup_file) {
			string out_sorted_dup_bam_name = options.prefix + ".dup.sorted.bam";
			string out_sorted_dup_sorted_num_bam_name = options.prefix + ".dup.num.bam";
			if(!options.working_dir.empty()) {
				out_sorted_dup_bam_name = options.working_prefix + ".dup.sorted.bam";
				out_sorted_dup_sorted_num_bam_name = options.working_prefix + ".dup.num.bam";
			}
//			if (boost::filesystem::exists(out_sorted_dup_bam_name) && boost::filesystem::exists(out_sorted_dup_sorted_num_bam_name)) {
//				return;
//			}
		}
//		else {
//			return;
//		}
	}

	if (options.verbose) {
		cout << "Memorizing clipped discordant reads...\n";
	}
	collect_clipped_discordants_serial();
	write_discordants_serial();
}
void ParallelDiscordExtractor::collect_boundaries(const int64_t size_block) {
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_boundaries");
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
	cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] # refs: %d\n") % n_refs).str();
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
		cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
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
				cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
		} else {
			if (a_ref_data.RefLength > size_block) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
				++boundary_id;
				if (verbose) {
					cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
				}
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			if (verbose) {
				cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
			++boundary_id;
		}
	}

	//for(uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
	//auto& a_ref_data = a_ref_vector[ref_id];
	//boundary_positions.push_back(make_pair(ref_id, a_ref_data.RefLength));
	//}

	int64_t calculated_n_blocks = boundary_positions.size();

	cout << "[ParallelDiscordExtractor.collect_boundaries] total ref. size: " << size_total_ref << "\n";

	fixed_size_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	cout << "[ParallelDiscordExtractor.collect_boundaries] collect boundary positions\n";
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
	tasks.push_back([&] {
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
	for (auto& an_offset : unmapped_offsets) {
		if (-1 == an_offset.ref_id) {
			BlockBoundary a_boundary;
			a_boundary.offset = an_offset.offset;
			a_boundary.ref_id = an_offset.ref_id;
			a_boundary.pos = an_offset.position;
			unmapped_included_blocks.push_back(a_boundary);
		}
	}

	cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] actual # blocks: %d\n") % fixed_size_blocks.size()).str();
	if (verbose) {
		for (uint64_t block_id = 0; block_id < fixed_size_blocks.size(); ++block_id) {
			auto& a_block = fixed_size_blocks[block_id];
			cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] BlockBoundary(%d): %d %s %d-%d(%d)\n") % block_id % a_block.offset % a_block.read_name % a_block.ref_id % a_block.pos % a_block.jump_pos).str();
		}
	}
	independent_blocks = fixed_size_blocks;
	cout << checker;
}
void ParallelDiscordExtractor::collect_boundaries_alt(vector<BlockOffset>& offset_blocks) {
	string a_path(options.infile_name);
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
	cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries_alt] # blocks: %d\n") % offset_blocks.size()).str();

}
void ParallelDiscordExtractor::collect_independent_boundaries() {
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_boundaries");
	if (silent) {
		checker.start_without_output();
	} else {
		checker.start();
	}
	BamTools::BamAlignment al;
	bool has_index = reader.HasIndex();
	if (!has_index) {
		cout << "[ParallelDiscordExtractor.collect_boundaries] no index is found and generating.";
		reader.CreateIndex();
	}

	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	if (!silent) {
		cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] # refs: %d\n") % n_refs).str();
	}
	int64_t size_total_ref = 0;
	for (int32_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		auto& a_ref_data = a_ref_vector[ref_id];
		//		cout << a_ref_data.RefName << ": " << a_ref_data.RefLength << "\n";
		size_total_ref += a_ref_data.RefLength;
	}

	int64_t estimated_n_blocks = size_total_ref / (double) size_block;
	++estimated_n_blocks;

	// calculate the positions of boundary entries.
	vector<pair<int32_t, int32_t>> boundary_positions;
	boundary_positions.push_back(make_pair(0, 0));
	if (!silent) {
		if (verbose) {
			cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
		}
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
				cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
		} else {
			if (a_ref_data.RefLength > size_block) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
				++boundary_id;
				if (verbose) {
					cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
				}
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			if (verbose) {
				cout << "[ParallelDiscordExtractor.collect_boundaries] Inserted-" << (boundary_positions.size() - 1) << ":" << boundary_positions.back().first << "/" << boundary_positions.back().second << "\n";
			}
			++boundary_id;
		}
	}

	//		for(uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
	//			auto& a_ref_data = a_ref_vector[ref_id];
	//			boundary_positions.push_back(make_pair(ref_id, a_ref_data.RefLength));
	//		}

	int64_t calculated_n_blocks = boundary_positions.size();

	if (!silent) {
		cout << "[ParallelDiscordExtractor.collect_boundaries] estimated # blocks: " << estimated_n_blocks << "\n";
		cout << "[ParallelDiscordExtractor.collect_boundaries] calculated # blocks: " << calculated_n_blocks << "\n";
		cout << "[ParallelDiscordExtractor.collect_boundaries] total ref. size: " << size_total_ref << "\n";
	}

	independent_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	if (!silent) {
		cout << "[ParallelDiscordExtractor.collect_boundaries] collect boundary positions\n";
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
						// ensure the boundary element does not have the same ref id and position with the previous element
						//									int64_t target_ref_id = local_alignment_entry.RefID;
						//									int64_t target_pos = local_alignment_entry.Position;
						//									int64_t current_ref_id = target_ref_id;
						//									int64_t current_pos = target_pos;
						//									do {
						//										success = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
						//										if(!success) {
						//											break;
						//										}
						//										current_ref_id = local_alignment_entry.RefID;
						//										current_pos = local_alignment_entry.Position;
						//									}while(target_ref_id == current_ref_id && target_pos == current_pos);
				independent_blocks[block_id].read_name = local_alignment_entry.Name;
				independent_blocks[block_id].ref_id = local_alignment_entry.RefID;
				independent_blocks[block_id].pos = local_alignment_entry.Position;
				independent_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
				independent_blocks[block_id].jump_pos = boundary_positions[block_id].second;
			}
			else {
				cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] Failed: %s %d:%d\n")
						% local_alignment_entry.Name
						% local_alignment_entry.RefID
						% local_alignment_entry.Position).str();
				break;
			}
			if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
				cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] unaligned: %s %d:%d\n")
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
	sort(independent_blocks.begin(), independent_blocks.end());
	independent_blocks.erase(unique(independent_blocks.begin(), independent_blocks.end()), independent_blocks.end());
	BlockBoundary a_block_boundary;
	a_block_boundary.read_name = "last";
	a_block_boundary.ref_id = n_refs;
	a_block_boundary.pos = 0;
	a_block_boundary.jump_pos = -1;
	independent_blocks.push_back(a_block_boundary);

	if (!silent) {
		cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] actual # blocks: %d\n") % independent_blocks.size()).str();
	}
	if (verbose) {
		for (uint64_t block_id = 0; block_id < independent_blocks.size(); ++block_id) {
			auto a_block = independent_blocks[block_id];
			cout << (boost::format("[ParallelDiscordExtractor.collect_boundaries] BlockBoundary(%d): %s %d-%d\n") % block_id % a_block.read_name % a_block.ref_id % a_block.pos).str();
		}
	}
	if (!silent) {
		cout << checker;
	}
}
void ParallelDiscordExtractor::collect_chromosome_wide_boundaries() {
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_chromosome_wide_boundaries");
	if (silent) {
		checker.start_without_output();
	} else {
		checker.start();
	}
	BamTools::BamAlignment al;
	bool has_index = reader.HasIndex();
	if (!has_index) {
		cout << "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] no index is found and generating.";
		reader.CreateIndex();
	}

	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	if (!silent) {
		cout << (boost::format("[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] # refs: %d\n") % n_refs).str();
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
//					<< "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] Inserted-"
//					<< (boundary_positions.size() - 1) << ":"
//					<< boundary_positions.back().first << "/"
//					<< boundary_positions.back().second << "\n";
//				}
//			}
	for (int64_t ref_id = 0; ref_id < n_refs; ++ref_id) {
		boundary_positions.push_back(make_pair(ref_id, 0));
	}
//			int64_t boundary_id = 0;
//			int32_t ref_id = 0;
//			int64_t n_remaining_bases = a_ref_vector[ref_id].RefLength;
//			while (n_remaining_bases >= 0) {
//				auto& a_ref_data = a_ref_vector[ref_id];
//				int64_t last_base_pos = boundary_positions[boundary_id].second;
//				n_remaining_bases = a_ref_data.RefLength - last_base_pos;
//				if (n_remaining_bases >= size_block) {
//					boundary_positions.push_back(
//							make_pair(ref_id, last_base_pos + size_block));
//					++boundary_id;
//					if (verbose) {
//						cout
//						<< "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] Inserted-"
//						<< (boundary_positions.size() - 1) << ":"
//						<< boundary_positions.back().first << "/"
//						<< boundary_positions.back().second << "\n";
//					}
//				} else {
//					if (a_ref_data.RefLength > size_block) {
//						boundary_positions.push_back(
//								make_pair(ref_id,
//										min(
//												static_cast<int32_t>(last_base_pos
//														+ size_block),
//														a_ref_data.RefLength)));
//						++boundary_id;
//						if (verbose) {
//							cout
//							<< "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] Inserted-"
//							<< (boundary_positions.size() - 1) << ":"
//							<< boundary_positions.back().first << "/"
//							<< boundary_positions.back().second << "\n";
//						}
//					}
//
//					++ref_id;
//					if (ref_id >= n_refs) {
//						break;
//					}
//					boundary_positions.push_back(make_pair(ref_id, 0));
//					if (verbose) {
//						cout
//						<< "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] Inserted-"
//						<< (boundary_positions.size() - 1) << ":"
//						<< boundary_positions.back().first << "/"
//						<< boundary_positions.back().second << "\n";
//					}
//					++boundary_id;
//				}
//			}

	//		for(uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
	//			auto& a_ref_data = a_ref_vector[ref_id];
	//			boundary_positions.push_back(make_pair(ref_id, a_ref_data.RefLength));
	//		}

	int64_t calculated_n_blocks = boundary_positions.size();

	if (!silent) {
//				cout
//				<< "[ParallelDiscordExtractor.collect_boundaries] estimated # blocks: "
//				<< estimated_n_blocks << "\n";
//				cout
//				<< "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] calculated # blocks: "
//				<< calculated_n_blocks << "\n";
		cout << "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] total ref. size: " << size_total_ref << "\n";
	}

	chromosome_wide_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	if (!silent) {
		cout << "[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] collect boundary positions\n";
	}
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			BamAlignment local_alignment_entry;
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
						cout << (boost::format("[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						cout << (boost::format("[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] unaligned: %s %d:%d\n")
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
		cout << (boost::format("[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] actual # blocks: %d\n") % chromosome_wide_blocks.size()).str();
	}
	if (verbose) {
		for (uint64_t block_id = 0; block_id < chromosome_wide_blocks.size(); ++block_id) {
			auto a_block = chromosome_wide_blocks[block_id];
			cout << (boost::format("[ParallelDiscordExtractor.collect_chromosome_wide_boundaries] BlockBoundary(%d): %s %d-%d\n") % block_id % a_block.read_name % a_block.ref_id % a_block.pos).str();
		}
	}
	if (!silent) {
		cout << checker;
	}
}

void ParallelDiscordExtractor::collect_clipped_discordants() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_clipped_discordants");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	vector<BlockBoundary> independent_blocks = unmapped_included_blocks;
	int64_t calculated_n_blocks = independent_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');
	vector<unordered_map<string, bool>> memo_clipped_disc_lists(calculated_n_blocks - 1);
	vector<int64_t> not_paired_lists(calculated_n_blocks - 1);
	vector<int64_t> not_both_mapped_lists(calculated_n_blocks - 1);
	vector<int64_t> same_pos_lists(calculated_n_blocks - 1);
	vector<int64_t> not_same_chr_lists(calculated_n_blocks - 1);
	vector<int64_t> not_diff_strands_lists(calculated_n_blocks - 1);
	vector<int64_t> weird_insert_size_lists(calculated_n_blocks - 1);
	vector<int64_t> big_insert_size_lists(calculated_n_blocks - 1);
	vector<int64_t> num_memorized_clipped_discs_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;
//			int32_t the_current_ref_id = independent_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = independent_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = independent_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = independent_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = independent_blocks[block_id + 1].aln_flag;
				string the_next_block_read_name = independent_blocks[block_id + 1].read_name;

				string str_block_id = boost::lexical_cast<string>(block_id);

//				bool jump_success = local_reader.Jump(independent_blocks[block_id].ref_id, independent_blocks[block_id].jump_pos);
//				if(!jump_success) {
//					cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] (Jump fail) Block-%d (%d/%d)-(%d/%d)\n")
//					% block_id % the_current_ref_id % the_current_ref_pos
//					% the_next_ref_id % the_next_ref_pos).str();
//					exit(1);
//				}
				int64_t the_current_ref_offset = independent_blocks[block_id].offset;
				int64_t the_next_ref_offset = independent_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;
//				if(cur_offset != the_current_ref_offset) {
//				cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] block-%d (start) offset: jump: %d, calc: %d\n")
//				% block_id % cur_offset % the_current_ref_offset).str();
//				}

//				if(verbose) {
//					cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] (start) Block-%d (%d/%d)-(%d/%d)\n")
//							% block_id % the_current_ref_id % the_current_ref_pos
//							% the_next_ref_id % the_next_ref_pos).str();
//				}
				string rg_name;
				BamAlignment local_alignment_entry;

				auto& local_memo_clipped_disc = memo_clipped_disc_lists[block_id];
				int64_t num_memorized_clipped_discs = 0;
				int64_t not_paired = 0;
				int64_t not_both_mapped = 0;
				int64_t same_pos = 0;
				int64_t not_same_chr = 0;
				int64_t not_diff_strands = 0;
				int64_t weird_insert_size = 0;
				int64_t big_insert_size = 0;
				/* Use the slow method to read tags so we can get the read group */
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(!silent) {
//						if(verbose && 0 == num_total) {
//							string a_block_boundary_str = (boost::format("%s %d-%d %d")
//									% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//									% local_alignment_entry.AlignmentFlag).str();
//							cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] (first) Block-%d %s\n")
//									% block_id % a_block_boundary_str).str();
//						}
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
					++num_total;
					if(!local_alignment_entry.GetReadGroup(rg_name)) {
						rg_name = "none";
					}
					if(black_listed.end() != black_listed.find(rg_name)) {
						continue;
					}
					/* Ignore reads flagged as duplicate if option is specified */
					if (options.filter_dups_by_flag && local_alignment_entry.IsDuplicate()) {
						continue;
					}
					auto& p = options.isinfo[rg_name];
					if (is_discordant(local_alignment_entry, p.first, p.second, options.nstdevs, not_paired,
									not_both_mapped, same_pos, not_same_chr, not_diff_strands, weird_insert_size,
									big_insert_size) && is_clipped(local_alignment_entry)) {
						++num_memorized_clipped_discs;
						local_memo_clipped_disc[local_alignment_entry.Name] = true;
					}
				}
				local_reader.Close();
				done_vector[block_id] = 'D';
				not_paired_lists[block_id] = not_paired;
				not_both_mapped_lists[block_id] = not_both_mapped;
				same_pos_lists[block_id] = same_pos;
				not_same_chr_lists[block_id] = not_same_chr;
				not_diff_strands_lists[block_id] = not_diff_strands;
				weird_insert_size_lists[block_id] = weird_insert_size;
				big_insert_size_lists[block_id] = big_insert_size;
				num_memorized_clipped_discs_lists[block_id] = num_memorized_clipped_discs;
//				if(prev_offset != the_next_ref_offset) {
//					cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n") % block_id % prev_offset % cur_offset % the_next_ref_offset).str();
//				}

				if(verbose) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] (last) Block-%d %s\n")
							% block_id % a_block_boundary_str).str();
				}
//				else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[ParallelDiscordExtractor.collect_clipped_discordants] gathers scattered information\n";
	int64_t memo_clipped_disc_size = 0;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		auto& local_memo_clipped_disc = memo_clipped_disc_lists[block_id];
		memo_clipped_disc_size += local_memo_clipped_disc.size();
	}
	//		memo_clipped_disc.resize(memo_clipped_disc_size);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		auto& local_memo_clipped_disc = memo_clipped_disc_lists[block_id];
		memo_clipped_disc.insert(local_memo_clipped_disc.begin(), local_memo_clipped_disc.end());
		//			local_memo_clipped_disc.clear();

		not_paired += not_paired_lists[block_id];
		not_both_mapped += not_both_mapped_lists[block_id];
		same_pos += same_pos_lists[block_id];
		not_same_chr += not_same_chr_lists[block_id];
		not_diff_strands += not_diff_strands_lists[block_id];
		weird_insert_size += weird_insert_size_lists[block_id];
		big_insert_size += big_insert_size_lists[block_id];
		num_memorized_clipped_discs += num_memorized_clipped_discs_lists[block_id];
	}
	if (options.verbose) {
		cout << (boost::format("recorded %d/%d clipped discordant read IDs in memory\n") % num_memorized_clipped_discs % memo_clipped_disc.size()).str();
	}
	cout << checker;
}

void ParallelDiscordExtractor::collect_clipped_discordants_serial() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_clipped_discordants");
	checker.start();
	vector<function<void()> > tasks;
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		return;
	}
	int64_t num_total = 0;
	string rg_name;
	BamAlignment local_alignment_entry;

	auto& local_memo_clipped_disc = memo_clipped_disc;
	/* Use the slow method to read tags so we can get the read group */
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		++num_total;
		if(!local_alignment_entry.GetReadGroup(rg_name)) {
			rg_name = "none";
		}
		if(black_listed.end() != black_listed.find(rg_name)) {
			continue;
		}
		/* Ignore reads flagged as duplicate if option is specified */
		if (options.filter_dups_by_flag && local_alignment_entry.IsDuplicate()) {
			continue;
		}
		auto& p = options.isinfo[rg_name];
		if (is_discordant(local_alignment_entry, p.first, p.second, options.nstdevs, not_paired,
						not_both_mapped, same_pos, not_same_chr, not_diff_strands, weird_insert_size,
						big_insert_size) && is_clipped(local_alignment_entry)) {
			++num_memorized_clipped_discs;
			local_memo_clipped_disc[local_alignment_entry.Name] = true;
		}
	}
	local_reader.Close();

	if (options.verbose) {
		cout << (boost::format("recorded %d/%d clipped discordant read IDs in memory\n") % num_memorized_clipped_discs % memo_clipped_disc.size()).str();
	}
	cout << checker;
}

void ParallelDiscordExtractor::reset_stats() {
	not_paired = 0;
	not_both_mapped = 0;
	not_same_chr = 0;
	not_diff_strands = 0;
	weird_insert_size = 0;
	big_insert_size = 0;
	same_pos = 0;
}


void ParallelDiscordExtractor::write_discordants_alt() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtrator.write_discordants_alt");
	checker.start();
	vector<function<void()> > tasks;

	vector<BlockBoundary> independent_blocks = unmapped_included_blocks;
	int64_t calculated_n_blocks = independent_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');
	vector<int64_t> num_in_lists(calculated_n_blocks - 1);
	vector<int64_t> num_blacklisted_lists(calculated_n_blocks - 1);
	vector<int64_t> num_out_lists(calculated_n_blocks - 1);
	vector<int64_t> num_pcr_marked_dups_lists(calculated_n_blocks - 1);
	vector<int64_t> not_paired_lists(calculated_n_blocks - 1);
	vector<int64_t> not_both_mapped_lists(calculated_n_blocks - 1);
	vector<int64_t> same_pos_lists(calculated_n_blocks - 1);
	vector<int64_t> not_same_chr_lists(calculated_n_blocks - 1);
	vector<int64_t> not_diff_strands_lists(calculated_n_blocks - 1);
	vector<int64_t> weird_insert_size_lists(calculated_n_blocks - 1);
	vector<int64_t> big_insert_size_lists(calculated_n_blocks - 1);
	vector<int64_t> num_memorized_clipped_discs_lists(calculated_n_blocks - 1);
	reset_stats(); /* reset counters associated with is_discordant */
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
//					continue;
			}
			const RefVector& refs = local_reader.GetReferenceData();
			int64_t num_in = 0;
			int64_t num_out = 0;
			int64_t num_blacklisted = 0;
			int64_t num_pcr_marked_dups = 0;
//			int32_t the_current_ref_id = independent_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = independent_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = independent_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = independent_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = independent_blocks[block_id + 1].aln_flag;
			string the_next_block_read_name = independent_blocks[block_id + 1].read_name;

			string str_block_id = boost::lexical_cast<string>(block_id);

//				bool jump_success = local_reader.Jump(independent_blocks[block_id].ref_id, independent_blocks[block_id].jump_pos);
//				if(!jump_success) {
//					cout << (boost::format("[ParallelDiscordExtrator.write_discordants_alt] (Jump fail) Block-%d (%d/%d)-(%d/%d)\n")
//					% block_id % the_current_ref_id % the_current_ref_pos
//					% the_next_ref_id % the_next_ref_pos).str();
//					local_reader.Close();
//					return;
//				}
			int64_t the_current_ref_offset = independent_blocks[block_id].offset;
			int64_t the_next_ref_offset = independent_blocks[block_id + 1].offset;

			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;

//			if(verbose) {
//				cout << (boost::format("[ParallelDiscordExtrator.write_discordants_alt] (start) Block-%d (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			string rg_name;
			BamAlignment local_alignment_entry;

			int64_t num_memorized_clipped_discs = 0;
			int64_t not_paired = 0;
			int64_t not_both_mapped = 0;
			int64_t same_pos = 0;
			int64_t not_same_chr = 0;
			int64_t not_diff_strands = 0;
			int64_t weird_insert_size = 0;
			int64_t big_insert_size = 0;

			/* Group stores up reads that start at the same genomic coordinate */
			/* Only used to detect PCR duplicates */
			Group g;
			/* Purpose of prev_dups is to deal with the fact that PCR duplicates
			 * are defined by having the same mapping positions on the forward
			 * AND reverse read.  That is, PCR duplicate blocks have two "ends":
			 * one with all of the forward reads and one with all of the reverse
			 * reads.  Depending on the direction of the original read, we don't
			 * know which of these blocks we will run into first (because we will
			 * see the one with the smallest genomic coordinate first).  In any
			 * case, when we do see the first block, we will have to arbitrarily
			 * choose which read to keep.  When we witness the other end of the
			 * duplicate block, we want to make sure to choose the mate that
			 * corresponds to the SAME read we chose in the initial block.  The
			 * prev_dups structure records which read was chosen, by its name,
			 * for every coordinate at which a PCR duplicate block was witnessed.
			 */
			record prev_dups;
			string local_outfile_name = options.outfile_name + "." + str_block_id;
			BamWriter writer;
			BamWriter dups;
			const bool dup_file_write = options.filter_dups;
			if(0 == block_id) {
				if(!writer.SAMOpen(local_outfile_name, local_reader.GetHeaderText(), local_reader.GetReferenceData())) {
					cout << "ERROR: could not open output BAM file '" << local_outfile_name << "' for writing\n";
					exit(1);
				}
				if (dup_file_write) {
					string local_dupfile_name = options.dupfile_name + "." + str_block_id;
					if (!dups.SAMOpen(local_dupfile_name, local_reader.GetHeaderText(), local_reader.GetReferenceData())) {
						cout << "ERROR: could not open dup file '" << local_dupfile_name << "' for writing\n";
						exit(1);
					}
				}
			} else {
				if(!writer.SAMOpenNoHeader(local_outfile_name, local_reader.GetReferenceData())) {
					cout << "ERROR: could not open output BAM file '" << local_outfile_name << "' for writing\n";
					exit(1);
				}
				if (dup_file_write) {
					string local_dupfile_name = options.dupfile_name + "." + str_block_id;
					if (!dups.SAMOpenNoHeader(local_dupfile_name, local_reader.GetReferenceData())) {
						cout << "ERROR: could not open dup file '" << local_dupfile_name << "' for writing\n";
						exit(1);
					}
				}
			}
//			if (!writer.Open(local_outfile_name, local_reader.GetHeaderText(), local_reader.GetReferenceData())) {
//				cout << "ERROR: could not open output BAM file '"
//				<< local_outfile_name << "' for writing\n";
//				exit(1);
//			}


			/* Use the slow method to read tags so we can get the read group */
			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//				if(verbose && 0 == num_in) {
//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//							% local_alignment_entry.AlignmentFlag).str();
//					cout << (boost::format("[ParallelDiscordExtrator.write_discordants_alt] (first) Block-%d %s\n")
//							% block_id % a_block_boundary_str).str();
//				}
//				if(local_alignment_entry.RefID == the_next_ref_id
//						&& local_alignment_entry.Position == the_next_ref_pos
//						&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//						&& local_alignment_entry.Name == the_next_block_read_name
//				) {
//					break;
//				}
				cur_offset = m_bgzf.Tell();
				if(prev_offset >= the_next_ref_offset) {
					break;
				}
				prev_offset = cur_offset;

				++num_in;
				/* Ignore reads flagged as duplicate if option is specified */
				if (options.filter_dups_by_flag && local_alignment_entry.IsDuplicate()) {
					++num_pcr_marked_dups;
					continue;
				}
				if(!local_alignment_entry.GetReadGroup(rg_name)) {
					rg_name = "none";
				}
				if(options.rg_blacklist.end() != options.rg_blacklist.find(rg_name)) {
					continue;
				}
				auto rg_name_itr = options.isinfo.find(rg_name);
				if(options.isinfo.end() == rg_name_itr) {
					continue;
				}
				auto& p = rg_name_itr->second;
				auto memit = memo_clipped_disc.find(local_alignment_entry.Name);
				/* not discordant and not in the memo buffer */
				if (!is_discordant(local_alignment_entry, p.first, p.second, options.nstdevs, not_paired,
								not_both_mapped, same_pos, not_same_chr, not_diff_strands, weird_insert_size,
								big_insert_size) && memo_clipped_disc.end() == memit) {
					continue;
				}

				if (is_blacklisted(local_alignment_entry, refs, options.blacklist)) {
					++num_blacklisted;
					continue;
				}

				/* Handle PCR dup filtering if specified, else write */
				if (dup_file_write) {
					if (g.should_write(local_alignment_entry)) {
						num_out += g.write(writer, prev_dups, dups, dup_file_write);
//						num_out += g.write_alt(writer, prev_dups, dups, dup_file_write);
						g.clear();
					}
					g.add(local_alignment_entry);
				} else {
//					writer.SaveAlignment(local_alignment_entry);
					writer.SaveSAMAlignment(local_alignment_entry);
					++num_out;
				}
			}
			if (dup_file_write) {
				num_out += g.write(writer, prev_dups, dups, dup_file_write);
//				num_out += g.write_alt(writer, prev_dups, dups, dup_file_write);
			}
			local_reader.Close();
			writer.Close();
			dups.Close();
			done_vector[block_id] = 'D';
			num_in_lists[block_id] = num_in;
			num_blacklisted_lists[block_id] = num_blacklisted;
			num_out_lists[block_id] = num_out;
			num_pcr_marked_dups_lists[block_id] = num_pcr_marked_dups;
			not_paired_lists[block_id] = not_paired;
			not_both_mapped_lists[block_id] = not_both_mapped;
			same_pos_lists[block_id] = same_pos;
			not_same_chr_lists[block_id] = not_same_chr;
			not_diff_strands_lists[block_id] = not_diff_strands;
			weird_insert_size_lists[block_id] = weird_insert_size;
			big_insert_size_lists[block_id] = big_insert_size;
			num_memorized_clipped_discs_lists[block_id] = num_memorized_clipped_discs;
			if(verbose) {
				string a_block_boundary_str = (boost::format("%s %d-%d %d")
						% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
						% local_alignment_entry.AlignmentFlag).str();
				cout << (boost::format("[ParallelDiscordExtrator.write_discordants] (last) Block-%d %s\n")
						% block_id % a_block_boundary_str).str();
			}
//			else {
//				size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//				double processed = n/(double)done_vector.size() * 100.0;
//				cout << (boost::format("%.2f %%\n") % processed).str();
//			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[ParallelDiscordExtrator.write_discordants] gathers scattered information\n";
	//		int64_t memo_clipped_disc_size = 0;
	int64_t num_in = 0;
	int64_t num_blacklisted = 0;
	int64_t num_out = 0;
	int64_t num_pcr_marked_dups = 0;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		num_in += num_in_lists[block_id];
		num_blacklisted += num_blacklisted_lists[block_id];
		num_out += num_out_lists[block_id];
		num_pcr_marked_dups += num_pcr_marked_dups_lists[block_id];
		not_paired += not_paired_lists[block_id];
		not_both_mapped += not_both_mapped_lists[block_id];
		same_pos += same_pos_lists[block_id];
		not_same_chr += not_same_chr_lists[block_id];
		not_diff_strands += not_diff_strands_lists[block_id];
		weird_insert_size += weird_insert_size_lists[block_id];
		big_insert_size += big_insert_size_lists[block_id];
		num_memorized_clipped_discs += num_memorized_clipped_discs_lists[block_id];
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			string local_outfile_name = options.outfile_name + "."
//					+ str_block_id;
//			string local_dupfile_name = options.dupfile_name + "."
//					+ str_block_id;
	}

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	if (options.verbose) {
		cout << "Finished file scan.\n";
	}

	cout << options.infile_name << ": " << num_in << " reads read\n";
	if (!options.bfile_name.empty()) {
		cout << options.bfile_name << ": rejected " << num_blacklisted << " discordant pairs\n";
	}
	cout << options.outfile_name << ": " << num_out << " reads written\n";
	cout << "   not paired:             " << not_paired //
			<< "\n   not both mapped:        " << not_both_mapped //
			<< "\n   not same chrom:         " << not_same_chr //
			<< "\n   not diff strand:        " << not_diff_strands //
			<< "\n   weird insert size:      " << weird_insert_size //
			<< "\n   big insert size:        " << big_insert_size //
			<< "\n   mates map to same posn: " << same_pos << "\n";
	cout << "discarded " << (not_same_chr + not_diff_strands + weird_insert_size + big_insert_size - num_out - num_blacklisted - same_pos) << " PCR duplicates\n";
	if (options.filter_dups_by_flag) {
		cout << "   -m PCR dups ignored:   " << num_pcr_marked_dups << "\n";
	}
	cout << checker;
}


void ParallelDiscordExtractor::write_discordants_serial() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtrator.write_discordants_serial");
	checker.start();
	vector<function<void()> > tasks;

	reset_stats(); /* reset counters associated with is_discordant */
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		return;
	}
	const RefVector& refs = local_reader.GetReferenceData();
	int64_t num_in = 0;
	int64_t num_out = 0;
	int64_t num_blacklisted = 0;
	int64_t num_pcr_marked_dups = 0;

	string rg_name;
	BamAlignment local_alignment_entry;

	/* Group stores up reads that start at the same genomic coordinate */
	/* Only used to detect PCR duplicates */
	Group g;
	/* Purpose of prev_dups is to deal with the fact that PCR duplicates
	 * are defined by having the same mapping positions on the forward
	 * AND reverse read.  That is, PCR duplicate blocks have two "ends":
	 * one with all of the forward reads and one with all of the reverse
	 * reads.  Depending on the direction of the original read, we don't
	 * know which of these blocks we will run into first (because we will
	 * see the one with the smallest genomic coordinate first).  In any
	 * case, when we do see the first block, we will have to arbitrarily
	 * choose which read to keep.  When we witness the other end of the
	 * duplicate block, we want to make sure to choose the mate that
	 * corresponds to the SAME read we chose in the initial block.  The
	 * prev_dups structure records which read was chosen, by its name,
	 * for every coordinate at which a PCR duplicate block was witnessed.
	 */
	record prev_dups;
	string local_outfile_name = options.outfile_name;
	BamWriter writer;
	BamWriter dups;
	const bool dup_file_write = options.filter_dups;
	if(!writer.SAMOpen(local_outfile_name, local_reader.GetHeaderText(), local_reader.GetReferenceData())) {
		cout << "ERROR: could not open output BAM file '" << local_outfile_name << "' for writing\n";
		exit(1);
	}
	if (dup_file_write) {
		string local_dupfile_name = options.dupfile_name;
		if (!dups.SAMOpen(local_dupfile_name, local_reader.GetHeaderText(), local_reader.GetReferenceData())) {
			cout << "ERROR: could not open dup file '" << local_dupfile_name << "' for writing\n";
			exit(1);
		}
	}


	/* Use the slow method to read tags so we can get the read group */
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		++num_in;
		/* Ignore reads flagged as duplicate if option is specified */
		if (options.filter_dups_by_flag && local_alignment_entry.IsDuplicate()) {
			++num_pcr_marked_dups;
			continue;
		}
		if(!local_alignment_entry.GetReadGroup(rg_name)) {
			rg_name = "none";
		}
		if(options.rg_blacklist.end() != options.rg_blacklist.find(rg_name)) {
			continue;
		}
		auto rg_name_itr = options.isinfo.find(rg_name);
		if(options.isinfo.end() == rg_name_itr) {
			continue;
		}
		auto& p = rg_name_itr->second;
		auto memit = memo_clipped_disc.find(local_alignment_entry.Name);
		/* not discordant and not in the memo buffer */
		if (!is_discordant(local_alignment_entry, p.first, p.second, options.nstdevs, not_paired,
						not_both_mapped, same_pos, not_same_chr, not_diff_strands, weird_insert_size,
						big_insert_size) && memo_clipped_disc.end() == memit) {
			continue;
		}
		if (is_blacklisted(local_alignment_entry, refs, options.blacklist)) {
			++num_blacklisted;
			continue;
		}

		/* Handle PCR dup filtering if specified, else write */
		if (dup_file_write) {
			if (g.should_write(local_alignment_entry)) {
				num_out += g.write(writer, prev_dups, dups, dup_file_write);
//						num_out += g.write_alt(writer, prev_dups, dups, dup_file_write);
				g.clear();
			}
			g.add(local_alignment_entry);
		} else {
//					writer.SaveAlignment(local_alignment_entry);
			writer.SaveSAMAlignment(local_alignment_entry);
			++num_out;
		}
	}
	if (dup_file_write) {
		num_out += g.write(writer, prev_dups, dups, dup_file_write);
//				num_out += g.write_alt(writer, prev_dups, dups, dup_file_write);
	}
	local_reader.Close();
	writer.Close();
	dups.Close();

	if (options.verbose) {
		cout << "Finished file scan.\n";
	}

	cout << options.infile_name << ": " << num_in << " reads read\n";
	if (!options.bfile_name.empty()) {
		cout << options.bfile_name << ": rejected " << num_blacklisted << " discordant pairs\n";
	}
	cout << options.outfile_name << ": " << num_out << " reads written\n";
	cout << "   not paired:             " << not_paired //
			<< "\n   not both mapped:        " << not_both_mapped //
			<< "\n   not same chrom:         " << not_same_chr //
			<< "\n   not diff strand:        " << not_diff_strands //
			<< "\n   weird insert size:      " << weird_insert_size //
			<< "\n   big insert size:        " << big_insert_size //
			<< "\n   mates map to same posn: " << same_pos << "\n";
	cout << "discarded " << (not_same_chr + not_diff_strands + weird_insert_size + big_insert_size - num_out - num_blacklisted - same_pos) << " PCR duplicates\n";
	if (options.filter_dups_by_flag) {
		cout << "   -m PCR dups ignored:   " << num_pcr_marked_dups << "\n";
	}
	cout << checker;
}


void ParallelDiscordExtractor::collect_discordants() {
	// opens the BAM files (by default without checking for indexes)
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtrator.collect_discordants");
	if (silent) {
		checker.start_without_output();
	} else {
		checker.start();
	}
	//		const bool verbose = false;
	vector<function<void()> > tasks;
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = independent_blocks.size();
	//		int64_t calculated_n_blocks = 833;
	int64_t disc_file_size = 0;
	int64_t dup_file_size = 0;
	vector<string> disc_files;
	vector<string> dup_files;
	vector<int64_t> disc_boundary_pos;
	vector<int64_t> dup_boundary_pos;
	disc_boundary_pos.push_back(0);
	dup_boundary_pos.push_back(0);
	for (int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string local_outfile_name = options.outfile_name + "." + str_block_id;
		string local_dupfile_name = options.dupfile_name + "." + str_block_id;
		disc_files.push_back(local_outfile_name);

		int64_t local_disc_file_size = boost::filesystem::file_size(local_outfile_name);
		disc_file_size += local_disc_file_size;
		disc_boundary_pos.push_back(disc_file_size);
		cout << disc_boundary_pos.back() << "\n";
		if (!options.dupfile_name.empty()) {
			dup_files.push_back(local_dupfile_name);

			int64_t local_dup_file_size = boost::filesystem::file_size(local_dupfile_name);
			dup_file_size += local_dup_file_size;
			dup_boundary_pos.push_back(dup_file_size);
		}
	}

	{
		ofstream out(options.outfile_name + ".sam", ios::binary | ios::ate);
		out.seekp(disc_file_size - 1, ios::beg);
		out << '\n';
	}
	{
		ofstream out(options.dupfile_name + ".sam", ios::binary | ios::ate);
		out.seekp(dup_file_size - 1, ios::beg);
		out << '\n';
	}
	{
		for (int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
			tasks.push_back([&, block_id] {
				string str_block_id = boost::lexical_cast<string>(block_id);
				string local_outfile_name = options.outfile_name + "."
				+ str_block_id;
				int64_t skip_pos = disc_boundary_pos[block_id];
				fstream out(options.outfile_name + ".sam", ios::in | ios::out | ios::binary|ios::ate);
				out.seekp(skip_pos, ios::beg);
				out << castle::IOUtils::read_fully(local_outfile_name);
			});
			if (!options.dupfile_name.empty()) {
				tasks.push_back([&, block_id] {
					string str_block_id = boost::lexical_cast<string>(block_id);
					string local_dupfile_name = options.dupfile_name + "."
					+ str_block_id;
					int64_t skip_pos = dup_boundary_pos[block_id];
					fstream out(options.dupfile_name + ".sam", ios::in | ios::out | ios::binary|ios::ate);
					out.seekp(skip_pos, ios::beg);
					out << castle::IOUtils::read_fully(local_dupfile_name);
				});
			}
			//				string str_block_id = boost::lexical_cast<string>(block_id);
			//				string local_outfile_name = options.outfile_name + "."
			//						+ str_block_id;
			//				string local_dupfile_name = options.dupfile_name + "."
			//						+ str_block_id;
			//				disc_files.push_back(local_outfile_name);
			//
			//				int64_t local_disc_file_size = boost::filesystem::file_size(local_outfile_name);
			//				disc_file_size += local_disc_file_size;
			//				disc_boundary_pos.push_back(disc_file_size);
			//				if(!options.dupfile_name.empty()) {
			//					dup_files.push_back(local_dupfile_name);
			//
			//					int64_t local_dup_file_size = boost::filesystem::file_size(local_dupfile_name);
			//					dup_file_size += local_dup_file_size;
			//					dup_boundary_pos.push_back(dup_file_size);
			//				}
		}
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (uint64_t f_id = 0; f_id < disc_files.size(); ++f_id) {
		tasks.push_back([&, f_id] {
			auto& a_file_name = disc_files[f_id];
			boost::filesystem::remove(a_file_name);
		});
	}
	if (!options.dupfile_name.empty()) {
		for (uint64_t f_id = 0; f_id < dup_files.size(); ++f_id) {
			tasks.push_back([&, f_id] {
				auto& a_file_name = dup_files[f_id];
				boost::filesystem::remove(a_file_name);
			});
		}
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void ParallelDiscordExtractor::collect_discordants_alt() {
	// opens the BAM files (by default without checking for indexes)
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtrator.collect_discordants_alt");
	if (silent) {
		checker.start_without_output();
	} else {
		checker.start();
	}
	//		const bool verbose = false;
	vector<function<void()> > tasks;
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = independent_blocks.size();

	//		int64_t calculated_n_blocks = 833;
	vector<string> disc_files;
	vector<string> dup_files;
	const bool write_dup_file = !options.dupfile_name.empty();
	for (int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string local_outfile_name = options.outfile_name + "." + str_block_id;
		string local_dupfile_name = options.dupfile_name + "." + str_block_id;
		if (castle::IOUtils::get_file_size(local_outfile_name) > 0) {
			disc_files.push_back(local_outfile_name);
		}

		if (write_dup_file) {
			if (castle::IOUtils::get_file_size(local_dupfile_name) > 0) {
				dup_files.push_back(local_dupfile_name);
			}
		}
	}
	tasks.push_back([&] {
		BamMultiReader reader;
		if (!reader.Open(disc_files, false, true)) {
			cerr << "ERROR: Could not open input BAM file(s)... Aborting.\n";
			exit(1);
		}

		// retrieve header & reference dictionary info
			std::string mergedHeader = reader.GetHeaderText();
			RefVector references = reader.GetReferenceData();

			// open writer
			BamWriter writer;
			if ( !writer.Open(options.outfile_name, mergedHeader, references, false) ) {
				cerr << "ERROR: Could not open BAM file " << options.outfile_name << " for writing... Aborting.\n";
				reader.Close();
				exit(1);
			}

			// if no region specified, store entire contents of file(s)
			BamAlignment al;
			while ( reader.GetNextAlignmentCore(al) ) {
				writer.SaveAlignment(al);
			}
		});

	tasks.push_back([&, write_dup_file] {
		if(!write_dup_file) {
			return;
		}
		BamMultiReader reader;
		if (!reader.Open(dup_files, false, true)) {
			cerr << "ERROR: Could not open input BAM file(s)... Aborting.\n";
			exit(1);
		}

		// retrieve header & reference dictionary info
			std::string mergedHeader = reader.GetHeaderText();
			RefVector references = reader.GetReferenceData();

			// open writer
			BamWriter writer;
			if ( !writer.Open(options.dupfile_name, mergedHeader, references) ) {
				cerr << "ERROR: Could not open BAM file " << options.outfile_name << " for writing... Aborting.\n";
				reader.Close();
				exit(1);
			}

			// if no region specified, store entire contents of file(s)
			BamAlignment al;
			while ( reader.GetNextAlignmentCore(al) ) {
				writer.SaveAlignment(al);
			}
		});

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (uint64_t f_id = 0; f_id < disc_files.size(); ++f_id) {
		tasks.push_back([&, f_id] {
			auto& a_file_name = disc_files[f_id];
			boost::filesystem::remove(a_file_name);
		});
	}
	if (!options.dupfile_name.empty()) {
		for (uint64_t f_id = 0; f_id < dup_files.size(); ++f_id) {
			tasks.push_back([&, f_id] {
				auto& a_file_name = dup_files[f_id];
				boost::filesystem::remove(a_file_name);
			});
		}
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void ParallelDiscordExtractor::collect_discordants_sambamba() {
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.collect_discordants_sambamba");
	checker.start();
	string a_path(options.infile_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	vector<BlockBoundary> independent_blocks = unmapped_included_blocks;
	int64_t calculated_n_blocks = independent_blocks.size();
	vector<string> disc_files;
	vector<string> dup_files;
	vector<string> remove_files;
	const bool write_dup_file = !options.dupfile_name.empty();
	for (int64_t block_id = 0; block_id < calculated_n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string local_outfile_name = options.outfile_name + "." + str_block_id;
		string local_dupfile_name = options.dupfile_name + "." + str_block_id;
		if (castle::IOUtils::get_file_size(local_outfile_name) > 0) {
			disc_files.push_back(local_outfile_name);
		} else {
			remove_files.push_back(local_outfile_name);
		}

		if (write_dup_file) {
			if (castle::IOUtils::get_file_size(local_dupfile_name) > 0) {
				dup_files.push_back(local_dupfile_name);
			} else {
				remove_files.push_back(local_dupfile_name);
			}
		}
	}
	string the_merge_out_sam_file = options.outfile_name + ".tmp.sam";
	string the_merge_out_bam_file = options.outfile_name + ".tmp.bam";
	string out_sorted_bam_name = options.prefix + ".disc.sorted.bam";
	string out_sorted_num_bam_name = options.prefix + ".disc.num.bam";

	if(!options.working_dir.empty()) {
		out_sorted_bam_name = options.working_prefix + ".disc.sorted.bam";
		out_sorted_num_bam_name = options.working_prefix + ".disc.num.bam";
	}
	castle::IOUtils::plain_file_merge(the_merge_out_sam_file, disc_files, n_cores, true);
	string samtools_cmd_1 = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % the_merge_out_bam_file % the_merge_out_sam_file).str();
	system(samtools_cmd_1.c_str());
//	string sambamba_cmd_1 = (boost::format("samtools merge -f -1 -@ %d %s %s") % n_cores % the_merge_outfile % (castle::StringUtils::join(disc_files, " "))).str();
//	system(sambamba_cmd_1.c_str());
	string sambamba_sort_cmd_1 = (boost::format("sambamba sort --tmpdir=%s -l 1 -t %d -o %s %s") % options.tmp_directory % n_cores % out_sorted_num_bam_name % the_merge_out_bam_file).str();
	system(sambamba_sort_cmd_1.c_str());

	string sambamba_sort_n_cmd_1 = (boost::format("sambamba sort --tmpdir=%s -N -l 1 -t %d -o %s %s") % options.tmp_directory % n_cores % out_sorted_bam_name % out_sorted_num_bam_name).str();
	system(sambamba_sort_n_cmd_1.c_str());

//	sambamba sort -N -l 1 -t $threads_bwa -o $disc_sort $disc_bam
//	disc_files.push_back(the_merge_out_bam_file);
	castle::IOUtils::remove_files(disc_files, n_cores);
	castle::IOUtils::remove_files(remove_files, n_cores);

	if (write_dup_file) {
		string the_merge_out_dup_sam_file = options.dupfile_name + ".tmp.sam";
		castle::IOUtils::plain_file_merge(the_merge_out_dup_sam_file, dup_files, n_cores, true);
		string the_merge_out_dup_bam_file = options.dupfile_name + ".tmp.bam";
//		string sambamba_cmd_2 = (boost::format("samtools merge -f -1 -@ %d %s %s") % n_cores % the_merge_dupfile % (castle::StringUtils::join(dup_files, " "))).str();
//		system(sambamba_cmd_2.c_str());
		string samtools_cmd_2 = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % the_merge_out_dup_bam_file % the_merge_out_dup_sam_file).str();
		system(samtools_cmd_2.c_str());
		string out_sorted_dup_bam_name = options.prefix + ".dup.sorted.bam";
		string out_sorted_dup_sorted_num_bam_name = options.prefix + ".dup.num.bam";
		if(!options.working_dir.empty()) {
			out_sorted_dup_bam_name = options.working_prefix + ".dup.sorted.bam";
			out_sorted_dup_sorted_num_bam_name = options.working_prefix + ".dup.num.bam";
		}

		string sambamba_sort_cmd_2 = (boost::format("sambamba sort --tmpdir=%s -l 1 -t %d -o %s %s") % options.tmp_directory % n_cores % out_sorted_dup_sorted_num_bam_name % the_merge_out_dup_bam_file).str();
		system(sambamba_sort_cmd_2.c_str());
		string sambamba_sort_n_cmd_2 = (boost::format("sambamba sort --tmpdir=%s -N -l 1 -t %d -o %s %s") % options.tmp_directory % n_cores % out_sorted_dup_bam_name % out_sorted_dup_sorted_num_bam_name).str();
		system(sambamba_sort_n_cmd_2.c_str());
//		dup_files.push_back(the_merge_out_dup_bam_file);
		castle::IOUtils::remove_files(dup_files, n_cores);
	}
	// remove temporary files
	if(options.is_cleaning) {
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
		vector<string> all_removal_names;

		for (auto& rgname : read_groups) {
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
		castle::IOUtils::remove_files(all_removal_names, n_cores);
	}
	cout << checker;
}

void ParallelDiscordExtractor::create_bni_even_index(const string& a_path, string bni_index_path) {
	get_bni_index_path(a_path, bni_index_path);
	if (boost::filesystem::exists(bni_index_path)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("ParallelDiscordExtractor.create_even_index");
	checker.start();


	string bai_index_path;
	get_bai_index_path(a_path, bai_index_path);

	cout << (boost::format("[ParallelDiscordExtractor.create_even_index] path: %s\n") % a_path).str();
	cout << (boost::format("[ParallelDiscordExtractor.create_even_index] bai: %s\n") % bai_index_path).str();
	cout << (boost::format("[ParallelDiscordExtractor.create_even_index] bni: %s\n") % bni_index_path).str();
	BamTools::BamReader serial_reader;
	if (!serial_reader.Open(a_path, bai_index_path)) {
		cout << checker;
		return;
	}

	vector<BamTools::BlockOffset> block_boundary;
	serial_reader.GetBlockOffsets(block_boundary);
	ofstream out(bni_index_path, ios::binary);
	for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
		out << block_boundary[block_id].offset << "\t" << block_boundary[block_id].ref_id << "\t" << block_boundary[block_id].position << "\n";
	}
	cout << (boost::format("[ParallelDiscordExtractor.create_even_index] %d entries were indexed\n") % block_boundary.size()).str();
	cout << checker;
}

} /* namespace meerkat */
