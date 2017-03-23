/*
 * AlternativeMapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#include "../cluster/AlternativeMapper.hpp"

namespace meerkat {

AlternativeMapper::AlternativeMapper() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

AlternativeMapper::~AlternativeMapper() {
}
void AlternativeMapper::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	if (options.blistname.empty()) {
		options.blistname = options.prefix + ".blacklist.gz.bak";
	}
	options.read_point_black_lists();
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
}

// The algorithm select entries in each chromosome mappings, hence the number of thread can be smaller than the number of physical cores.
void AlternativeMapper::select_alternative_mapping() {
	string mapping_raw_file = options.prefix + ".mapping.raw";
	if (!options.working_dir.empty()) {
		mapping_raw_file = options.working_prefix + ".mapping.raw";
	}
//	if(boost::filesystem::exists(mapping_raw_file)) {
//		return;
//	}
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.select_alternative_mapping");
	checker.start();
	string input_file_name = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		input_file_name = options.working_prefix + ".alt.map.sort";
	}
//cout << input_file_name << "\n";
// calculate the positions of boundary entries.
	vector<function<void()> > tasks;
// prepare for the boundary elements

	deque<int64_t> current_boundary_positions;
	{
// all the small cluster with ref id within each 4096000-byte block won't be separated.
// For instance, unlike the large chromosome 1, small chromosomes such as GL000191.1, GL000192.1, ... and hs37d5 will be included in a single large block.
		int64_t size_total_ref = castle::IOUtils::get_file_size(input_file_name);
		int64_t n_remaining_bases = size_total_ref;
		const int64_t size_block = 4096000;
		int64_t boundary_id = 0;
		current_boundary_positions.push_back(0);
		while (n_remaining_bases >= 0) {
			int64_t last_base_pos = current_boundary_positions[boundary_id];
			n_remaining_bases = size_total_ref - last_base_pos;
			if (n_remaining_bases >= size_block) {
				current_boundary_positions.push_back(last_base_pos + size_block);
				++boundary_id;
			} else {
				current_boundary_positions.push_back(size_total_ref);
				break;
			}
		}

		cout << (boost::format("[AlternativeMapper.select_alternative_mapping] collect boundary positions of file size %d\n") % size_total_ref).str();
		int64_t n_boundaries = current_boundary_positions.size();
		for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				ifstream local_reader(input_file_name, ios::binary);
				int64_t skip_bytes = current_boundary_positions[block_id];
				int64_t next_skip_bytes = current_boundary_positions[block_id + 1];
				local_reader.seekg(skip_bytes, ios::beg);
				string line;
				int64_t boundary_bytes = skip_bytes;
				int64_t previous_boundary_bytes = skip_bytes;
				const char* delims = "\t";
				vector<string> data;
				string prev_ref_id;
				bool found_skip_point = false;
				getline(local_reader, line, '\n');
				boundary_bytes += line.size() + 1;
				while(getline(local_reader, line, '\n')) {
					previous_boundary_bytes = boundary_bytes;
					boundary_bytes += line.size() + 1;
					if(boundary_bytes > next_skip_bytes) {
						break;
					}
					castle::StringUtils::c_string_multi_split(line, delims, data);
					if(!prev_ref_id.empty() && data[1] != prev_ref_id) {
						found_skip_point = true;
						break;
					}
					prev_ref_id = data[1];
				}
				if(found_skip_point) {
					current_boundary_positions[block_id] = previous_boundary_bytes;
				} else {
					current_boundary_positions[block_id] = -1;
				}
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		current_boundary_positions.erase(remove_if(current_boundary_positions.begin(), current_boundary_positions.end(), [](int64_t value) {return -1 == value;}), current_boundary_positions.end());
	}

	if (0 != current_boundary_positions[0]) {
		current_boundary_positions.push_front(0);
	}

	int64_t n_boundaries = current_boundary_positions.size();

	cout << (boost::format("[AlternativeMapper.select_alternative_mapping] # blocks: %d\n") % n_boundaries).str();
//	for (int64_t block_id = 0; block_id < n_boundaries; ++block_id) {
//		cout
//				<< (boost::format("[AlternativeMapper.select_alternative_mapping] block-%d: %d\n") % block_id % current_boundary_positions[block_id]).str();
//	}

//vector<int64_t> n_entries(n_boundaries);

	/* input data format
	 *read name       = data[0];
	 * ref. id         = data[1];
	 * pair start      = data[2];
	 * pair strand     = data[3];
	 * mate ref. id    = data[4];
	 * mate start      = data[5];
	 * mate strand     = data[6];
	 * insert size     = data[7];
	 * nm              = data[8];
	 * mate nm         = data[9];
	 * read group      = data[10];
	 * aln. len        = data[11];
	 * mate aln. len   = data[12];
	 * weight          = data[13];
	 */

	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
//if(0 != block_id) {
//return;
//}
				int64_t skip_bytes = current_boundary_positions[block_id];
				int64_t next_skip_bytes = current_boundary_positions[block_id + 1];

				int64_t count_alt = 0;
				string last_ref_id;
				int64_t boundary_bytes = skip_bytes;

				map<int64_t, map<int64_t, AlternativeMap>> cluster_alt;
				string str_block_id = boost::lexical_cast<string>(block_id);
//ofstream MAPPING(options.prefix + "." + str_block_id + ".mapping.raw", ios::binary);
//stringstream buffer;
//int64_t un_flushed_size = 0;
//const int64_t BLOCK_SIZE = castle::IOUtils::MEDIUM_WRITE_BUFFER_SIZE;

				vector<string> data;
				AlternativeMap lhs;
				vector<AlternativeMap> the_alt_lists;
				string line;
				ifstream local_reader(input_file_name, ios::binary);
				local_reader.seekg(skip_bytes, ios::beg);

				while (getline(local_reader, line, '\n')) {
					boundary_bytes += line.size() + 1;
					if(boundary_bytes > next_skip_bytes) {
						break;
					}
					if(!lhs.read_from_string(line, data)) {
						continue;
					}
					the_alt_lists.push_back(lhs);
				}

				sort(the_alt_lists.begin(), the_alt_lists.end());

				auto& ref_is = options.is;

				string local_mapping_raw_file = options.prefix + "." + str_block_id + ".mapping.raw";
				if(!options.working_dir.empty()) {
					local_mapping_raw_file = options.working_prefix + "." + str_block_id + ".mapping.raw";
				}
				ofstream MAPPING(local_mapping_raw_file, ios::binary);

				int64_t the_size_entry_id = the_alt_lists.size();
				vector<bool> participated_in_already(the_alt_lists.size());
// ST-E00104:502:HFJN5CCXX:1:1104:16995:27521
				for(uint64_t current_head_id = 0; current_head_id < the_alt_lists.size(); ++current_head_id) {
					auto& the_head = the_alt_lists[current_head_id];
//cout << "head:\t" << the_head.str_cl() << "\n";
					int64_t the_last_successful_id = current_head_id;
					int64_t the_last_entry_id = current_head_id + 1;
// the status indicates if the entry belongs to a cluster head
					vector<bool> the_cluster_status;
					while(the_last_entry_id < the_size_entry_id) {
						if(participated_in_already[the_last_entry_id]) {
							the_cluster_status.push_back(false);
							++the_last_entry_id;
							continue;
						}
// lhs represents the last_successful
						auto& lhs = the_alt_lists[the_last_successful_id];
						auto& rhs = the_alt_lists[the_last_entry_id];
//cout << "cur:\t" << rhs.str_cl() << "\n";
						if (lhs.ref_id != rhs.ref_id) {
							break;
						}
// the read_start is the most important factor whether the entry will be assigned to a cluster or not.
						if (abs(rhs.read_start - lhs.read_start) > ref_is["isu"]["selected"]) {
							break;
						}
						auto& rg = lhs.read_group;
						auto& rgt = rhs.read_group;
						double devi = sqrt(pow(ref_is[rg]["sd"], 2.0) + pow(ref_is[rgt]["sd"], 2.0));
						int64_t isize_cutoff_u = ref_is[rg]["median"]
						- ref_is[rgt]["median"] + options.sd_cutoff_cl * devi;
						int64_t isize_cutoff_d = ref_is[rg]["median"]
						- ref_is[rgt]["median"] - options.sd_cutoff_cl * devi;

						int64_t map_cutoff_u = (ref_is[rgt]["isu"] > ref_is[rg]["isu"]) ? ref_is[rgt]["isu"] : ref_is[rg]["isu"];
// it's around 1000
//cout << "Isize u: " << isize_cutoff_u << ", I size d: " << isize_cutoff_d << ", map_cutoff_u: " << map_cutoff_u << "/" << ref_is["isu"]["selected"]<< "\n";

						if(rhs.mate_ref_id == lhs.mate_ref_id && rhs.read_strand == lhs.read_strand
								&& rhs.mate_strand == lhs.mate_strand
								&& abs(rhs.read_start - lhs.read_start) <= map_cutoff_u
								&& abs(rhs.mate_start - lhs.mate_start) <= map_cutoff_u) {
							if (lhs.read_strand == lhs.mate_strand) {
//cout << "type 1:\t" << rhs.str_cl() << "\n";
								the_cluster_status.push_back(true);
								the_last_successful_id = the_last_entry_id;
							} else {
								int64_t delta_isize = -1;
								if(-1 != rhs.insert_size && -1 != lhs.insert_size) {
									delta_isize = abs(rhs.insert_size - lhs.insert_size);
								}
								if ((delta_isize <= isize_cutoff_u) && (delta_isize >= isize_cutoff_d)) {
//cout << "type 2:\t" << rhs.str_cl() << "\n";
									the_cluster_status.push_back(true);
									the_last_successful_id = the_last_entry_id;
								} else {
									the_cluster_status.push_back(false);
								}
							}
						} else {
							the_cluster_status.push_back(false);
						}
						++the_last_entry_id;
					}

					if(!participated_in_already[current_head_id]) {
//cout << "H-" << count_alt << ":" << the_head.str_cl() << "\n";
						MAPPING << count_alt << "\t" << the_head.str_cl() << "\n";
					}

					int64_t the_unvisited = current_head_id + 1;
					bool any_unvisted = false;
					bool the_unvisited_changed = false;
					for(uint64_t cl_id = 0; cl_id < the_cluster_status.size(); ++cl_id) {
						int64_t current_global_cl_id = current_head_id + 1 + cl_id;
						if(the_cluster_status[cl_id]) {
// cl_id = 0 indicates the the_last_entry_id in the previous code
							auto& rhs = the_alt_lists[current_global_cl_id];
//cout << "type 3:\t" << rhs.str_cl() << "\n";
							if(!participated_in_already[current_global_cl_id]) {
//cout << "C-" << count_alt << ":" << rhs.str_cl() << "\n";
								MAPPING << count_alt << "\t" << rhs.str_cl() << "\n";
							}
							participated_in_already[current_global_cl_id] = true;
						} else {
//auto& rhs = the_alt_lists[current_global_cl_id];
//cout << "type 4:\t" << rhs.str_cl() << "\n";
// this ensures that only the smallest cl_id will be used for the next round of clustering
							if(!the_unvisited_changed) {
								the_unvisited = current_global_cl_id - 1;
								the_unvisited_changed = true;
							}
							any_unvisted = true;
						}
					}

//if(the_last_entry_id < static_cast<int64_t>(the_alt_lists.size())) {
//auto& rhs = the_alt_lists[the_last_entry_id];
//cout << "tail:\t" << rhs.str_cl() << "\n";
//}

					if(!any_unvisted) {
						the_unvisited = current_head_id + the_cluster_status.size();
					}
					if(!participated_in_already[current_head_id]) {
						++count_alt;
					}
					participated_in_already[current_head_id] = true;

					current_head_id = the_unvisited;
				}

			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	ofstream MAPPING(mapping_raw_file, ios::binary);
	int64_t acc_count = 0;
	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		string line;
		string str_block_id = boost::lexical_cast<string>(block_id);
		vector<string> data;
		string input_mapping_raw_file = options.prefix + "." + str_block_id + ".mapping.raw";
		if (!options.working_dir.empty()) {
			input_mapping_raw_file = options.working_prefix + "." + str_block_id + ".mapping.raw";
		}
		ifstream in(input_mapping_raw_file, ios::binary);
		int64_t the_last_cluster_id = 0;
		while (getline(in, line, '\n')) {
			int64_t first_tab = line.find('\t');
			string cluster_id = line.substr(0, first_tab);
			string line_remain = line.substr(first_tab);
			the_last_cluster_id = boost::lexical_cast<int64_t>(cluster_id);
			string n_cluster_id = boost::lexical_cast<string>(the_last_cluster_id + acc_count);
			MAPPING << n_cluster_id << line_remain << "\n";
		}
		acc_count += the_last_cluster_id + 1;
	}
	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string a_file_name = options.prefix + "." + str_block_id + ".mapping.raw";
			if(!options.working_dir.empty()) {
				a_file_name = options.working_prefix + "." + str_block_id + ".mapping.raw";
			}
			boost::filesystem::remove(a_file_name);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void AlternativeMapper::select_alternative_mapping_alt() {
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.select_alternative_mapping_alt");
	checker.start();
	string input_file_name = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		input_file_name = options.working_prefix + ".alt.map.sort";
	}
	string mapping_raw_file = options.prefix + ".mapping.raw";
	if (!options.working_dir.empty()) {
		mapping_raw_file = options.working_prefix + ".mapping.raw";
	}
//	if(boost::filesystem::exists(mapping_raw_file)) {
//		return;
//	}
//cout << input_file_name << "\n";
// calculate the positions of boundary entries.
	vector<function<void()> > tasks;
// prepare for the boundary elements

	deque<int64_t> current_boundary_positions;
	int64_t size_total_ref = castle::IOUtils::get_file_size(input_file_name);
	if(size_total_ref >= 4096000) {
// all the small cluster with ref id within each 4096000-byte block won't be separated.
// For instance, unlike the large chromosome 1, small chromosomes such as GL000191.1, GL000192.1, ... and hs37d5 will be included in a single large block.
		int64_t n_remaining_bases = size_total_ref;
		const int64_t size_block = 4096000;
		int64_t boundary_id = 0;
		current_boundary_positions.push_back(0);
		while (n_remaining_bases >= 0) {
			int64_t last_base_pos = current_boundary_positions[boundary_id];
			n_remaining_bases = size_total_ref - last_base_pos;
			if (n_remaining_bases >= size_block) {
				current_boundary_positions.push_back(last_base_pos + size_block);
				++boundary_id;
			} else {
				current_boundary_positions.push_back(size_total_ref);
				break;
			}
		}

		cout << (boost::format("[AlternativeMapper.select_alternative_mapping_alt] collect boundary positions of file size %d\n") % size_total_ref).str();
		int64_t n_boundaries = current_boundary_positions.size();
		for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				ifstream local_reader(input_file_name, ios::binary);
				int64_t skip_bytes = current_boundary_positions[block_id];
				int64_t next_skip_bytes = current_boundary_positions[block_id + 1];
				local_reader.seekg(skip_bytes, ios::beg);
				string line;
				int64_t boundary_bytes = skip_bytes;
				int64_t previous_boundary_bytes = skip_bytes;
				const char* delims = "\t";
				vector<string> data;
				string prev_ref_id;
				bool found_skip_point = false;
				getline(local_reader, line, '\n');
				boundary_bytes += line.size() + 1;
				while(getline(local_reader, line, '\n')) {
					previous_boundary_bytes = boundary_bytes;
					boundary_bytes += line.size() + 1;
					if(boundary_bytes > next_skip_bytes) {
						break;
					}
					castle::StringUtils::c_string_multi_split(line, delims, data);
					if(!prev_ref_id.empty() && data[1] != prev_ref_id) {
						found_skip_point = true;
						break;
					}
					prev_ref_id = data[1];
				}
				if(found_skip_point) {
					current_boundary_positions[block_id] = previous_boundary_bytes;
				} else {
					current_boundary_positions[block_id] = -1;
				}
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		current_boundary_positions.erase(remove_if(current_boundary_positions.begin(), current_boundary_positions.end(), [](int64_t value) {return -1 == value;}), current_boundary_positions.end());
	} else {
		current_boundary_positions.push_back(0);
		current_boundary_positions.push_back(size_total_ref);
	}

	if (0 != current_boundary_positions[0]) {
		current_boundary_positions.push_front(0);
	}

	int64_t n_boundaries = current_boundary_positions.size();

	cout << (boost::format("[AlternativeMapper.select_alternative_mapping_alt] # blocks: %d\n") % n_boundaries).str();
//	for (int64_t block_id = 0; block_id < n_boundaries; ++block_id) {
//		cout
//				<< (boost::format("[AlternativeMapper.select_alternative_mapping] block-%d: %d\n") % block_id % current_boundary_positions[block_id]).str();
//	}

//vector<int64_t> n_entries(n_boundaries);

	/* input data format
	 *read name       = data[0];
	 * ref. id         = data[1];
	 * pair start      = data[2];
	 * pair strand     = data[3];
	 * mate ref. id    = data[4];
	 * mate start      = data[5];
	 * mate strand     = data[6];
	 * insert size     = data[7];
	 * nm              = data[8];
	 * mate nm         = data[9];
	 * read group      = data[10];
	 * aln. len        = data[11];
	 * mate aln. len   = data[12];
	 * weight          = data[13];
	 */

	vector<int64_t> max_count_alt_lists(n_boundaries - 1);
	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			int64_t cur_pos = current_boundary_positions[block_id];
			int64_t next_pos = current_boundary_positions[block_id + 1];

			string last_ref_id;

//			unordered_map<int64_t, unordered_map<int64_t, vector<string>>> cluster_alt;
			boost::unordered_map<int64_t, boost::unordered_map<int64_t, AlternativeMap>> cluster_alt;
			string str_block_id = boost::lexical_cast<string>(block_id);

			const char* delims = "\t";
			int64_t count_alt = 0;
			int64_t max_count_alt = 0;
			string lastseqid;
			vector<string> data;
			AlternativeMap lhs;
//			vector<AlternativeMap> the_alt_lists;
			string line;

			ifstream local_reader(input_file_name, ios::binary);
			local_reader.seekg(cur_pos, ios::beg);

			const int64_t v_max = numeric_limits<int64_t>::max();
			string out_file_name = mapping_raw_file + "." + str_block_id;
			ofstream MAPPING(out_file_name, ios::binary);
			while(getline(local_reader, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delims, data);
				if (data.size() < 13) {
					continue;
				}
				try {
				string& seqid = data[1];
//				const bool debug = ("ST-E00104:502:HFJN5CCXX:1:1203:24119:10117" == seqid);
				lhs.read_name = data[0];
				lhs.ref_id = data[1];
				lhs.read_start = boost::lexical_cast<int64_t>(data[2]);
				lhs.read_strand = boost::lexical_cast<int32_t>(data[3]);
				lhs.mate_ref_id = data[4];
				lhs.mate_start = boost::lexical_cast<int64_t>(data[5]);
				lhs.mate_strand = boost::lexical_cast<int32_t>(data[6]);
				lhs.insert_size = v_max;
				lhs.nm = -1;
				lhs.mate_nm = -1;
				lhs.aln_len = -1;
				lhs.mate_aln_len = -1;
				lhs.weight = -1;
//				if (13 == data.size()) {
//					lhs.nm = boost::lexical_cast<int64_t>(data[7]);
//					lhs.mate_nm = boost::lexical_cast<int64_t>(data[8]);
//					lhs.read_group = data[9];
//					lhs.aln_len = boost::lexical_cast<int64_t>(data[10]);
//					lhs.mate_aln_len = boost::lexical_cast<int64_t>(data[11]);
//					lhs.weight = boost::lexical_cast<int64_t>(data[12]);
//				} else if (data.size() > 13) {
				if(data[7].empty()) {
					lhs.insert_size = 0;
				} else {
					lhs.insert_size = boost::lexical_cast<int64_t>(data[7]);
				}
					lhs.nm = boost::lexical_cast<int64_t>(data[8]);
					lhs.mate_nm = boost::lexical_cast<int64_t>(data[9]);
					lhs.read_group = data[10];
					lhs.aln_len = boost::lexical_cast<int64_t>(data[11]);
					lhs.mate_aln_len = boost::lexical_cast<int64_t>(data[12]);
					lhs.weight = boost::lexical_cast<double>(data[13]);
//				}

				bool incluster = false;
				int64_t topush = -1;
				bool is_terminating_cl5 = false;
				bool is_continuing_cl4 = false;
				auto& ref_is = options.is;
				if (seqid != lastseqid) {
					cluster_alt.clear();
				}
				lastseqid = seqid;
				for (int64_t i = count_alt; i >= 0; --i) {
					bool outofrange = false;
					if (cluster_alt.end() == cluster_alt.find(i)) {
						continue;
					}
					for (auto& datat : cluster_alt[i]) {
						auto& d = datat.second;
						const string& readnamet = d.read_name;
						const string& seqidt = d.ref_id;
						int64_t startt = d.read_start;
						int32_t strandt = d.read_strand;
						string& mseqidt = d.mate_ref_id;
						int64_t mstartt = d.mate_start;
						int32_t mstrandt = d.mate_strand;
						int64_t isizet = d.insert_size;
						string rgt = d.read_group;

						double devi = sqrt(pow(ref_is[lhs.read_group]["sd"], 2.0) + pow(ref_is[rgt]["sd"], 2.0));
						int64_t isize_cutoff_u = ref_is[lhs.read_group]["median"] - ref_is[rgt]["median"] + options.sd_cutoff_cl * devi;
						int64_t isize_cutoff_d = ref_is[lhs.read_group]["median"] - ref_is[rgt]["median"] - options.sd_cutoff_cl * devi;

						int64_t map_cutoff_u = (ref_is[rgt]["isu"] > ref_is[lhs.read_group]["isu"]) ? ref_is[rgt]["isu"] : ref_is[lhs.read_group]["isu"];
						if (seqid == seqidt) {
							if (readnamet == lhs.read_name && lhs.mate_ref_id == mseqidt &&
									startt == lhs.read_start && mstartt == lhs.mate_start && lhs.read_strand == strandt && lhs.mate_strand == mstrandt) {
								is_continuing_cl4 = true;
								break;
							} else if (lhs.mate_ref_id == mseqidt && lhs.read_strand == strandt && lhs.mate_strand == mstrandt && abs(lhs.read_start - startt)
							<= map_cutoff_u && abs(lhs.mate_start - mstartt) <= map_cutoff_u) {
								if (lhs.read_strand == lhs.mate_strand) {
									topush = i;
									incluster = true;
									is_terminating_cl5 = true;
									break;
								} else {
									int64_t delta_isize = v_max;
//									if (v_max != lhs.insert_size && v_max != isizet) {
										delta_isize = lhs.insert_size - isizet;
//									}
									if ((delta_isize <= isize_cutoff_u) && (delta_isize >= isize_cutoff_d)) {
										topush = i;
										incluster = true;
										is_terminating_cl5 = true;
										break;
									}
								}
							}
							if (abs(lhs.read_start - startt) > ref_is["isu"]["selected"]) {
								outofrange = true;
								break;
							}
						} else {
							is_terminating_cl5 = true;
							break;
						}
					}
					if (outofrange || is_terminating_cl5 || is_continuing_cl4) {
						break;
					}
				}
				if (is_continuing_cl4) {
					continue;
				}
				if (!incluster) {
//					cluster_alt[count_alt][0].clear();
//					cluster_alt[count_alt][0].insert(cluster_alt[count_alt][0].end(), data.begin(), data.end());
					cluster_alt[count_alt][0] = lhs;
					if(seqid != lhs.mate_ref_id) {
						lhs.insert_size = 0;
					}
					string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % count_alt % lhs.read_name % lhs.read_group % seqid % lhs.read_strand %
							lhs.read_start % lhs.aln_len % lhs.mate_ref_id % lhs.mate_strand % lhs.mate_start % lhs.mate_aln_len % (v_max == lhs.insert_size ? "" : boost::lexical_cast<string>(lhs.insert_size)) % lhs.nm % lhs.mate_nm % lhs.weight).str();
					max_count_alt = max(max_count_alt, count_alt);
					MAPPING << an_entry;
					++count_alt;
				}
				if (-1 != topush) {
					int64_t k = cluster_alt[topush].size();
//					cluster_alt[topush][k].insert(cluster_alt[topush][k].end(), data.begin(), data.end());
					cluster_alt[topush][k] = lhs;
					max_count_alt = max(max_count_alt, topush);
					if(seqid != lhs.mate_ref_id) {
						lhs.insert_size = 0;
					}
					string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % topush % lhs.read_name % lhs.read_group % seqid % lhs.read_strand %
							lhs.read_start % lhs.aln_len % lhs.mate_ref_id % lhs.mate_strand % lhs.mate_start % lhs.mate_aln_len % (v_max == lhs.insert_size ? "" : boost::lexical_cast<string>(lhs.insert_size)) % lhs.nm % lhs.mate_nm % lhs.weight).str();
					MAPPING << an_entry;
				}

				if(cur_pos >= next_pos) {
					break;
				}
				} catch(exception& ex) {
					cout << line << "\n";
					cout << data.size() << "\n";
					cout << ex.what() << "\n";
				}
			}
			max_count_alt_lists[block_id] = max_count_alt;
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	ofstream MAPPING(mapping_raw_file, ios::binary);
	int64_t acc_count = 0;
	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		string line;
		string str_block_id = boost::lexical_cast<string>(block_id);
		vector<string> data;
		string input_mapping_raw_file = mapping_raw_file + "." + str_block_id;
		ifstream in(input_mapping_raw_file, ios::binary);
		int64_t the_last_cluster_id = 0;
		while (getline(in, line, '\n')) {
			int64_t first_tab = line.find('\t');
			string cluster_id = line.substr(0, first_tab);
			string line_remain = line.substr(first_tab);
			the_last_cluster_id = boost::lexical_cast<int64_t>(cluster_id);
			// the last + 1 makes the entry id starts from 1.
			string n_cluster_id = boost::lexical_cast<string>(the_last_cluster_id + acc_count + 1);
			MAPPING << n_cluster_id << line_remain << "\n";
		}
		acc_count += max_count_alt_lists[block_id] + 1;
	}
	for (int64_t block_id = 0; block_id < n_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string input_mapping_raw_file = mapping_raw_file + "." + str_block_id;
			boost::filesystem::remove(input_mapping_raw_file);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void AlternativeMapper::select_alternative_mapping_serial() {
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.select_alternative_mapping_serial");
	checker.start();
	string input_file_name = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		input_file_name = options.working_prefix + ".alt.map.sort";
	}
	ifstream local_reader(input_file_name, ios::binary);
	string line;
	vector<string> data;
//std::string delims("\t");
	const char* delims = "\t";
	/* input data format
	 *read name       = data[0];
	 * ref. id         = data[1];
	 * pair start      = data[2];
	 * pair strand     = data[3];
	 * mate ref. id    = data[4];
	 * mate start      = data[5];
	 * mate strand     = data[6];
	 * insert size     = data[7];
	 * nm              = data[8];
	 * mate nm         = data[9];
	 * read group      = data[10];
	 * aln. len        = data[11];
	 * mate aln. len   = data[12];
	 * weight          = data[13];
	 */

	int64_t count_alt = 0;
	string lastseqid;

//	unordered_map<int64_t, unordered_map<int64_t, vector<string>>> cluster_alt;
	boost::unordered_map<int64_t, boost::unordered_map<int64_t, AlternativeMap>> cluster_alt;
	string mapping_raw_file = options.prefix + ".mapping.raw";
	if (!options.working_dir.empty()) {
		mapping_raw_file = options.working_prefix + ".mapping.raw";
	}
	ofstream MAPPING(mapping_raw_file, ios::binary);
	AlternativeMap lhs;
	stringstream buffer;
	const int64_t v_max = numeric_limits<int64_t>::max();
//int64_t un_flushed_size = 0;
//const int64_t BLOCK_SIZE = castle::IOUtils::MEDIUM_WRITE_BUFFER_SIZE;
	while (getline(local_reader, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delims, data);
		if (data.size() < 13) {
			continue;
		}
		string& seqid = data[1];
		lhs.read_name = data[0];
		lhs.ref_id = data[1];
		lhs.read_start = boost::lexical_cast<int64_t>(data[2]);
		lhs.read_strand = boost::lexical_cast<int32_t>(data[3]);
		lhs.mate_ref_id = data[4];
		lhs.mate_start = boost::lexical_cast<int64_t>(data[5]);
		lhs.mate_strand = boost::lexical_cast<int32_t>(data[6]);
		lhs.insert_size = v_max;
		lhs.nm = -1;
		lhs.mate_nm = -1;
		lhs.aln_len = -1;
		lhs.mate_aln_len = -1;
		lhs.weight = -1;

		const bool debug = ("ST-E00104:502:HFJN5CCXX:1:1203:24119:10117" == lhs.read_name);
		if (13 == data.size()) {
			lhs.nm = boost::lexical_cast<int64_t>(data[7]);
			lhs.mate_nm = boost::lexical_cast<int64_t>(data[8]);
			lhs.read_group = data[9];
			lhs.aln_len = boost::lexical_cast<int64_t>(data[10]);
			lhs.mate_aln_len = boost::lexical_cast<int64_t>(data[11]);
			lhs.weight = boost::lexical_cast<int64_t>(data[12]);
		} else if (data.size() > 13) {
			lhs.insert_size = boost::lexical_cast<int64_t>(data[7]);
			lhs.nm = boost::lexical_cast<int64_t>(data[8]);
			lhs.mate_nm = boost::lexical_cast<int64_t>(data[9]);
			lhs.read_group = data[10];
			lhs.aln_len = boost::lexical_cast<int64_t>(data[11]);
			lhs.mate_aln_len = boost::lexical_cast<int64_t>(data[12]);
			lhs.weight = boost::lexical_cast<int64_t>(data[13]);
		}

		bool incluster = false;
		int64_t topush = -1;
		bool is_terminating_cl5 = false;
		bool is_continuing_cl4 = false;
		auto& ref_is = options.is;
		if (seqid != lastseqid) {
			cluster_alt.clear();
		}
		lastseqid = seqid;
		if(debug) {
			cout << "[AlternativeMapper.select_alternative_mapping_serial] count_alt before: " << count_alt << "\n";
		}
		for (int64_t i = count_alt; i >= 0; --i) {
			bool outofrange = false;
			if (cluster_alt.end() == cluster_alt.find(i)) {
				continue;
			}
			for (auto& datat : cluster_alt[i]) {
				auto& d = datat.second;
				const string& readnamet = d.read_name;
				const string& seqidt = d.ref_id;
				int64_t startt = d.read_start;
				int32_t strandt = d.read_strand;
				string& mseqidt = d.mate_ref_id;
				int64_t mstartt = d.mate_start;
				int32_t mstrandt = d.mate_strand;
				int64_t isizet = d.insert_size;
				string rgt = d.read_group;

				double devi = sqrt(pow(ref_is[lhs.read_group]["sd"], 2.0) + pow(ref_is[rgt]["sd"], 2.0));
				int64_t isize_cutoff_u = ref_is[lhs.read_group]["median"] - ref_is[rgt]["median"] + options.sd_cutoff_cl * devi;
				int64_t isize_cutoff_d = ref_is[lhs.read_group]["median"] - ref_is[rgt]["median"] - options.sd_cutoff_cl * devi;

				int64_t map_cutoff_u = (ref_is[rgt]["isu"] > ref_is[lhs.read_group]["isu"]) ? ref_is[rgt]["isu"] : ref_is[lhs.read_group]["isu"];
				if (seqid == seqidt) {
					if (readnamet == lhs.read_name && lhs.mate_ref_id == mseqidt &&
							startt == lhs.read_start && mstartt == lhs.mate_start && lhs.read_strand == strandt && lhs.mate_strand == mstrandt) {
						if(debug) {
							cout << "[AlternativeMapper.select_alternative_mapping_serial] is_continuing_cl4: " << i << "\n";
						}
						is_continuing_cl4 = true;
						break;
					} else if (lhs.mate_ref_id == mseqidt && lhs.read_strand == strandt && lhs.mate_strand == mstrandt && abs(lhs.read_start - startt)
					<= map_cutoff_u && abs(lhs.mate_start - mstartt) <= map_cutoff_u) {
						if(debug) {
							cout << "[AlternativeMapper.select_alternative_mapping_serial] here-1: " << i << "\n";
						}
						if (lhs.read_strand == lhs.mate_strand) {
							if(debug) {
								cout << "[AlternativeMapper.select_alternative_mapping_serial] here-2: " << i << "\n";
							}
							topush = i;
							incluster = true;
							is_terminating_cl5 = true;
							break;
						} else {
							if(debug) {
								cout << "[AlternativeMapper.select_alternative_mapping_serial] here-3: " << i << ", " << lhs.insert_size << "/" << isizet << "\n";
							}
							int64_t delta_isize = v_max;
//							if (v_max != lhs.insert_size && v_max != isizet) {
								delta_isize = lhs.insert_size - isizet;
								if(debug) {
									cout << "[AlternativeMapper.select_alternative_mapping_serial] isize: " << delta_isize << "\n";
								}
//							}
							if ((delta_isize <= isize_cutoff_u) && (delta_isize >= isize_cutoff_d)) {
								if(debug) {
									cout << "[AlternativeMapper.select_alternative_mapping_serial] here-4: " << i << "\n";
								}
								topush = i;
								incluster = true;
								is_terminating_cl5 = true;
								break;
							}
						}
					}
					if (abs(lhs.read_start - startt) > ref_is["isu"]["selected"]) {
						if(debug) {
							cout << "[AlternativeMapper.select_alternative_mapping_serial] here-5: " << i << "\n";
						}
						outofrange = true;
						break;
					}
				} else {
					if(debug) {
						cout << "[AlternativeMapper.select_alternative_mapping_serial] here-6: " << i << "\n";
					}
					is_terminating_cl5 = true;
					break;
				}
			}
			if (outofrange || is_terminating_cl5 || is_continuing_cl4) {
				break;
			}
		}
		if (is_continuing_cl4) {
			continue;
		}
		if (!incluster) {
			if(debug) {
				cout << "[AlternativeMapper.select_alternative_mapping_serial] count_alt : " << count_alt << "\n";
			}
	//					cluster_alt[count_alt][0].clear();
	//					cluster_alt[count_alt][0].insert(cluster_alt[count_alt][0].end(), data.begin(), data.end());
			cluster_alt[count_alt][0] = lhs;
			string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % (count_alt + 1) % lhs.read_name % lhs.read_group % seqid % lhs.read_strand %
					lhs.read_start % lhs.aln_len % lhs.mate_ref_id % lhs.mate_strand % lhs.mate_start % lhs.mate_aln_len % (v_max == lhs.insert_size ? "" : boost::lexical_cast<string>(lhs.insert_size)) % lhs.nm % lhs.mate_nm % lhs.weight).str();
			MAPPING << an_entry;
			++count_alt;
		}
		if (-1 != topush) {
			if(debug) {
				cout << "[AlternativeMapper.select_alternative_mapping_serial] topush : " << topush << "\n";
			}
			int64_t k = cluster_alt[topush].size();
	//					cluster_alt[topush][k].insert(cluster_alt[topush][k].end(), data.begin(), data.end());
			cluster_alt[topush][k] = lhs;
			string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % (topush + 1) % lhs.read_name % lhs.read_group % seqid % lhs.read_strand %
					lhs.read_start % lhs.aln_len % lhs.mate_ref_id % lhs.mate_strand % lhs.mate_start % lhs.mate_aln_len % (v_max == lhs.insert_size ? "" : boost::lexical_cast<string>(lhs.insert_size)) % lhs.nm % lhs.mate_nm % lhs.weight).str();
			MAPPING << an_entry;
		}
	}
	cout << checker;
}

// The current BAM tools API is not able to create index for BAM file sorted by read name
void AlternativeMapper::create_raw_alternative_mapping_serial() {
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.create_raw_alternative_mapping_serial");
	checker.start();
	string a_path(options.prefix + ".disc.sorted.bam");
	if (!options.working_dir.empty()) {
		a_path = options.working_prefix + ".disc.sorted.bam";
	}
//string a_path("/dev/shm/example.bam");

//	cout << a_path << "\n";
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path)) {
		cout << checker;
		return;
	}
//	cout << "here-0\n";
	BamTools::BamAlignment al;
	vector<BamTools::BamAlignment> recent_als;
	const RefVector& ref_vec = local_reader.GetReferenceData();
	boost::unordered_map<string, int64_t> ref_reverse_index;
	for (uint64_t ref_id = 0; ref_id < ref_vec.size(); ++ref_id) {
		auto& a_ref = ref_vec[ref_id];
		ref_reverse_index[a_ref.RefName] = ref_id;
	}
	auto& ref_blacklist = options.point_black_lists;
//string a_tag;
//uint32_t nm;
//stringstream out;
////ofstream out("/dev/shm/test.alt.map.raw.intest", ios::binary);
	string out_path(options.prefix + ".alt.map.raw");
	if (!options.working_dir.empty()) {
		out_path = options.working_prefix + ".alt.map.raw";
	}
	int64_t n_entries = 0;
	const bool debug = false;
	ofstream out(out_path, ios::binary);

	while (local_reader.LoadNextAlignmentCore(al)) {
		if (debug) {
			cout << "here-1\n";
			cout << al.Name << "\n";
		}
		++n_entries;
		if (0 == (n_entries & 1048575)) {
			cout << (boost::format("[AlternativeMapper.create_raw_alternative_mapping_serial] Processed: %d\n") % n_entries).str();
		}
		if (recent_als.empty()) {
			recent_als.push_back(al);
			continue;
		}
		if (recent_als.back().Name == al.Name) {
			recent_als.push_back(al);
		} else {

			if (debug) {
				cout << "in processing\n";
			}
			if (3 <= recent_als.size()) {
				if (debug) {
					cout << "# recent aln.: " << recent_als.size() << "\n";
					for (uint64_t e_id = 0; e_id < recent_als.size(); ++e_id) {
						cout << recent_als[e_id].Name << "\n";
						for (auto c : recent_als[e_id].CigarData) {
							cout << c.Length << c.Type;
						}
						cout << "\n";
					}
					cout << "\n";
				}
				select_most_likely_alignments(recent_als);
				if (debug) {
					cout << "# remains: " << recent_als.size() << "\n";
				}
				if (2 == recent_als.size()) {
					recent_als[0].Name += "sc";
					recent_als[1].Name += "sc";
				}
				if (debug) {
					cout << "here-2\n";
					for (auto aln : recent_als) {
						cout << aln.Name << "\n";
					}
				}
//for(auto aln : recent_als) {
//aln.Name += "sc";
//if(debug) {
//cout << aln.Name << "\n";
//}
//}
//}
			}
			if (debug) {
				cout << "here-3\n";
			}
// There are now only two entries
			if (2 == recent_als.size()) {
				bool found_valid_entry = false;
				if (debug) {
					cout << "here-4\n";
				}
				if ((recent_als[0].RefID < recent_als[1].RefID) || (recent_als[0].RefID == recent_als[1].RefID && recent_als[0].Position > recent_als[1].Position)) {
					if (debug) {
						cout << "here-5\n";
					}
					iter_swap(recent_als.begin(), recent_als.begin() + 1);
				}
				if (debug) {
					cout << "here-6\n";
				}
				auto& data1 = recent_als[0];
				auto& data2 = recent_als[1];
				if (-1 == data1.RefID || -1 == data2.RefID) {
					recent_als.clear();
					recent_als.push_back(al);
					continue;
				}
				if (options.min_mapq) {
					if (debug) {
						cout << "here-7\n";
					}
					if ((data1.MapQuality < options.min_mapq) || data2.MapQuality < options.min_mapq) {
						if (debug) {
							cout << "here-8\n";
						}
						recent_als.clear();
						recent_als.push_back(al);
						continue;
					}
				}
// the XT tag won't be filled by bwa mem, but I just adopted codes from the original meerkat Perl script.

				char xt = 'U';
				char mate_xt = 'U';
				string xt_tag_str;
				string mate_xt_tag_str;
				if (debug) {
					cout << "here-9\n";
				}
				if (data1.GetTag("XT", xt_tag_str)) {
					if (debug) {
						cout << "here-10\n";
					}
					xt = xt_tag_str[0];
				}
				if (debug) {
					cout << "here-11\n";
				}
				if (data2.GetTag("XT", mate_xt_tag_str)) {
					if (debug) {
						cout << "here-12\n";
					}
					mate_xt = mate_xt_tag_str[0];
				}
				if (debug) {
					cout << "here-13\n";
				}
				if (!data1.AlignmentFlag) {
					if (debug) {
						cout << "here-14\n";
					}
					xt = 'R';
				}
				if (debug) {
					cout << "here-15\n";
				}
				if (!data2.AlignmentFlag) {
					if (debug) {
						cout << "here-16\n";
					}
					mate_xt = 'R';
				}
				if (debug) {
					cout << "here-17\n";
				}
				if (options.ad_align) {
					if (('R' == xt && 'R' == mate_xt) || 'M' == xt || 'M' == mate_xt) {
						if (debug) {
							cout << "here-18\n";
						}
						recent_als.clear();
						recent_als.push_back(al);
						continue;
					} else {
						if (!options.use_all_align) {
							if (!('U' == xt && 'U' == mate_xt)) {
								if (debug) {
									cout << "here-19\n";
								}
								recent_als.clear();
								recent_als.push_back(al);
								continue;
							}
						}
					}
				}
				if (debug) {
					cout << "here-20\n";
					cout << data1.RefID << "/" << data2.RefID << "/" << ref_vec.size() << "\n";
				}
// check coverage based black list
				if (-1 == data1.RefID || -1 == data2.RefID) {
					recent_als.clear();
					recent_als.push_back(al);
					continue;
				}
				auto a_key_1 = make_pair(ref_vec[data1.RefID].RefName, data1.Position);
				auto a_key_2 = make_pair(ref_vec[data2.RefID].RefName, data2.Position);
				auto a_key_3 = make_pair(ref_vec[data1.RefID].RefName, static_cast<int64_t>(data1.Position + (options.cut_sr << 1)));
				auto a_key_4 = make_pair(ref_vec[data2.RefID].RefName, static_cast<int64_t>(data2.Position + (options.cut_sr << 1)));
				auto the_black_itr_1 = ref_blacklist.find(a_key_1);
				auto the_black_itr_2 = ref_blacklist.find(a_key_2);
				auto the_black_itr_3 = ref_blacklist.find(a_key_3);
				auto the_black_itr_4 = ref_blacklist.find(a_key_4);
				if (ref_blacklist.end() != the_black_itr_1 || ref_blacklist.end() != the_black_itr_2 || ref_blacklist.end() != the_black_itr_3 || ref_blacklist.end() != the_black_itr_4) {
					if (debug) {
						cout << "here-21\n";
					}
					recent_als.clear();
					recent_als.push_back(al);
					continue;
				}
				if (!has_valid_insertion_size(out, ref_vec, ref_reverse_index, recent_als, found_valid_entry)) {
					if (debug) {
						cout << "here-22\n";
					}
					recent_als.clear();
					recent_als.push_back(al);
					continue;
				}
			}
			recent_als.clear();
			recent_als.push_back(al);
		}
	}

	cout << checker;
}
void AlternativeMapper::create_raw_alternative_mapping() {
//	if (boost::filesystem::exists(options.output_filename)) {
//		return;
//	}
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.create_raw_alternative_mapping");
	checker.start();
	cout << (boost::format("[AlternativeMapper.create_raw_alternative_mapping] input: %s\n") % options.input_filename).str();
	string a_path(options.input_filename);
	BamReader serial_reader;
	string a_path_index(a_path + ".bfi");
	if (!serial_reader.Open(a_path)) {
		cout << checker;
		return;
	}
//	const string debug_read_name = "ST-E00104:502:HFJN5CCXX:7:1123:12266:39510";
	vector<function<void()> > tasks;
	vector<int64_t> block_boundary;
	if (boost::filesystem::exists(a_path_index)) {
		string line;
		ifstream in(a_path_index);
		while (getline(in, line, '\n')) {
			block_boundary.push_back(boost::lexical_cast<int64_t>(line));
		}
	} else {
		const int64_t N_BLOCK_ENTRIES = 262144;
		auto& data = serial_reader.GetBGZF();
		BamAlignment al;
		int64_t n_entries = 0;
		string previous_name;
		block_boundary.push_back(0);
		int64_t global_n_entries = 0;
		cout << "[AlternativeMapper.create_raw_alternative_mapping] start indexing the BAM file\n";
		int64_t previous_pos = -1;
		vector<char> tmp_buf;
		while (serial_reader.ReadNextAlignment(tmp_buf)) {
			++global_n_entries;
			++n_entries;
			if (n_entries > N_BLOCK_ENTRIES) {
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
		block_boundary.push_back(numeric_limits<int64_t>::max());
		ofstream out(a_path_index, ios::binary);
		for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
			out << block_boundary[block_id] << "\n";
		}
		cout << (boost::format("[AlternativeMapper.create_raw_alternative_mapping] %d entries were indexed\n") % global_n_entries).str();
	}
	int64_t n_block_boundaries = block_boundary.size();
	vector<string> output_files(n_block_boundaries - 1);

//	boost::unordered_set<string> the_keys;
//	{
//		string line;
//		ifstream in("/home/el174/vincent_lim/preserved_exp/alt_map_debug.txt", ios::binary);
//		while(getline(in, line, '\n')) {
//			the_keys.insert(line);
//		}
//	}
	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path)) {
				return;
			}
			const RefVector& ref_vec = local_reader.GetReferenceData();
			boost::unordered_map<string, int64_t> ref_reverse_index;
			for (uint64_t ref_id = 0; ref_id < ref_vec.size(); ++ref_id) {
				auto& a_ref = ref_vec[ref_id];
				ref_reverse_index[a_ref.RefName] = ref_id;
			}
//			auto& ref_blacklist = options.point_black_lists;
				auto& local_data = local_reader.GetBGZF();
				if(0 != block_boundary[block_id]) {
					local_data.Seek(block_boundary[block_id]);
				}
				int64_t next_block_pos = block_boundary[block_id + 1];
				int64_t cur_pos = -1;
				string str_block_id = boost::lexical_cast<string>(block_id);

				vector<BamAlignment> recent_als;
				BamAlignment al;

				string out_path(options.output_filename + "." + str_block_id);
				output_files[block_id] = out_path;
				ofstream out(out_path, ios::binary);
				bool is_terminating = false;
				while (local_reader.LoadNextAlignmentCore(al)) {
					if (recent_als.empty()) {
						recent_als.push_back(al);
						continue;
					}
					if (recent_als.back().Name == al.Name) {
						if(-1 != al.RefID) {
							if (options.min_mapq) {
								if (al.MapQuality >= options.min_mapq) {
									recent_als.push_back(al);
								}
							} else {
								recent_als.push_back(al);
							}
						}
					} else {
//					const bool debug = (recent_als.size() > 0) && (string::npos != recent_als[0].Name.find(debug_read_name));
//					if(debug) {
//						cout << "[AlternativeMapper.create_raw_alternative_mapping] alignments: " << recent_als.size() << "\n";
//					}

//						if(string::npos != recent_als[0].Name.find("ST-E00104:502:HFJN5CCXX:6:2215:14671:61960")) {
//							for(auto& lal: recent_als) {
//								cout << BamWriter::GetSAMAlignment(lal, local_reader.GetReferenceData()) << "\n";
//							}
//						}
						vector<BamAlignment> original_als(recent_als);
						if (3 <= recent_als.size()) {
//						if(debug) {
//							cout << "[AlternativeMapper.create_raw_alternative_mapping] select most-likely alignments: " << recent_als.size() << "\n";
//							for(auto& local_al : recent_als) {
//								cout << "[AlternativeMapper.create_raw_alternative_mapping] the al: " << BamWriter::GetSAMAlignment(local_al, ref_vec) << "\n";
//							}
//						}
							select_most_likely_alignments(recent_als);
							recent_als[0].Name += "sc";
							recent_als[1].Name += "sc";
//						if(debug) {
//							cout << "[AlternativeMapper.create_raw_alternative_mapping] select most-likely alignments: " << recent_als.size() << "\n";
//						}
						}
//					if(debug) {
//						cout << "[AlternativeMapper.create_raw_alternative_mapping] start selecting maximum-likely alignments? " << recent_als.size() << "\n";
//					}
//					bool found_valid_entry = false;
//						write_alt_map(out, ref_vec, ref_reverse_index, recent_als);
						if(!write_alt_map(out, ref_vec, ref_reverse_index, recent_als)) {
//						if(debug) {
//							cout << "[AlternativeMapper.create_raw_alternative_mapping] select maximum-likely alignments: " << recent_als.size() << "\n";
//						}
							vector<BamAlignment> copied_als(original_als);
							select_maximum_likely_alignments(copied_als);
							if (3 <= copied_als.size()) {
								copied_als[0].Name += "sc";
								copied_als[1].Name += "sc";
							}
//							write_alt_map(out, ref_vec, ref_reverse_index, copied_als);
							if(!write_alt_map(out, ref_vec, ref_reverse_index, copied_als)) {
								boost::unordered_map<int8_t, vector<BamAlignment>> paired_aln_map;
								for(auto& al: original_als) {
									if(al.IsFirstMate()) {
										paired_aln_map[1].push_back(al);
									} else if(al.IsSecondMate()) {
										paired_aln_map[2].push_back(al);
									}
								}
//							if (debug) {
//								cout << "[AlternativeMapper.create_raw_alternative_mapping] here-0: " << paired_aln_map[1].size() << "/" << paired_aln_map[2].size() << "\n";
//							}
								bool found_pair = false;
								for(auto& first_pair : paired_aln_map[1]) {
									for(auto& second_pair : paired_aln_map[2]) {
										vector<BamAlignment> temp_als;
										temp_als.push_back(first_pair);
										temp_als.push_back(second_pair);
										if (3 <= copied_als.size()) {
											temp_als[0].Name += "sc";
											temp_als[1].Name += "sc";
										}
//									if (debug) {
//										cout << "[AlternativeMapper.create_raw_alternative_mapping] here-0-a: " << BamWriter::GetSAMAlignment(temp_als[0], ref_vec) << "\n";
//										cout << "[AlternativeMapper.create_raw_alternative_mapping] here-0-b: " << BamWriter::GetSAMAlignment(temp_als[1], ref_vec) << "\n";
//									}
										if(write_alt_map(out, ref_vec, ref_reverse_index, temp_als)) {
											found_pair = true;
											break;
										}
									}
									if(found_pair) {
										break;
									}
								}
							}
						}
						recent_als.clear();
						recent_als.push_back(al);
					}
					cur_pos = local_data.Tell();
					if(cur_pos >= next_block_pos) {
						if(!is_terminating) {
							is_terminating = true;
							continue;
						}
						if(is_terminating) {
							break;
						}
					}
				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge(options.output_filename, output_files, n_cores, true);
	string sort_alt_map_file = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		sort_alt_map_file = options.working_prefix + ".alt.map.sort";
	}
	string sort_cmd = (boost::format("sort -k 2,2 -k 3,3n %s > %s") % options.output_filename % sort_alt_map_file).str();
	system(sort_cmd.c_str());
	cout << checker;
}

void AlternativeMapper::create_raw_alternative_mapping_alt() {
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.create_raw_alternative_mapping_alt");
	checker.start();
	string a_path(options.input_filename);
	BamReader SAM;
	if (!SAM.Open(a_path)) {
		return;
	}
	const RefVector& ref_vec = SAM.GetReferenceData();
//	boost::unordered_map<string, int64_t> ref_reverse_index;
//	for (uint64_t ref_id = 0; ref_id < ref_vec.size(); ++ref_id) {
//		auto& a_ref = ref_vec[ref_id];
//		ref_reverse_index[a_ref.RefName] = ref_id;
//	}
//	auto& local_data = SAM.GetBGZF();

	vector<BamAlignment> disc_alg;
	string last_readname;
	BamAlignment al;

	string out_path(options.output_filename);
	ofstream rawalt_fh(out_path, ios::binary);
	while (SAM.LoadNextAlignmentCore(al)) {
		if (al.Name != last_readname) {
			alt_map(disc_alg, rawalt_fh, ref_vec);
			disc_alg.clear();
			disc_alg.push_back(al);
			last_readname = al.Name;
		} else {
			disc_alg.push_back(al);
		}
	}
	SAM.Close();
	alt_map(disc_alg, rawalt_fh, ref_vec);

	if (options.clip) {
		string a_path(options.prefix + ".cl.sorted.bam");
		if(!options.working_dir.empty()) {
			a_path = options.working_prefix + ".cl.sorted.bam";
		}
		BamReader CLSAM;
		if (!CLSAM.Open(a_path)) {
			return;
		}
		const RefVector& ref_vec = CLSAM.GetReferenceData();
		while (CLSAM.LoadNextAlignmentCore(al)) {
			if (al.Name != last_readname) {
				alt_map(disc_alg, rawalt_fh, ref_vec);
				disc_alg.clear();
				disc_alg.push_back(al);
				last_readname = al.Name;
			} else {
				disc_alg.push_back(al);
			}
		}
		alt_map(disc_alg, rawalt_fh, ref_vec);
		CLSAM.Close();
	}
	rawalt_fh.close();
	cerr << "local constructed all mappings\n";
	string sort_alt_map_file = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		sort_alt_map_file = options.working_prefix + ".alt.map.sort";
	}
	string sort_cmd = (boost::format("sort -k 2,2 -k 3,3n %s > %s") % options.output_filename % sort_alt_map_file).str();
	system(sort_cmd.c_str());
	cout << checker;
}

void AlternativeMapper::create_raw_alternative_mapping_par() {
	string sort_alt_map_file = options.prefix + ".alt.map.sort";
	if (!options.working_dir.empty()) {
		sort_alt_map_file = options.working_prefix + ".alt.map.sort";
	}
//	if(boost::filesystem::exists(sort_alt_map_file) && boost::filesystem::exists(options.output_filename)) {
//		return;
//	}
	string disc_sorted_bam = options.prefix + ".disc.sorted.bam";
	string cl_sorted_disc_sorted_bam = options.prefix + ".cl.sorted.disc.sorted.bam";
	if(!options.working_dir.empty()) {
		disc_sorted_bam = options.working_prefix + ".disc.sorted.bam";
		cl_sorted_disc_sorted_bam = options.working_prefix + ".cl.sorted.disc.sorted.bam";
	}
	vector<string> out_file_names;
	out_file_names.push_back(options.output_filename + "_1");
	out_file_names.push_back(options.output_filename + "_2");
	_create_raw_alternative_mapping_par(out_file_names[0], disc_sorted_bam);
	if(options.clip && boost::filesystem::exists(cl_sorted_disc_sorted_bam)) {
		_create_raw_alternative_mapping_par(out_file_names[1], cl_sorted_disc_sorted_bam);
	}

	castle::IOUtils::plain_file_merge(options.output_filename, out_file_names, n_cores, true);
	cerr << "local constructed all mappings\n";

	string sort_cmd = (boost::format("sort -k 2,2 -k 3,3n %s > %s") % options.output_filename % sort_alt_map_file).str();
	system(sort_cmd.c_str());
}

void AlternativeMapper::_create_raw_alternative_mapping_par(const string& output_path, const string& input_path) {
	castle::TimeChecker checker;
	checker.setTarget("AlternativeMapper.create_raw_alternative_mapping_par");
	checker.start();
	cout << (boost::format("[AlternativeMapper.create_raw_alternative_mapping_par] input: %s\n") % input_path).str();
	cout << (boost::format("[AlternativeMapper.create_raw_alternative_mapping_par] output: %s\n") % output_path).str();
	string a_path_index(input_path + ".bfi");
	if (!boost::filesystem::exists(input_path)) {
		cout << checker;
		return;
	}
//	const string debug_read_name = "ST-E00104:502:HFJN5CCXX:7:1123:12266:39510";
	vector<function<void()> > tasks;
	vector<int64_t> block_boundary;
	if (boost::filesystem::exists(a_path_index)) {
		castle::ParallelRunner::load_bfi_index(block_boundary, a_path_index);
	} else {
		const int64_t N_BLOCK_ENTRIES = 262144;
		castle::ParallelRunner::create_bfi_index(block_boundary, input_path, a_path_index, N_BLOCK_ENTRIES);
	}
	int64_t n_block_boundaries = block_boundary.size();
	vector<string> output_files(n_block_boundaries - 1);

	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader SAM;
			if (!SAM.Open(input_path)) {
				return;
			}
			const RefVector& ref_vec = SAM.GetReferenceData();
			boost::unordered_map<string, int64_t> ref_reverse_index;
			for (uint64_t ref_id = 0; ref_id < ref_vec.size(); ++ref_id) {
				auto& a_ref = ref_vec[ref_id];
				ref_reverse_index[a_ref.RefName] = ref_id;
			}
//			auto& ref_blacklist = options.point_black_lists;
				auto& local_data = SAM.GetBGZF();
				if(0 != block_boundary[block_id]) {
					local_data.Seek(block_boundary[block_id]);
				}
				int64_t next_block_pos = block_boundary[block_id + 1];
				int64_t cur_pos = -1;
				string str_block_id = boost::lexical_cast<string>(block_id);

				vector<BamAlignment> disc_alg;
				BamAlignment al;
				string last_readname;

				string out_path(output_path + "." + str_block_id);
				output_files[block_id] = out_path;
				ofstream rawalt_fh(out_path, ios::binary);
//				const bool debug = string::npos != al.Name.find("ar2278164");
				const bool debug = false;
				while (SAM.LoadNextAlignmentCore(al)) {

					if(debug) {
						cout << BamWriter::GetSAMAlignment(al, ref_vec) << "\n";
					}
					if (al.Name != last_readname) {
						alt_map(disc_alg, rawalt_fh, ref_vec);
						disc_alg.clear();
						disc_alg.push_back(al);
						last_readname = al.Name;
					} else {
						disc_alg.push_back(al);
					}
					cur_pos = local_data.Tell();
					if(cur_pos >= next_block_pos) {
						break;
					}
				}
				SAM.Close();
				alt_map(disc_alg, rawalt_fh, ref_vec);
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	castle::IOUtils::plain_file_merge(output_path, output_files, n_cores, true);
	cout << checker;
}

bool AlternativeMapper::alt_map(vector<BamAlignment>& ref_disc_alg, ofstream& rawalt_fh, const RefVector& ref_vec) {
	const bool debug = ref_disc_alg.size() > 0 && string::npos != ref_disc_alg[0].Name.find("ar2278164");
//	const bool debug = false;
	auto& ref_blacklist = options.point_black_lists;
	auto& ref_is = options.is;


	if (ref_disc_alg.size() > 2) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] ref_disc_alg size > 2\n";
		}
		//# 4 entries per pair
		if (ref_disc_alg.size() > 3) {
			if(debug) {
				cout << "[AlternativeMapper.alt_map] ref_disc_alg size > 4\n";
			}
			return false;
		}
		//# 3 entries per pair
		else {
			auto& data1 = ref_disc_alg[0];
			auto& data2 = ref_disc_alg[1];
			auto& data3 = ref_disc_alg[2];
			uint32_t maxmatch_value2 = 0;
			uint32_t maxmatch_value3 = 0;
			//# entry 2 and 3 from the same read
			if ((data2.IsFirstMate() && data3.IsFirstMate()) || (data2.IsSecondMate() && data3.IsSecondMate())) {
				auto cigar2 = data2.CigarData;
				auto cigar3 = data3.CigarData;
				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length > maxmatch_value2) {
						maxmatch_value2 = a_cigar.Length;
					}
				}

				for (uint64_t i = 0; i < cigar3.size(); ++i) {
					auto& a_cigar = cigar3[i];
					if ('M' == a_cigar.Type && a_cigar.Length > maxmatch_value3) {
						maxmatch_value3 = a_cigar.Length;
					}
				}

				int64_t matchup2 = 0;
				int64_t matchup3 = 0;
				int64_t matchdown2 = 0;
				int64_t matchdown3 = 0;

				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value2) {
						break;
					}
					matchup2 += a_cigar.Length;
				}
				for (uint64_t i = 0; i < cigar3.size(); ++i) {
					auto& a_cigar = cigar3[i];
//			#print "$data2[5]\t$i\t$cigar2[$i][0]\t$cigar2[$i][1]\n";
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value3) {
						break;
					}
					matchup3 += a_cigar.Length;
				}

				int64_t tocount = 0;
				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value2) {
						tocount = 1;
						continue;
					}
					if (1 == tocount) {
						matchdown2 += a_cigar.Length;
					}
				}

				tocount = 0;
				for (uint64_t i = 0; i < cigar3.size(); ++i) {
					auto& a_cigar = cigar3[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value3) {
						tocount = 1;
						continue;
					}
					if (1 == tocount) {
						matchdown3 += a_cigar.Length;
					}
				}

				if (data2.IsReverseStrand() && data3.IsReverseStrand()) {
					if (matchup2 < matchup3) {
						data2 = data3;
					}
				} else if (data2.IsReverseStrand() && !data3.IsReverseStrand()) {
					if (matchup2 < matchdown3) {
						data2 = data3;
					}
				} else if (!data2.IsReverseStrand() && data3.IsReverseStrand()) {
					if (matchdown2 < matchup3) {
						data2 = data3;
					}
				} else if (!data2.IsReverseStrand() && !data3.IsReverseStrand()) {
					if (matchdown2 < matchdown3) {
						data2 = data3;
					}
				}
			}
			//# entry 1 and 2 from the same read
			else if ((data1.IsFirstMate() && data2.IsFirstMate()) || (data1.IsSecondMate() && data2.IsSecondMate())) {
				auto& cigar1 = data1.CigarData;
				auto& cigar2 = data2.CigarData;
				uint32_t maxmatch_value1 = 0;
				uint32_t maxmatch_value2 = 0;

				for (uint64_t i = 0; i < cigar1.size(); ++i) {
					auto& a_cigar = cigar1[i];
					if ('M' == a_cigar.Type && a_cigar.Length > maxmatch_value1) {
						maxmatch_value1 = a_cigar.Length;
					}
				}

				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length > maxmatch_value2) {
						maxmatch_value2 = a_cigar.Length;
					}
				}

				int64_t matchup1 = 0;
				int64_t matchup2 = 0;
				int64_t matchdown1 = 0;
				int64_t matchdown2 = 0;

				for (uint64_t i = 0; i < cigar1.size(); ++i) {
					auto& a_cigar = cigar1[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value1) {
						break;
					}
					matchup1 += a_cigar.Length;
				}
				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value2) {
						break;
					}
					matchup2 += a_cigar.Length;
				}

				int64_t tocount = 0;
				for (uint64_t i = 0; i < cigar1.size(); ++i) {
					auto& a_cigar = cigar1[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value1) {
						tocount = 1;
						continue;
					}
					if (1 == tocount) {
						matchdown1 += a_cigar.Length;
					}
				}
				tocount = 0;
				for (uint64_t i = 0; i < cigar2.size(); ++i) {
					auto& a_cigar = cigar2[i];
					if ('M' == a_cigar.Type && a_cigar.Length == maxmatch_value2) {
						tocount = 1;
						continue;
					}
					if (1 == tocount) {
						matchdown2 += a_cigar.Length;
					}
				}

				if (data1.IsReverseStrand() && data2.IsReverseStrand()) {
					if (matchup1 < matchup2) {
						data1 = data2;
						data2 = data3;
					} else {
						data2 = data3;
					}

				} else if (data1.IsReverseStrand() && !data2.IsReverseStrand()) {
					if (matchup1 < matchdown2) {
						data1 = data2;
						data2 = data3;
					} else {
						data2 = data3;
					}
				} else if (!data1.IsReverseStrand() && data2.IsReverseStrand()) {
					if (matchdown1 < matchup2) {
						data1 = data2;
						data2 = data3;
					} else {
						data2 = data3;
					}
				} else if (!data1.IsReverseStrand() && !data2.IsReverseStrand()) {
					if (matchdown1 < matchdown2) {
						data1 = data2;
						data2 = data3;
					} else {
						data2 = data3;
					}
				}
			}
			data1.Name += "sc";
			data2.Name += "sc";
		}
	}

	if (ref_disc_alg.size() < 2) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] ref_disc_alg size < 2\n";
		}
		return false;
	}

//	return 0 unless ( defined( $data2[0] ) );    # skip unpaired reads

	auto data1 = ref_disc_alg[0];
	auto data2 = ref_disc_alg[1];
	auto data1_ref_id = ref_vec[data1.RefID].RefName;
	auto data2_ref_id = ref_vec[data2.RefID].RefName;

	if ((data1_ref_id > data2_ref_id) || (data1.RefID == data2.RefID && data1.Position > data2.Position)) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] swap data\n";
		}
		swap(data1, data2);
		swap(data1_ref_id, data2_ref_id);
//		auto& temp = data2;
//		data2 = data1;
//		data1 = temp;
	}
	if (options.min_mapq) {
		if (data1.MapQuality < options.min_mapq || data2.MapQuality < options.min_mapq) {
			return false;
		}
	}
	auto readname = data1.Name;

//# a read is non-uniq mapped if XT tag with R, or has XA tag or mapping qual is 0

	char xt = 'U';
	char mate_xt = 'U';
	string xt_tag_str;
	string mate_xt_tag_str;
	string xa_tag_str;
	string mate_xa_tag_str;

	if (data1.GetTag("XT", xt_tag_str)) {
		xt = xt_tag_str[0];
	}
	if (data1.GetTag("XA", xa_tag_str)) {
		xt = 'R';
	}
	if (data2.GetTag("XT", mate_xt_tag_str)) {
		mate_xt = mate_xt_tag_str[0];
	}
	if (data2.GetTag("XA", mate_xa_tag_str)) {
		mate_xt = 'R';
	}
	if (!data1.AlignmentFlag) {
		xt = 'R';
	}
	if (!data2.AlignmentFlag) {
		mate_xt = 'R';
	}
	if (options.ad_align) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] ad_align1\n";
		}
		if (('R' == xt && 'R' == mate_xt) || 'M' == xt || 'M' == mate_xt) {
			if(debug) {
				cout << "[AlternativeMapper.alt_map] R == xt && R == mate_xt || M == xt || M == mate_xt\n";
			}
			return false;
		}
	} else {
		if (!options.use_all_align) {
			if(debug) {
				cout << "[AlternativeMapper.alt_map] !use_all_align\n";
			}
			if (!('U' == xt && 'U' == mate_xt)) {
				if(debug) {
					cout << "[AlternativeMapper.alt_map] !('U' == xt && 'U' == mate_xt)\n";
				}
				return false;
			}
		}
	}

	if (-1 == data1.RefID || -1 == data2.RefID) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] unmapped read pair\n";
		}
		return false;
	}
	string seqid = ref_vec[data1.RefID].RefName;
	string mseqid = ref_vec[data2.RefID].RefName;
	auto a_key_1 = make_pair(seqid, data1.Position);
	auto a_key_2 = make_pair(mseqid, data2.Position);
	auto a_key_3 = make_pair(seqid, static_cast<int64_t>(data1.Position + (options.cut_sr << 1)));
	auto a_key_4 = make_pair(mseqid, static_cast<int64_t>(data2.Position + (options.cut_sr << 1)));
	auto the_black_itr_1 = ref_blacklist.find(a_key_1);
	auto the_black_itr_2 = ref_blacklist.find(a_key_2);
	auto the_black_itr_3 = ref_blacklist.find(a_key_3);
	auto the_black_itr_4 = ref_blacklist.find(a_key_4);
	if (ref_blacklist.end() != the_black_itr_1 || ref_blacklist.end() != the_black_itr_2 || ref_blacklist.end() != the_black_itr_3 || ref_blacklist.end() != the_black_itr_4) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] blacklisted\n";
		}
		return false;
	}

	int64_t start = data1.Position;
	int64_t mstart = data2.Position;
	int64_t strand = 1;
	if (data1.IsReverseStrand()) {
		strand = -1;
	}

	int64_t mstrand = 1;
	if (data2.IsReverseStrand()) {
		mstrand = -1;
	}

	int64_t len = 0;
	int64_t mlen = 0;

	auto& cigar1 = data1.CigarData;
	auto& cigar2 = data2.CigarData;

	for (uint64_t i = 0; i < cigar1.size(); ++i) {
		auto& a_cigar = cigar1[i];
		if ('M' == a_cigar.Type && a_cigar.Length > len) {
			len = a_cigar.Length;
		}
	}

	for (uint64_t i = 0; i < cigar2.size(); ++i) {
		auto& a_cigar = cigar2[i];
		if ('M' == a_cigar.Type && a_cigar.Length > mlen) {
			mlen = a_cigar.Length;
		}
	}

	const int64_t MIN_VALUE = numeric_limits<int64_t>::min();
	int64_t isize = MIN_VALUE;
	if (seqid == mseqid) {
		if (strand == mstrand) {
			isize = abs(start - mstart) + 1;
		} else if (1 == strand && -1 == mstrand) {
			isize = abs(start - mstart) + mlen;
		} else if (strand == -1 && mstrand == 1) {
			isize = abs(start - mstart) - mlen + 2;
		}
	}

	string rg;
//	my $rg = 'none';
	if (!data1.GetReadGroup(rg)) {
		rg = "none";
	}
	if (seqid == mseqid && 1 == strand && -1 == mstrand && isize <= ref_is[rg]["isu"]) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] seqid == mseqid && 1 == strand && -1 == mstrand && isize <= ref_is[rg][isu]\n";
		}
		return false;
	}

	uint32_t nm = 0;
	uint32_t mnm = 0;
	int64_t tr = 0;
	int64_t mtr = 0;
	vector<string> alt_mapping;
	vector<string> an_alt_mapping;
	vector<string> malt_mapping;
	vector<string> an_malt_mapping;
	int64_t n_alt_mapping = 0;
	int64_t n_malt_mapping = 0;
	string mxa_tag_str;
	bool xa_found = false;
	bool mxa_found = false;
	const char* delim_semi_colon = ";";
	const char* delim_comma = ",";
	if (options.ad_align) {
		if(debug) {
			cout << "[AlternativeMapper.alt_map] ad_align2\n";
		}
		data1.GetEditDistance(nm);
		data2.GetEditDistance(mnm);
		if ('R' == xt) {
			xa_found = data1.GetTag("XA", xa_tag_str);
			if (!xa_found) {
				tr = 1;
			}
			castle::StringUtils::c_string_multi_split(xa_tag_str, delim_semi_colon, alt_mapping);
			n_alt_mapping = alt_mapping.size();
			if (-1 != options.alt_map_max && n_alt_mapping > options.alt_map_max) {
				tr = 1;
			}
		}
		if ('R' == mate_xt) {
			mxa_found = data2.GetTag("XA", mxa_tag_str);
			if (!mxa_found) {
				mtr = 1;
			}
			castle::StringUtils::c_string_multi_split(mxa_tag_str, delim_semi_colon, malt_mapping);
			n_malt_mapping = malt_mapping.size();
			if (-1 != options.alt_map_max && n_malt_mapping > options.alt_map_max) {
				mtr = 1;
			}
		}
		if ('R' == xt || 'R' == mate_xt) {
			if (1 == tr || 1 == mtr) {
				if(debug) {
					cout << "[AlternativeMapper.alt_map] 1 == tr || 1 == mtr\n";
				}
				return false;
			}
		}
		double weight = 1;
		if ('R' == xt || 'R' == mate_xt) {
			weight = (n_alt_mapping > n_malt_mapping) ? n_alt_mapping + 1 : n_malt_mapping + 1;
			weight = (double) 1 / weight;
//			$weight = sprintf( "%.9f", 1 / $weight );
		}
		vector<string> toprint;
		string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % readname % seqid % (start + 1) % strand % mseqid
				% (mstart + 1) % mstrand % ((isize == MIN_VALUE) ? "" : boost::lexical_cast<string>(isize)) % nm % mnm % rg % len % mlen % weight).str();
		if(debug) {
			cout << "[AlternativeMapper.alt_map] here-0\n";
			cout << an_entry;
		}
		toprint.push_back(an_entry);
		for (auto& an_alt : alt_mapping) {
			castle::StringUtils::c_string_multi_split(an_alt, delim_comma, an_alt_mapping);
			string seqid_alt = an_alt_mapping[0];
			string start_alt_str = an_alt_mapping[1];
			int64_t start_alt = boost::lexical_cast<int64_t>(an_alt_mapping[1]);
			string mapping = an_alt_mapping[2];
			uint32_t nm_alt = boost::lexical_cast<uint32_t>(an_alt_mapping[3]);
			int64_t nosign_start_alt = boost::lexical_cast<int64_t>(start_alt_str.substr(1));
			auto a_key_1 = make_pair(seqid_alt, nosign_start_alt);
			auto a_key_2 = make_pair(seqid_alt, nosign_start_alt + (options.cut_sr << 1));
			auto the_black_itr_1 = ref_blacklist.find(a_key_1);
			auto the_black_itr_2 = ref_blacklist.find(a_key_2);
			if (ref_blacklist.end() != the_black_itr_1 || ref_blacklist.end() != the_black_itr_2) {
				continue;
			}

			string mseqid_alt = mseqid;
			int64_t mstart_alt = mstart;
			int64_t mstrand_alt = mstrand;
			uint32_t mnm_alt = mnm;
			int64_t len_alt = len;
			int64_t mlen_alt = mlen;
			if (seqid == seqid_alt && start == start_alt) {
				continue;
			}
			int64_t strand_alt = 1;
			if (start_alt < 0) {
				strand_alt = -1;
			}
			start_alt = nosign_start_alt;
			int64_t isize_alt = MIN_VALUE;
			if (seqid_alt == mseqid_alt) {
				if (start_alt <= mstart_alt) {
					if (strand_alt == mstrand_alt) {
						isize_alt = abs(mstart_alt - start_alt) + 1;
					} else {
						if (-1 == strand_alt && 1 == mstrand_alt) {
							isize_alt = abs(mstart_alt - start_alt) - mlen_alt + 2;
						} else {
							isize_alt = abs(mstart_alt - start_alt) + mlen_alt;
						}
					}
				} else {
					if (strand_alt == mstrand_alt) {
						isize_alt = abs(start_alt - mstart_alt) + 1;
					} else {
						if (-1 == strand_alt && 1 == mstrand_alt) {
							isize_alt = abs(start_alt - mstart_alt) + mlen_alt;
						} else {
							isize_alt = abs(start_alt - mstart_alt) - mlen_alt + 2;
						}
					}
					swap(seqid_alt, mseqid_alt);
					swap(start_alt, mstart_alt);
					swap(strand_alt, mstrand_alt);
					swap(nm_alt, mnm_alt);
					swap(len_alt, mlen_alt);
				}
			} else if (seqid_alt > mseqid_alt) {
				swap(seqid_alt, mseqid_alt);
				swap(start_alt, mstart_alt);
				swap(strand_alt, mstrand_alt);
				swap(nm_alt, mnm_alt);
				swap(len_alt, mlen_alt);
			}
			if (seqid_alt == mseqid_alt && 1 == strand_alt && -1 == mstrand_alt && isize_alt <= ref_is[rg]["isu"]) {
				continue;
			}
			string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
			% readname % seqid_alt % (start_alt + 1) % strand_alt % mseqid_alt % (mstart_alt + 1) % mstrand_alt % ((isize_alt == MIN_VALUE) ? "" : boost::lexical_cast<string>(isize_alt)) % nm_alt % mnm_alt % rg % len_alt % mlen_alt % weight).str();
			if(debug) {
				cout << "[AlternativeMapper.alt_map] here-1\n";
				cout << an_entry;
			}
			toprint.push_back(an_entry);
		}
		for (auto& an_alt : malt_mapping) {
			castle::StringUtils::c_string_multi_split(an_alt, delim_comma, an_alt_mapping);
			string mseqid_alt = an_alt_mapping[0];
			string mstart_alt_str = an_alt_mapping[1];
			int64_t mstart_alt = boost::lexical_cast<int64_t>(an_alt_mapping[1]);
			string mmapping = an_alt_mapping[2];
			uint32_t mnm_alt = boost::lexical_cast<uint32_t>(an_alt_mapping[3]);
			int64_t nosign_mstart_alt = boost::lexical_cast<int64_t>(mstart_alt_str.substr(1));
			//		( $seqid_alt, $start_alt, $mapping, $nm_alt ) = split( /,/, $_ );
			//		my $nosign_start_alt = substr( $start_alt, 1 );
			auto a_key_1 = make_pair(mseqid_alt, nosign_mstart_alt);
			auto a_key_2 = make_pair(mseqid_alt, nosign_mstart_alt + (options.cut_sr << 1));
			auto the_black_itr_1 = ref_blacklist.find(a_key_1);
			auto the_black_itr_2 = ref_blacklist.find(a_key_2);
			if (ref_blacklist.end() != the_black_itr_1 || ref_blacklist.end() != the_black_itr_2) {
				continue;
			}

			string seqid_alt = seqid;
			int64_t start_alt = start;
			int64_t strand_alt = strand;
			uint32_t nm_alt = nm;
			int64_t len_alt = len;
			int64_t mlen_alt = mlen;
			if (mseqid == mseqid_alt && mstart == mstart_alt) {
				continue;
			}
			int64_t mstrand_alt = 1;
			if (mstart_alt < 0) {
				mstrand_alt = -1;
			}
			mstart_alt = nosign_mstart_alt;
			int64_t isize_alt = MIN_VALUE;

			if (seqid_alt == mseqid_alt) {
				if (start_alt <= mstart_alt) {
					if (strand_alt == mstrand_alt) {
						isize_alt = abs(mstart_alt - start_alt) + 1;
					} else {
						if (-1 == strand_alt && 1 == mstrand_alt) {
							isize_alt = abs(mstart_alt - start_alt) - mlen_alt + 2;
						} else {
							isize_alt = abs(mstart_alt - start_alt) + mlen_alt;
						}
					}
				} else {
					if (strand_alt == mstrand_alt) {
						isize_alt = abs(start_alt - mstart_alt) + 1;
					} else {
						if (-1 == strand_alt && 1 == mstrand_alt) {
							isize_alt = abs(start_alt - mstart_alt) + mlen_alt;
						} else {
							isize_alt = abs(start_alt - mstart_alt) - mlen_alt + 2;
						}
					}
					swap(seqid_alt, mseqid_alt);
					swap(start_alt, mstart_alt);
					swap(strand_alt, mstrand_alt);
					swap(nm_alt, mnm_alt);
					swap(len_alt, mlen_alt);
				}
			} else if (seqid_alt > mseqid_alt) {
				swap(seqid_alt, mseqid_alt);
				swap(start_alt, mstart_alt);
				swap(strand_alt, mstrand_alt);
				swap(nm_alt, mnm_alt);
				swap(len_alt, mlen_alt);
			}

			if (seqid_alt == mseqid_alt && 1 == strand_alt && -1 == mstrand_alt && isize_alt <= ref_is[rg]["isu"]) {
				continue;
			}
			string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
			% readname % seqid_alt % (start_alt + 1) % strand_alt % mseqid_alt % (mstart_alt + 1) % mstrand_alt % ((isize_alt == MIN_VALUE) ? "" : boost::lexical_cast<string>(isize_alt)) % nm_alt % mnm_alt % rg % len_alt % mlen_alt % weight).str();
			if(debug) {
				cout << "[AlternativeMapper.alt_map] here-2\n";
				cout << an_entry;
			}
			toprint.push_back(an_entry);
		}
		for (auto& an_entry : toprint) {
			if(debug) {
				cout << "[AlternativeMapper.alt_map] here-3\n";
				cout << an_entry;
			}
			rawalt_fh << an_entry;
		}
	} else {
		string an_entry = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% readname % seqid % (start + 1) % strand % mseqid % (mstart + 1) % mstrand % ((isize == MIN_VALUE) ? "" : boost::lexical_cast<string>(isize)) % nm % mnm % rg % len % mlen).str();
		if(debug) {
			cout << "[AlternativeMapper.alt_map] here-4\n";
			cout << an_entry;
		}
		rawalt_fh << an_entry;
	}

	return true;
}

bool AlternativeMapper::write_alt_map(ostream& out, const RefVector& ref_vec, const boost::unordered_map<string, int64_t>& ref_reverse_index, vector<BamAlignment>& recent_als) {
	if (2 != recent_als.size()) {
		return false;
	}

	bool found_valid_entry = false;
//	auto& ref_blacklist = options.point_black_lists;
	if ((recent_als[0].RefID > recent_als[1].RefID) || (recent_als[0].RefID == recent_als[1].RefID && recent_als[0].Position > recent_als[1].Position)) {
		iter_swap(recent_als.begin(), recent_als.begin() + 1);
	}

	auto& data1 = recent_als[0];
	auto& data2 = recent_als[1];
//	const bool debug = string::npos != recent_als[0].Name.find("ST-E00104:502:HFJN5CCXX:3:1123:31000:16270");
//	if (debug) {
//		cout << "[AlternativeMapper.write_alt_map] here-0-1-a: " << BamWriter::GetSAMAlignment(data1, ref_vec) << "\n";
//		cout << "[AlternativeMapper.write_alt_map] here-0-1-b: " << BamWriter::GetSAMAlignment(data2, ref_vec) << "\n";
//	}
//					if(-1 == data1.RefID || -1 == data2.RefID) {
//						recent_als.clear();
//						recent_als.push_back(al);
//						continue;
//					}
//					if (options.min_mapq) {
//						if ((data1.MapQuality < options.min_mapq)
//								|| data2.MapQuality < options.min_mapq) {
//							recent_als.clear();
//							recent_als.push_back(al);
//							continue;
//						}
//					}
// the XT tag won't be filled by bwa mem, but I just adopted codes from the original meerkat Perl script.

	char xt = 'U';
	char mate_xt = 'U';
	string xt_tag_str;
	string mate_xt_tag_str;
	if (data1.GetTag("XT", xt_tag_str)) {
		xt = xt_tag_str[0];
	}
	if (data2.GetTag("XT", mate_xt_tag_str)) {
		mate_xt = mate_xt_tag_str[0];
	}
	if (!data1.AlignmentFlag) {
		xt = 'R';
	}
	if (!data2.AlignmentFlag) {
		mate_xt = 'R';
	}
	if (options.ad_align) {
		if (('R' == xt && 'R' == mate_xt) || 'M' == xt || 'M' == mate_xt) {
//			if(debug) {
//				cout << "[AlternativeMapper.write_alt_map] here-0-2: xt: " << xt_tag_str << ", mate_xt: " << mate_xt_tag_str << "\n";
//			}
			return false;
		} else {
			if (!options.use_all_align) {
				if (!('U' == xt && 'U' == mate_xt)) {
//					if(debug) {
//						cout << "[AlternativeMapper.write_alt_map] here-0-2: xt: " << xt_tag_str << ", mate_xt: " << mate_xt_tag_str << "\n";
//					}
					return false;
				}
			}
		}
	}
// check coverage based black list
////	auto a_key_1 = make_pair(ref_vec[data1.RefID].RefName,
////			data1.Position);
////	auto a_key_2 = make_pair(ref_vec[data2.RefID].RefName,
////			data2.Position);
////	auto a_key_3 = make_pair(ref_vec[data1.RefID].RefName,
////			static_cast<int64_t>(data1.Position
////					+ (options.cut_sr << 1)));
////	auto a_key_4 = make_pair(ref_vec[data2.RefID].RefName,
////			static_cast<int64_t>(data2.Position
////					+ (options.cut_sr << 1)));
////	auto the_black_itr_1 = ref_blacklist.find(a_key_1);
////	auto the_black_itr_2 = ref_blacklist.find(a_key_2);
////	auto the_black_itr_3 = ref_blacklist.find(a_key_3);
////	auto the_black_itr_4 = ref_blacklist.find(a_key_4);
//	if (ref_blacklist.end() == the_black_itr_1
//			&& ref_blacklist.end() == the_black_itr_2
//			&& ref_blacklist.end() == the_black_itr_3
//			&& ref_blacklist.end() == the_black_itr_4) {
	if (!has_valid_insertion_size(out, ref_vec, ref_reverse_index, recent_als, found_valid_entry)) {
		if (found_valid_entry) {
//				if(debug) {
//					cout << "[AlternativeMapper.write_alt_map] found valid entry\n";
//				}
			return true;
		} else {
//				if(debug) {
//					cout << "[AlternativeMapper.write_alt_map] not found valid entry\n";
//				}
			return false;
		}
	}
//	}
//	if(debug) {
//		cout << "[AlternativeMapper::write_alt_map] black listed\n";
//	}
	return found_valid_entry;
}
} /* namespace meerkat */
