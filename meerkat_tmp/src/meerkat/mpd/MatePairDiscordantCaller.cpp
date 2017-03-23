/*
 * MatePairDiscordantCaller.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

/*
 * MatePairDiscordantCaller.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: el174
 */

#include "MatePairDiscordantCaller.hpp"

namespace meerkat {

MatePairDiscordantCaller::MatePairDiscordantCaller() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

MatePairDiscordantCaller::~MatePairDiscordantCaller() {
}
void MatePairDiscordantCaller::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}

void MatePairDiscordantCaller::select_candidate_events_in_clusters() {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.select_candidate_events_in_clusters");
	checker.start();
	select_discordant_clusters();
	call_initial_events();
	cout << checker;
}

void MatePairDiscordantCaller::select_candidate_events_in_clusters_alt() {
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
//	if(0 < castle::IOUtils::get_file_size(mpintra_outfile) && 0 < castle::IOUtils::get_file_size(mpinter_outfile)) {
//		return;
//	}
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.select_candidate_events_in_clusters_alt");
	checker.start();
	select_discordant_clusters();
	call_initial_events_serial();
	cout << checker;
}

void MatePairDiscordantCaller::select_discordant_clusters() {
	string clusterfile = options.prefix + ".clusters";
	string discclfile = options.prefix + ".discord";
	if(!options.working_dir.empty()) {
		clusterfile = options.working_prefix + ".clusters";
		discclfile = options.working_prefix + ".discord";
	}
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.select_discordant_clusters");
	checker.start();
	vector<function<void()> > tasks;

	const int64_t BLOCK_SIZE = 4 * 1024 * 1024;
	const int64_t ref_file_size = castle::IOUtils::get_file_size(clusterfile);
	cout << (boost::format("[MatePairDiscordantCaller.select_discordant_clusters] Ref. File size: %d\n") % ref_file_size).str();
	int64_t n_blocks = (ref_file_size / (double) BLOCK_SIZE) + 1;
	vector<int64_t> block_boundary;
	if(ref_file_size > BLOCK_SIZE) {
		block_boundary.resize(n_blocks);
		block_boundary[0] = 0;
		block_boundary[n_blocks - 1] = numeric_limits<int64_t>::max();
		for (int64_t block_id = 1; block_id < n_blocks - 1; ++block_id) {
			tasks.push_back([&, ref_file_size, block_id, BLOCK_SIZE] {
				int64_t cur_boundary_pos = block_id * BLOCK_SIZE;
				if(cur_boundary_pos >= ref_file_size) {
					return;
				}
				string line;
				ifstream in(clusterfile, ios::binary);
				in.seekg(cur_boundary_pos, ios::beg);
				getline(in, line, '\n');
				cur_boundary_pos += line.size() + 1;
				string previous_pid;
				while(getline(in, line, '\n')) {
					cout << line << "\n";
					size_t pos = line.find_first_of('\t');
					string cur_pid = line.substr(0, pos);
					if(previous_pid.empty()) {
						previous_pid = cur_pid;
						cur_boundary_pos += line.size() + 1;
						continue;
					}

					if(previous_pid != cur_pid) {
						break;
					}
					cur_boundary_pos += line.size() + 1;
				}
				block_boundary[block_id] = cur_boundary_pos;
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	} else {
		n_blocks = 2;
		block_boundary.resize(n_blocks);
		block_boundary[0] = 0;
		block_boundary[n_blocks - 1] = numeric_limits<int64_t>::max();
	}
	cout << (boost::format("[MatePairDiscordantCaller.select_discordant_clusters] # blocks: %d\n") % n_blocks).str();
	{
		vector<string> patterns;
		patterns.push_back("mu1");
		patterns.push_back("mu2");
		patterns.push_back("sc");
		vector<map<string, int64_t>> support_lists(n_blocks - 1);
		vector<map<string, int64_t>> support_f_lists(n_blocks - 1);
		vector<set<string>> read_names_lists(n_blocks - 1);
		vector<vector<ClusterEntry>> cluster_lists(n_blocks - 1);

		vector<map<string, map<string, OrientationEntry>>> orientation_lists(n_blocks - 1);
		vector<map<string, map<string, pair<int64_t, int64_t>>> > starts_map_lists(n_blocks - 1);
		vector<map<string, map<string, pair<int64_t, int64_t>>> > mstarts_map_lists(n_blocks - 1);
		vector<map<string, map<string, vector<int64_t>>> > cbps_map_lists(n_blocks - 1);
		vector<map<string, map<string, map<int64_t, vector<string>>> >> cluster_type_lists(n_blocks - 1);
		for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				int64_t cur_boundary_pos = block_boundary[block_id];
				int64_t cur_pos = cur_boundary_pos;
				int64_t next_boundary_pos = block_boundary[block_id + 1];

				vector<string> data;
				string lastpid;
				string lastsid;
				int64_t k = 0;		// number of supporting entries
				int64_t l = 0;		// number of supporting read pairs
				int64_t j = 0;		// number of supporting full length read pairs
				string line;
				const char* delims = "\t";
				auto& local_cluster = cluster_lists[block_id];
				auto& read_names = read_names_lists[block_id];
				auto& local_support = support_lists[block_id];
				auto& local_support_f = support_f_lists[block_id];
				auto& local_orientation = orientation_lists[block_id];
				auto& local_starts_map = starts_map_lists[block_id];
				auto& local_mstarts_map = mstarts_map_lists[block_id];
				auto& local_cbps_map = cbps_map_lists[block_id];
				auto& local_cluster_type = cluster_type_lists[block_id];

				string str_block_id = boost::lexical_cast<string>(block_id);

				ofstream out(discclfile + "." + str_block_id, ios::binary);
				ifstream in(clusterfile, ios::binary);
				in.seekg(cur_boundary_pos, ios::beg);
				while (getline(in, line, '\n')) {
					cur_pos += line.size() + 1;
					castle::StringUtils::tokenize(line, delims, data);
					if (data.size() < 16) {
						continue;
					}

					ClusterEntry ce;
					ce.readname = data[5];
//const bool debug = "ST-E00104:502:HFJN5CCXX:1:1102:11525:23073" == ce.readname;
//const bool debug = "8" == data[0];
					ce.rg = data[6];
					ce.seqid = data[7];
					ce.strand = boost::lexical_cast<int32_t>(data[8]);
					ce.start = boost::lexical_cast<int64_t>(data[9]);
					ce.len = boost::lexical_cast<int64_t>(data[10]);
					ce.end = ce.start + ce.len - 1;
					ce.mseqid = data[11];
					ce.mstrand = boost::lexical_cast<int32_t>(data[12]);
					ce.mstart = boost::lexical_cast<int64_t>(data[13]);
					ce.mlen = boost::lexical_cast<int64_t>(data[14]);
					ce.mend = ce.mstart + ce.mlen - 1;
					ce.isize = boost::lexical_cast<int64_t>(data[15]);
					if (local_cluster.empty() || (lastpid == data[0] && lastsid == data[1])) {
						local_cluster.push_back(ce);
						if (read_names.end() == read_names.find(ce.readname)) {
							++k;
							read_names.insert(ce.readname);
						}
						size_t pos = castle::StringUtils::rfind_one_of_multiple_patterns(
								ce.readname, patterns);
						if (string::npos != pos) {
							string original_read_name = ce.readname.substr(0, pos);
//if (read_names.end() == read_names.find(original_read_name)) {
							++l;
//}
							read_names.insert(original_read_name);
//cout << original_read_name << "/" << ce.readname << "\n";
						} else {
							++j;
//if (read_names.end() == read_names.find(ce.readname)) {
							++l;
//}
							read_names.insert(ce.readname);
						}
					} else {
						local_support[lastpid] = l;
						local_support_f[lastpid] = j;
//						if(lastpid == "11946") {
//							cout << "lastpid: " << l << "/" << j << "\n";
//						}
						determine_a_cluster(out, local_cluster, local_orientation, local_starts_map, local_mstarts_map, local_cbps_map, local_cluster_type, local_support, local_support_f, lastpid, lastsid, l, j);
						local_cluster.clear();
						read_names.clear();
						k = 0;
						l = 0;
						j = 0;

						local_cluster.push_back(ce);
						if (read_names.end() == read_names.find(ce.readname)) {
							++k;
							read_names.insert(ce.readname);
						}
						size_t pos = castle::StringUtils::rfind_one_of_multiple_patterns(
								ce.readname, patterns);
						if (string::npos != pos) {
							string original_read_name = ce.readname.substr(0, pos);
//if (read_names.end() == read_names.find(original_read_name)) {
							++l;
//}
							read_names.insert(original_read_name);
//cout << original_read_name << "/" << ce.readname << "\n";
						} else {
							++j;
//if (read_names.end() == read_names.find(ce.readname)) {
							++l;
//}
							read_names.insert(ce.readname);
						}
					}
					lastpid = data[0];
					lastsid = data[1];
					if(cur_pos >= next_boundary_pos) {
						break;
					}
				}
				{
// process the last cluster
					local_support[lastpid] = l;
					local_support_f[lastpid] = j;
					determine_a_cluster(out, local_cluster, local_orientation, local_starts_map, local_mstarts_map, local_cbps_map, local_cluster_type, local_support, local_support_f, lastpid, lastsid, l, j);
					local_cluster.clear();
					read_names.clear();
				}
				if(block_id == n_blocks - 2) {
					the_last_pid = lastpid;
				}
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		{
			ofstream out(discclfile, ios::binary);
			for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				string str_block_id = boost::lexical_cast<string>(block_id);
				out << castle::IOUtils::read_fully(discclfile + "." + str_block_id);
			}
		}
		{
			for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				tasks.push_back([&, block_id] {
					string str_block_id = boost::lexical_cast<string>(block_id);
					boost::filesystem::remove(discclfile + "." + str_block_id);
				});
			}
			castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		}
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_orientation = orientation_lists[block_id];
				for(auto& pid_entry: local_orientation) {
					for(auto& sid_entry : pid_entry.second) {
//						bool debug = (pid_entry.first == "6710");
//						if(debug) {
//							cout << "orientation: " << sid_entry.second.strand << "/" << sid_entry.second.mate_strand << "\n";
//						}
						orientation[pid_entry.first][sid_entry.first] = sid_entry.second;
					}
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_starts_map = starts_map_lists[block_id];
				for(auto& pid_entry : local_starts_map) {
					for(auto& sid_entry : pid_entry.second) {
						if(sid_entry.second.first > sid_entry.second.second) {
							swap(sid_entry.second.first, sid_entry.second.second);
						}
//						if("22" == pid_entry.first) {
//							cout << "pid_entry_second: " << sid_entry.first << ": " << sid_entry.second.first << "/" << sid_entry.second.second << "\n";
//						}
						starts_map[pid_entry.first][sid_entry.first] = sid_entry.second;
					}
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_mstarts_map = mstarts_map_lists[block_id];
				for(auto& pid_entry : local_mstarts_map) {
					for(auto& sid_entry : pid_entry.second) {
						if(sid_entry.second.first > sid_entry.second.second) {
							swap(sid_entry.second.first, sid_entry.second.second);
						}
						mstarts_map[pid_entry.first][sid_entry.first] = sid_entry.second;
					}
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_support = support_lists[block_id];
				for(auto& pid_entry : local_support) {
//					bool debug = pid_entry.first == "2449";
//					if(debug) {
//						cout << "select discordant cluster: " << pid_entry.first << "/" << pid_entry.second << "\n";
//					}
					support[pid_entry.first] = pid_entry.second;
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_support_f = support_f_lists[block_id];
				for(auto& pid_entry : local_support_f) {
					support_f[pid_entry.first] = pid_entry.second;
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_cbps_map = cbps_map_lists[block_id];
				for(auto& pid_entry : local_cbps_map) {
					for(auto& sid_entry : pid_entry.second) {
						auto& tmp_vec = cbps_map[pid_entry.first][sid_entry.first];
						tmp_vec.insert(tmp_vec.end(), sid_entry.second.begin(), sid_entry.second.end());
					}
				}
			}
		});
		tasks.push_back([&, n_blocks] {
			for(int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
				auto& local_cluster_type = cluster_type_lists[block_id];
				for(auto& pid_entry : local_cluster_type) {
					for(auto& sid_entry : pid_entry.second) {
						for(auto& the_type_of_strand : sid_entry.second) {
							auto& tmp_vec = cluster_type[pid_entry.first][sid_entry.first][the_type_of_strand.first];
//							bool debug = ("1" == pid_entry.first && "hs37d5" == sid_entry.first && 3 == the_type_of_strand.first);
//							if(debug) {
//								for(auto& the_last_entry : the_type_of_strand.second) {
//									cout << "the_last_entry: " << the_last_entry << "\n";
//								}
//							}
							tmp_vec.insert(tmp_vec.end(), the_type_of_strand.second.begin(), the_type_of_strand.second.end());
						}
					}
				}
			}
		});
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	}
	{
		const int32_t N_TYPES_OF_STRANDS = 4;
		the_intervals.resize(N_TYPES_OF_STRANDS);
		for (int8_t type_id = 0; type_id < N_TYPES_OF_STRANDS; ++type_id) {
//			tasks.push_back([&, type_id] {
				vector<string> cols;
				for(auto& first_chr_entry : cluster_type) {
					for(auto& second_chr_entry : first_chr_entry.second) {
						for(auto& cluster_entry : second_chr_entry.second) {
							auto the_type = cluster_entry.first;
							if(the_type != type_id) {
								continue;
							}
							auto& the_local_interval = the_intervals[the_type];
							for(auto the_id : cluster_entry.second) {
								castle::StringUtils::c_string_multi_split(the_id, "_", cols);
								auto& p_id = cols[0];
								auto& s_id = cols[1];
								IntervalEntryType additional_entry;
								additional_entry.p_id = p_id;
								additional_entry.s_id = s_id;
								auto start_map_start = starts_map[p_id][s_id].first;
								auto start_map_end = starts_map[p_id][s_id].second;
								if(start_map_start < start_map_end) {
									swap(start_map_start, start_map_end);
								}
								auto mstart_map_start = mstarts_map[p_id][s_id].first;
								auto mstart_map_end = mstarts_map[p_id][s_id].second;
								if(mstart_map_start < mstart_map_end) {
									swap(mstart_map_start, mstart_map_end);
								}
//								if(starts_map[p_id][s_id].first > starts_map[p_id][s_id].second) {
//									swap(starts_map[p_id][s_id].first, starts_map[p_id][s_id].second);
//								}
//								if(mstarts_map[p_id][s_id].first > mstarts_map[p_id][s_id].second) {
//									swap(mstarts_map[p_id][s_id].first, mstarts_map[p_id][s_id].second);
//								}
//IntervalEntry an_entry(starts_map[p_id][s_id].first, starts_map[p_id][s_id].second, additional_entry);
//IntervalEntry mate_entry(mstarts_map[p_id][s_id].first, mstarts_map[p_id][s_id].second, additional_entry);
//the_local_interval.push_back(an_entry);
//the_local_interval.push_back(mate_entry);
					IntervalEntry an_entry(start_map_start, mstart_map_end, additional_entry);
					the_local_interval.push_back(an_entry);
				}
			}
		}
	}
//});
		}
//		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	}

	cout << checker;
}

void MatePairDiscordantCaller::select_discordant_clusters_serial() {
	string clusterfile = options.prefix + ".clusters";
	string discclfile = options.prefix + ".discord";
	if(!options.working_dir.empty()) {
		clusterfile = options.working_prefix + ".clusters";
		discclfile = options.working_prefix + ".discord";
	}
	string line;
	const char* delims = "\t";
	vector<string> data;
	string lastpid;
	string lastsid;
	int64_t k = 0;    // number of supporting entries
	int64_t l = 0;    // number of supporting read pairs
	int64_t j = 0;    // number of supporting full length read pairs
	map<string, int64_t> support;
	map<string, int64_t> support_f;
	set<string> read_names;
	vector<ClusterEntry> cluster;
	vector<string> patterns;
	patterns.push_back("mu1");
	patterns.push_back("mu2");
	patterns.push_back("sc");

	map<string, map<string, OrientationEntry>> orientation;
	map<string, map<string, pair<int64_t, int64_t>>> starts_map;
	map<string, map<string, pair<int64_t, int64_t>>> mstarts_map;
	map<string, map<string, vector<int64_t>>> cbps_map;
	map<string, map<string, map<int64_t, vector<string>>> > cluster_type;
	ofstream DISCCL(discclfile, ios::binary);
	ifstream CLMPD(clusterfile, ios::binary);

	while (getline(CLMPD, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delims, data);
		if (data.size() < 16) {
			continue;
		}

		ClusterEntry ce;
		ce.readname = data[5];
//const bool debug = "ST-E00104:502:HFJN5CCXX:1:1102:11525:23073" == ce.readname;
//const bool debug = "8" == data[0];
		ce.rg = data[6];
		ce.seqid = data[7];
		ce.strand = boost::lexical_cast<int32_t>(data[8]);
		ce.start = boost::lexical_cast<int64_t>(data[9]);
		ce.len = boost::lexical_cast<int64_t>(data[10]);
		ce.end = ce.start + ce.len - 1;
		ce.mseqid = data[11];
		ce.mstrand = boost::lexical_cast<int32_t>(data[12]);
		ce.mstart = boost::lexical_cast<int64_t>(data[13]);
		ce.mlen = boost::lexical_cast<int64_t>(data[14]);
		ce.mend = ce.mstart + ce.mlen - 1;
		ce.isize = boost::lexical_cast<int64_t>(data[15]);
		if (cluster.empty() || (lastpid == data[0] && lastsid == data[1])) {
			cluster.push_back(ce);
			if (read_names.end() == read_names.find(ce.readname)) {
				++k;
				read_names.insert(ce.readname);
			}
			size_t pos = castle::StringUtils::rfind_one_of_multiple_patterns(ce.readname, patterns);
			if (string::npos != pos) {
				string original_read_name = ce.readname.substr(0, pos);
//if (read_names.end() == read_names.find(original_read_name)) {
				++l;
//}
				read_names.insert(original_read_name);
//cout << original_read_name << "/" << ce.readname << "\n";
			} else {
				++j;
//if (read_names.end() == read_names.find(ce.readname)) {
				++l;
//}
				read_names.insert(ce.readname);
			}
		} else {
			support[lastpid] = l;
			support_f[lastpid] = j;
			determine_a_cluster(DISCCL, cluster, orientation, starts_map, mstarts_map, cbps_map, cluster_type, support, support_f, lastpid, lastsid,
					l, j);
			cluster.clear();
			read_names.clear();
			k = 0;
			l = 0;
			j = 0;

			cluster.push_back(ce);
			if (read_names.end() == read_names.find(ce.readname)) {
				++k;
				read_names.insert(ce.readname);
			}
			size_t pos = castle::StringUtils::rfind_one_of_multiple_patterns(ce.readname, patterns);
			if (string::npos != pos) {
				string original_read_name = ce.readname.substr(0, pos);
//if (read_names.end() == read_names.find(original_read_name)) {
				++l;
//}
				read_names.insert(original_read_name);
//cout << original_read_name << "/" << ce.readname << "\n";
			} else {
				++j;
//if (read_names.end() == read_names.find(ce.readname)) {
				++l;
//}
				read_names.insert(ce.readname);
			}
		}
		lastpid = data[0];
		lastsid = data[1];
	}
	{
// process the last cluster
		support[lastpid] = l;
		support_f[lastpid] = j;
		determine_a_cluster(DISCCL, cluster, orientation, starts_map, mstarts_map, cbps_map, cluster_type, support, support_f, lastpid, lastsid, l,
				j);
		cluster.clear();
		read_names.clear();
	}
}

void MatePairDiscordantCaller::call_initial_events() {
	vector<EventEntry> result_mp;    // store all results
//cout << "[MatePairDiscordantCaller.call_initial_events_serial] call intra-chromosomal events\n";
	call_intra_chromosomal_events(result_mp);
	call_inter_chromosomal_events(result_mp);
	call_rest_events(result_mp);
	call_smaller_events_among_intra_events(result_mp);
	call_smaller_events_among_inter_events(result_mp);
}

void MatePairDiscordantCaller::call_initial_events_serial() {
	vector<EventEntry> result_mp;    // store all results
	call_intra_chromosomal_events(result_mp);
	call_inter_chromosomal_events(result_mp);
	call_rest_events_alt(result_mp);
	call_smaller_events_among_intra_events_alt(result_mp);
	call_smaller_events_among_inter_events_alt(result_mp);
}


// this function consumes more amount of memory but faster if you have a multi-core machine
void MatePairDiscordantCaller::call_intra_chromosomal_events(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_intra_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	auto& ref_is = options.is;
	cout << "[MatePairDiscordantCaller.call_intra_chromosomal_events] collect candidate entries\n";

	vector<pair<string, string>> ii_lists_id3;
	vector<pair<string, string>> ii_lists_id1;
	tasks.push_back([&] {
		const char* underscore = "_";
		vector<string> cols;
		for (auto& pid_entry : cluster_type) {
			auto& chr = pid_entry.first;
			if (cluster_type[chr][chr].end() == cluster_type[chr][chr].find(0) || cluster_type[chr][chr].end() == cluster_type[chr][chr].find(3)) {
				continue;
			}
			for (auto& id3 : cluster_type[chr][chr][3]) {
				castle::StringUtils::c_string_multi_split(id3, underscore, cols);
				auto i = cols[0];
				auto ii = cols[1];
//				cout << "[MatePairDiscordantCaller.call_intra_chromosomal_events] 0: " << id3 << "\n";
				if (support[i] < options.support_mps || support_f[i] < options.support_mpf) {
					continue;
				}
//				cout << "[MatePairDiscordantCaller.call_intra_chromosomal_events] 1: " << id3 << "\n";
				if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
					continue;
				}
//				cout << "[MatePairDiscordantCaller.call_intra_chromosomal_events] 2: " << id3 << "\n";
				ii_lists_id3.push_back(make_pair(i, ii));
			}
		}
	});
	tasks.push_back([&] {
		const char* underscore = "_";
		vector<string> cols;
		for (auto& pid_entry : cluster_type) {
			auto& chr = pid_entry.first;
			if (cluster_type[chr][chr].end() == cluster_type[chr][chr].find(1) || cluster_type[chr][chr].end() == cluster_type[chr][chr].find(2)) {
				continue;
			}
			for (auto& id1 : cluster_type[chr][chr][1]) {
				castle::StringUtils::c_string_multi_split(id1, underscore, cols);
				auto i = cols[0];
				auto ii = cols[1];
				if (support[i] < options.support_mps || support_f[i] < options.support_mpf) {
					continue;
				}
				if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
					continue;
				}
				ii_lists_id1.push_back(make_pair(i, ii));
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
// each block looks at 10,000 entries
	int64_t n_total_entries_id3 = ii_lists_id3.size();
	cout << (boost::format("[MatePairDiscordantCaller.call_intra_chromosomal_events] start finding intra chromosomal events: id3(%d)\n") % n_total_entries_id3).str();
	int64_t block_size_id3 = n_total_entries_id3 / (double) n_cores;
	int64_t n_blocks_id3 = n_total_entries_id3 / (double) block_size_id3;
	if(n_total_entries_id3 < n_cores) {
		block_size_id3 = 1;
		n_blocks_id3 = n_total_entries_id3;
	}

	vector<vector<EventEntry>> result_mp_lists_id3(n_blocks_id3);
//IntervalClusterTree local_interval_tree_type_zero(the_intervals[0]);
	for (int64_t block_id = 0; block_id < n_blocks_id3; ++block_id) {
		tasks.push_back([&, block_id] {
			const int64_t n_total_entries_id3 = ii_lists_id3.size();
			int64_t block_size_id3 = n_total_entries_id3 / (double)n_cores;
			int64_t n_blocks_id3 = n_total_entries_id3 / (double)block_size_id3;
			if(n_total_entries_id3 < n_cores) {
				block_size_id3 = 1;
				n_blocks_id3 = n_total_entries_id3;
			}
			const int64_t start_id = castle::ParallelRunner::get_next_start_pos(block_id, block_size_id3, 0, 0);
			const int64_t end_id = castle::ParallelRunner::get_next_end_pos(block_id, n_blocks_id3 - 1, block_size_id3, n_total_entries_id3, 0);
			IntervalEntryVector results;
			results.reserve(10000);
// must use the copied version

			auto& local_result_mp = result_mp_lists_id3[block_id];
			IntervalEntryVector local_intervals = the_intervals[0];
			IntervalClusterTree local_interval_tree_type_zero(local_intervals);

			for(int64_t entry_id = start_id; entry_id < end_id; ++entry_id) {
				auto i = ii_lists_id3[entry_id].first;
				auto ii = ii_lists_id3[entry_id].second;
//const bool debug = "23903" == i;
				const bool debug = false;
				results.clear();
				map<int64_t, EventEntry> current_intra_event;
				int64_t_value_sortedset ca_event_size;
				int64_t cei = 0;// index for current event
//IntervalClusterTree* left = local_interval_tree_type_zero.left;
//IntervalClusterTree* right = local_interval_tree_type_zero.right;
//int64_t center = local_interval_tree_type_zero.center;

//local_interval_tree_type_zero.find_overlap_mt(starts_map[i][ii].first, starts_map[i][ii].second, center, results, local_intervals, left, right);
				if(starts_map.end() == starts_map.find(i) || starts_map[i].end() == starts_map[i].find(ii) ||
						mstarts_map.end() == mstarts_map.find(i) || mstarts_map[i].end() == mstarts_map[i].find(ii)) {
					continue;
				}
				local_interval_tree_type_zero.find_overlap(starts_map[i][ii].first, mstarts_map[i][ii].second, results);
				if(debug) {
					cout << "the start: " << starts_map[i][ii].first << "/" << starts_map[i][ii].second << "\n";
				}
//local_interval_tree_type_zero.find_point(starts_map[i][ii].first, results);
//local_interval_tree_type_zero.find_point(starts_map[i][ii].second, results);
//local_interval_tree_type_zero.find_point(mstarts_map[i][ii].first, results);
//local_interval_tree_type_zero.find_point(mstarts_map[i][ii].second, results);
//local_interval_tree_type_zero.find_overlap(starts_map[i][ii].first, starts_map[i][ii].second, results);
//local_interval_tree_type_zero.find_overlap(mstarts_map[i][ii].first, mstarts_map[i][ii].second, results);
				for(uint64_t item_id = 0; item_id < results.size(); ++item_id) {
					auto& an_entry = results[item_id].value;
					auto& k = an_entry.p_id;
					auto& kk = an_entry.s_id;
					if(debug) {
						cout << "first: " << k << "_" << kk << "\n";
					}
					if(support.end() == support.find(k) || support_f.end() == support_f.find(k) ||
							orientation.end() == orientation.find(i) || orientation[i].end() == orientation[i].find(ii) ||
							orientation.end() == orientation.find(k) || orientation[k].end() == orientation[k].find(kk)) {
						continue;
					}
					if (support[k] < options.support_mps
							|| support_f[k] < options.support_mpf || orientation[i][ii].ref_id != orientation[i][ii].mate_ref_id
							|| orientation[k][kk].ref_id != orientation[k][kk].mate_ref_id
							|| orientation[i][ii].ref_id != orientation[k][kk].ref_id) {
						continue;
					}
					if(starts_map.end() == starts_map.find(k) || starts_map[k].end() == starts_map[k].find(kk) ||
							mstarts_map.end() == mstarts_map.find(k) || mstarts_map[k].end() == mstarts_map[k].find(kk)) {
						continue;
					}
					if(cbps_map.end() == cbps_map.find(k) || cbps_map[k].end() == cbps_map[k].find(kk)) {
						continue;
					}
					if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1])
							|| (0 == cbps_map[k][kk][2]
									&& 0 == cbps_map[k][kk][3])) {
						continue;
					}
					if (i == k) {
						continue;
					}
//if(debug) {
//cout << "second: " << k << "_" << kk << "\n";
//}
// deletion with insertion in the same orientation
// the first two seqid ensures that the entry is an intra-chromosomal event.
// i: reverse-forward strand, k: forward-reverse strand
					if (-1 == orientation[i][ii].strand
							&& 1 == orientation[i][ii].mate_strand
							&& 1 == orientation[k][kk].strand
							&& -1 == orientation[k][kk].mate_strand
					) {
						if(debug) {
							cout << "confuse-3\n";
						}
// k is of a smaller start pos.
// starts_map.first: cluster head pos.
// starts_map.second: cluster tail pos.
// mstart_map.first: mate cluster head pos.
// mstart_map.second: mate cluster tail pos.
// i: reverse-forward strand, k: forward-reverse strand
						if (starts_map[k][kk].first < starts_map[i][ii].first &&
// i: cluster head-cluster tail (reverse reads), k: cluster tail-mate cluster tail(reverse-forward reads)
								covered(starts_map[i][ii].first, starts_map[i][ii].second,
										starts_map[k][kk].second, mstarts_map[k][kk].first) &&
// k: mate head-mate tail(reverse reads), i: cluster head-mate cluster tail(reverse-forward reads)
								covered(mstarts_map[k][kk].first, mstarts_map[k][kk].second,
										starts_map[i][ii].first, mstarts_map[i][ii].second)
						) {
							int64_t del_size = starts_map[i][ii].first
							- starts_map[k][kk].second - 1;
							int64_t ins_size = mstarts_map[i][ii].second
							- mstarts_map[k][kk].first + 1;
							int64_t distance = mstarts_map[k][kk].first
							- starts_map[i][ii].first;
							int64_t event_size = abs(del_size) + abs(ins_size);
							if(debug) {
								cout << "confuse-4\n";
							}

							if (del_size < options.sv_size_cutoff
									&& ins_size < options.sv_size_cutoff && 0 != event_size) {
//if(debug) {
//cout << "here-3\n";
//}
								if (del_size < ref_is["rlu"]["selected"]) {
									current_intra_event[cei].type = "inssd";
								} else {
									current_intra_event[cei].type = "del_inssd";
								}
								current_intra_event[cei].cluster_id = i + "_"
								+ ii;
								current_intra_event[cei].mate_cluster_id = k + "_"
								+ kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id =
								orientation[i][ii].ref_id;
								current_intra_event[cei].event_start =
								starts_map[k][kk].second;
								current_intra_event[cei].event_end =
								starts_map[i][ii].first;
								current_intra_event[cei].event_size_1 =
								del_size;

								current_intra_event[cei].mate_ref_id =
								orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start =
								mstarts_map[k][kk].first;
								current_intra_event[cei].mate_event_end =
								mstarts_map[i][ii].second;
								current_intra_event[cei].event_size_2 =
								ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] =
								cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] =
								cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] =
								cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] =
								cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] =
								cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] =
								cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] =
								cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] =
								cbps_map[k][kk][3];
								Int64Pair a_pair(cei, event_size);
								ca_event_size.insert(a_pair);
//ca_event_size[cei] = event_size;
								++cei;
							}
						}
// k is of a smaller start pos.
// starts_map.first: cluster head pos.
// starts_map.second: cluster tail pos.
// mstart_map.first: mate cluster head pos.
// mstart_map.second: mate cluster tail pos.

// i: reverse-forward strand, k: forward-reverse strand
						if (mstarts_map[i][ii].second < mstarts_map[k][kk].second &&
// k: cluster head-cluster tail (forward reads), i: cluster tail-mate cluster tail(forward-reverse reads)
								covered(starts_map[k][kk].first, starts_map[k][kk].second,
										starts_map[i][ii].first, mstarts_map[i][ii].second)
// i: mate cluster head-mate cluster tail (forward reads), k: cluster tail-mate cluster head(forward-reverse reads)
								&& covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second,
										starts_map[k][kk].second, mstarts_map[k][kk].first)
						) {
							int64_t del_size = mstarts_map[k][kk].first
							- mstarts_map[i][ii].second - 1;
							int64_t ins_size = starts_map[k][kk].second
							- starts_map[i][ii].first + 1;
							int64_t distance = mstarts_map[i][ii].second
							- starts_map[k][kk].second;
							int64_t event_size = abs(del_size) + abs(ins_size);
							if(debug) {
								cout << "confuse-5\n";
							}

							if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff
									&& 0 != event_size) {
								if(debug) {
									cout << "confuse-6: " << event_size << "\n";
								}
								if (del_size < ref_is["rlu"]["selected"]) {
									current_intra_event[cei].type = "inssu";
								} else {
									current_intra_event[cei].type = "del_inssu";
								}
								current_intra_event[cei].cluster_id = i + "_"
								+ ii;
								current_intra_event[cei].mate_cluster_id = k + "_"
								+ kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id =
								orientation[i][ii].ref_id;
								current_intra_event[cei].event_start =
								mstarts_map[i][ii].second;
								current_intra_event[cei].event_end =
								mstarts_map[k][kk].first;
								current_intra_event[cei].event_size_1 =
								del_size;

								current_intra_event[cei].mate_ref_id =
								orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start =
								starts_map[i][ii].first;
								current_intra_event[cei].mate_event_end =
								starts_map[k][kk].second;
								current_intra_event[cei].event_size_2 =
								ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] =
								cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] =
								cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] =
								cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] =
								cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] =
								cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] =
								cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] =
								cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] =
								cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
//if(debug) {
//cout << current_intra_event[cei].str() << "\n";
//}
								Int64Pair a_pair(cei, event_size);
								ca_event_size.insert(a_pair);
								++cei;
							}
						}
//if(debug) {
//cout << "here-6\n";
//}
					}
//if(debug) {
//cout << "here-7\n";
//exit(1);
//}
				}
// call smallest event for cluster i
				int64_t top1 = -1;
				int64_t top2 = -1;
				int64_t top3 = -1;
				int64_t chi = 0;
// ca_event_size is a set sorted by value, thus the key is already sorted.
				for (auto cei : ca_event_size) {
					if(debug) {
						cout << "confuse event size: " << cei.first << ":" << cei.second << "\n";
					}
					if (0 == chi) {
						top1 = cei.first;
					} else if (1 == chi) {
						top2 = cei.first;
					} else if (2 == chi) {
						top3 = cei.first;
						break;
					}
					++chi;
				}
				if(-1 != top1) {
//TODO: the following code is adopted from the original meerkat.
// However, there are no inverse, insod, and insoud calling in the above codes.
// hence it should be more simplified.
					if ((string::npos
									!= current_intra_event[top1].type.find("insod")
									|| (string::npos
											!= current_intra_event[top1].type.find("insou")))
							&& ("invers" == current_intra_event[top2].type
									|| "invers" == current_intra_event[top3].type)) {
						if (current_intra_event[top1].event_size_1 < 20
								&& current_intra_event[top1].event_size_2 < 20
								&& current_intra_event[top1].distance > 20) {
							if (-1 != top2 && "invers" == current_intra_event[top2].type) {
//if(current_intra_event[top2].type.empty()) {
//cout << (boost::format("top2-1: %s (%s)\n") % current_intra_event[top2].str() % i).str();
//}
								local_result_mp.push_back(current_intra_event[top2]);
							} else if (-1 != top3 && "invers" == current_intra_event[top3].type) {
//if(current_intra_event[top3].type.empty()) {
//cout << (boost::format("top3-2: %s (%s)\n") % current_intra_event[top3].str() % i).str();
//}
								local_result_mp.push_back(current_intra_event[top3]);
							}
						} else {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-3: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
							local_result_mp.push_back(current_intra_event[top1]);
						}
					} else {
						if (current_intra_event.end()
								!= current_intra_event.find(top1)) {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-4: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
							local_result_mp.push_back(current_intra_event[top1]);
						}
					}
				}
			}
		});
	}
//castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	int64_t n_total_entries_id1 = ii_lists_id1.size();
	cout << (boost::format("[MatePairDiscordantCaller.call_intra_chromosomal_events] start finding intra chromosomal events: id1(%d)\n") % n_total_entries_id1).str();
	int64_t block_size_id1 = n_total_entries_id1 / (double) n_cores;
	int64_t n_blocks_id1 = n_total_entries_id1 / (double) block_size_id1;
	if(n_total_entries_id1 < n_cores) {
		block_size_id1 = 1;
		n_blocks_id1 = n_total_entries_id1;
	}
	vector<vector<EventEntry>> result_mp_lists_id1(n_blocks_id1);
//cout << n_blocks_id1 << "/" << n_blocks_id3 << "\n";
//IntervalClusterTree local_interval_tree_type_two(the_intervals[2]);
	for (int64_t block_id = 0; block_id < n_blocks_id1; ++block_id) {
		tasks.push_back([&, block_id] {
			const int64_t n_total_entries_id1 = ii_lists_id1.size();
			int64_t block_size_id1 = n_total_entries_id1 / (double)n_cores;
			int64_t n_blocks_id1 = n_total_entries_id1 / (double)block_size_id1;
			if(n_total_entries_id1 < n_cores) {
				block_size_id1 = 1;
				n_blocks_id1 = n_total_entries_id1;
			}
			const int64_t start_id = castle::ParallelRunner::get_next_start_pos(block_id, block_size_id1, 0, 0);
			const int64_t end_id = castle::ParallelRunner::get_next_end_pos(block_id, n_blocks_id1 - 1, block_size_id1, n_total_entries_id1, 0);
			IntervalEntryVector results;
			results.reserve(10000);
			auto& local_result_mp = result_mp_lists_id1[block_id];
// should use a copied version
				IntervalEntryVector local_interval = the_intervals[2];
				IntervalClusterTree local_interval_tree_type_two(local_interval);
//IntervalEntryVector copied = local_interval_tree_type_two.intervals;
//IntervalEntryVector local_intervals = local_interval_tree_type_two.intervals;

				for(int64_t entry_id = start_id; entry_id < end_id; ++entry_id) {
					auto i = ii_lists_id1[entry_id].first;
					auto ii = ii_lists_id1[entry_id].second;
					results.clear();
//					const bool debug = "71237" == i;
					const bool debug = false;
					map<int64_t, EventEntry> current_intra_event;
					int64_t_value_sortedset ca_event_size;
					int64_t cei = 0;// index for current event
//IntervalClusterTree* left = local_interval_tree_type_two.left;
//IntervalClusterTree* right = local_interval_tree_type_two.right;
//int64_t center = local_interval_tree_type_zero.center;
//local_interval_tree_type_two.find_overlap_mt(starts_map[i][ii].first, starts_map[i][ii].second, center, results, local_intervals, left, right);
//local_interval_tree_type_two.find_overlap(starts_map[i][ii].first, starts_map[i][ii].second, results);
					local_interval_tree_type_two.find_overlap(starts_map[i][ii].first, mstarts_map[i][ii].second, results);
//local_interval_tree_type_two.find_overlap(mstarts_map[i][ii].first, mstarts_map[i][ii].second, results);
					for(uint64_t item_id = 0; item_id < results.size(); ++item_id) {
						auto& an_entry = results[item_id].value;
						auto& k = an_entry.p_id;
						auto& kk = an_entry.s_id;
						if(debug) {
							cout << "inv-1: " << k << "_" << kk << "\n";
						}
						if(support.end() == support.find(k) || support_f.end() == support_f.find(k) ||
								orientation.end() == orientation.find(i) || orientation[i].end() == orientation[i].find(ii) ||
								orientation.end() == orientation.find(k) || orientation[k].end() == orientation[k].find(kk)) {
							continue;
						}
						if (support[k] < options.support_mps
								|| support_f[k] < options.support_mpf
								|| orientation[i][ii].ref_id != orientation[i][ii].mate_ref_id
								|| orientation[k][kk].ref_id != orientation[k][kk].mate_ref_id
								|| orientation[i][ii].ref_id != orientation[k][kk].ref_id) {
							continue;
						}
						if(cbps_map.end() == cbps_map.find(k) || cbps_map[k].end() == cbps_map[k].find(kk)) {
							continue;
						}
						if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1])
								|| (0 == cbps_map[k][kk][2]
										&& 0 == cbps_map[k][kk][3])) {
							continue;
						}
						if (i == k) {
							continue;
						}
						if(debug) {
							cout << "inv-2\n";
						}
						if(starts_map.end() == starts_map.find(k) || starts_map[k].end() == starts_map[k].find(kk) ||
								mstarts_map.end() == mstarts_map.find(k) || mstarts_map[k].end() == mstarts_map[k].find(kk)) {
							continue;
						}
//# inversions
						if (starts_map[i][ii].first < starts_map[k][kk].first
								&& starts_map[i][ii].second < starts_map[k][kk].second
								&& mstarts_map[i][ii].first < mstarts_map[k][kk].first
								&& mstarts_map[i][ii].second < mstarts_map[k][kk].second &&
								covered(starts_map[k][kk].first, starts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second)
								&& covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first)) {
//#print "14: inversion\n";
							int64_t inv_size = mstarts_map[i][ii].second - starts_map[k][kk].first;
							if(debug) {
								cout << "inv-3: " << inv_size << "\n";
							}
							if(inv_size < options.sv_size_cutoff ) {
//#print "15: inversion\n";
								if(covered(starts_map[i][ii].second - 10, starts_map[i][ii].second + ref_is["rlu"]["selected"],
												starts_map[k][kk].first - ref_is["rlu"]["selected"], starts_map[k][kk].first + 10) &&
										covered(mstarts_map[i][ii].second - 10, mstarts_map[i][ii].second + ref_is["rlu"]["selected"],
												mstarts_map[k][kk].first - ref_is["rlu"]["selected"], mstarts_map[k][kk].first + 10)) {
									if(debug) {
										cout << "inv-4\n";
									}
									current_intra_event[cei].type = "invers";
								} else {
									if(debug) {
										cout << "inv-5\n";
									}
									current_intra_event[cei].type = "del_invers";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = starts_map[i][ii].second;
								current_intra_event[cei].event_end = starts_map[k][kk].first;
								current_intra_event[cei].event_size_1 = inv_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = mstarts_map[i][ii].second;
								current_intra_event[cei].mate_event_end = mstarts_map[k][kk].first;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
								Int64Pair a_pair(cei, abs(inv_size));
								ca_event_size.insert(a_pair);
								++cei;
							}
						}

//# deletion with insertion in opposite orientation
						if (starts_map[i][ii].first < starts_map[k][kk].first &&
								covered(starts_map[k][kk].first, starts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second) &&
								covered(mstarts_map[k][kk].first, mstarts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second)) {
//#print "16: del with insertion opposite direction\n";
							int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
							int64_t ins_size = mstarts_map[i][ii].second - mstarts_map[k][kk].first + 1;
							int64_t distance = mstarts_map[k][kk].first - starts_map[k][kk].first;
							if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
//#print "17: del with insertion opposite direction\n";
								if (del_size < ref_is["rlu"]["selected"]) {
									current_intra_event[cei].type = "insod";
								} else {
									current_intra_event[cei].type = "del_insod";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = starts_map[i][ii].second;
								current_intra_event[cei].event_end = starts_map[k][kk].first;
								current_intra_event[cei].event_size_1 = del_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = mstarts_map[k][kk].first;
								current_intra_event[cei].mate_event_end = mstarts_map[i][ii].second;
								current_intra_event[cei].event_size_2 = ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
								Int64Pair a_pair(cei, abs(del_size) + abs(ins_size));
								ca_event_size.insert(a_pair);
								++cei;
							}
						}

//# deletion with insertion in opposite orientation
						if (mstarts_map[i][ii].second < mstarts_map[k][kk].second &&
								covered(starts_map[i][ii].first, starts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first) &&
								covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first)) {
//#print "18: del with insertion opposite direction\n";
							int64_t del_size = mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1;
							int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
							int64_t distance = mstarts_map[i][ii].second - starts_map[i][ii].second;
							if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
								if (del_size < ref_is["rlu"]["selected"] ) {
									current_intra_event[cei].type = "insou";
								} else {
									current_intra_event[cei].type = "del_insou";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = mstarts_map[i][ii].second;
								current_intra_event[cei].event_end = mstarts_map[k][kk].first;
								current_intra_event[cei].event_size_1 = del_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = starts_map[k][kk].first;
								current_intra_event[cei].mate_event_end = starts_map[i][ii].second;
								current_intra_event[cei].event_size_2 = ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
								Int64Pair a_pair(cei, abs(del_size) + abs(ins_size));
								ca_event_size.insert(a_pair);
								++cei;
							}
						}
					}
// call smallest event for cluster i
					int64_t top1 = -1;
					int64_t top2 = -1;
					int64_t top3 = -1;
					int64_t chi = 0;
// ca_event_size is a set sorted by value, thus the key is already sorted.
					for (auto cei : ca_event_size) {
						if (debug) {
							cout << "inv_size: " << cei.first << ":" << cei.second << "\n";
						}
						if (0 == chi) {
							top1 = cei.first;
						} else if (1 == chi) {
							top2 = cei.first;
						} else if (2 == chi) {
							top3 = cei.first;
							break;
						}
						++chi;
					}
					if(-1 != top1) {
						if ((string::npos != current_intra_event[top1].type.find("insod")
										|| (string::npos != current_intra_event[top1].type.find("insou")))
								&& ("invers" == current_intra_event[top2].type
										|| "invers" == current_intra_event[top3].type)) {
							if (current_intra_event[top1].event_size_1 < 20
									&& current_intra_event[top1].event_size_2 < 20
									&& current_intra_event[top1].distance > 20) {
								if (-1 != top2 && "invers" == current_intra_event[top2].type) {
//if(current_intra_event[top2].type.empty()) {
//cout << (boost::format("top2-5: %s (%s)\n") % current_intra_event[top2].str() % i).str();
//}
									local_result_mp.push_back(current_intra_event[top2]);
								} else if (-1 != top3 && "invers" == current_intra_event[top3].type) {
//if(current_intra_event[top3].type.empty()) {
//cout << (boost::format("top3-6: %s (%s)\n") % current_intra_event[top3].str() % i).str();
//}
									local_result_mp.push_back(current_intra_event[top3]);
								}
							} else {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top3-7: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
								local_result_mp.push_back(current_intra_event[top1]);
							}
						} else {
							if (current_intra_event.end()
									!= current_intra_event.find(top1)) {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-8: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
								local_result_mp.push_back(current_intra_event[top1]);
							}
						}
					}
				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (int64_t block_id = 0; block_id < n_blocks_id3; ++block_id) {
		auto& local_result_mp = result_mp_lists_id3[block_id];
		result_mp.insert(result_mp.end(), local_result_mp.begin(), local_result_mp.end());
//for(auto& an_entry : local_result_mp) {
//cout << an_entry.str() << "\n";
//}
	}
	for (int64_t block_id = 0; block_id < n_blocks_id1; ++block_id) {
		auto& local_result_mp = result_mp_lists_id1[block_id];
		result_mp.insert(result_mp.end(), local_result_mp.begin(), local_result_mp.end());
//for(auto& an_entry : local_result_mp) {
//cout << an_entry.str() << "\n";
//}
	}
	cout << checker;
}

void MatePairDiscordantCaller::call_intra_chromosomal_events_alt(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_intra_chromosomal_events_alt");
	checker.start();
	vector<function<void()> > tasks;
	auto& ref_is = options.is;
	cout << "[MatePairDiscordantCaller.call_intra_chromosomal_events_alt] collect candidate entries\n";
	const char* delim_underscore = "_";
	vector<string> cols;
//# intra-chr events
	IntervalClusterTree local_interval_tree_type_zero(the_intervals[0]);
	IntervalClusterTree local_interval_tree_type_two(the_intervals[2]);
	IntervalEntryVector results;
	results.reserve(10000);

	for(auto& chr_entry : cluster_type) {
		auto& chr = chr_entry.first;
		auto& cl_map = cluster_type[chr][chr];
		if (cl_map.end() != cl_map.find(0) && cl_map.end() != cl_map.find(3)){
			auto id3_list_itr = cl_map.find(3);
			auto& id3_list = id3_list_itr->second;
			for(auto& id3 : id3_list) {
				castle::StringUtils::c_string_multi_split(id3, delim_underscore, cols);
				auto i = cols[0];
				auto ii = cols[1];
				if(support[i] < options.support_mps || support_f[i] < options.support_mpf) {
					continue;
				}

//				#print "$chr\t$id3\n";
				if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
					continue;
				}
//				     # @current_intra_event: all possible events from cluster i
				map<int64_t, EventEntry> current_intra_event;
//				     # %ca_event_size: all event sizes for @current_intra_event
				int64_t_value_sortedset ca_event_size;
				//# index for current event
				int64_t cei = 0;

				results.clear();
				local_interval_tree_type_zero.find_overlap(starts_map[i][ii].first, mstarts_map[i][ii].second, results);
				for(uint64_t item_id = 0; item_id < results.size(); ++item_id) {
					auto& an_entry = results[item_id].value;
					auto& k = an_entry.p_id;
					auto& kk = an_entry.s_id;
					if (support[k] < options.support_mps
							|| support_f[k] < options.support_mpf || orientation[i][ii].ref_id != orientation[i][ii].mate_ref_id
							|| orientation[k][kk].ref_id != orientation[k][kk].mate_ref_id
							|| orientation[i][ii].ref_id != orientation[k][kk].ref_id) {
						continue;
					}
					if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1])
							|| (0 == cbps_map[k][kk][2]
									&& 0 == cbps_map[k][kk][3])) {
						continue;
					}
					if (i == k) {
						continue;
					}
//					# deletion with insertion in the same orientation
					if (orientation[i][ii].ref_id == orientation[i][ii].mate_ref_id
							&& orientation[k][kk].ref_id == orientation[k][kk].mate_ref_id
							&& -1 == orientation[i][ii].strand
							&& 1 == orientation[i][ii].mate_strand
							&& 1 == orientation[k][kk].strand
							&& -1 == orientation[k][kk].mate_strand
					) {

						if (covered(starts_map[i][ii].first, starts_map[i][ii].second, starts_map[k][kk].second, mstarts_map[k][kk].first)
							&& covered(mstarts_map[k][kk].first, mstarts_map[k][kk].second, starts_map[i][ii].first,  mstarts_map[i][ii].second)
							&& starts_map[k][kk].first < starts_map[i][ii].first) {
							int64_t del_size = starts_map[i][ii].first - starts_map[k][kk].second - 1;
							int64_t ins_size = mstarts_map[i][ii].second - mstarts_map[k][kk].first + 1;
							int64_t distance = mstarts_map[k][kk].first - starts_map[i][ii].first;
							if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
								if (del_size < ref_is["rlu"]["selected"] ) {
									current_intra_event[cei].type = "inssd";
								} else {
									current_intra_event[cei].type = "del_inssd";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];
								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = starts_map[k][kk].second;
								current_intra_event[cei].event_end = starts_map[i][ii].first;
								current_intra_event[cei].event_size_1 = del_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = mstarts_map[k][kk].first;
								current_intra_event[cei].mate_event_end = mstarts_map[i][ii].second;
								current_intra_event[cei].event_size_2 = ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
								int64_t event_size = abs(del_size) + abs(ins_size);
								Int64Pair a_pair(cei, event_size);
								ca_event_size.insert(a_pair);
//ca_event_size[cei] = event_size;
								++cei;
							}
						}

						if (covered(starts_map[k][kk].first, starts_map[k][kk].second, starts_map[i][ii].first, mstarts_map[i][ii].second) &&
								covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second,
								starts_map[k][kk].second, mstarts_map[k][kk].first)
							&& mstarts_map[i][ii].second < mstarts_map[k][kk].second) {
							int64_t del_size = mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1;
							int64_t ins_size = starts_map[k][kk].second - starts_map[i][ii].first + 1;
							int64_t distance = mstarts_map[i][ii].second - starts_map[k][kk].second;
							if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
								if (del_size < ref_is["rlu"]["selected"] ) {
									current_intra_event[cei].type = "inssu";
								} else {
									current_intra_event[cei].type = "del_inssu";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = mstarts_map[i][ii].second;
								current_intra_event[cei].event_end = mstarts_map[k][kk].first;
								current_intra_event[cei].event_size_1 = del_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = starts_map[i][ii].first;
								current_intra_event[cei].mate_event_end = starts_map[k][kk].second;
								current_intra_event[cei].event_size_2 = ins_size;
								current_intra_event[cei].distance = distance;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
//if(debug) {
//cout << current_intra_event[cei].str() << "\n";
//}
								int64_t event_size = abs(del_size) + abs(ins_size);
								Int64Pair a_pair(cei, event_size);
								ca_event_size.insert(a_pair);
								++cei;
							}
						}
					}
				}

//				# call smallest event for cluster i
				int64_t top1 = -1;
				int64_t top2 = -1;
				int64_t top3 = -1;
				int64_t chi = 0;
// ca_event_size is a set sorted by value, thus the key is already sorted.
				for (auto cei : ca_event_size) {
					if (0 == chi) {
						top1 = cei.first;
					} else if (1 == chi) {
						top2 = cei.first;
					} else if (2 == chi) {
						top3 = cei.first;
						break;
					}
					++chi;
				}
				if(-1 != top1) {
//TODO: the following code is adopted from the original meerkat.
// However, there are no inverse, insod, and insoud calling in the above codes.
// hence it should be more simplified.
					if ((string::npos != current_intra_event[top1].type.find("insod") ||
							(string::npos != current_intra_event[top1].type.find("insou")))
							&& ("invers" == current_intra_event[top2].type || "invers" == current_intra_event[top3].type)) {
						if (current_intra_event[top1].event_size_1 < 20
								&& current_intra_event[top1].event_size_2 < 20
								&& current_intra_event[top1].distance > 20) {
							if (-1 != top2 && "invers" == current_intra_event[top2].type) {
//if(current_intra_event[top2].type.empty()) {
//cout << (boost::format("top2-1: %s (%s)\n") % current_intra_event[top2].str() % i).str();
//}
								result_mp.push_back(current_intra_event[top2]);
							} else if (-1 != top3 && "invers" == current_intra_event[top3].type) {
//if(current_intra_event[top3].type.empty()) {
//cout << (boost::format("top3-2: %s (%s)\n") % current_intra_event[top3].str() % i).str();
//}
								result_mp.push_back(current_intra_event[top3]);
							}
						} else {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-3: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
							result_mp.push_back(current_intra_event[top1]);
						}
					} else {
						if (current_intra_event.end()
								!= current_intra_event.find(top1)) {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-4: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
							result_mp.push_back(current_intra_event[top1]);
						}
					}
				}
			}
		}

		if (cl_map.end() != cl_map.find(1) && cl_map.end() != cl_map.find(2)){
			auto id1_list_itr = cl_map.find(1);
			auto& id1_list = id1_list_itr->second;
			for(auto& id1 : id1_list) {
				castle::StringUtils::c_string_multi_split(id1, delim_underscore, cols);
				auto i = cols[0];
				auto ii = cols[1];
				if (support[i] < options.support_mps || support_f[i] < options.support_mpf) {
					continue;
				}
				if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
					continue;
				}
//				# @current_intra_event: all possible events from cluster i
				map<int64_t, EventEntry> current_intra_event;
//				# %ca_event_size: all event sizes for @current_intra_event
				int64_t_value_sortedset ca_event_size;
//				# index for current event
				int64_t cei = 0;
				results.clear();
				local_interval_tree_type_two.find_overlap(starts_map[i][ii].first, mstarts_map[i][ii].second, results);
				for(uint64_t item_id = 0; item_id < results.size(); ++item_id) {
					auto& an_entry = results[item_id].value;
					auto& k = an_entry.p_id;
					auto& kk = an_entry.s_id;
					if (support[k] < options.support_mps
							|| support_f[k] < options.support_mpf
							|| orientation[i][ii].ref_id != orientation[i][ii].mate_ref_id
							|| orientation[k][kk].ref_id != orientation[k][kk].mate_ref_id
							|| orientation[i][ii].ref_id != orientation[k][kk].ref_id) {
						continue;
					}
					if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1])
							|| (0 == cbps_map[k][kk][2] && 0 == cbps_map[k][kk][3])) {
						continue;
					}
					if (i == k) {
						continue;
					}
//					# inversions
					if (starts_map[i][ii].first < starts_map[k][kk].first && starts_map[i][ii].second < starts_map[k][kk].second
						&& mstarts_map[i][ii].first < mstarts_map[k][kk].first && mstarts_map[i][ii].second < mstarts_map[k][kk].second &&
						covered(starts_map[k][kk].first, starts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second)
						&& covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first)) {
//#print "14: inversion\n";
							int64_t inv_size = mstarts_map[i][ii].second - starts_map[k][kk].first;
							if(inv_size < options.sv_size_cutoff ) {
//#print "15: inversion\n";
								if(covered(starts_map[i][ii].second - 10, starts_map[i][ii].second + ref_is["rlu"]["selected"],
												starts_map[k][kk].first - ref_is["rlu"]["selected"], starts_map[k][kk].first + 10) &&
										covered(mstarts_map[i][ii].second - 10, mstarts_map[i][ii].second + ref_is["rlu"]["selected"],
												mstarts_map[k][kk].first - ref_is["rlu"]["selected"], mstarts_map[k][kk].first + 10)) {
									current_intra_event[cei].type = "invers";
								} else {
									current_intra_event[cei].type = "del_invers";
								}
								current_intra_event[cei].cluster_id = i + "_" + ii;
								current_intra_event[cei].mate_cluster_id = k + "_" + kk;
								current_intra_event[cei].n_supports = support[i];
								current_intra_event[cei].n_mate_support = support[k];

								current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].event_start = starts_map[i][ii].second;
								current_intra_event[cei].event_end = starts_map[k][kk].first;
								current_intra_event[cei].event_size_1 = inv_size;

								current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
								current_intra_event[cei].mate_event_start = mstarts_map[i][ii].second;
								current_intra_event[cei].mate_event_end = mstarts_map[k][kk].first;

								current_intra_event[cei].cbps.resize(8);
								current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
								current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
								current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
								current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
								current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
								current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
								current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
								current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
								Int64Pair a_pair(cei, abs(inv_size));
								ca_event_size.insert(a_pair);
								++cei;
							}
						}

					//# deletion with insertion in opposite orientation
					if (starts_map[i][ii].first < starts_map[k][kk].first &&
							covered(starts_map[k][kk].first, starts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second) &&
							covered(mstarts_map[k][kk].first, mstarts_map[k][kk].second, starts_map[i][ii].second, mstarts_map[i][ii].second)) {
//#print "16: del with insertion opposite direction\n";
						int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
						int64_t ins_size = mstarts_map[i][ii].second - mstarts_map[k][kk].first + 1;
						int64_t distance = mstarts_map[k][kk].first - starts_map[k][kk].first;
						if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
//#print "17: del with insertion opposite direction\n";
							if (del_size < ref_is["rlu"]["selected"]) {
								current_intra_event[cei].type = "insod";
							} else {
								current_intra_event[cei].type = "del_insod";
							}
							current_intra_event[cei].cluster_id = i + "_" + ii;
							current_intra_event[cei].mate_cluster_id = k + "_" + kk;
							current_intra_event[cei].n_supports = support[i];
							current_intra_event[cei].n_mate_support = support[k];

							current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
							current_intra_event[cei].event_start = starts_map[i][ii].second;
							current_intra_event[cei].event_end = starts_map[k][kk].first;
							current_intra_event[cei].event_size_1 = del_size;

							current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
							current_intra_event[cei].mate_event_start = mstarts_map[k][kk].first;
							current_intra_event[cei].mate_event_end = mstarts_map[i][ii].second;
							current_intra_event[cei].event_size_2 = ins_size;
							current_intra_event[cei].distance = distance;

							current_intra_event[cei].cbps.resize(8);
							current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
							current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
							current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
							current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
							current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
							current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
							current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
							current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
//ca_event_size[cei] = event_size;
							Int64Pair a_pair(cei, abs(del_size) + abs(ins_size));
							ca_event_size.insert(a_pair);
							++cei;
						}
					}

					//# deletion with insertion in opposite orientation
					if (mstarts_map[i][ii].second < mstarts_map[k][kk].second &&
							covered(starts_map[i][ii].first, starts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first) &&
							covered(mstarts_map[i][ii].first, mstarts_map[i][ii].second, starts_map[k][kk].first, mstarts_map[k][kk].first)) {
//#print "18: del with insertion opposite direction\n";
						int64_t del_size = mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1;
						int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
						int64_t distance = mstarts_map[i][ii].second - starts_map[i][ii].second;
						if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
							if (del_size < ref_is["rlu"]["selected"] ) {
								current_intra_event[cei].type = "insou";
							} else {
								current_intra_event[cei].type = "del_insou";
							}
							current_intra_event[cei].cluster_id = i + "_" + ii;
							current_intra_event[cei].mate_cluster_id = k + "_" + kk;
							current_intra_event[cei].n_supports = support[i];
							current_intra_event[cei].n_mate_support = support[k];

							current_intra_event[cei].ref_id = orientation[i][ii].ref_id;
							current_intra_event[cei].event_start = mstarts_map[i][ii].second;
							current_intra_event[cei].event_end = mstarts_map[k][kk].first;
							current_intra_event[cei].event_size_1 = del_size;

							current_intra_event[cei].mate_ref_id = orientation[i][ii].ref_id;
							current_intra_event[cei].mate_event_start = starts_map[k][kk].first;
							current_intra_event[cei].mate_event_end = starts_map[i][ii].second;
							current_intra_event[cei].event_size_2 = ins_size;
							current_intra_event[cei].distance = distance;

							current_intra_event[cei].cbps.resize(8);
							current_intra_event[cei].cbps[0] = cbps_map[i][ii][0];
							current_intra_event[cei].cbps[1] = cbps_map[i][ii][1];
							current_intra_event[cei].cbps[2] = cbps_map[i][ii][2];
							current_intra_event[cei].cbps[3] = cbps_map[i][ii][3];
							current_intra_event[cei].cbps[4] = cbps_map[k][kk][0];
							current_intra_event[cei].cbps[5] = cbps_map[k][kk][1];
							current_intra_event[cei].cbps[6] = cbps_map[k][kk][2];
							current_intra_event[cei].cbps[7] = cbps_map[k][kk][3];
							Int64Pair a_pair(cei, abs(del_size) + abs(ins_size));
							ca_event_size.insert(a_pair);
							++cei;
						}
					}

					// call smallest event for cluster i
					int64_t top1 = -1;
					int64_t top2 = -1;
					int64_t top3 = -1;
					int64_t chi = 0;
// ca_event_size is a set sorted by value, thus the key is already sorted.
					for (auto cei : ca_event_size) {
						if (0 == chi) {
							top1 = cei.first;
						} else if (1 == chi) {
							top2 = cei.first;
						} else if (2 == chi) {
							top3 = cei.first;
							break;
						}
						++chi;
					}
					if(-1 != top1) {
						if ((string::npos != current_intra_event[top1].type.find("insod")
										|| (string::npos != current_intra_event[top1].type.find("insou")))
								&& ("invers" == current_intra_event[top2].type
										|| "invers" == current_intra_event[top3].type)) {
							if (current_intra_event[top1].event_size_1 < 20
									&& current_intra_event[top1].event_size_2 < 20
									&& current_intra_event[top1].distance > 20) {
								if (-1 != top2 && "invers" == current_intra_event[top2].type) {
//if(current_intra_event[top2].type.empty()) {
//cout << (boost::format("top2-5: %s (%s)\n") % current_intra_event[top2].str() % i).str();
//}
									result_mp.push_back(current_intra_event[top2]);
								} else if (-1 != top3 && "invers" == current_intra_event[top3].type) {
//if(current_intra_event[top3].type.empty()) {
//cout << (boost::format("top3-6: %s (%s)\n") % current_intra_event[top3].str() % i).str();
//}
									result_mp.push_back(current_intra_event[top3]);
								}
							} else {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top3-7: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
								result_mp.push_back(current_intra_event[top1]);
							}
						} else {
							if (current_intra_event.end()
									!= current_intra_event.find(top1)) {
//if(current_intra_event[top1].type.empty()) {
//cout << (boost::format("top1-8: %s (%s)\n") % current_intra_event[top1].str() % i).str();
//}
								result_mp.push_back(current_intra_event[top1]);
							}
						}
					}
				}
			}
		}

	}

	cout << checker;
}

void MatePairDiscordantCaller::call_inter_chromosomal_events(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_inter_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	cout << "[MatePairDiscordantCaller.call_inter_chromosomal_events] collect candidate entries\n";
	auto& ref_is = options.is;
	const char* underscore = "_";
	vector<string> cols;
	vector<pair<string, string>> p_ids;
	for (auto& chr1_entry : cluster_type) {
		auto& chr1 = chr1_entry.first;
		for (auto& chr2_entry : chr1_entry.second) {
			auto& chr2 = chr2_entry.first;
			if (chr1 == chr2) {
				continue;
			}
			for (auto& the_type_entry : chr2_entry.second) {
//auto the_type = the_type_entry.first;
				for (auto& id : the_type_entry.second) {
					castle::StringUtils::c_string_multi_split(id, underscore, cols);
					if(cols.size() >= 2) {
						auto& i = cols[0];
						auto& ii = cols[1];
						bool is_supporting = support[i] < options.support_mps || support_f[i] < options.support_mpf;
						if(cbps_map[i][ii].size() < 4) {
							cbps_map[i][ii].resize(4);
						}
						bool is_valid_cbps = (0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3]);
						if (is_supporting) {
							continue;
						}
						if (is_valid_cbps) {
							continue;
						}
						p_ids.push_back(make_pair(i, ii));
					} else if(1 == cols.size()){
						auto& i = cols[0];
						auto ii = "0";
						if(cbps_map[i][ii].size() < 4) {
							cbps_map[i][ii].resize(4);
						}
						bool is_supporting = support[i] < options.support_mps || support_f[i] < options.support_mpf;
						bool is_valid_cbps = (0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3]);
						if (is_supporting) {
							continue;
						}
						if (is_valid_cbps) {
							continue;
						}
						p_ids.push_back(make_pair(i, ii));
					}
				}
			}
		}
	}
	int64_t n_total_entries = p_ids.size();
	int64_t block_size = (n_total_entries + n_cores - 1) / (double) n_cores;
	int64_t n_blocks = (n_total_entries + block_size - 1) / (double) block_size;
	vector<vector<EventEntry>> result_mp_lists(n_blocks);
	cout << "[MatePairDiscordantCaller.call_inter_chromosomal_events] start finding inter chromosomal events\n";
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			auto& local_result_mp = result_mp_lists[block_id];
			const char* underscore = "_";
			vector<string> cols;
			const int64_t n_total_entries = p_ids.size();
			const int64_t block_size = (n_total_entries + n_cores - 1) / (double)n_cores;
			const int64_t n_blocks = (n_total_entries + block_size - 1) / (double)block_size;
			const int64_t start_id = castle::ParallelRunner::get_next_start_pos(block_id, block_size, 0, 0);
			const int64_t end_id = castle::ParallelRunner::get_next_end_pos(block_id, n_blocks - 1, block_size, n_total_entries, 0);
			for(int64_t entry_id = start_id; entry_id < end_id; ++entry_id) {
				auto i = p_ids[entry_id].first;
				auto ii = p_ids[entry_id].second;

//				const bool debug = (string::npos != i.find("8883"));
//				const bool debug = (string::npos != i.find("8834"));
//				const bool debug = "22" == i;
				const bool debug = false;
// distance1, distance2: event size of possible paired clusters to query cluster
// m, n: cluster id of all possible paired clusters

				if (support[i] < options.support_mps || support_f[i] < options.support_mpf) {
					if(debug) {
						cout << "call_inter here-0: " << i << "\n";
					}
					continue;
				}
				if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
					if(debug) {
						cout << "call_inter here-1\n";
					}
					continue;
				}
// distance1, distance2: event size of possible paired clusters to query cluster
// m, n: cluster id of all possible paired clusters
				string_int64_t_value_sortedset distance1;
				string_int64_t_value_sortedset distance2;

				for (int64_t local_ii = 0;; ++local_ii) {
					string str_entry_id = boost::lexical_cast<string>(local_ii);
					if (starts_map.end() == starts_map.find(i) || starts_map[i].end() == starts_map[i].find(str_entry_id)) {
						if(debug) {
							cout << "call_inter here-2\n";
						}
						break;
					}
					if ((0 == cbps_map[i][str_entry_id][0] && 0 == cbps_map[i][str_entry_id][1])
							|| (0 == cbps_map[i][str_entry_id][2] && 0 == cbps_map[i][str_entry_id][3])) {
						if(debug) {
							cout << "call_inter here-3\n";
						}
						continue;
					}
					if (orientation.end() == orientation.find(i) || orientation[i].end() == orientation[i].find(str_entry_id)) {
						if(debug) {
							cout << "call_inter here-4\n";
						}
						continue;
					}

					if(debug) {
						cout << "cal_inter here-4-b: " << local_ii << "\n";
					}

					if (1 == orientation[i][str_entry_id].strand && -1 == orientation[i][str_entry_id].mate_strand) {
						if(debug) {
							cout << "call_inter here-5\n";
						}
						for (auto id3 : cluster_type[orientation[i][str_entry_id].ref_id][orientation[i][str_entry_id].mate_ref_id][3]) {
							castle::StringUtils::c_string_multi_split(id3, underscore, cols);
							auto k = cols[0];
							auto kk = cols[1];
							if(debug) {
								cout << "call_inter here-5b: " << k << "/" << orientation[i][str_entry_id].ref_id << "/" << orientation[i][str_entry_id].mate_ref_id << "\n";
							}
							if (support[k] < options.support_mps || support_f[k] < options.support_mpf) {
								if(debug) {
									cout << "call_inter here-6: " << k << "/" << orientation[i][str_entry_id].ref_id << "/" << orientation[i][str_entry_id].mate_ref_id << "\n";
								}
								continue;
							}
							if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1]) || (0 == cbps_map[k][kk][2] && 0 == cbps_map[k][kk][3])) {
								if(debug) {
									cout << "call_inter here-7\n";
								}
								continue;
							}
							if (i == k) {
								if(debug) {
									cout << "call_inter here-8\n";
								}
								continue;
							}

							string dis_id = i + "_" + str_entry_id + "_" + k + "_" + kk;
							StringInt64Pair a_pair_1(dis_id, abs(starts_map[k][kk].first - starts_map[i][str_entry_id].first));
							StringInt64Pair a_pair_2(dis_id, abs(mstarts_map[k][kk].first - mstarts_map[i][str_entry_id].first));
							distance1.insert(a_pair_1);
							distance2.insert(a_pair_2);

//   #print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
						}
					} else if (1 == orientation[i][str_entry_id].strand && 1 == orientation[i][str_entry_id].mate_strand) {
						if(debug) {
							cout << "call_inter here-9\n";
						}
						for (auto id2 : cluster_type[orientation[i][str_entry_id].ref_id][orientation[i][str_entry_id].mate_ref_id][2]) {
							castle::StringUtils::c_string_multi_split(id2, underscore, cols);
							auto k = cols[0];
							auto kk = cols[1];
							if (support[k] < options.support_mps || support_f[k] < options.support_mpf) {
								if(debug) {
									cout << "call_inter here-10\n";
								}
								continue;
							}
							if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1]) || (0 == cbps_map[k][kk][2] && 0 == cbps_map[k][kk][3])) {
								if(debug) {
									cout << "call_inter here-11\n";
								}
								continue;
							}
							if (i == k) {
								if(debug) {
									cout << "call_inter here-12\n";
								}
								continue;
							}
							string dis_id = i + "_" + str_entry_id + "_" + k + "_" + kk;
							StringInt64Pair a_pair_1(dis_id, abs(starts_map[k][kk].first - starts_map[i][str_entry_id].first));
							StringInt64Pair a_pair_2(dis_id, abs(mstarts_map[k][kk].first - mstarts_map[i][str_entry_id].first));
							distance1.insert(a_pair_1);
							distance2.insert(a_pair_2);
//   #print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
						}
					}
				}
//				int64_t smallestm = 0;
//				int64_t smallestn = 0;
				vector<string> mm;
				vector<string> nn;
				int64_t sei = 0;
				auto& tmp_1 = distance1.get<0>();
				if(debug) {
					for (auto key : tmp_1) {
						cout << "tmp1: " << key.first << ":" << key.second << "\n";
					}
				}
				for (auto key : tmp_1) {
					if (0 == sei) {
						mm.push_back(key.first);
//						smallestm = key.second;
						sei++;
						continue;
					}
//					if (smallestm == key.second) {
						mm.push_back(key.first);
//					}
					if(mm.size() > 3) {
						break;
					}
//					if (key.second > smallestm) {
//						break;
//					}
				}
				sei = 0;
				auto& tmp_2 = distance2.get<0>();
				if(debug) {
					for (auto key : tmp_2) {
						cout << "tmp2: " << key.first << ":" << key.second << "\n";
					}
				}
				for (auto key : tmp_2) {
					if (0 == sei) {
						nn.push_back(key.first);
//						smallestn = key.second;
						++sei;
						continue;
					}
//					if (smallestn == key.second) {
						nn.push_back(key.first);
//					}
					if(nn.size() > 3) {
						break;
					}
//					if (key.second > smallestn) {
//						break;
//					}
				}

//#print "$i\t@m\t@n\n";
				string found_mn;
				for (auto local_m : mm) {
					for (auto local_n : nn) {
						if (local_m == local_n) {
							found_mn = local_m;
							break;
						}
					}
					if (!found_mn.empty()) {
						break;
					}
				}
				if(debug) {
					cout << "found_mn: " << found_mn << "\n";
				}

				if (!found_mn.empty()) {
					castle::StringUtils::c_string_multi_split(found_mn, underscore, cols);
					auto ii = cols[1];
					auto k = cols[2];
					auto kk = cols[3];
					if(debug) {
						cout << "found_mn: here-0\n";
					}
//my ( $trash, $ii, $k, $kk ) = split( /_/, $m );
					if (k != i && orientation[i][ii].ref_id < orientation[i][ii].mate_ref_id) {
						if(debug) {
							cout << "found_mn: here-1\n";
						}
						if (orientation[i][ii].strand == orientation[i][ii].mate_strand) {
							if(debug) {
								cout << "found_mn: here-2\n";
							}
//# acceptor on smaller chr
//							if (starts_map[k][kk].first - starts_map[i][ii].second - 1 > -ref_is["rlu"]["selected"] && mstarts_map[i][ii].second > mstarts_map[k][kk].first) {
							if (mstarts_map[i][ii].second > mstarts_map[k][kk].first) {
								if(debug) {
									cout << "found_mn: here-3\n";
								}
								int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
								int64_t ins_size = mstarts_map[i][ii].second - mstarts_map[k][kk].first + 1;
								if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
									if(debug) {
										cout << "found_mn: here-4\n";
									}
									EventEntry an_entry;
									if (del_size < ref_is["rlu"]["selected"]) {
										an_entry.type = "inso";
									} else {
										an_entry.type = "del_inso";
									}
									an_entry.cluster_id = i + "_" + ii;
									an_entry.mate_cluster_id = k + "_" + kk;
									an_entry.n_supports = support[i];
									an_entry.n_mate_support = support[k];
									an_entry.ref_id = orientation[i][ii].ref_id;
									an_entry.event_start = starts_map[i][ii].second;
									an_entry.event_end = starts_map[k][kk].first;
									an_entry.event_size_1 = del_size;

									an_entry.mate_ref_id = orientation[i][ii].mate_ref_id;
									an_entry.mate_event_start = mstarts_map[k][kk].first;
									an_entry.mate_event_end = mstarts_map[i][ii].second;
									an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
									an_entry.distance = -1;
									an_entry.cbps.resize(8);
									an_entry.cbps[0] = cbps_map[i][ii][0];
									an_entry.cbps[1] = cbps_map[i][ii][1];
									an_entry.cbps[2] = cbps_map[i][ii][2];
									an_entry.cbps[3] = cbps_map[i][ii][3];
									an_entry.cbps[4] = cbps_map[k][kk][0];
									an_entry.cbps[5] = cbps_map[k][kk][1];
									an_entry.cbps[6] = cbps_map[k][kk][2];
									an_entry.cbps[7] = cbps_map[k][kk][3];
									local_result_mp.push_back(an_entry);
								}
							}
//# donor on smaller chr
							else if (starts_map[i][ii].second > starts_map[k][kk].first) {
//							else if (mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1 > -ref_is["rlu"]["selected"] && starts_map[i][ii].second > starts_map[k][kk].first) {
								if(debug) {
									cout << "found_mn: here-5\n";
								}
								int64_t del_size = mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1;
								int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
								if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
									if(debug) {
										cout << "found_mn: here-6\n";
									}
//push @used_cluster, $i;
//push @used_cluster, $k;
//$used_inter{$i} = 1;
//$used_inter{$k} = 1;
									EventEntry an_entry;
									if (del_size < ref_is["rlu"]["selected"]) {
										an_entry.type = "inso";
									} else {
										an_entry.type = "del_inso";
									}
									an_entry.cluster_id = i + "_" + ii;
									an_entry.mate_cluster_id = k + "_" + kk;
									an_entry.n_supports = support[i];
									an_entry.n_mate_support = support[k];
									an_entry.ref_id = orientation[i][ii].mate_ref_id;
									an_entry.event_start = mstarts_map[i][ii].second;
									an_entry.event_end = mstarts_map[k][kk].first;
									an_entry.event_size_1 = del_size;

									an_entry.mate_ref_id = orientation[i][ii].ref_id;
									an_entry.mate_event_start = starts_map[k][kk].first;
									an_entry.mate_event_end = starts_map[i][ii].second;
									an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
									an_entry.distance = -1;
									an_entry.cbps.resize(8);
									an_entry.cbps[0] = cbps_map[i][ii][2];
									an_entry.cbps[1] = cbps_map[i][ii][3];
									an_entry.cbps[2] = cbps_map[i][ii][0];
									an_entry.cbps[3] = cbps_map[i][ii][1];
									an_entry.cbps[4] = cbps_map[k][kk][2];
									an_entry.cbps[5] = cbps_map[k][kk][3];
									an_entry.cbps[6] = cbps_map[k][kk][0];
									an_entry.cbps[7] = cbps_map[k][kk][1];
									local_result_mp.push_back(an_entry);
								}
							}
						} else {
//# acceptor on smaller chr
							if(debug) {
								cout << "found_mn: here-7-a: starts_map[k][kk].first: " << starts_map[k][kk].first << ", starts_map[i][ii].second: " << starts_map[i][ii].second <<
								", mstarts_map[k][kk].second: " << mstarts_map[k][kk].second << ", mstarts_map[i][ii].first: " << mstarts_map[i][ii].first << "\n";
								cout << "found_mn: here-7-b: starts_map calc: " << (starts_map[k][kk].first - starts_map[i][ii].second - 1) <<
										", rlu: " << (-ref_is["rlu"]["selected"]) << "\n";
							}
//							bool first_condition = false;
//							int64_t calculated_event_size = (starts_map[k][kk].first - starts_map[i][ii].second - 1);
//							if(calculated_event_size >= 0) {
//								first_condition = calculated_event_size < ref_is["rlu"]["selected"];
//							} else {
//								first_condition = -calculated_event_size > ref_is["rlu"]["selected"];
//							}
//							if (first_condition && mstarts_map[k][kk].second > mstarts_map[i][ii].first) {
							if (mstarts_map[k][kk].second > mstarts_map[i][ii].first) {
								if(debug) {
									cout << "found_mn: here-8\n";
								}
								int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
								int64_t ins_size = mstarts_map[k][kk].second - mstarts_map[i][ii].first + 1;
								if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
									if(debug) {
										cout << "found_mn: here-9\n";
									}
									EventEntry an_entry;
									if (del_size < ref_is["rlu"]["selected"]) {
										an_entry.type = "inss";
									} else {
										an_entry.type = "del_inss";
									}
									an_entry.cluster_id = i + "_" + ii;
									an_entry.mate_cluster_id = k + "_" + kk;
									an_entry.n_supports = support[i];
									an_entry.n_mate_support = support[k];
									an_entry.ref_id = orientation[i][ii].ref_id;
									an_entry.event_start = starts_map[i][ii].second;
									an_entry.event_end = starts_map[k][kk].first;
									an_entry.event_size_1 = del_size;

									an_entry.mate_ref_id = orientation[i][ii].mate_ref_id;
									an_entry.mate_event_start = mstarts_map[i][ii].first;
									an_entry.mate_event_end = mstarts_map[k][kk].second;
									an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
									an_entry.distance = -1;
									an_entry.cbps.resize(8);
									an_entry.cbps[0] = cbps_map[i][ii][0];
									an_entry.cbps[1] = cbps_map[i][ii][1];
									an_entry.cbps[2] = cbps_map[i][ii][2];
									an_entry.cbps[3] = cbps_map[i][ii][3];
									an_entry.cbps[4] = cbps_map[k][kk][0];
									an_entry.cbps[5] = cbps_map[k][kk][1];
									an_entry.cbps[6] = cbps_map[k][kk][2];
									an_entry.cbps[7] = cbps_map[k][kk][3];
									local_result_mp.push_back(an_entry);
								}
							}
//# donor on smaller chr
							else if (starts_map[i][ii].second > starts_map[k][kk].first) {
//							else if (mstarts_map[i][ii].first - mstarts_map[k][kk].second - 1 > -ref_is["rlu"]["selected"]
//									&& starts_map[i][ii].second > starts_map[k][kk].first) {
								if(debug) {
									cout << "found_mn: here-10\n";
								}
								int64_t del_size = mstarts_map[i][ii].first - mstarts_map[k][kk].second - 1;
								int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
								if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
									if(debug) {
										cout << "found_mn: here-11\n";
									}
									EventEntry an_entry;
									if (del_size < ref_is["rlu"]["selected"]) {
										an_entry.type = "inss";
									} else {
										an_entry.type = "del_inss";
									}
									an_entry.cluster_id = i + "_" + ii;
									an_entry.mate_cluster_id = k + "_" + kk;
									an_entry.n_supports = support[i];
									an_entry.n_mate_support = support[k];
									an_entry.ref_id = orientation[i][ii].mate_ref_id;
									an_entry.event_start = mstarts_map[k][kk].second;
									an_entry.event_end = mstarts_map[i][ii].first;
									an_entry.event_size_1 = del_size;

									an_entry.mate_ref_id = orientation[i][ii].ref_id;
									an_entry.mate_event_start = starts_map[k][kk].first;
									an_entry.mate_event_end = starts_map[i][ii].second;
									an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
									an_entry.distance = -1;
									an_entry.cbps.resize(8);
									an_entry.cbps[0] = cbps_map[i][ii][2];
									an_entry.cbps[1] = cbps_map[i][ii][3];
									an_entry.cbps[2] = cbps_map[i][ii][0];
									an_entry.cbps[3] = cbps_map[i][ii][1];
									an_entry.cbps[4] = cbps_map[k][kk][2];
									an_entry.cbps[5] = cbps_map[k][kk][3];
									an_entry.cbps[6] = cbps_map[k][kk][0];
									an_entry.cbps[7] = cbps_map[k][kk][1];
									local_result_mp.push_back(an_entry);
								}
							}
						}
					}
				}

			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		auto& local_result_mp = result_mp_lists[block_id];
		result_mp.insert(result_mp.end(), local_result_mp.begin(), local_result_mp.end());
	}
	cout << checker;
}

//inter-chr events
void MatePairDiscordantCaller::call_inter_chromosomal_events_serial(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_inter_chromosomal_events_serial");
	checker.start();
//set<string> used_cluster;
	const char* underscore = "_";
	vector<string> cols;
	auto& ref_is = options.is;

	cout << "[MatePairDiscordantCaller.call_inter_chromosomal_events_serial] start finding inter chromosomal events\n";

	for (auto& chr1_entry : cluster_type) {
		auto& chr1 = chr1_entry.first;
		for (auto chr2_entry : chr1_entry.second) {
			auto& chr2 = chr2_entry.first;
			if (chr1 == chr2) {
				continue;
			}
//cout << chr1 << "/" << chr2 << "\n";
			for (auto& the_type_entry : chr2_entry.second) {
//auto the_type = the_type_entry.first;
				for (auto& id : the_type_entry.second) {
					castle::StringUtils::c_string_multi_split(id, underscore, cols);
					auto i = cols[0];
					auto ii = cols[1];
					if (support[i] < options.support_mps || support_f[i] < options.support_mpf) {
						continue;
					}
					if ((0 == cbps_map[i][ii][0] && 0 == cbps_map[i][ii][1]) || (0 == cbps_map[i][ii][2] && 0 == cbps_map[i][ii][3])) {
						continue;
					}
// distance1, distance2: event size of possible paired clusters to query cluster
// m, n: cluster id of all possible paired clusters
					string_int64_t_value_sortedset distance1;
					string_int64_t_value_sortedset distance2;

					for (int64_t local_ii = 0;; ++local_ii) {
						string str_entry_id = boost::lexical_cast<string>(local_ii);
						if (starts_map.end() == starts_map.find(i) || starts_map[i].end() == starts_map[i].find(str_entry_id)) {
							break;
						}
						if ((0 == cbps_map[i][str_entry_id][0] && 0 == cbps_map[i][str_entry_id][1])
								|| (0 == cbps_map[i][str_entry_id][2] && 0 == cbps_map[i][str_entry_id][3])) {
							continue;
						}
						if (orientation.end() == orientation.find(i) || orientation[i].end() == orientation[i].find(str_entry_id)) {
							continue;
						}

						if (1 == orientation[i][str_entry_id].strand && -1 == orientation[i][str_entry_id].mate_strand) {
							for (auto id3 : cluster_type[orientation[i][str_entry_id].ref_id][orientation[i][str_entry_id].mate_ref_id][3]) {
								castle::StringUtils::c_string_multi_split(id3, underscore, cols);
								auto k = cols[0];
								auto kk = cols[1];
								if (support[k] < options.support_mps || support_f[k] < options.support_mpf) {
									continue;
								}
								if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1]) || (0 == cbps_map[k][kk][2] && 0 == cbps_map[k][kk][3])) {
									continue;
								}
								if (i == k) {
									continue;
								}

								string dis_id = i + "_" + str_entry_id + "_" + k + "_" + kk;
								StringInt64Pair a_pair_1(dis_id, abs(starts_map[k][kk].first - starts_map[i][str_entry_id].first));
								StringInt64Pair a_pair_2(dis_id, abs(mstarts_map[k][kk].first - mstarts_map[i][str_entry_id].first));
								distance1.insert(a_pair_1);
								distance2.insert(a_pair_2);

//   #print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
							}
						} else if (1 == orientation[i][str_entry_id].strand && 1 == orientation[i][str_entry_id].mate_strand) {
							for (auto id2 : cluster_type[orientation[i][str_entry_id].ref_id][orientation[i][str_entry_id].mate_ref_id][2]) {
								castle::StringUtils::c_string_multi_split(id2, underscore, cols);
								auto k = cols[0];
								auto kk = cols[1];
								if (support[k] < options.support_mps || support_f[k] < options.support_mpf) {
									continue;
								}
								if ((0 == cbps_map[k][kk][0] && 0 == cbps_map[k][kk][1]) || (0 == cbps_map[k][kk][2] && 0 == cbps_map[k][kk][3])) {
									continue;
								}
								if (i == k) {
									continue;
								}
								string dis_id = i + "_" + str_entry_id + "_" + k + "_" + kk;
								StringInt64Pair a_pair_1(dis_id, abs(starts_map[k][kk].first - starts_map[i][str_entry_id].first));
								StringInt64Pair a_pair_2(dis_id, abs(mstarts_map[k][kk].first - mstarts_map[i][str_entry_id].first));
								distance1.insert(a_pair_1);
								distance2.insert(a_pair_2);
//   #print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
							}
						}
					}
//					int64_t smallestm = 0;
//					int64_t smallestn = 0;
					vector<string> mm;
					vector<string> nn;
					int64_t sei = 0;
					bool debug = i == "25080";
					auto& tmp_1 = distance1.get<0>();
					for (auto key : tmp_1) {
						if (debug) {
							cout << "dist_1: " << key.first << ":" << key.second << "\n";
						}
						if (0 == sei) {
							mm.push_back(key.first);
//							smallestm = key.second;
							sei++;
							continue;
						}
//						if (smallestm == key.second) {
							mm.push_back(key.first);
							if(mm.size() > 3) {
								break;
							}
//						}
//						if (key.second > smallestm) {
//							break;
//						}
					}
					sei = 0;
					auto& tmp_2 = distance2.get<0>();
					for (auto key : tmp_2) {
						if (debug) {
							cout << "dist_2: " << key.first << ":" << key.second << "\n";
						}
						if (0 == sei) {
							nn.push_back(key.first);
//							smallestn = key.second;
							++sei;
							continue;
						}
//						if (smallestn == key.second) {
							nn.push_back(key.first);
//						}
//						if (key.second > smallestn) {
//							break;
//						}
						if(nn.size() > 3) {
							break;
						}
					}

//#print "$i\t@m\t@n\n";
					string found_mn;
					for (auto local_m : mm) {
						for (auto local_n : nn) {
							if (local_m == local_n) {
								found_mn = local_m;
								break;
							}
						}
						if (!found_mn.empty()) {
							break;
						}
					}
					if (debug) {
						cout << "found mn: " << found_mn << "\n";
					}
//#print "$i\t$m\t$n\t$distance1{$m}\t$distance2{$m}\n";
					if (!found_mn.empty()) {
						castle::StringUtils::c_string_multi_split(found_mn, underscore, cols);
						auto ii = cols[1];
						auto k = cols[2];
						auto kk = cols[3];
//my ( $trash, $ii, $k, $kk ) = split( /_/, $m );
						if (debug) {
							cout << "inter-1\n";
						}
						if (k != i && orientation[i][ii].ref_id < orientation[i][ii].mate_ref_id) {
							if (debug) {
								cout << "inter-2\n";
							}
							if (orientation[i][ii].strand == orientation[i][ii].mate_strand) {
								if (debug) {
									cout << "inter-3\n";
								}
//# acceptor on smaller chr
								if (starts_map[k][kk].first - starts_map[i][ii].second - 1 > -ref_is["rlu"]["selected"]
										&& mstarts_map[i][ii].second > mstarts_map[k][kk].first) {
									if (debug) {
										cout << "inter-4\n";
									}
									int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
									int64_t ins_size = mstarts_map[i][ii].second - mstarts_map[k][kk].first + 1;
									if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
//push @used_cluster, $i;
//push @used_cluster, $k;
//$used_inter{$i} = 1;
//$used_inter{$k} = 1;
										EventEntry an_entry;
										if (del_size < ref_is["rlu"]["selected"]) {
											an_entry.type = "inso";
										} else {
											an_entry.type = "del_inso";
										}
										an_entry.cluster_id = i + "_" + ii;
										an_entry.mate_cluster_id = k + "_" + kk;
										an_entry.n_supports = support[i];
										an_entry.n_mate_support = support[k];
										an_entry.ref_id = orientation[i][ii].ref_id;
										an_entry.event_start = starts_map[i][ii].second;
										an_entry.event_end = starts_map[k][kk].first;
										an_entry.event_size_1 = del_size;

										an_entry.mate_ref_id = orientation[i][ii].mate_ref_id;
										an_entry.mate_event_start = mstarts_map[k][kk].first;
										an_entry.mate_event_end = mstarts_map[i][ii].second;
										an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
										an_entry.distance = -1;
										an_entry.cbps.resize(8);
										an_entry.cbps[0] = cbps_map[i][ii][0];
										an_entry.cbps[1] = cbps_map[i][ii][1];
										an_entry.cbps[2] = cbps_map[i][ii][2];
										an_entry.cbps[3] = cbps_map[i][ii][3];
										an_entry.cbps[4] = cbps_map[k][kk][0];
										an_entry.cbps[5] = cbps_map[k][kk][1];
										an_entry.cbps[6] = cbps_map[k][kk][2];
										an_entry.cbps[7] = cbps_map[k][kk][3];
										result_mp.push_back(an_entry);
									}
								}
//# donor on smaller chr
								else if (mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1 > -ref_is["rlu"]["selected"]
										&& starts_map[i][ii].second > starts_map[k][kk].first) {
									if (debug) {
										cout << "inter-5\n";
									}
									int64_t del_size = mstarts_map[k][kk].first - mstarts_map[i][ii].second - 1;
									int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
									if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
										if (debug) {
											cout << "inter-6\n";
										}
//push @used_cluster, $i;
//push @used_cluster, $k;
//$used_inter{$i} = 1;
//$used_inter{$k} = 1;
										EventEntry an_entry;
										if (del_size < ref_is["rlu"]["selected"]) {
											an_entry.type = "inso";
										} else {
											an_entry.type = "del_inso";
										}
										an_entry.cluster_id = i + "_" + ii;
										an_entry.mate_cluster_id = k + "_" + kk;
										an_entry.n_supports = support[i];
										an_entry.n_mate_support = support[k];
										an_entry.ref_id = orientation[i][ii].mate_ref_id;
										an_entry.event_start = mstarts_map[i][ii].second;
										an_entry.event_end = mstarts_map[k][kk].first;
										an_entry.event_size_1 = del_size;

										an_entry.mate_ref_id = orientation[i][ii].ref_id;
										an_entry.mate_event_start = starts_map[k][kk].first;
										an_entry.mate_event_end = starts_map[i][ii].second;
										an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
										an_entry.distance = -1;
										an_entry.cbps.resize(8);
										an_entry.cbps[0] = cbps_map[i][ii][2];
										an_entry.cbps[1] = cbps_map[i][ii][3];
										an_entry.cbps[2] = cbps_map[i][ii][0];
										an_entry.cbps[3] = cbps_map[i][ii][1];
										an_entry.cbps[4] = cbps_map[k][kk][2];
										an_entry.cbps[5] = cbps_map[k][kk][3];
										an_entry.cbps[6] = cbps_map[k][kk][0];
										an_entry.cbps[7] = cbps_map[k][kk][1];
										result_mp.push_back(an_entry);
									}
								}
							} else {
								if (debug) {
									cout << "inter-7\n";
									cout << starts_map[k][kk].first << "/" << starts_map[i][ii].second << "\n";
									cout << mstarts_map[k][kk].second << "/" << mstarts_map[i][ii].first << "\n";
									cout << ref_is["rlu"]["selected"] << "\n";
								}
//# acceptor on smaller chr
								if (starts_map[k][kk].first - starts_map[i][ii].second - 1 > -ref_is["rlu"]["selected"]
										&& mstarts_map[k][kk].second > mstarts_map[i][ii].first) {
									if (debug) {
										cout << "inter-7-2\n";
									}
									int64_t del_size = starts_map[k][kk].first - starts_map[i][ii].second - 1;
									int64_t ins_size = mstarts_map[k][kk].second - mstarts_map[i][ii].first + 1;
									if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
//push @used_cluster, $i;
//push @used_cluster, $k;
//$used_inter{$i} = 1;
//$used_inter{$k} = 1;
										EventEntry an_entry;
										if (del_size < ref_is["rlu"]["selected"]) {
											if (debug) {
												cout << "inter-8\n";
											}
											an_entry.type = "inss";
										} else {
											if (debug) {
												cout << "inter-9\n";
											}
											an_entry.type = "del_inss";
										}
										an_entry.cluster_id = i + "_" + ii;
										an_entry.mate_cluster_id = k + "_" + kk;
										an_entry.n_supports = support[i];
										an_entry.n_mate_support = support[k];
										an_entry.ref_id = orientation[i][ii].ref_id;
										an_entry.event_start = starts_map[i][ii].second;
										an_entry.event_end = starts_map[k][kk].first;
										an_entry.event_size_1 = del_size;

										an_entry.mate_ref_id = orientation[i][ii].mate_ref_id;
										an_entry.mate_event_start = mstarts_map[i][ii].first;
										an_entry.mate_event_end = mstarts_map[k][kk].second;
										an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
										an_entry.distance = -1;
										an_entry.cbps.resize(8);
										an_entry.cbps[0] = cbps_map[i][ii][0];
										an_entry.cbps[1] = cbps_map[i][ii][1];
										an_entry.cbps[2] = cbps_map[i][ii][2];
										an_entry.cbps[3] = cbps_map[i][ii][3];
										an_entry.cbps[4] = cbps_map[k][kk][0];
										an_entry.cbps[5] = cbps_map[k][kk][1];
										an_entry.cbps[6] = cbps_map[k][kk][2];
										an_entry.cbps[7] = cbps_map[k][kk][3];
										result_mp.push_back(an_entry);
									}
								}
//# donor on smaller chr
								else if (mstarts_map[i][ii].first - mstarts_map[k][kk].second - 1 > -ref_is["rlu"]["selected"]
										&& starts_map[i][ii].second > starts_map[k][kk].first) {
									int64_t del_size = mstarts_map[i][ii].first - mstarts_map[k][kk].second - 1;
									int64_t ins_size = starts_map[i][ii].second - starts_map[k][kk].first + 1;
									if (debug) {
										cout << "inter-10\n";
									}
									if (del_size < options.sv_size_cutoff && ins_size < options.sv_size_cutoff) {
										if (debug) {
											cout << "inter-11\n";
										}
//push @used_cluster, $i;
//push @used_cluster, $k;
//$used_inter{$i} = 1;
//$used_inter{$k} = 1;
										EventEntry an_entry;
										if (del_size < ref_is["rlu"]["selected"]) {
											if (debug) {
												cout << "inter-12\n";
											}
											an_entry.type = "inss";
										} else {
											an_entry.type = "del_inss";
											if (debug) {
												cout << "inter-13\n";
											}
										}
										an_entry.cluster_id = i + "_" + ii;
										an_entry.mate_cluster_id = k + "_" + kk;
										an_entry.n_supports = support[i];
										an_entry.n_mate_support = support[k];
										an_entry.ref_id = orientation[i][ii].mate_ref_id;
										an_entry.event_start = mstarts_map[k][kk].second;
										an_entry.event_end = mstarts_map[i][ii].first;
										an_entry.event_size_1 = del_size;

										an_entry.mate_ref_id = orientation[i][ii].ref_id;
										an_entry.mate_event_start = starts_map[k][kk].first;
										an_entry.mate_event_end = starts_map[i][ii].second;
										an_entry.event_size_2 = ins_size;
// inter-chromosomal event should not have a distance value, which only can be defined in a chromosome.
										an_entry.distance = -1;
										an_entry.cbps.resize(8);
										an_entry.cbps[0] = cbps_map[i][ii][2];
										an_entry.cbps[1] = cbps_map[i][ii][3];
										an_entry.cbps[2] = cbps_map[i][ii][0];
										an_entry.cbps[3] = cbps_map[i][ii][1];
										an_entry.cbps[4] = cbps_map[k][kk][2];
										an_entry.cbps[5] = cbps_map[k][kk][3];
										an_entry.cbps[6] = cbps_map[k][kk][0];
										an_entry.cbps[7] = cbps_map[k][kk][1];
										result_mp.push_back(an_entry);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cout << checker;
}
void MatePairDiscordantCaller::call_rest_events(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_rest_events");
	checker.start();
//set<string> used_cluster;
//	const char* underscore = "_";
	vector<string> cols;
//auto& ref_is = options.is;
//	set<string> visited;
	set<string> del_cluster;
//	for (auto& entry : result_mp) {
//		if (entry.type.empty()) {
//			continue;
//		}
//		castle::StringUtils::c_string_multi_split(entry.cluster_id, underscore, cols);
//		if (!cols.empty()) {
//			visited.insert(cols[0]);
//		}
//		castle::StringUtils::c_string_multi_split(entry.mate_cluster_id, underscore, cols);
//		if (!cols.empty()) {
//			visited.insert(cols[0]);
//		}
//	}
	cout << "[MatePairDiscordantCaller.call_rest_events] start finding rest events\n";
	cout << (boost::format("[MatePairDiscordantCaller.call_rest_events] # items: %d\n") % the_last_pid).str();

	string str_id;
	int64_t num_the_last_pid = boost::lexical_cast<int64_t>(the_last_pid);

	for (int64_t i = 0; i < num_the_last_pid; ++i) {
//		const bool debug = i == 6710;
		const bool debug = false;
		str_id = boost::lexical_cast<string>(i);
		if(debug) {
			cout << "rest here-0\n";
		}
		if (support[str_id] < options.support_mps || support_f[str_id] < options.support_mpf) {
			continue;
		}
		if(debug) {
			cout << "rest here-1\n";
		}
		if (cbps_map.end() == cbps_map.find(str_id) || cbps_map[str_id].end() == cbps_map[str_id].find("0")) {
			continue;
		}
		if(debug) {
			cout << "rest here-2\n";
		}
		if ((0 == cbps_map[str_id]["0"][0] && 0 == cbps_map[str_id]["0"][1]) || (0 == cbps_map[str_id]["0"][2] && 0 == cbps_map[str_id]["0"][3])) {
			continue;
		}
		if(debug) {
			cout << "rest here-3\n";
		}
//		if (visited.end() != visited.find(str_id)) {
//			continue;
//		}
//#print "$i\tstarts_map[i][0].first\tstarts_map[i][0].second\torientation[i][0][2]\torientation[i][0][3]\n";
		if (0 == support[str_id]) {
			continue;
		}
		if(debug) {
			cout << "rest here-4\n";
		}

		if (orientation[str_id]["0"].ref_id == orientation[str_id]["0"].mate_ref_id) {
//# simple deletion
			if (1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				bool printed = del_cluster.end() != del_cluster.find(str_id);
				if (printed) {
					continue;
				}
				int64_t del_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].second - 1;
				if (del_size > 0 && del_size < options.sv_size_cutoff) {
					del_cluster.insert(str_id);
					EventEntry an_entry;
					an_entry.type = "del";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].second;
					an_entry.event_end = mstarts_map[str_id]["0"].first;
					an_entry.event_size_1 = del_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				int64_t invers_size = mstarts_map[str_id]["0"].second - starts_map[str_id]["0"].second;
				if (invers_size < options.sv_size_cutoff) {
					EventEntry an_entry;
					an_entry.type = "invers_f";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].second;
					an_entry.event_end = mstarts_map[str_id]["0"].second;
					an_entry.event_size_1 = invers_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (-1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				int64_t invers_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].first;
				if (invers_size < options.sv_size_cutoff) {
					EventEntry an_entry;
					an_entry.type = "invers_r";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].first;
					an_entry.event_end = mstarts_map[str_id]["0"].first;
					an_entry.event_size_1 = invers_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (-1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {

				int64_t transl_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].first;
				if (transl_size < options.sv_size_cutoff) {
					if (debug) {
						cout << "tandem_dup-0\n";
					}
					EventEntry an_entry;
					an_entry.type = "tandem_dup";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].first;
					an_entry.event_end = mstarts_map[str_id]["0"].second;
					an_entry.event_size_1 = transl_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
		}

//# inter chromosomal events
		else {
			if (1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				if(debug) {
					cout << "transl_inter here-0\n";
				}
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].second;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].first;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (-1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				if(debug) {
					cout << "transl_inter here-1\n";
				}
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].first;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].second;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				if(debug) {
					cout << "transl_inter here-2\n";
				}
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].second;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].second;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (-1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				if(debug) {
					cout << "transl_inter here-3\n";
				}
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].first;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].first;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			}
		}
	}
	cout << checker;
}

void MatePairDiscordantCaller::call_rest_events_alt(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_rest_events_alt");
	checker.start();

	set<string> used_cluster;
	set<string> del_cluster;
//#	left over clusters
	cout << "[MatePairDiscordantCaller.call_rest_events_alt] start finding rest events\n";
	cout << (boost::format("[MatePairDiscordantCaller.call_rest_events_alt] # items: %d\n") % the_last_pid).str();

	string str_id;
	int64_t num_the_last_pid = boost::lexical_cast<int64_t>(the_last_pid);

	for (int64_t i = 0; i < num_the_last_pid; ++i) {
		str_id = boost::lexical_cast<string>(i);
		if(support.end() == support.find(str_id)) {
			continue;
		}
		if (support[str_id] < options.support_mps || support_f[str_id] < options.support_mpf) {
			continue;
		}
		if (cbps_map.end() == cbps_map.find(str_id) || cbps_map[str_id].end() == cbps_map[str_id].find("0")) {
			continue;
		}
		if (used_cluster.end() != used_cluster.find(str_id)) {
			continue;
		}

//#print "$i\t$start{$i}{0}[0]\t$start{$i}{0}[1]\t$orientation{$i}{0}[2]\t$orientation{$i}{0}[3]\n";
		if (orientation[str_id]["0"].ref_id == orientation[str_id]["0"].mate_ref_id) {
//			# simple deletion
			if (1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				bool printed = del_cluster.end() != del_cluster.find(str_id);
				if (printed) {
					continue;
				}
				int64_t del_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].second - 1;
				if (del_size > 0 && del_size < options.sv_size_cutoff) {
					del_cluster.insert(str_id);
					EventEntry an_entry;
					an_entry.type = "del";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].second;
					an_entry.event_end = mstarts_map[str_id]["0"].first;
					an_entry.event_size_1 = del_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				int64_t invers_size = mstarts_map[str_id]["0"].second - starts_map[str_id]["0"].second;
				if (invers_size < options.sv_size_cutoff) {
					EventEntry an_entry;
					an_entry.type = "invers_f";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].second;
					an_entry.event_end = mstarts_map[str_id]["0"].second;
					an_entry.event_size_1 = invers_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (-1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				int64_t invers_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].first;
				if (invers_size < options.sv_size_cutoff) {
					EventEntry an_entry;
					an_entry.type = "invers_r";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].first;
					an_entry.event_end = mstarts_map[str_id]["0"].first;
					an_entry.event_size_1 = invers_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
			if (-1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				int64_t transl_size = mstarts_map[str_id]["0"].first - starts_map[str_id]["0"].first;
				if (transl_size < options.sv_size_cutoff) {
					EventEntry an_entry;
					an_entry.type = "tandem_dup";
					an_entry.cluster_id = str_id;
					an_entry.n_supports = support[str_id];
					an_entry.ref_id = orientation[str_id]["0"].ref_id;
					an_entry.event_start = starts_map[str_id]["0"].first;
					an_entry.event_end = mstarts_map[str_id]["0"].second;
					an_entry.event_size_1 = transl_size;
					an_entry.cbps.resize(8);
					an_entry.cbps[0] = cbps_map[str_id]["0"][0];
					an_entry.cbps[1] = cbps_map[str_id]["0"][1];
					an_entry.cbps[2] = cbps_map[str_id]["0"][2];
					an_entry.cbps[3] = cbps_map[str_id]["0"][3];
					result_mp.push_back(an_entry);
				}
			}
		}

//		# inter chromosomal events
		else {
			if (1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].second;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].first;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (-1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].first;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].second;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (1 == orientation[str_id]["0"].strand && 1 == orientation[str_id]["0"].mate_strand) {
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].second;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].second;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			} else if (-1 == orientation[str_id]["0"].strand && -1 == orientation[str_id]["0"].mate_strand) {
				EventEntry an_entry;
				an_entry.type = "transl_inter";
				an_entry.cluster_id = str_id;
				an_entry.n_supports = support[str_id];
				an_entry.ref_id = orientation[str_id]["0"].ref_id;
				an_entry.event_start = starts_map[str_id]["0"].first;
				an_entry.strand = orientation[str_id]["0"].strand;

				an_entry.mate_ref_id = orientation[str_id]["0"].mate_ref_id;
				an_entry.mate_event_start = mstarts_map[str_id]["0"].first;
				an_entry.mate_strand = orientation[str_id]["0"].mate_strand;

				an_entry.cbps.resize(8);
				an_entry.cbps[0] = cbps_map[str_id]["0"][0];
				an_entry.cbps[1] = cbps_map[str_id]["0"][1];
				an_entry.cbps[2] = cbps_map[str_id]["0"][2];
				an_entry.cbps[3] = cbps_map[str_id]["0"][3];
				if(an_entry.cbps[2] > an_entry.mate_event_start || an_entry.cbps[3] < an_entry.mate_event_end) {
					swap(an_entry.mate_event_start, an_entry.mate_event_end);
				}
				result_mp.push_back(an_entry);
			}
		}
	}
	cout << checker;
}

//TODO: should remove multiple calls for some clusters.
// call smaller event with same primary cluster id for intra-chr events
void MatePairDiscordantCaller::call_smaller_events_among_intra_events(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_smaller_events_among_intra_events");
	checker.start();
	const string debug_str = "10";
//map<int64_t, EventEntry> current_intra_event;
//int64_t_value_sortedset ca_event_size;
//# $data{event id} = complete event string
//# $cluster_event_map{$cluster_id} = array of event ids
//# $eventsize{event id} = event size

	map<string, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;

	int64_t mpai = 0;
	vector<string> cols;
	for (auto& data : result_mp) {
		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type) {
			++mpai;
			continue;
		}
		cout << "[MatePairDiscordantCaller.call_smaller_events_among_intra_events] " << data.sr_str() << "\n";
//		bool debug = string::npos != data.cluster_id.find(debug_str) || data.mate_cluster_id.find(debug_str);
		const bool debug = false;

		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string i_pid = data.cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
			string k_pid = data.mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			if(debug) {
				cout << "[MatePairDiscordantCaller.call_smaller_events_among_intra_events] here-0: " << data.event_size_1 << ", " << data.event_size_2 << "\n";
			}
			cluster_event_map[i_pid].push_back(mpai);
			cluster_event_map[k_pid].push_back(mpai);
			eventsize[mpai] += abs(data.event_size_1);
			eventsize[mpai] += abs(data.event_size_2);
		} else {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string a_pid = data.cluster_id;
			if (!cols.empty()) {
				a_pid = cols[0];
			}
			cluster_event_map[a_pid].push_back(mpai);
			eventsize[mpai] += abs(data.event_size_1);
			if(debug) {
				cout << "[MatePairDiscordantCaller.call_smaller_events_among_intra_events] here-1: " << data.event_size_1 << "\n";
			}
		}
		++mpai;
	}

	set<int64_t> ignoring_ids;
	for (auto cluster_id : cluster_event_map) {
		if (cluster_id.second.size() < 2) {
//			if (!cluster_id.second.empty()) {
//				ignoring_ids.insert(cluster_id.second.begin(), cluster_id.second.end());
//			}
			continue;
		}
		set<int64_t> remove_candidate_ids;
		int64_t select_id = cluster_event_map[cluster_id.first][0];
		for (uint64_t c_id = 1; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (eventsize[cur_id] < eventsize[select_id]) {
				select_id = cur_id;
			}
		}
//auto& selected_event = result_mp[select_id];
//cout << "Selected: " << cluster_id.first << ":" << selected_event.str() << "\n";
		for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (cur_id != select_id) {
//auto& remove_candidate_event = result_mp[cur_id];
//cout << "Candidate: " << cluster_id.first << ":" << remove_candidate_event.str() << "\n";
				remove_candidate_ids.insert(cur_id);
				ignoring_ids.insert(cur_id);
			}
		}
		for (auto r_id : remove_candidate_ids) {
			auto& current_event = result_mp[r_id];
			castle::StringUtils::c_string_multi_split(current_event.cluster_id, "_", cols);
			bool debug = string::npos != current_event.cluster_id.find(debug_str) || current_event.mate_cluster_id.find(debug_str);
//			const bool debug = false;

			if (current_event.mate_cluster_id.empty()) {
				if(debug) {
					cout << "mate_empty: " << current_event.mpd_str() << "\n";
				}
				string i_pid = current_event.cluster_id;
				if (!cols.empty()) {
					i_pid = cols[0];
				}
				ignoring_ids.insert(boost::lexical_cast<int64_t>(i_pid));
				continue;
			}
			if(debug) {
				cout << "removal candidate: " << current_event.mpd_str() << "\n";
			}

			if (string::npos != current_event.type.find("inssd")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;

				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-0\n";
					}
					current_event.type = "del";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-1\n";
					}
					current_event.type = "tandem_dup";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
//					if(debug) {
//						cout << "call_smaller_events_among_intra_events: here-1\n";
//						cout << "cluster_id.first: " << cluster_id.first << ", cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//						string a_str = current_event.mpd_str();
//						cout << "call_smaller_events_among_intra_events result: " << a_str << "\n";
//					}
				}
			} else if (string::npos != current_event.type.find("inssu")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-2\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-2\n";
					}
					current_event.type = "del";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-3\n";
					}
					current_event.type = "tandem_dup";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.event_start;
					current_event.event_start = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("insod")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-3\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-4\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-5\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("insou")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-4\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-6\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-7\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.event_start;
					current_event.event_start = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("del_invers")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-5\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-8\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-9\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			}
		}
	}

//	for (auto& data : result_mp) {
//		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
//				|| "inso" == data.type) {
//			++mpai;
//			continue;
//		}
//		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
//			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
//			string i_pid = data.cluster_id;
//			if (!cols.empty()) {
//				i_pid = cols[0];
//			}
//			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
//			string k_pid = data.mate_cluster_id;
//			if (!cols.empty()) {
//				k_pid = cols[0];
//			}
//			bool has_other_than_del = false;
//			for(auto r_id : cluster_event_map[i_pid]) {
//				auto& current_event = result_mp[r_id];
//				if("del" != current_event.type) {
//					has_other_than_del = true;
//					break;
//				}
//			}
//			if(!has_other_than_del) {
//				for(auto r_id : cluster_event_map[k_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						has_other_than_del = true;
//						break;
//					}
//				}
//			}
//			if(has_other_than_del) {
//				for(auto r_id : cluster_event_map[i_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//				for(auto r_id : cluster_event_map[k_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//			}
//		} else {
//			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
//			string a_pid = data.cluster_id;
//			if (!cols.empty()) {
//				a_pid = cols[0];
//			}
//			bool has_other_than_del = false;
//			for(auto r_id : cluster_event_map[a_pid]) {
//				auto& current_event = result_mp[r_id];
//				if("del" != current_event.type) {
//					has_other_than_del = true;
//					break;
//				}
//			}
//			if(has_other_than_del) {
//				for(auto r_id : cluster_event_map[a_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//			}
//		}
//	}

	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}
	ofstream MPINTRA(mpintra_outfile, ios::binary);
	for (uint64_t i = 0; i < result_mp.size(); ++i) {
//if(ignoring_ids.end() != ignoring_ids.find(i)) {
//continue;
//}
		auto& data = result_mp[i];
// the following statements are for inter-chromosomal events
		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type) {
			continue;
		}
		string a_str = result_mp[i].mpd_str();
		castle::StringUtils::trim(a_str);
		if (a_str.empty()) {
			continue;
		}
		MPINTRA << a_str << "\n";
	}
	cout << checker;
}

void MatePairDiscordantCaller::call_smaller_events_among_intra_events_alt(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_smaller_events_among_intra_events_alt");
	checker.start();
//map<int64_t, EventEntry> current_intra_event;
//int64_t_value_sortedset ca_event_size;
//# $data{event id} = complete event string
//# $cluster_event_map{$cluster_id} = array of event ids
//# $eventsize{event id} = event size

	set<string> all_id_events;
//	for (auto& data : result_mp) {
//		bool found_zero_values = false;
//		for(auto& a_bp : data.cbps) {
//			if(0 == a_bp) {
//				found_zero_values = true;
//				break;
//			}
//		}
//		if(found_zero_values) {
//			continue;
//		}
//		all_id_events.insert(data.cluster_id);
//	}
//	vector<EventEntry> temp_mp;
//	for (auto& data : result_mp) {
//		bool found_zero_values = false;
//		for(auto& a_bp : data.cbps) {
//			if(0 == a_bp) {
//				found_zero_values = true;
//				break;
//			}
//		}
//		if(found_zero_values) {
//			continue;
//		}
//		if(string::npos != data.cluster_id.rfind("_0")) {
//			temp_mp.push_back(data);
//			continue;
//		}
//		string the_cluster_id = data.cluster_id + "_0";
//		if(all_id_events.end() == all_id_events.find(the_cluster_id)) {
//			temp_mp.push_back(data);
//		}
//	}
//	temp_mp.swap(result_mp);

	map<string, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;

	int64_t mpai = 0;
	vector<string> cols;
	for (auto& data : result_mp) {
		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type) {
			++mpai;
			continue;
		}
		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string i_pid = data.cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
			string k_pid = data.mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			cluster_event_map[i_pid].push_back(mpai);
			cluster_event_map[k_pid].push_back(mpai);
			eventsize[mpai] += abs(data.event_size_1);
			eventsize[mpai] += abs(data.event_size_2);
		} else {
//			cout << "[MatePairDiscordantCaller.call_smaller_events_among_intra_events_alt] " << data.mpd_pure_str() << "/" << data.sr_str() << "\n";
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string a_pid = data.cluster_id;
			if (!cols.empty()) {
				a_pid = cols[0];
			}
			cluster_event_map[a_pid].push_back(mpai);
			eventsize[mpai] += abs(data.event_size_1);
		}
		++mpai;
	}
//	const string debug_str = "7310";

	set<int64_t> ignoring_ids;
	set<string> survived_ids;
	set<string> terminate_ids;
	for (auto cluster_id : cluster_event_map) {
		if (cluster_id.second.size() < 2) {
//			if (!cluster_id.second.empty()) {
//				ignoring_ids.insert(cluster_id.second.begin(), cluster_id.second.end());
//			}
			continue;
		}
//		const bool debug = "9825" == cluster_id.first || "10476" == cluster_id.first;
		const bool debug = false;
//		const bool debug = "10" == cluster_id.first;
		set<int64_t> remove_candidate_ids;
		int64_t select_id = cluster_event_map[cluster_id.first][0];
		for (uint64_t c_id = 1; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (eventsize[cur_id] < eventsize[select_id]) {
				select_id = cur_id;
			}
		}
		if(cluster_id.second.size() > 1 && result_mp[select_id].mate_cluster_id.empty()) {
			if(debug) {
				cout << "rivals: here-0\n";
				for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
					int64_t cur_id = cluster_id.second[c_id];
					cout << result_mp[cur_id].mpd_str();
				}
			}
			for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
				int64_t cur_id = cluster_id.second[c_id];
				if(result_mp[cur_id].mate_cluster_id.empty()) {
					continue;
				}
//				if (eventsize[cur_id] < eventsize[select_id]) {
					select_id = cur_id;
					if(debug) {
						cout << "here-1: " << result_mp[select_id].mpd_str();
					}
					break;
//				}
			}
			for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
				int64_t cur_id = cluster_id.second[c_id];
				if(result_mp[cur_id].mate_cluster_id.empty()) {
					continue;
				}
				if (eventsize[cur_id] < eventsize[select_id]) {
					select_id = cur_id;
					if(debug) {
						cout << "here-2: " << result_mp[select_id].mpd_str();
					}
				}
			}
		}

		survived_ids.insert(result_mp[select_id].cluster_id);
		{
			castle::StringUtils::c_string_multi_split(result_mp[select_id].cluster_id, "_", cols);
			string i_pid = result_mp[select_id].cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(result_mp[select_id].mate_cluster_id, "_", cols);
			string k_pid = result_mp[select_id].mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			terminate_ids.insert(i_pid);
			terminate_ids.insert(k_pid);
		}
//auto& selected_event = result_mp[select_id];
//cout << "Selected: " << cluster_id.first << ":" << selected_event.str() << "\n";
		for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if(debug) {
				cout << "rival size: cur_id: " << eventsize[cur_id] << ", " << eventsize[select_id] << "\n";
				cout << "rival: " << result_mp[cur_id].mpd_str();
			}
			if (cur_id != select_id && survived_ids.end() == survived_ids.find(result_mp[cur_id].cluster_id)) {
//auto& remove_candidate_event = result_mp[cur_id];
//cout << "Candidate: " << cluster_id.first << ":" << remove_candidate_event.str() << "\n";
				remove_candidate_ids.insert(cur_id);
				ignoring_ids.insert(cur_id);
			}
		}
		for (auto r_id : remove_candidate_ids) {
			auto& current_event = result_mp[r_id];
			castle::StringUtils::c_string_multi_split(current_event.cluster_id, "_", cols);
//			const bool debug = string::npos != current_event.cluster_id.find(debug_str) || string::npos != current_event.mate_cluster_id.find(debug_str);
//			const bool debug = false;

			if (current_event.mate_cluster_id.empty()) {
				if(debug) {
					cout << "mate_empty: " << current_event.mpd_str() << "\n";
				}
				string i_pid = current_event.cluster_id;
				if (!cols.empty()) {
					i_pid = cols[0];
				}
				if(terminate_ids.end() != terminate_ids.find(i_pid)) {
					current_event.type = "";
				}
				ignoring_ids.insert(boost::lexical_cast<int64_t>(i_pid));
				continue;
			}
			if(debug) {
				cout << "removal candidate: " << current_event.mpd_str() << "\n";
			}

			if (string::npos != current_event.type.find("inssd")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;

				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-0\n";
					}
					current_event.type = "del";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-1\n";
					}
					current_event.type = "tandem_dup";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
//					if(debug) {
//						cout << "call_smaller_events_among_intra_events: here-1\n";
//						cout << "cluster_id.first: " << cluster_id.first << ", cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//						string a_str = current_event.mpd_str();
//						cout << "call_smaller_events_among_intra_events result: " << a_str << "\n";
//					}
				}
			} else if (string::npos != current_event.type.find("inssu")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-2\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-2\n";
					}
					current_event.type = "del";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-3\n";
					}
					current_event.type = "tandem_dup";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.event_start;
					current_event.event_start = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("insod")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-3\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-4\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-5\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("insou")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-4\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-6\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-7\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.event_start;
					current_event.event_start = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			} else if (string::npos != current_event.type.find("del_invers")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
				auto& mpd2 = current_event.n_mate_support;
//				if(debug) {
//					cout << "call_smaller_events_among_intra_events: here-5\n";
//					cout << "cl_id1: " << cl_id1 << ", cl_id2: " << cl_id2 << ", mpd1: " << mpd1 << ", mpd2: " << mpd2 << "\n";
//				}
				if (string::npos != cl_id1.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-8\n";
					}
					current_event.type = "invers_r";
					current_event.cluster_id = cl_id2;
					current_event.n_supports = mpd2;
					current_event.event_start = current_event.event_end;
					current_event.event_end = current_event.mate_event_end;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				} else if (string::npos != cl_id2.find(cluster_id.first)) {
					if(debug) {
						cout << "smaller intra: here-9\n";
					}
					current_event.type = "invers_f";
					current_event.cluster_id = cl_id1;
					current_event.n_supports = mpd1;
					current_event.event_end = current_event.mate_event_start;
					current_event.event_size_1 = current_event.event_end - current_event.event_start;
				}
			}
		}
	}

//	for (auto& data : result_mp) {
//		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
//				|| "inso" == data.type) {
//			++mpai;
//			continue;
//		}
//		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
//			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
//			string i_pid = data.cluster_id;
//			if (!cols.empty()) {
//				i_pid = cols[0];
//			}
//			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
//			string k_pid = data.mate_cluster_id;
//			if (!cols.empty()) {
//				k_pid = cols[0];
//			}
//			bool has_other_than_del = false;
//			for(auto r_id : cluster_event_map[i_pid]) {
//				auto& current_event = result_mp[r_id];
//				if("del" != current_event.type) {
//					has_other_than_del = true;
//					break;
//				}
//			}
//			if(!has_other_than_del) {
//				for(auto r_id : cluster_event_map[k_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						has_other_than_del = true;
//						break;
//					}
//				}
//			}
//			if(has_other_than_del) {
//				for(auto r_id : cluster_event_map[i_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//				for(auto r_id : cluster_event_map[k_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//			}
//		} else {
//			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
//			string a_pid = data.cluster_id;
//			if (!cols.empty()) {
//				a_pid = cols[0];
//			}
//			bool has_other_than_del = false;
//			for(auto r_id : cluster_event_map[a_pid]) {
//				auto& current_event = result_mp[r_id];
//				if("del" != current_event.type) {
//					has_other_than_del = true;
//					break;
//				}
//			}
//			if(has_other_than_del) {
//				for(auto r_id : cluster_event_map[a_pid]) {
//					auto& current_event = result_mp[r_id];
//					if("del" != current_event.type) {
//						continue;
//					}
//					current_event.type = "";
//				}
//			}
//		}
//	}

	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}
	ofstream MPINTRA(mpintra_outfile, ios::binary);
	for (uint64_t i = 0; i < result_mp.size(); ++i) {
//if(ignoring_ids.end() != ignoring_ids.find(i)) {
//continue;
//}
		auto& data = result_mp[i];
// the following statements are for inter-chromosomal events
		if ("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type) {
			continue;
		}
		string a_str = result_mp[i].mpd_str();
		castle::StringUtils::trim(a_str);
		if (a_str.empty()) {
			continue;
		}
		MPINTRA << a_str << "\n";
	}
	cout << checker;
}

//call smaller event if same cluster is used for inter-chr events
void MatePairDiscordantCaller::call_smaller_events_among_inter_events(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_smaller_events_among_inter_events");
	checker.start();
	map<string, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;
	vector<string> cols;
	int64_t mpei = 0;
	for (auto& data : result_mp) {
		if (!("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type)) {
			++mpei;
			continue;
		}
		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string i_pid = data.cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
			string k_pid = data.mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			cluster_event_map[i_pid].push_back(mpei);
			cluster_event_map[k_pid].push_back(mpei);
			eventsize[mpei] += abs(data.event_size_1);
			eventsize[mpei] += abs(data.event_size_2);
		} else {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string a_pid = data.cluster_id;
			if (!cols.empty()) {
				a_pid = cols[0];
			}
			cluster_event_map[a_pid].push_back(mpei);
			eventsize[mpei] += 1000000000;
		}
		++mpei;
	}
	set<int64_t> ignoring_ids;
	for (auto cluster_id : cluster_event_map) {
		if (cluster_id.second.size() < 2) {
//if (!cluster_id.second.empty()) {
//ignoring_ids.insert(cluster_id.second.begin(), cluster_id.second.end());
//}
			continue;
		}
		set<int64_t> remove_candidate_ids;
		int64_t select_id = cluster_event_map[cluster_id.first][0];
		for (uint64_t c_id = 1; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (eventsize[cur_id] < eventsize[select_id]) {
				select_id = cur_id;
			}
		}
//auto& selected_event = result_mp[select_id];
//cout << "Selected: " << cluster_id.first << ":" << selected_event.str() << "\n";
		for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (cur_id != select_id) {
//auto& remove_candidate_event = result_mp[cur_id];
//cout << "Candidate: " << cluster_id.first << ":" << remove_candidate_event.str() << "\n";
//				remove_candidate_ids.insert(cur_id);
				ignoring_ids.insert(cur_id);
			}
		}
		for (auto r_id : remove_candidate_ids) {
//#print "$cluster_id\t$select_id\t$_\n";
			auto& current_event = result_mp[r_id];
//my @current_event = split( /\t/, $data{$_} );
			if (current_event.mate_cluster_id.empty()) {
				castle::StringUtils::c_string_multi_split(current_event.cluster_id, "_", cols);
				string i_pid = current_event.cluster_id;
				if (!cols.empty()) {
					i_pid = cols[0];
				}
				ignoring_ids.insert(boost::lexical_cast<int64_t>(i_pid));
				continue;
			}
//			const bool debug = string::npos != current_event.cluster_id.find("8834");
			const bool debug = false;
			if (string::npos != current_event.type.find("inss")) {
				auto& cl_id1 = current_event.cluster_id;
//auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
//auto& mpd2 = current_event.n_mate_support;

				current_event.type = "transl_inter";
				current_event.cluster_id = cl_id1;
				current_event.mate_cluster_id = "";
				current_event.n_supports = mpd1;
				current_event.strand = 1;
				current_event.mate_strand = -1;
				if (current_event.ref_id > current_event.mate_ref_id) {
					swap(current_event.ref_id, current_event.mate_ref_id);
					swap(current_event.event_start, current_event.mate_event_start);
					swap(current_event.event_end, current_event.mate_event_end);
					swap(current_event.event_start, current_event.event_end);
					swap(current_event.cbps[0], current_event.cbps[2]);
					swap(current_event.cbps[1], current_event.cbps[3]);
					if(current_event.cbps[2] > current_event.mate_event_start || current_event.cbps[3] < current_event.mate_event_end) {
						if(debug) {
							cout << "transl_inter-1: " << current_event.mpd_str() << "\n";
						}
						swap(current_event.mate_event_start, current_event.mate_event_end);
					}
				}
			} else if (string::npos != current_event.type.find("inso")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& mpd1 = current_event.n_supports;
				current_event.type = "transl_inter";
				current_event.cluster_id = cl_id1;
				current_event.mate_cluster_id = "";
				current_event.n_supports = mpd1;
				current_event.strand = 1;
				current_event.mate_strand = 1;
				if (current_event.ref_id > current_event.mate_ref_id) {
					swap(current_event.ref_id, current_event.mate_ref_id);
					swap(current_event.event_start, current_event.mate_event_start);
					swap(current_event.event_end, current_event.mate_event_end);
					swap(current_event.event_start, current_event.event_end);
					swap(current_event.cbps[0], current_event.cbps[2]);
					swap(current_event.cbps[1], current_event.cbps[3]);
					if(current_event.cbps[2] > current_event.mate_event_start || current_event.cbps[3] < current_event.mate_event_end) {
						if(debug) {
							cout << "transl_inter-2: " << current_event.mpd_str() << "\n";
						}
						swap(current_event.mate_event_start, current_event.mate_event_end);
					}
//					swap(current_event.mate_event_start, current_event.mate_event_end);
				}

			}
		}
	}
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
	ofstream MPINTER(mpinter_outfile, ios::binary);
	for (uint64_t i = 0; i < result_mp.size(); ++i) {
//if(ignoring_ids.end() != ignoring_ids.find(i)) {
//continue;
//}
		auto& data = result_mp[i];
// the following statements are for inter-chromosomal events
		if (!("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
						|| "inso" == data.type)) {
			continue;
		}
		string a_str = result_mp[i].mpd_str();
		castle::StringUtils::trim(a_str);
		if (a_str.empty()) {
			continue;
		}
		MPINTER << a_str << "\n";
	}
	cout << checker;
}

void MatePairDiscordantCaller::call_smaller_events_among_inter_events_alt(vector<EventEntry>& result_mp) {
	castle::TimeChecker checker;
	checker.setTarget("MatePairDiscordantCaller.call_smaller_events_among_inter_events_alt");
	checker.start();
//	set<string> all_id_events;
//	for (auto& data : result_mp) {
//		bool found_zero_values = false;
//		for(auto& a_bp : data.cbps) {
//			if(0 == a_bp) {
//				found_zero_values = true;
//				break;
//			}
//		}
//		if(found_zero_values) {
//			continue;
//		}
//		all_id_events.insert(data.cluster_id);
//	}
//	vector<EventEntry> temp_mp;
//	for (auto& data : result_mp) {
//		bool found_zero_values = false;
//		for(auto& a_bp : data.cbps) {
//			if(0 == a_bp) {
//				found_zero_values = true;
//				break;
//			}
//		}
//		if(found_zero_values) {
//			continue;
//		}
//		if(string::npos != data.cluster_id.rfind("_0")) {
//			temp_mp.push_back(data);
//			continue;
//		}
//		string the_cluster_id = data.cluster_id + "_0";
//		if(all_id_events.end() == all_id_events.find(the_cluster_id)) {
//			temp_mp.push_back(data);
//		}
//	}
//	temp_mp.swap(result_mp);

	map<string, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;
	vector<string> cols;
	int64_t mpei = 0;
	for (auto& data : result_mp) {
		if (!("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
				|| "inso" == data.type)) {
			++mpei;
			continue;
		}
		if (!data.cluster_id.empty() && !data.mate_cluster_id.empty()) {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string i_pid = data.cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(data.mate_cluster_id, "_", cols);
			string k_pid = data.mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			cluster_event_map[i_pid].push_back(mpei);
			cluster_event_map[k_pid].push_back(mpei);
			eventsize[mpei] += abs(data.event_size_1);
			eventsize[mpei] += abs(data.event_size_2);
		} else {
			castle::StringUtils::c_string_multi_split(data.cluster_id, "_", cols);
			string a_pid = data.cluster_id;
			if (!cols.empty()) {
				a_pid = cols[0];
			}
			cluster_event_map[a_pid].push_back(mpei);
			eventsize[mpei] += 1000000000;
		}
		++mpei;
	}
	set<int64_t> ignoring_ids;
	set<string> survived_ids;
	set<string> terminate_ids;
	for (auto cluster_id : cluster_event_map) {
		if (cluster_id.second.size() < 2) {
//if (!cluster_id.second.empty()) {
//ignoring_ids.insert(cluster_id.second.begin(), cluster_id.second.end());
//}
			continue;
		}
		set<int64_t> remove_candidate_ids;
		int64_t select_id = cluster_event_map[cluster_id.first][0];
		for (uint64_t c_id = 1; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (eventsize[cur_id] < eventsize[select_id]) {
				select_id = cur_id;
			}
		}

		if(cluster_id.second.size() > 1 && result_mp[select_id].mate_cluster_id.empty()) {
			for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
				int64_t cur_id = cluster_id.second[c_id];
				if(result_mp[cur_id].mate_cluster_id.empty()) {
					continue;
				}
//				if (eventsize[cur_id] < eventsize[select_id]) {
					select_id = cur_id;
					break;
//				}
			}
			for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
				int64_t cur_id = cluster_id.second[c_id];
				if(result_mp[cur_id].mate_cluster_id.empty()) {
					continue;
				}
				if (eventsize[cur_id] < eventsize[select_id]) {
					select_id = cur_id;
				}
			}
		}
		survived_ids.insert(result_mp[select_id].cluster_id);
		{
			castle::StringUtils::c_string_multi_split(result_mp[select_id].cluster_id, "_", cols);
			string i_pid = result_mp[select_id].cluster_id;
			if (!cols.empty()) {
				i_pid = cols[0];
			}
			castle::StringUtils::c_string_multi_split(result_mp[select_id].mate_cluster_id, "_", cols);
			string k_pid = result_mp[select_id].mate_cluster_id;
			if (!cols.empty()) {
				k_pid = cols[0];
			}
			terminate_ids.insert(i_pid);
			terminate_ids.insert(k_pid);
		}
		for (uint64_t c_id = 0; c_id < cluster_id.second.size(); ++c_id) {
			int64_t cur_id = cluster_id.second[c_id];
			if (cur_id != select_id) {
//auto& remove_candidate_event = result_mp[cur_id];
//cout << "Candidate: " << cluster_id.first << ":" << remove_candidate_event.str() << "\n";
				remove_candidate_ids.insert(cur_id);
				ignoring_ids.insert(cur_id);
			}
		}
		for (auto r_id : remove_candidate_ids) {
//#print "$cluster_id\t$select_id\t$_\n";
			auto& current_event = result_mp[r_id];
//my @current_event = split( /\t/, $data{$_} );
			if (current_event.mate_cluster_id.empty()) {
				castle::StringUtils::c_string_multi_split(current_event.cluster_id, "_", cols);
				string i_pid = current_event.cluster_id;
				if (!cols.empty()) {
					i_pid = cols[0];
				}
				if(terminate_ids.end() != terminate_ids.find(i_pid)) {
					current_event.type = "";
				}
				ignoring_ids.insert(boost::lexical_cast<int64_t>(i_pid));
				continue;
			}
//			const bool debug = string::npos != current_event.cluster_id.find("8834");
			const bool debug = false;
			if (string::npos != current_event.type.find("inss")) {
				auto& cl_id1 = current_event.cluster_id;
//auto& cl_id2 = current_event.mate_cluster_id;
				auto& mpd1 = current_event.n_supports;
//auto& mpd2 = current_event.n_mate_support;

				current_event.type = "transl_inter";
				current_event.cluster_id = cl_id1;
				current_event.mate_cluster_id = "";
				current_event.n_supports = mpd1;
				current_event.strand = 1;
				current_event.mate_strand = -1;
				if (current_event.ref_id > current_event.mate_ref_id) {
					swap(current_event.ref_id, current_event.mate_ref_id);
					swap(current_event.event_start, current_event.mate_event_start);
					swap(current_event.event_end, current_event.mate_event_end);
					swap(current_event.event_start, current_event.event_end);
					swap(current_event.cbps[0], current_event.cbps[2]);
					swap(current_event.cbps[1], current_event.cbps[3]);
					if(current_event.cbps[2] > current_event.mate_event_start || current_event.cbps[3] < current_event.mate_event_end) {
						if(debug) {
							cout << "transl_inter-1: " << current_event.mpd_str() << "\n";
						}
						swap(current_event.mate_event_start, current_event.mate_event_end);
					}
				}
			} else if (string::npos != current_event.type.find("inso")) {
				auto& cl_id1 = current_event.cluster_id;
				auto& mpd1 = current_event.n_supports;
				current_event.type = "transl_inter";
				current_event.cluster_id = cl_id1;
				current_event.mate_cluster_id = "";
				current_event.n_supports = mpd1;
				current_event.strand = 1;
				current_event.mate_strand = 1;
				if (current_event.ref_id > current_event.mate_ref_id) {
					swap(current_event.ref_id, current_event.mate_ref_id);
					swap(current_event.event_start, current_event.mate_event_start);
					swap(current_event.event_end, current_event.mate_event_end);
					swap(current_event.event_start, current_event.event_end);
					swap(current_event.cbps[0], current_event.cbps[2]);
					swap(current_event.cbps[1], current_event.cbps[3]);
					if(current_event.cbps[2] > current_event.mate_event_start || current_event.cbps[3] < current_event.mate_event_end) {
						if(debug) {
							cout << "transl_inter-2: " << current_event.mpd_str() << "\n";
						}
						swap(current_event.mate_event_start, current_event.mate_event_end);
					}
//					swap(current_event.mate_event_start, current_event.mate_event_end);
				}

			}
		}
	}
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
	ofstream MPINTER(mpinter_outfile, ios::binary);
	for (uint64_t i = 0; i < result_mp.size(); ++i) {
//if(ignoring_ids.end() != ignoring_ids.find(i)) {
//continue;
//}
		auto& data = result_mp[i];
// the following statements are for inter-chromosomal events
		if (!("transl_inter" == data.type || "del_inss" == data.type || "del_inso" == data.type || "inss" == data.type
						|| "inso" == data.type)) {
			continue;
		}
		string a_str = result_mp[i].mpd_str();
		castle::StringUtils::trim(a_str);
		if (a_str.empty()) {
			continue;
		}
		MPINTER << a_str << "\n";
	}
	cout << checker;

}


} /* namespace meerkat */
