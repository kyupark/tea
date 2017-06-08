/*
 * CNVOverlapper.cpp
 *
 *  Created on: Oct 14, 2016
 *      Author: el174
 */

#include "CNVOverlapper.hpp"

namespace meerkat {

CNVOverlapper::CNVOverlapper() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

CNVOverlapper::~CNVOverlapper() {
}

void CNVOverlapper::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	n_cores = options.n_cores;
}

void CNVOverlapper::find_overlaps() {
	boost::unordered_map<string, StringIntervalEntryVector> cnv_interval_vec_map;
	string line;
	const char* delim_tab = "\t";
	vector<string> data;
	{
		ifstream in_cnv(options.input_cnv_filename, ios::binary);
		while(getline(in_cnv, line, '\n')) {
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			boost::replace_all(data[0], "\"", "");
			if("chrom" == data[0]) {
				continue;
			}
			boost::replace_all(data[0], "chr", "");
//			cout << castle::StringUtils::join(data, "\t") << "\n";
			int64_t start = boost::lexical_cast<int64_t>(data[1]);
			int64_t end = boost::lexical_cast<int64_t>(data[2]);
			auto& chr = data[0];
			if(start > end) {
				swap(start, end);
			}
			StringIntervalEntry an_entry(start, end, castle::StringUtils::join(data, "\t"));
			cnv_interval_vec_map[chr].push_back(an_entry);
		}
	}
	boost::unordered_map<string, StringIntervalClusterTree> cnv_interval_tree_map;
	for(auto& key_entry: cnv_interval_vec_map) {
		cnv_interval_tree_map.insert(make_pair(key_entry.first, StringIntervalClusterTree(key_entry.second)));
	}
	{
		const uint64_t n_cutoff = options.n_cutoff;
		StringIntervalEntryVector results;
		ifstream in_sv(options.input_sv_filename, ios::binary);
		ofstream out_sv(options.outfile_name, ios::binary);
		while(getline(in_sv, line, '\n')) {
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			results.clear();
			if("del" == data[0] || "del_ins" == data[0] || "invers" == data[0] || "invers_f" == data[0] || "invers_r" == data[0] || "tandem_dup" == data[0]) {
				auto& chr = data[5];
				auto tree_itr = cnv_interval_tree_map.find(chr);
				if(cnv_interval_tree_map.end() == tree_itr) {
					continue;
				}
				int64_t start = boost::lexical_cast<int64_t>(data[6]);
				int64_t end = boost::lexical_cast<int64_t>(data[7]);
				if(start > end) {
					swap(start, end);
				}
				tree_itr->second.find_overlap(start, end, results);
				if(results.size() > n_cutoff) {
					continue;
				}
				for(auto& a_result : results) {
					out_sv << a_result.value << "<=>" << line << "\n";
				}
			} else if("transl_inter" == data[0]) {
				auto& chr1 = data[5];
				auto& chr2 = data[8];
				auto tree_itr_1 = cnv_interval_tree_map.find(chr1);
				auto tree_itr_2 = cnv_interval_tree_map.find(chr2);

				if(cnv_interval_tree_map.end() == tree_itr_1 || cnv_interval_tree_map.end() == tree_itr_2) {
					continue;
				}
				int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);

				tree_itr_1->second.find_overlap(start_1, start_1, results);
				set<string> valid_results;
				for(auto& a_result : results) {
					valid_results.insert(a_result.value);
				}
				results.clear();

				int64_t start_2 = boost::lexical_cast<int64_t>(data[9]);
				tree_itr_2->second.find_overlap(start_2, start_2, results);
				for(auto& a_result : results) {
					valid_results.insert(a_result.value);
				}
				if(valid_results.size() > n_cutoff) {
					continue;
				}
				for(auto& a_result : valid_results) {
					out_sv << a_result << "<=>" << line << "\n";
				}
			} else if(string::npos != data[0].find("inso") || string::npos != data[0].find("inss") || "del_invers" == data[0]) {
				auto& chr1 = data[5];
				auto& chr2 = data[9];
				auto tree_itr_1 = cnv_interval_tree_map.find(chr1);
				auto tree_itr_2 = cnv_interval_tree_map.find(chr2);

				if(cnv_interval_tree_map.end() == tree_itr_1 || cnv_interval_tree_map.end() == tree_itr_2) {
					continue;
				}
				int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);
				int64_t end_1 = boost::lexical_cast<int64_t>(data[7]);
				if(start_1 > end_1) {
					swap(start_1, end_1);
				}
				tree_itr_1->second.find_overlap(start_1, end_1, results);

				set<string> valid_results;
				for(auto& a_result : results) {
					valid_results.insert(a_result.value);
				}
				results.clear();
				int64_t start_2 = boost::lexical_cast<int64_t>(data[10]);
				int64_t end_2 = boost::lexical_cast<int64_t>(data[11]);
				if(start_2 > end_2) {
					swap(start_2, end_2);
				}
				tree_itr_2->second.find_overlap(start_2, end_2, results);
				for(auto& a_result : results) {
					valid_results.insert(a_result.value);
				}
				if(valid_results.size() > n_cutoff) {
					continue;
				}
				for(auto& a_result : valid_results) {
					out_sv << a_result << "<=>" << line << "\n";
				}
			}
//			else {
//				cout << line << "\n";
//			}
		}
	}
}

void CNVOverlapper::find_sv_overlaps() {
	castle::TimeChecker checker;
	checker.setTarget("CNVOverlapper.find_sv_overlaps");
	checker.start();
	string control_path = options.input_sv1_filename;
	options.expand_path(control_path);
	string treatment_path(options.input_sv2_filename);
	options.expand_path(treatment_path);
	string output_path(options.outfile_name);
	options.expand_path(output_path);

	cout << "[CNVOverlapper.find_sv_overlaps] Control Group file: " << control_path << "\n";
	cout << "[CNVOverlapper.find_sv_overlaps] Treatment Group file: " << treatment_path << "\n";
	cout << "[CNVOverlapper.find_sv_overlaps] Output file: " << options.outfile_name << "\n";
	boost::unordered_map<string, ColoredStringIntervalEntryVector> svs_interval_vec_map;
	{
		const char* delim_tab = "\t";
		vector<string> data;
		int64_t n_group_id = 0;
		string line;
		string entry;
		ifstream control_list_in(control_path, ios::binary);
		while(getline(control_list_in, line, '\n')) {
			options.expand_path(line);
			ifstream in(line, ios::binary);
			while(getline(in, entry, '\n')) {
				castle::StringUtils::c_string_multi_split(entry, delim_tab, data);
				if("del" == data[0] || "del_ins" == data[0] || "invers" == data[0] || "invers_f" == data[0] || "invers_r" == data[0] || "tandem_dup" == data[0]) {
					auto& chr = data[5];
					int64_t start = boost::lexical_cast<int64_t>(data[6]);
					int64_t end = boost::lexical_cast<int64_t>(data[7]);
					if(start > end) {
						swap(start, end);
					}
					ColoredString a_colored_string;
					a_colored_string.group_id = n_group_id;
					a_colored_string.type = data[0];
					a_colored_string.value = entry;
					a_colored_string.the_size = boost::lexical_cast<int64_t>(8);
					ColoredStringIntervalEntry an_entry(start, end, a_colored_string);
					svs_interval_vec_map[chr].push_back(an_entry);
				} else if("transl_inter" == data[0]) {
					auto& chr1 = data[5];
					auto& chr2 = data[8];
					int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);
					int64_t end_1 = start_1 + 1000;
					{
						ColoredString a_colored_string;
						a_colored_string.group_id = n_group_id;
						a_colored_string.type = data[0];
						a_colored_string.value = entry;
						a_colored_string.the_size = -1;
						ColoredStringIntervalEntry an_entry(start_1, end_1, a_colored_string);
						svs_interval_vec_map[chr1].push_back(an_entry);
					}

					int64_t start_2 = boost::lexical_cast<int64_t>(data[9]);
					int64_t end_2 = start_2 + 1000;
					{
						ColoredString a_colored_string;
						a_colored_string.group_id = n_group_id;
						a_colored_string.type = data[0];
						a_colored_string.value = entry;
						a_colored_string.the_size = -1;
						ColoredStringIntervalEntry an_entry(start_2, end_2, a_colored_string);
						svs_interval_vec_map[chr2].push_back(an_entry);
					}
				} else if(string::npos != data[0].find("inso") || string::npos != data[0].find("inss") || "del_invers" == data[0]) {
					auto& chr1 = data[5];
					auto& chr2 = data[9];
					int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);
					int64_t end_1 = boost::lexical_cast<int64_t>(data[7]);
					{
						ColoredString a_colored_string;
						a_colored_string.group_id = n_group_id;
						a_colored_string.type = data[0];
						a_colored_string.value = entry;
						a_colored_string.the_size = boost::lexical_cast<int64_t>(12);
						ColoredStringIntervalEntry an_entry(start_1, end_1, a_colored_string);
						svs_interval_vec_map[chr1].push_back(an_entry);
					}
					if(start_1 > end_1) {
						swap(start_1, end_1);
					}
					int64_t start_2 = boost::lexical_cast<int64_t>(data[10]);
					int64_t end_2 = boost::lexical_cast<int64_t>(data[11]);
					{
						ColoredString a_colored_string;
						a_colored_string.group_id = n_group_id;
						a_colored_string.type = data[0];
						a_colored_string.value = entry;
						ColoredStringIntervalEntry an_entry(start_2, end_2, a_colored_string);
						svs_interval_vec_map[chr2].push_back(an_entry);
					}
				}
			}
			++n_group_id;
		}
	}
	boost::unordered_map<string, ColoredStringIntervalClusterTree> svs_interval_tree_map;
	int64_t n_total_entries = 0;
	for(auto& key_entry: svs_interval_vec_map) {
		svs_interval_tree_map.insert(make_pair(key_entry.first, ColoredStringIntervalClusterTree(key_entry.second)));
		n_total_entries += key_entry.second.size();
	}
	cout << (boost::format("[CNVOverlapper.find_sv_overlaps] # chromosomes: %d\n") % svs_interval_vec_map.size()).str();
	cout << (boost::format("[CNVOverlapper.find_sv_overlaps] # entries: %d\n") % n_total_entries).str();
	{
		ColoredStringIntervalEntryVector results;
		const char* delim_tab = "\t";
		vector<string> data;
		string line;
		string entry;
		ifstream in_sv(treatment_path, ios::binary);
		ofstream out_sv(output_path, ios::binary);
		while(getline(in_sv, line, '\n')) {
			options.expand_path(line);
			ifstream in(line, ios::binary);
			while(getline(in, entry, '\n')) {
				castle::StringUtils::c_string_multi_split(entry, delim_tab, data);
				results.clear();
				if("del" == data[0] || "del_ins" == data[0] || "invers" == data[0] || "invers_f" == data[0] || "invers_r" == data[0] || "tandem_dup" == data[0]) {
					auto& chr = data[5];
					auto tree_itr = svs_interval_tree_map.find(chr);
					if(svs_interval_tree_map.end() == tree_itr) {
						out_sv << "D\t" << entry << "\n";
						continue;
					}

					int64_t start = boost::lexical_cast<int64_t>(data[6]);
					int64_t end = boost::lexical_cast<int64_t>(data[7]);
					if(start > end) {
						swap(start, end);
					}
					int64_t cur_size = boost::lexical_cast<int64_t>(8);
					tree_itr->second.find_overlap(start, end, results);
					cout << "del chr: " << chr << ":" << results.size() << "\n";
					if(0 == results.size()) {
						out_sv << "D\t" << entry << "\n";
					} else {
						int64_t n_matches = 0;
						for(auto& a_result : results) {
							if(a_result.value.type == data[0] && abs(a_result.value.the_size - cur_size) < 10) {
								out_sv << "S\t" << a_result.value.group_id << "\t" << a_result.value.value << "<=>" << entry << "\n";
								++n_matches;
							}
						}
						if(0 == n_matches) {
							out_sv << "D\t" << entry << "\n";
						}
					}
				} else if("transl_inter" == data[0]) {
					auto& chr1 = data[5];
					auto& chr2 = data[8];

					auto tree_itr_1 = svs_interval_tree_map.find(chr1);
					auto tree_itr_2 = svs_interval_tree_map.find(chr2);

					if(svs_interval_tree_map.end() == tree_itr_1 || svs_interval_tree_map.end() == tree_itr_2) {
						out_sv << "D\t" << entry << "\n";
						continue;
					}
					int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);

					tree_itr_1->second.find_overlap(start_1, start_1, results);
					set<ColoredString> valid_results;
					for(auto& a_result : results) {
						valid_results.insert(a_result.value);
					}
					results.clear();

					int64_t start_2 = boost::lexical_cast<int64_t>(data[9]);
					tree_itr_2->second.find_overlap(start_2, start_2, results);
					for(auto& a_result : results) {
						valid_results.insert(a_result.value);
					}
					if(0 == valid_results.size()) {
						out_sv << "D\t" << entry << "\n";
					} else {
						int64_t n_matches = 0;
						for(auto& a_result : valid_results) {
							if(a_result.type == data[0]) {
								out_sv << "S\t" << a_result.group_id << "\t" << a_result.value << "<=>" << entry << "\n";
								++n_matches;
							}
						}
						if(0 == n_matches) {
							out_sv << "D\t" << entry << "\n";
						}
					}
				} else if(string::npos != data[0].find("inso") || string::npos != data[0].find("inss") || "del_invers" == data[0]) {
					auto& chr1 = data[5];
					auto& chr2 = data[9];
					auto tree_itr_1 = svs_interval_tree_map.find(chr1);
					auto tree_itr_2 = svs_interval_tree_map.find(chr2);

					if(svs_interval_tree_map.end() == tree_itr_1 || svs_interval_tree_map.end() == tree_itr_2) {
						out_sv << "D\t" << entry << "\n";
						continue;
					}
					cout << "transl chr: " << chr1 << "/" << chr2 << "\n";
					int64_t start_1 = boost::lexical_cast<int64_t>(data[6]);
					int64_t end_1 = boost::lexical_cast<int64_t>(data[7]);
					if(start_1 > end_1) {
						swap(start_1, end_1);
					}
					tree_itr_1->second.find_overlap(start_1, end_1, results);

					set<ColoredString> valid_results;
					for(auto& a_result : results) {
						valid_results.insert(a_result.value);
					}
					results.clear();
					int64_t start_2 = boost::lexical_cast<int64_t>(data[10]);
					int64_t end_2 = boost::lexical_cast<int64_t>(data[11]);
					if(start_2 > end_2) {
						swap(start_2, end_2);
					}
					tree_itr_2->second.find_overlap(start_2, end_2, results);
					for(auto& a_result : results) {
						valid_results.insert(a_result.value);
					}
					if(0 == valid_results.size()) {
						out_sv << "D\t" << entry << "\n";
					} else {
						int64_t n_matches = 0;
						for(auto& a_result : valid_results) {
							if(a_result.type == data[0]) {
								out_sv << "S\t" << a_result.group_id << "\t" << a_result.value << "<=>" << entry << "\n";
								++n_matches;
							}
						}
						if(0 == n_matches) {
							out_sv << "D\t" << entry << "\n";
						}
					}
				}
			}
		}
	}
	cout << checker;
}

} /* namespace meerkat */
