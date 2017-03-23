/*
 * ResultCompartor.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: el174
 */

#include "ResultCompartor.hpp"

namespace tea {

ResultCompartor::ResultCompartor() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

ResultCompartor::~ResultCompartor() {
}
void ResultCompartor::set_option_parser(const TEAOptionParser& the_options) {
	options = the_options;
	n_cores = options.n_cores;
}
void ResultCompartor::find_overlaps(const string& comp_suffix) {
	castle::TimeChecker checker;
	checker.setTarget("ResultCompartor.find_overlaps");
	checker.start();
	string control_path = options.l1;
	options.expand_path(control_path);
	string treatment_path(options.l2);
	options.expand_path(treatment_path);
	string output_path(options.comp_out + comp_suffix);
	options.expand_path(output_path);

	cout << "[ResultCompartor.find_overlaps] Control Group file: " << control_path << "\n";
	cout << "[ResultCompartor.find_overlaps] Treatment Group file: " << treatment_path << "\n";
	cout << "[ResultCompartor.find_overlaps] Output file: " << output_path << "\n";
	boost::unordered_map<string, GermlineStringIntervalEntryVector> germline_interval_vec_map;
	{
		const char* delim_tab = "\t";
		vector<string> data;
		string entry;
		ifstream in(control_path, ios::binary);
		while(getline(in, entry, '\n')) {
			castle::StringUtils::c_string_multi_split(entry, delim_tab, data);
			if(data[data.size() -1] != "5") {
				continue;
			}
			auto& chr = data[1];
			auto& family = data[9];
			int64_t pbp = boost::lexical_cast<int64_t>(data[6]);
			int64_t nbp = boost::lexical_cast<int64_t>(data[7]);
			if(pbp > nbp) {
				swap(pbp, nbp);
			}
			if(-1 == pbp) {
				pbp = nbp;
			}
			string a_key = chr + ":" + family;
			GermlineString a_colored_string;
			a_colored_string.value = entry;
			a_colored_string.pbp = pbp;
			a_colored_string.nbp = nbp;
			GermlineStringIntervalEntry an_entry(pbp - 100, nbp + 100, a_colored_string);
			germline_interval_vec_map[a_key].push_back(an_entry);
		}
	}
	boost::unordered_map<string, GermlineStringIntervalClusterTree> germline_interval_tree_map;
	int64_t n_total_entries = 0;
	for(auto& key_entry: germline_interval_vec_map) {
		germline_interval_tree_map.insert(make_pair(key_entry.first, GermlineStringIntervalClusterTree(key_entry.second)));
		n_total_entries += key_entry.second.size();
	}
	cout << (boost::format("[ResultCompartor.find_overlaps] # chromosomes: %d\n") % germline_interval_vec_map.size()).str();
	cout << (boost::format("[ResultCompartor.find_overlaps] # entries: %d\n") % n_total_entries).str();
	{
		GermlineStringIntervalEntryVector results;
		const char* delim_tab = "\t";
		vector<string> data;
		string entry;
		ifstream in(treatment_path, ios::binary);
		ofstream out_germline(output_path, ios::binary);
		getline(in, entry, '\n');
		out_germline << entry << "\n";
		while(getline(in, entry, '\n')) {
			castle::StringUtils::c_string_multi_split(entry, delim_tab, data);
			if(data[data.size() -1] != "5") {
				continue;
			}
			results.clear();
			auto& chr = data[1];
			auto& family = data[9];
			string a_key = chr + ":" + family;

			auto tree_itr = germline_interval_tree_map.find(a_key);
			if(germline_interval_tree_map.end() == tree_itr) {
				out_germline << "D\t" << entry << "\n";
				continue;
			}
			int64_t pbp = boost::lexical_cast<int64_t>(data[6]);
			int64_t nbp = boost::lexical_cast<int64_t>(data[7]);

			if(pbp > nbp) {
				swap(pbp, nbp);
			}
			tree_itr->second.find_overlap(pbp, nbp, results);
			if(0 == results.size()) {
				out_germline << "D\t" << entry << "\n";
			} else {
				int64_t min_delta = numeric_limits<int64_t>::max();
				bool has_match = false;
				string match;

				for(auto& a_result : results) {
					int64_t pbp_delta = abs(a_result.value.pbp - pbp);
					int64_t nbp_delta = abs(a_result.value.nbp - nbp);
					int64_t cur_delta = pbp_delta + nbp_delta;
					if(min_delta > cur_delta) {
						match = a_result.value.value;
						has_match = true;

						min_delta = cur_delta;
					}
				}
				if(has_match) {
					out_germline << "S\t" << match << "<=>" << entry << "\n";
				} else {
					out_germline << "D\t" << entry << "\n";
				}
			}
		}
	}
	cout << checker;

}
} /* namespace tea */
