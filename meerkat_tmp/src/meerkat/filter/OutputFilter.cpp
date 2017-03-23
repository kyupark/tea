/*
 * OutputFilter.cpp
 *
 *  Created on: Jul 18, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#include "OutputFilter.hpp"

namespace meerkat {

OutputFilter::OutputFilter() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

OutputFilter::~OutputFilter() {
}
void OutputFilter::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}

void OutputFilter::filter() {
	filter_intra_chromosomal_events();
	filter_inter_chromosomal_events();
}
void OutputFilter::filter_intra_chromosomal_events() {
	castle::TimeChecker checker;
	checker.setTarget("OutputFilter.filter_intra_chromosomal_events");
	checker.start();
	string& prefix = options.prefix;
	string srintra_out = prefix + ".sr.intra.out";
	string srintra_filter = prefix + ".sr.intra.filtered";
	if(!options.working_dir.empty()) {
		srintra_out = options.working_prefix + ".sr.intra.out";
		srintra_filter = options.working_prefix + ".sr.intra.filtered";
	}

	map<int64_t, string> the_events;
	vector<string> data;
	vector<string> cl_ids;
	vector<string> current_event;
	vector<string> cluster_ids;
	vector<string> mpds;
	vector<string> srds;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";

	map<int64_t, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;

//	# data{event id} = complete event string
//	# $cluster_event_map{$cluster_id} = array of event ids
//	# $eventsize{event id} = event size
	int64_t ftai = 0;
	string line;
	ifstream SRINTRAOUTFT(srintra_out, ios::binary);

	while (getline(SRINTRAOUTFT, line, '\n')) {
//		cout << line << "\n";
		castle::StringUtils::tokenize(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl_ids);
		int64_t local_event_sizes = 0;
		int64_t max_id = min(static_cast<int64_t>(cl_ids.size()), static_cast<int64_t>(2));
		for (int64_t c_id = 0; c_id < max_id; ++c_id) {
			auto& a_cl = cl_ids[c_id];
			auto a_pos = a_cl.rfind("_");
			if (string::npos != a_pos) {
				a_cl = a_cl.substr(0, a_pos);
			}
			int64_t a = boost::lexical_cast<int64_t>(a_cl);
			cluster_event_map[a].push_back(ftai);
			auto& type = data[0];
			if("tandem_dup" == type || "del" == type) {
				if (0 == c_id) {
					local_event_sizes += abs(boost::lexical_cast<int64_t>(data[7]));
				} else {
					if(data.size() > 11) {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[11]));
					}
				}
			} else {
				if (0 == c_id) {
					local_event_sizes += abs(boost::lexical_cast<int64_t>(data[7]));
				} else {
					if(data.size() > 11) {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[11]));
					}
				}
			}

		}
		eventsize[ftai] = local_event_sizes;
		the_events[ftai] = line;
		++ftai;
	}
	for (auto& an_event_entry : cluster_event_map) {
		auto& cluster_events = an_event_entry.second;

		if (cluster_events.size() > 1) {
			int64_t select_id = cluster_events[0];
			const auto an_event_entry = eventsize.find(select_id);
			if(eventsize.end() == an_event_entry) {
				continue;
			}

			vector<int64_t> delete_ids;
			for (auto& the_cluster_id : cluster_events) {
				const auto& first_event_entry = eventsize.find(the_cluster_id);
				if(eventsize.end() == first_event_entry) {
					continue;
				}
				const auto& second_event_entry = eventsize.find(select_id);
				if(eventsize.end() == second_event_entry) {
					continue;
				}
				if (first_event_entry->second < second_event_entry->second) {
					select_id = the_cluster_id;
				}
			}
			for (auto& the_cluster_id : cluster_events) {
				const auto& first_event_entry = eventsize.find(the_cluster_id);
				if(eventsize.end() == first_event_entry) {
					continue;
				}
				if (the_cluster_id != select_id) {
					delete_ids.push_back(the_cluster_id);
				}
			}
			for (auto& a_cluster_id : delete_ids) {
				castle::StringUtils::tokenize(the_events[a_cluster_id], delim_tab, current_event);
				if (string::npos != current_event[1].find("/")) {
					if (string::npos != current_event[0].find("inssd")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);
						int64_t end_pos = boost::lexical_cast<int64_t>(current_event[6]);
						int64_t start_pos = boost::lexical_cast<int64_t>(current_event[5]);
						int64_t length = end_pos - start_pos;

						vector<string> temp_event;
						temp_event.push_back("tandem_dup");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);
//						temp_event.push_back("0");
						// ref_id
						temp_event.push_back(current_event[4]);
						// event end
						temp_event.push_back(current_event[6]);
						// mate event end
						temp_event.push_back(current_event[10]);
						temp_event.push_back(boost::lexical_cast<string>(length));
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					} else if (string::npos != current_event[0].find("inssu")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);
						int64_t end_pos = boost::lexical_cast<int64_t>(current_event[6]);
						int64_t start_pos = boost::lexical_cast<int64_t>(current_event[5]);
						int64_t length = end_pos - start_pos;
						vector<string> temp_event;
						temp_event.push_back("tandem_dup");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);
//						temp_event.push_back("0");
						temp_event.push_back(current_event[4]);
						// mate start
						temp_event.push_back(current_event[9]);
						// start
						temp_event.push_back(current_event[5]);
						temp_event.push_back(boost::lexical_cast<string>(length));
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					} else if (string::npos != current_event[0].find("insod")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);
						int64_t end_pos = boost::lexical_cast<int64_t>(current_event[6]);
						int64_t start_pos = boost::lexical_cast<int64_t>(current_event[5]);
						int64_t length = end_pos - start_pos;
						vector<string> temp_event;
						temp_event.push_back("invers_f");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);
						temp_event.push_back(current_event[4]);
						// start
						temp_event.push_back(current_event[5]);
						// mate end
						temp_event.push_back(current_event[10]);
						temp_event.push_back(boost::lexical_cast<string>(length));
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					} else if (string::npos != current_event[0].find("insou")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);
						int64_t end_pos = boost::lexical_cast<int64_t>(current_event[6]);
						int64_t start_pos = boost::lexical_cast<int64_t>(current_event[5]);
						int64_t length = end_pos - start_pos;
						vector<string> temp_event;
						temp_event.push_back("invers_f");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);
						temp_event.push_back(current_event[4]);
						// mate end
						temp_event.push_back(current_event[10]);
						// start
						temp_event.push_back(current_event[5]);
						temp_event.push_back(boost::lexical_cast<string>(length));
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					} else if ("del_invers" == current_event[0]) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);
						int64_t end_pos = boost::lexical_cast<int64_t>(current_event[6]);
						int64_t start_pos = boost::lexical_cast<int64_t>(current_event[5]);
						int64_t length = end_pos - start_pos;
						vector<string> temp_event;
						temp_event.push_back("invers_f");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);
						temp_event.push_back(current_event[4]);
						// start
						temp_event.push_back(current_event[5]);
						// mate end
						temp_event.push_back(current_event[10]);
						temp_event.push_back(boost::lexical_cast<string>(length));
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					}
				} else {
					the_events.erase(a_cluster_id);
				}
			}
		}
	}
	vector<string> dataa;
	vector<string> datab;
	map<int64_t, string> naive_events = the_events;
	for (int64_t i = 0; i <= ftai; ++i) {
		if(the_events.end() == the_events.find(i)) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(the_events[i], delim_tab, dataa);
		dataa.erase(dataa.begin() + 1, dataa.begin() + 3);
		naive_events[i] = castle::StringUtils::join(dataa, "__");
	}
	for (int64_t i = 0; i <= ftai; ++i) {
		if(the_events.end() == the_events.find(i)) {
			continue;
		}
//		castle::StringUtils::c_string_multi_split(the_events[i], delim_tab, dataa);
//		dataa.erase(dataa.begin() + 1, dataa.begin() + 3);
//		string stringa = castle::StringUtils::join(dataa, "__");
		string& stringa = naive_events[i];
		for (int64_t j = 0; j <= ftai; ++j) {
			if (i == j || the_events.end() == the_events.find(j)) {
				continue;
			}

//			castle::StringUtils::c_string_multi_split(the_events[j], delim_tab, datab);
//			datab.erase(datab.begin() + 1, datab.begin() + 3);
//			string stringb = castle::StringUtils::join(datab, "__");
			string& stringb = naive_events[j];
			if (stringa == stringb) {
				the_events.erase(j);
			}
		}
	}
	ofstream SRINTRAFILTER(srintra_filter, ios::binary);
	for (int64_t i = 0; i <= ftai; ++i) {
		if (the_events.end() == the_events.find(i)) {
			continue;
		}
		SRINTRAFILTER << the_events[i] << "\n";
	}
	cout << checker;
}

void OutputFilter::filter_inter_chromosomal_events() {
	castle::TimeChecker checker;
	checker.setTarget("OutputFilter.filter_inter_chromosomal_events");
	checker.start();
	string& prefix = options.prefix;
	string srinter_out = prefix + ".sr.inter.out";
	string srinter_filter = prefix + ".sr.inter.filtered";
	if(!options.working_dir.empty()) {
		srinter_out = options.working_prefix + ".sr.inter.out";
		srinter_filter = options.working_prefix + ".sr.inter.filtered";
	}
	map<int64_t, string> the_events;
	vector<string> data;
	vector<string> cl_ids;
	vector<string> current_event;
	vector<string> cluster_ids;
	vector<string> mpds;
	vector<string> srds;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";

//	my (%data, %cluster_event_map, %eventsize);

//	# $cluster_event_map{$cluster_id} = array of event ids
//	# $eventsize{event id} = event size
	map<int64_t, vector<int64_t>> cluster_event_map;
	map<int64_t, int64_t> eventsize;
	int64_t ftei = 0;
	string line;
	ifstream SRINTEROUTFT(srinter_out, ios::binary);
	while (getline(SRINTEROUTFT, line, '\n')) {
//		cout << line << "\n";
		castle::StringUtils::tokenize(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl_ids);
		int64_t local_event_sizes = 0;
		if (cl_ids.size() > 1) {
			for (int64_t c_id = 0; c_id < 2; ++c_id) {
				auto& a_cl = cl_ids[c_id];
				auto a_pos = a_cl.rfind("_");
				if (string::npos != a_pos) {
					a_cl = a_cl.substr(0, a_pos);
				}
				int64_t a = boost::lexical_cast<int64_t>(a_cl);
				cluster_event_map[a].push_back(ftei);
				auto& type = data[0];
				if("tandem_dup" == type || "del" == type) {
					if (0 == c_id) {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[8]));
					} else {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[12]));
					}
				} else {
					if (0 == c_id) {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[7]));
					} else {
						local_event_sizes += abs(boost::lexical_cast<int64_t>(data[11]));
					}
				}
			}
			eventsize[ftei] = local_event_sizes;

		} else if (1 == cl_ids.size()) {
			auto& a_cl = cl_ids[0];
			auto a_pos = a_cl.rfind("_");
			if (string::npos != a_pos) {
				a_cl = a_cl.substr(0, a_pos);
			}
			int64_t a = boost::lexical_cast<int64_t>(a_cl);
			cluster_event_map[a].push_back(ftei);
			eventsize[ftei] = 1000000000;
		}

		the_events[ftei] = line;
		++ftei;
	}
	for (auto& an_event_entry : cluster_event_map) {
		auto& cluster_events = an_event_entry.second;
		if (cluster_events.size() > 1) {
//				#print "$cluster_id\t@{$cluster_event_map{$cluster_id}}\n";
			vector<int64_t> delete_ids;
			int64_t select_id = cluster_events[0];
			for (auto& the_cluster_id : cluster_events) {
				const auto& first_event_entry = eventsize.find(the_cluster_id);
				const auto& second_event_entry = eventsize.find(select_id);
				if (first_event_entry->second < second_event_entry->second) {
					select_id = the_cluster_id;
				}
			}
			for (auto& the_cluster_id : cluster_events) {
				if (the_cluster_id != select_id) {
					delete_ids.push_back(the_cluster_id);
				}
			}
			for (auto& a_cluster_id : delete_ids) {
				castle::StringUtils::tokenize(the_events[a_cluster_id], delim_tab, current_event);
				if (string::npos != current_event[1].find("/")) {
					if (string::npos != current_event[0].find("inss")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);

						vector<string> temp_event;
						temp_event.push_back("transl_inter");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);

						if (current_event[4] < current_event[8]) {
							temp_event.push_back(current_event[4]);
							temp_event.push_back(current_event[5]);
							temp_event.push_back("1");
							temp_event.push_back(current_event[8]);
							temp_event.push_back(current_event[9]);
							temp_event.push_back("-1");
						} else {
							temp_event.push_back(current_event[8]);
							temp_event.push_back(current_event[10]);
							temp_event.push_back("1");
							temp_event.push_back(current_event[4]);
							temp_event.push_back(current_event[6]);
							temp_event.push_back("-1");
						}
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					} else if (string::npos != current_event[0].find("inso")) {
						castle::StringUtils::c_string_multi_split(current_event[1], delim_slash, cluster_ids);
						castle::StringUtils::c_string_multi_split(current_event[2], delim_slash, mpds);
						castle::StringUtils::c_string_multi_split(current_event[3], delim_slash, srds);

						vector<string> temp_event;
						temp_event.push_back("transl_inter");
						temp_event.push_back(cluster_ids[0]);
						temp_event.push_back(mpds[0]);
						temp_event.push_back(srds[0]);

						if (current_event[5] < current_event[9]) {
							temp_event.push_back(current_event[5]);
							temp_event.push_back(current_event[6]);
							temp_event.push_back("1");
							temp_event.push_back(current_event[9]);
							temp_event.push_back(current_event[11]);
							temp_event.push_back("1");
						} else {
							temp_event.push_back(current_event[9]);
							temp_event.push_back(current_event[11]);
							temp_event.push_back("1");
							temp_event.push_back(current_event[5]);
							temp_event.push_back(current_event[6]);
							temp_event.push_back("1");
						}
						string modified_event = castle::StringUtils::join(temp_event, "\t");
						the_events[a_cluster_id] = modified_event;
					}
				} else {
					the_events.erase(a_cluster_id);
				}
			}
		}
	}

	vector<string> dataa;
	vector<string> datab;
	map<int64_t, string> naive_events = the_events;
	for (int64_t i = 0; i <= ftei; ++i) {
		if(the_events.end() == the_events.find(i)) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(the_events[i], delim_tab, dataa);
		dataa.erase(dataa.begin() + 1, dataa.begin() + 3);
		naive_events[i] = castle::StringUtils::join(dataa, "__");
	}
	for (int64_t i = 0; i <= ftei; ++i) {
		if(the_events.end() == the_events.find(i)) {
			continue;
		}
		string& stringa = naive_events[i];
		for (int64_t j = 0; j <= ftei; ++j) {
			if (i == j || the_events.end() == the_events.find(j)) {
				continue;
			}
			string& stringb = naive_events[j];
			if (stringa == stringb) {
				the_events.erase(j);
			}
		}
	}

	ofstream SRINTERFILTER(srinter_filter, ios::binary);
	for (int64_t i = 0; i <= ftei; ++i) {
		if (the_events.end() == the_events.find(i)) {
			continue;
		}
		SRINTERFILTER << the_events[i] << "\n";
	}
//	my $srintra_coord = prefix + ".sr.intra.filtered.crd.sorted";
//	my $srinter_coord = prefix + ".sr.inter.filtered.crd.sorted";

//	#system "sort $srintra_filter -k 5,5 -k 6,6n > $srintra_coord";
//	#system "sort $srinter_filter -k 5,5 -k 6,6n > $srinter_coord";
	cout << checker;
}

} /* namespace meerkat */
