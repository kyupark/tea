/*
 * Mechanism.cpp
 *
 *  Created on: Jul 25, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#include "Mechanism.hpp"

namespace meerkat {

Mechanism::Mechanism() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

Mechanism::~Mechanism() {
}

void Mechanism::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}

void Mechanism::call_variants() {
	map<string, vector<RepeatEntry>> te;
	map<string, vector<RepeatEntry>> sr;
	collect_repeat_data(te, sr);
	vector<deque<string>> variants;
	collect_variants(variants);
	call_mechanism(te, sr, variants);
//	call_mechanism_alt(te, sr, variants);
//	clean_temporary_files();
//	call_mechanism_alt(te, sr, variants);
}

void Mechanism::call_variants_alt() {
	map<string, vector<RepeatEntry>> te;
	map<string, vector<RepeatEntry>> sr;
	collect_repeat_data(te, sr);
	vector<deque<string>> variants;
	collect_variants(variants);
	call_mechanism_alt(te, sr, variants);
}

void Mechanism::collect_repeat_data(map<string, vector<RepeatEntry>>& te, map<string, vector<RepeatEntry>>& sr) {
	castle::TimeChecker checker;
	checker.setTarget("Mechanism.collect_repeat_data");
	checker.start();

	vector<function<void()> > tasks;

	int64_t file_size = castle::IOUtils::get_file_size(options.rmskfile);
	const int64_t BLOCK_SIZE = file_size /(double)n_cores;
	int64_t n_blocks = file_size / (double)BLOCK_SIZE;
	vector<uint64_t> skip_points;
	castle::IOUtils::find_skip_points(skip_points, options.rmskfile, BLOCK_SIZE, file_size, n_blocks, n_cores);
	vector<map<string, vector<RepeatEntry>>> te_lists(n_blocks);
	vector<map<string, vector<RepeatEntry>>> sr_lists(n_blocks);
	for(int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id, file_size]{
			const char* delim_tab = "\t";
			vector<string> data;
			string line;

			string the_first_line;
			const int64_t beg_pos = skip_points[block_id];
			int64_t cur_pos = beg_pos;
			const int64_t next_pos = skip_points[block_id + 1];
			auto& local_te = te_lists[block_id];
			auto& local_sr = sr_lists[block_id];
			ifstream RMSK(options.rmskfile, ios::binary);
			if(0 != cur_pos) {
				RMSK.seekg(cur_pos, ios::beg);
			}
			while (getline(RMSK, line, '\n')) {
				if(cur_pos == beg_pos) {
					the_first_line = line;
				}
				cur_pos += line.size() + 1;
				castle::StringUtils::trim(line);
				if('#' == line[0]) {
					continue;
				}
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
				try {
				auto& repeat_class = data[11];
				auto& ref_id = data[5];
				int64_t start = boost::lexical_cast<int64_t>(data[6]);
				int64_t end = boost::lexical_cast<int64_t>(data[7]);

				RepeatEntry an_entry;
				an_entry.start = start;
				an_entry.end = end;
				an_entry.repeat_class = repeat_class;
				if ("LINE" == repeat_class || "SINE" == repeat_class || "LTR" == repeat_class || "DNA" == repeat_class) {
					auto& name = data[10];
					an_entry.name = name;
					local_te[ref_id].push_back(an_entry);
		//			#te[strand][tei[strand]][2] = data[9]; // # strand
		//			#print "te[strand][tei[strand]][0]\tte[strand][tei[strand]][1]\tte[strand][tei[strand]][2]\tte[strand][tei[strand]][3]\tte[strand][tei[strand]][4]\tte[strand][tei[strand]][5]\n";
				} else if (options.include_other && "Other" == repeat_class) {
					auto& name = data[10];
					an_entry.name = name;
					local_te[ref_id].push_back(an_entry);
		//			te[ref_id][tei[ref_id]][3] = data[10]; //# name
		//			#te[strand][tei[strand]][2] = data[9]; //# strand
		//			#print "te[strand][tei[strand]][0]\tte[strand][tei[strand]][1]\tte[strand][tei[strand]][2]\tte[strand][tei[strand]][3]\tte[strand][tei[strand]][4]\tte[strand][tei[strand]][5]\n";
				} else if ("Satellite" == repeat_class || "Simple_repeat" == repeat_class || "Low_complexity" == repeat_class) {
					local_sr[ref_id].push_back(an_entry);
		//			#sr[strand][sri[strand]][3] = data[12];# family
		//			#sr[strand][sri[strand]][4] = data[9]; # strand
		//			#sr[strand][sri[strand]][5] = data[10];# name
		//			#print "sr[strand][sri[strand]][0]\tsr[strand][sri[strand]][1]\tsr[strand][sri[strand]][2]\tsr[strand][sri[strand]][3]\tsr[strand][sri[strand]][4]\tsr[strand][sri[strand]][5]\n";
				}
				} catch(exception& ex) {
					cout << line << "\n";
					cout << ex.what() << "\n";
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
//			cout << (boost::format("%d: %d-%d\n") % block_id % beg_pos % next_pos).str();
//			cout << (boost::format("%d: %s\n") % block_id % the_first_line).str();
//			cout << (boost::format("%d: %s\n") % block_id % line).str();
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	tasks.push_back([&]{
		for(int64_t block_id = 0; block_id < n_blocks; ++block_id) {
			auto& local_te = te_lists[block_id];
			te.insert(local_te.begin(), local_te.end());
		}
	});
	tasks.push_back([&]{
		for(int64_t block_id = 0; block_id < n_blocks; ++block_id) {
			auto& local_sr = sr_lists[block_id];
			sr.insert(local_sr.begin(), local_sr.end());
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}
void Mechanism::collect_variants(vector<deque<string>>& variants) {
	string intra_refine_type = options.prefix + ".intra.refined.typ.sorted";
	string inter_refine_type = options.prefix + ".inter.refined.typ.sorted";
	if(!options.working_dir.empty()) {
		intra_refine_type = options.working_prefix + ".intra.refined.typ.sorted";
		inter_refine_type = options.working_prefix + ".inter.refined.typ.sorted";
	}
	int64_t intra_file_size = castle::IOUtils::get_file_size(intra_refine_type);
	int64_t intra_chunk_size = intra_file_size / (double) n_cores;
	int64_t intra_n_blocks = intra_file_size / (double) intra_chunk_size;
	vector<uint64_t> intra_skip_points;
	castle::IOUtils::find_skip_points(intra_skip_points, intra_refine_type, intra_chunk_size, intra_file_size, intra_n_blocks, n_cores);
	vector<vector<deque<string>>> intra_blocks(intra_n_blocks);
	vector<function<void()> > tasks;
	for(int64_t block_id = 0; block_id < intra_n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
		int64_t cur_pos = intra_skip_points[block_id];
		int64_t next_pos = intra_skip_points[block_id + 1];
		const char* delim_tab = "\t";
		deque<string> data;
		string line;
		auto& local_variants = intra_blocks[block_id];
		ifstream FILE(intra_refine_type, ios::binary);
		if(0 != block_id) {
			FILE.seekg(cur_pos, ios::beg);
		}
		while (getline(FILE, line, '\n')) {
			cur_pos += line.size() + 1;
			castle::StringUtils::tokenize(line, delim_tab, data);
			local_variants.push_back(data);
			if(cur_pos >= next_pos) {
				break;
			}
		}
//	}
	});
	}

	int64_t inter_file_size = castle::IOUtils::get_file_size(inter_refine_type);
	int64_t inter_chunk_size = inter_file_size / (double) n_cores;
	int64_t inter_n_blocks = inter_file_size / (double) intra_chunk_size;
	vector<uint64_t> inter_skip_points;
	castle::IOUtils::find_skip_points(inter_skip_points, inter_refine_type, inter_chunk_size, inter_file_size, inter_n_blocks, n_cores);
	vector<vector<deque<string>>> inter_blocks(inter_n_blocks);

	for(int64_t block_id = 0; block_id < inter_n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			int64_t cur_pos = inter_skip_points[block_id];
			int64_t next_pos = inter_skip_points[block_id + 1];
			const char* delim_tab = "\t";
			deque<string> data;
			string line;
			auto& local_variants = inter_blocks[block_id];
			ifstream FILE(inter_refine_type, ios::binary);
			if(0 != block_id) {
				FILE.seekg(cur_pos, ios::beg);
			}
			while (getline(FILE, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delim_tab, data);
				local_variants.push_back(data);
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	int64_t n_elements = 0;
	for(int64_t block_id = 0; block_id < intra_n_blocks; ++block_id) {
		auto& local_variants = intra_blocks[block_id];
		n_elements += local_variants.size();
	}
	for(int64_t block_id = 0; block_id < inter_n_blocks; ++block_id) {
		auto& local_variants = inter_blocks[block_id];
		n_elements += local_variants.size();
	}
	variants.reserve(n_elements);
	for(int64_t block_id = 0; block_id < intra_n_blocks; ++block_id) {
		auto& local_variants = intra_blocks[block_id];
		move(begin(local_variants), end(local_variants), back_inserter(variants));
//		variants.insert(variants.end(), local_variants.begin(), local_variants.end());
	}
	for(int64_t block_id = 0; block_id < inter_n_blocks; ++block_id) {
		auto& local_variants = inter_blocks[block_id];
//		variants.insert(variants.end(), local_variants.begin(), local_variants.end());
		move(begin(local_variants), end(local_variants), back_inserter(variants));
	}

}

void Mechanism::call_mechanism(const map<string, vector<RepeatEntry>>& te, const map<string, vector<RepeatEntry>>& sr, vector<deque<string>>& variants) {
	castle::TimeChecker checker;
	checker.setTarget("Mechanism.call_mechanism");
	checker.start();
	vector<function<void()> > tasks;
//	LOOP:
	const int64_t max_id = variants.size();
//	cout << (boost::format("[Mechanism.call_mechanism] # elems: %d\n") % max_id).str();
	int64_t BLOCK_SIZE = max_id / (double)n_cores;
	int64_t n_blocks = max_id / (double)BLOCK_SIZE;
	vector<string> output_file_names(n_blocks);

	for(int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id]{
			const int64_t max_id = variants.size();
			int64_t BLOCK_SIZE = max_id / (double)n_cores;
			int64_t n_blocks = max_id / (double)BLOCK_SIZE;
			const int64_t start_id = castle::ParallelRunner::get_next_start_pos(block_id, BLOCK_SIZE, 0, 0);
			const int64_t end_id = castle::ParallelRunner::get_next_end_pos(block_id, n_blocks - 1, BLOCK_SIZE, max_id, 0);

//			cout << (boost::format("[Mechanism.call_mechanism] Block-%d: %d-%d\n") % block_id % start_id % end_id).str();
			string str_block_id = boost::lexical_cast<string>(block_id);
			string outfile = options.prefix + ".variants." + str_block_id;
			if(!options.working_dir.empty()) {
				outfile = options.working_prefix + ".variants." + str_block_id;
			}
			output_file_names[block_id] = outfile;
			ofstream OUT(outfile, ios::binary);
			for(int64_t d_id = start_id; d_id < end_id; ++d_id) {
				auto& data = variants[d_id];
				auto type = data[0];
				string mechanism;
		//		my (mechanism, @bp_annotation);
		//		#print "@data\n";
		//		# call mechanism
				if ("del" == type) {
		//			# large events
		//			if (data[7] > options.sv_size_cutoff) {
		//				;
		//			} else
					int64_t sv_size = boost::lexical_cast<int64_t>(data[7]);
					if (sv_size <= options.sv_size_cutoff) {
						string& ref_id = data[4];
						int64_t start = boost::lexical_cast<int64_t>(data[5]);
						int64_t end = boost::lexical_cast<int64_t>(data[6]);
						int64_t indicator = boost::lexical_cast<int64_t>(data[8]);
						string te_class;
						string te_name;
						call_tei(te_class, te_name, te, ref_id, start, end);
						if ("complex" == te_class) {
							mechanism = "TEI_complex";
						} else if (!te_class.empty()) {
							mechanism = "TEI_" + te_class + "_" + te_name;
						}

						if (mechanism.empty()) {
							bool vntr = call_vntr(sr, ref_id, start, end);
							if (vntr) {
								mechanism = "VNTR";
							}
							if (mechanism.empty()) {
								if (indicator > 100) {
									mechanism = "NAHR";
								}
								if (mechanism.empty()) {
									if (indicator >= 2 && indicator <= 100) {
										mechanism = "alt-EJ";
									}
									if (mechanism.empty()) {
										mechanism = "NHEJ";
									}
								}
							}
		//				if (mechanism) {
		//					goto LEND;
		//				}
						}
					}
				} else if (string::npos != type.find("ins")) {
		//			# large events
		//			if (data[7] > options.sv_size_cutoff || data[11] > options.sv_size_cutoff) {
		//				;
		//			} else {
					int64_t sv_size_1 = boost::lexical_cast<int64_t>(data[7]);
					int64_t sv_size_2 = boost::lexical_cast<int64_t>(data[11]);
					if (!(sv_size_1 > options.sv_size_cutoff || sv_size_2 > options.sv_size_cutoff)) {
						string& ref_id = data[4];
						int64_t start = boost::lexical_cast<int64_t>(data[5]);
						int64_t end = boost::lexical_cast<int64_t>(data[6]);
						string& mate_ref_id = data[8];
						int64_t mate_start = 0;
						int64_t mate_end = 0;
						if ("-" != data[8]) {
							mate_start = boost::lexical_cast<int64_t>(data[9]);
							mate_end = boost::lexical_cast<int64_t>(data[10]);
						}
						int64_t indicator = boost::lexical_cast<int64_t>(data[11]);

		//				# deletion with insertion inside the break points
						if (string::npos != type.find("del")) {
		//					# known source of insertion
							if ("-" != data[8]) {
								string te_class_d;
								string te_name_d;
								string te_class_i;
								string te_name_i;
								call_tei(te_class_d, te_name_d, te, ref_id, start, end);
								call_tei(te_class_i, te_name_i, te, mate_ref_id, mate_start, mate_end);
								if (!te_class_d.empty() && !te_class_i.empty()) {
									mechanism = "TEA_" + te_class_d + "_" + te_class_i;
								}
								if (mechanism.empty()) {
									if (te_class_i == "complex") {
										mechanism = "TEI_complex";
									} else if (!te_class_i.empty()) {
										mechanism = "TEI_" + te_class_i + "_" + te_name_i;
									}

									if (mechanism.empty()) {
										bool vntr_d = call_vntr(sr, ref_id, start, end);
										bool vntr_i = call_vntr(sr, mate_ref_id, mate_start, mate_end);
										if (vntr_d && vntr_i) {
											mechanism = "VNTR";
										}
										if (mechanism.empty()) {
											if (indicator <= 10) {
												mechanism = "NHEJ";
											}
										}

		//						if (mechanism)
		//							goto LEND;
									}
								}
							}
		//					# unknown source of insertion
							else {
								string te_class;
								string te_name;
								call_tei(te_class, te_name, te, ref_id, start, end);
								if ("complex" == te_class) {
									mechanism = "TEI_complex";
								} else if (!te_class.empty()) {
									mechanism = "TEI_" + te_class + "_" + te_name;
								}

								if (mechanism.empty()) {
									if (call_vntr(sr, ref_id, start, end)) {
										mechanism = "VNTR";
									}
									if (mechanism.empty()) {
										if (indicator <= 10) {
											mechanism = "NHEJ";
										}
									}
								}
		//						if (mechanism)
		//							goto LEND;
							}
							if (mechanism.empty() && indicator > 10) {
								mechanism = "FoSTeS";
							}
		//					if (mechanism) {
		//						goto LEND;
		//					}
						}
		//				# insertion
						else {
							string te_class_i;
							string te_name_i;
							call_tei(te_class_i, te_name_i, te, mate_ref_id, mate_start, mate_end);
							if ("complex" == te_class_i) {
								mechanism = "TEI_complex";
							} else if (!te_class_i.empty()) {
								mechanism = "TEI_" + te_class_i + "_" + te_name_i;
							}

							if (mechanism.empty()) {
								if (call_vntr(sr, mate_ref_id, mate_start, mate_end)) {
									mechanism = "VNTR";
								}
							}
		//					if (mechanism) {
		//						goto LEND;
		//					}
						}
					}
				} else if (string::npos != type.find("invers")) {
		//			# del_invers
					if (string::npos != type.find("del")) {
						mechanism = "FoSTeS";
					}
		//			# invers_f, invers_r and invers
		//			else
		//			{
		//				;
		//			}
				} else if ("tandem_dup" == type) {
					string& ref_id = data[4];
					int64_t start = boost::lexical_cast<int64_t>(data[5]);
					int64_t end = boost::lexical_cast<int64_t>(data[6]);
					if (call_vntr(sr, ref_id, start, end)) {
						mechanism = "VNTR";
					}
		//			if (mechanism) {
		//				goto LEND;
		//			}
				}
		//		else if ("transl_inter" == type) {
		//			;
		//		}
		//		LEND:
				if (mechanism.empty()) {
					mechanism = "NA";
				}
		//		cout << data[0] << ":" << data[1] << ":" << data[2] << "\n";
				// the data is shifted by 1 at this point!
				data.insert(data.begin() + 1, mechanism);
		//		data.pop_front();
		//		data.push_front(mechanism);
		//		data.push_front(type);
		//		cout << data[0] << ":" << data[1] << ":" << data[2] << "\n";
		//		shift(data);
		//		unshift(data, mechanism);
		//		unshift(data, type);

				vector<string> bp_annotation;

		//		# annotate break points
				if ("del" == type || "tandem_dup" == type || "invers" == type || "invers_f" == type || "invers_r" == type) {
					bp_annotation.resize(2);
					string& ref_id = data[5];
					int64_t start = boost::lexical_cast<int64_t>(data[6]);
					int64_t end = boost::lexical_cast<int64_t>(data[7]);
					if (sr_te_overlap(sr, ref_id, start)) {
						bp_annotation[0] = "SR";
					} else if (sr_te_overlap(te, ref_id, start)) {
						bp_annotation[0] = "TE";
					} else {
						bp_annotation[0] = "";
					}
					if (sr_te_overlap(sr, ref_id, end)) {
						bp_annotation[1] = "SR";
					} else if (sr_te_overlap(te, ref_id, end)) {
						bp_annotation[1] = "TE";
					} else {
						bp_annotation[1] = "";
					}
				} else if (string::npos != type.find("ins") || "del_invers" == type) {
					string& ref_id = data[5];
					int64_t start = boost::lexical_cast<int64_t>(data[6]);
					int64_t end = boost::lexical_cast<int64_t>(data[7]);
		//			# known source of insertion
					if ("-" != data[9]) {
						bp_annotation.resize(4);
						string& mate_ref_id = data[9];
						int64_t mate_start = boost::lexical_cast<int64_t>(data[10]);
						int64_t mate_end = boost::lexical_cast<int64_t>(data[11]);
						if (sr_te_overlap(sr, ref_id, start)) {
							bp_annotation[0] = "SR";
						} else if (sr_te_overlap(te, ref_id, start)) {
							bp_annotation[0] = "TE";
						} else {
							bp_annotation[0] = "";
						}
						if (sr_te_overlap(sr, ref_id, end)) {
							bp_annotation[1] = "SR";
						} else if (sr_te_overlap(te, ref_id, end)) {
							bp_annotation[1] = "TE";
						} else {
							bp_annotation[1] = "";
						}
						if (sr_te_overlap(sr, mate_ref_id, mate_start)) {
							bp_annotation[2] = "SR";
						} else if (sr_te_overlap(te, mate_ref_id, mate_start)) {
							bp_annotation[2] = "TE";
						} else {
							bp_annotation[2] = "";
						}
						if (sr_te_overlap(sr, mate_ref_id, mate_end)) {
							bp_annotation[3] = "SR";
						} else if (sr_te_overlap(te, mate_ref_id, mate_end)) {
							bp_annotation[3] = "TE";
						} else {
							bp_annotation[3] = "";
						}
					}
		//			# unknown source of insertion
					else {
						bp_annotation.resize(2);
						if (sr_te_overlap(sr, ref_id, start)) {
							bp_annotation[0] = "SR";
						} else if (sr_te_overlap(te, ref_id, start)) {
							bp_annotation[0] = "TE";
						} else {
							bp_annotation[0] = "";
						}
						if (sr_te_overlap(sr, ref_id, end)) {
							bp_annotation[1] = "SR";
						} else if (sr_te_overlap(te, ref_id, end)) {
							bp_annotation[1] = "TE";
						} else {
							bp_annotation[1] = "";
						}
					}
				} else if ("transl_inter" == type) {
					bp_annotation.resize(2);
					auto& ref_id = data[5];
					int64_t start = boost::lexical_cast<int64_t>(data[6]);
					auto& mate_ref_id = data[8];
					int64_t mate_start = boost::lexical_cast<int64_t>(data[9]);
					if (sr_te_overlap(sr, ref_id, start)) {
						bp_annotation[0] = "SR";
					} else if (sr_te_overlap(te, ref_id, start)) {
						bp_annotation[0] = "TE";
					} else {
						bp_annotation[0] = "";
					}
					if (sr_te_overlap(sr, mate_ref_id, mate_start)) {
						bp_annotation[1] = "SR";
					} else if (sr_te_overlap(te, mate_ref_id, mate_start)) {
						bp_annotation[1] = "TE";
					} else {
						bp_annotation[1] = "";
					}
		//			data[4] += "/" + data[5];
		//			data.erase(data.begin() + 5);
				}
				string bp_annotation_str = castle::StringUtils::join(bp_annotation, "_");
				bp_annotation_str = "BP:" + bp_annotation_str;
				data.push_back(bp_annotation_str);
				string toprint = castle::StringUtils::join(data, "\t");
				OUT << toprint << "\n";
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string outfile = options.prefix + ".variants";
	if(!options.working_dir.empty()) {
		outfile = options.working_prefix + ".variants";
	}
	castle::IOUtils::plain_file_merge(outfile, output_file_names, n_cores, true);
	cout << checker;
}

void Mechanism::call_mechanism_alt(const map<string, vector<RepeatEntry>>& te, const map<string, vector<RepeatEntry>>& sr, vector<deque<string>>& variants) {
	castle::TimeChecker checker;
	checker.setTarget("Mechanism.call_mechanism_alt");
	checker.start();
	const int64_t max_id = variants.size();

	string outfile = options.prefix + ".variants";
	if(!options.working_dir.empty()) {
		outfile = options.working_prefix + ".variants";
	}
	bool debug = false;
	ofstream OUT(outfile, ios::binary);
	for(int64_t d_id = 0; d_id < max_id; ++d_id) {
		auto& data = variants[d_id];
		auto type = data[0];
		string a_line = castle::StringUtils::join(data, ",");
		cout << a_line << "\n";
		if("del_ins" == type && "6107" == data[1]) {
			debug = true;
		}
		if(!debug) {
			continue;
		}
		string mechanism;
//		# call mechanism
		if ("del" == type) {
//			# large events
//			if (data[7] > options.sv_size_cutoff) {
//				;
//			} else
			if(debug) {
				cout << "del: here-0\n";
			}
			int64_t sv_size = boost::lexical_cast<int64_t>(data[7]);
			if (sv_size <= options.sv_size_cutoff) {
				if(debug) {
					cout << "del: here-1\n";
				}
				string& ref_id = data[4];
				int64_t start = boost::lexical_cast<int64_t>(data[5]);
				int64_t end = boost::lexical_cast<int64_t>(data[6]);
				int64_t indicator = boost::lexical_cast<int64_t>(data[8]);
				string te_class;
				string te_name;
				call_tei(te_class, te_name, te, ref_id, start, end);
				if(debug) {
					cout << "del: here-2\n";
				}
				if ("complex" == te_class) {
					mechanism = "TEI_complex";
				} else if (!te_class.empty()) {
					mechanism = "TEI_" + te_class + "_" + te_name;
				}

				if (mechanism.empty()) {
					bool vntr = call_vntr(sr, ref_id, start, end);
					if (vntr) {
						mechanism = "VNTR";
					}
					if (mechanism.empty()) {
						if (indicator > 100) {
							mechanism = "NAHR";
						}
						if (mechanism.empty()) {
							if (indicator >= 2 && indicator <= 100) {
								mechanism = "alt-EJ";
							}
							if (mechanism.empty()) {
								mechanism = "NHEJ";
							}
						}
					}
//				if (mechanism) {
//					goto LEND;
//				}
				}
			}
		} else if (string::npos != type.find("ins")) {
//			# large events
//			if (data[7] > options.sv_size_cutoff || data[11] > options.sv_size_cutoff) {
//				;
//			} else {
			if(debug) {
				cout << "ins: here-0\n";
			}
			int64_t sv_size_1 = boost::lexical_cast<int64_t>(data[7]);
			int64_t sv_size_2 = boost::lexical_cast<int64_t>(data[11]);
			if (!(sv_size_1 > options.sv_size_cutoff || sv_size_2 > options.sv_size_cutoff)) {
				if(debug) {
					cout << "ins: here-1\n";
				}
				string& ref_id = data[4];
				int64_t start = boost::lexical_cast<int64_t>(data[5]);
				int64_t end = boost::lexical_cast<int64_t>(data[6]);
				string& mate_ref_id = data[8];
				int64_t mate_start = 0;
				int64_t mate_end = 0;
				if ("-" != data[8]) {
					mate_start = boost::lexical_cast<int64_t>(data[9]);
					mate_end = boost::lexical_cast<int64_t>(data[10]);
				}
				int64_t indicator = boost::lexical_cast<int64_t>(data[11]);

//				# deletion with insertion inside the break points
				if (string::npos != type.find("del")) {
					if(debug) {
						cout << "del_ins: here-0\n";
					}
//					# known source of insertion
					if ("-" != data[8]) {
						if(debug) {
							cout << "del_ins: here-1\n";
						}
						string te_class_d;
						string te_name_d;
						string te_class_i;
						string te_name_i;
						call_tei(te_class_d, te_name_d, te, ref_id, start, end);
						call_tei(te_class_i, te_name_i, te, mate_ref_id, mate_start, mate_end);
						if (!te_class_d.empty() && !te_class_i.empty()) {
							mechanism = "TEA_" + te_class_d + "_" + te_class_i;
						}
						if (mechanism.empty()) {
							if (te_class_i == "complex") {
								mechanism = "TEI_complex";
							} else if (!te_class_i.empty()) {
								mechanism = "TEI_" + te_class_i + "_" + te_name_i;
							}

							if (mechanism.empty()) {
								bool vntr_d = call_vntr(sr, ref_id, start, end);
								bool vntr_i = call_vntr(sr, mate_ref_id, mate_start, mate_end);
								if (vntr_d && vntr_i) {
									mechanism = "VNTR";
								}
								if (mechanism.empty()) {
									if (indicator <= 10) {
										mechanism = "NHEJ";
									}
								}

//						if (mechanism)
//							goto LEND;
							}
						}
					}
//					# unknown source of insertion
					else {
						if(debug) {
							cout << "del_ins: here-2\n";
						}
						string te_class;
						string te_name;
						call_tei(te_class, te_name, te, ref_id, start, end);
						if ("complex" == te_class) {
							mechanism = "TEI_complex";
						} else if (!te_class.empty()) {
							mechanism = "TEI_" + te_class + "_" + te_name;
						}

						if (mechanism.empty()) {
							if (call_vntr(sr, ref_id, start, end)) {
								mechanism = "VNTR";
							}
							if (mechanism.empty()) {
								if (indicator <= 10) {
									mechanism = "NHEJ";
								}
							}
						}
//						if (mechanism)
//							goto LEND;
					}
					if (mechanism.empty() && indicator > 10) {
						mechanism = "FoSTeS";
					}
//					if (mechanism) {
//						goto LEND;
//					}
				}
//				# insertion
				else {
					if(debug) {
						cout << "del_ins: here-3\n";
					}
					string te_class_i;
					string te_name_i;
					call_tei(te_class_i, te_name_i, te, mate_ref_id, mate_start, mate_end);
					if ("complex" == te_class_i) {
						mechanism = "TEI_complex";
					} else if (!te_class_i.empty()) {
						mechanism = "TEI_" + te_class_i + "_" + te_name_i;
					}

					if (mechanism.empty()) {
						if (call_vntr(sr, mate_ref_id, mate_start, mate_end)) {
							mechanism = "VNTR";
						}
					}
//					if (mechanism) {
//						goto LEND;
//					}
				}
			}
		} else if (string::npos != type.find("invers")) {
//			# del_invers
			if (string::npos != type.find("del")) {
				mechanism = "FoSTeS";
			}
//			# invers_f, invers_r and invers
//			else
//			{
//				;
//			}
		} else if ("tandem_dup" == type) {
			string& ref_id = data[4];
			int64_t start = boost::lexical_cast<int64_t>(data[5]);
			int64_t end = boost::lexical_cast<int64_t>(data[6]);
			if (call_vntr(sr, ref_id, start, end)) {
				mechanism = "VNTR";
			}
//			if (mechanism) {
//				goto LEND;
//			}
		}
//		else if ("transl_inter" == type) {
//			;
//		}
//		LEND:
		if (mechanism.empty()) {
			mechanism = "NA";
		}
		if(debug) {
			cout << "here-10-a: " << castle::StringUtils::join(data, ",") << "\n";;
		}
//		cout << data[0] << ":" << data[1] << ":" << data[2] << "\n";
		// the data is shifted by 1 at this point!
		data.insert(data.begin() + 1, mechanism);
//		data.pop_front();
//		data.push_front(mechanism);
//		data.push_front(type);
//		cout << data[0] << ":" << data[1] << ":" << data[2] << "\n";
//		shift(data);
//		unshift(data, mechanism);
//		unshift(data, type);

		vector<string> bp_annotation;
		if(debug) {
			cout << "here-10-b: " << castle::StringUtils::join(data, ",") << "\n";;
		}

//		# annotate break points
		if ("del" == type || "tandem_dup" == type || "invers" == type || "invers_f" == type || "invers_r" == type) {
			if(debug) {
				cout << "here-11\n";
			}
			bp_annotation.resize(2);
			string& ref_id = data[5];
			int64_t start = boost::lexical_cast<int64_t>(data[6]);
			int64_t end = boost::lexical_cast<int64_t>(data[7]);
			if (sr_te_overlap(sr, ref_id, start)) {
				if(debug) {
					cout << "here-12\n";
				}
				bp_annotation[0] = "SR";
			} else if (sr_te_overlap(te, ref_id, start)) {
				if(debug) {
					cout << "here-13\n";
				}
				bp_annotation[0] = "TE";
			} else {
				if(debug) {
					cout << "here-14\n";
				}
				bp_annotation[0] = "";
			}
			if (sr_te_overlap(sr, ref_id, end)) {
				if(debug) {
					cout << "here-15\n";
				}
				bp_annotation[1] = "SR";
			} else if (sr_te_overlap(te, ref_id, end)) {
				if(debug) {
					cout << "here-16\n";
				}
				bp_annotation[1] = "TE";
			} else {
				if(debug) {
					cout << "here-17\n";
				}
				bp_annotation[1] = "";
			}
		} else if (string::npos != type.find("ins") || "del_invers" == type) {
			if(debug) {
				cout << "del_ins: here-5: " << castle::StringUtils::join(data, ",") << "\n";
			}
			string& ref_id = data[5];
			int64_t start = boost::lexical_cast<int64_t>(data[6]);
			int64_t end = boost::lexical_cast<int64_t>(data[7]);
//			# known source of insertion
			if ("-" != data[9]) {
				if(debug) {
					cout << "del_ins: here-6\n";
				}
				bp_annotation.resize(4);
				string& mate_ref_id = data[9];
				int64_t mate_start = boost::lexical_cast<int64_t>(data[10]);
				int64_t mate_end = boost::lexical_cast<int64_t>(data[11]);
				if (sr_te_overlap(sr, ref_id, start)) {
					bp_annotation[0] = "SR";
				} else if (sr_te_overlap(te, ref_id, start)) {
					bp_annotation[0] = "TE";
				} else {
					bp_annotation[0] = "";
				}
				if (sr_te_overlap(sr, ref_id, end)) {
					bp_annotation[1] = "SR";
				} else if (sr_te_overlap(te, ref_id, end)) {
					bp_annotation[1] = "TE";
				} else {
					bp_annotation[1] = "";
				}
				if (sr_te_overlap(sr, mate_ref_id, mate_start)) {
					bp_annotation[2] = "SR";
				} else if (sr_te_overlap(te, mate_ref_id, mate_start)) {
					bp_annotation[2] = "TE";
				} else {
					bp_annotation[2] = "";
				}
				if (sr_te_overlap(sr, mate_ref_id, mate_end)) {
					bp_annotation[3] = "SR";
				} else if (sr_te_overlap(te, mate_ref_id, mate_end)) {
					bp_annotation[3] = "TE";
				} else {
					bp_annotation[3] = "";
				}
			}
//			# unknown source of insertion
			else {
				bp_annotation.resize(2);
				if (sr_te_overlap(sr, ref_id, start)) {
					bp_annotation[0] = "SR";
				} else if (sr_te_overlap(te, ref_id, start)) {
					bp_annotation[0] = "TE";
				} else {
					bp_annotation[0] = "";
				}
				if (sr_te_overlap(sr, ref_id, end)) {
					bp_annotation[1] = "SR";
				} else if (sr_te_overlap(te, ref_id, end)) {
					bp_annotation[1] = "TE";
				} else {
					bp_annotation[1] = "";
				}
			}
		} else if ("transl_inter" == type) {
			bp_annotation.resize(2);
			auto& ref_id = data[5];
			int64_t start = boost::lexical_cast<int64_t>(data[6]);
			auto& mate_ref_id = data[8];
			int64_t mate_start = boost::lexical_cast<int64_t>(data[9]);
			if (sr_te_overlap(sr, ref_id, start)) {
				bp_annotation[0] = "SR";
			} else if (sr_te_overlap(te, ref_id, start)) {
				bp_annotation[0] = "TE";
			} else {
				bp_annotation[0] = "";
			}
			if (sr_te_overlap(sr, mate_ref_id, mate_start)) {
				bp_annotation[1] = "SR";
			} else if (sr_te_overlap(te, mate_ref_id, mate_start)) {
				bp_annotation[1] = "TE";
			} else {
				bp_annotation[1] = "";
			}
//			data[4] += "/" + data[5];
//			data.erase(data.begin() + 5);
		}
		if(debug) {
			cout << "here-12\n";
		}
		string bp_annotation_str = castle::StringUtils::join(bp_annotation, "_");
		bp_annotation_str = "BP:" + bp_annotation_str;
		if(debug) {
			cout << "here-13: " << bp_annotation_str << "\n";
		}
		data.push_back(bp_annotation_str);
		string toprint = castle::StringUtils::join(data, "\t");
		OUT << toprint << "\n";
	}
	cout << checker;

}

void Mechanism::clean_temporary_files() {
	vector<function<void()> > tasks;
	vector<vector<string>> file_name_lists(4);
	tasks.push_back([&]{
		auto& local_file_names = file_name_lists[0];
		int64_t n_blocks = 0;
		while(true) {
			string str_block_id = boost::lexical_cast<string>(n_blocks);
			string file_name = options.prefix + ".softclips.fq.gz."	+ str_block_id;
			if(!options.working_dir.empty()) {
				file_name = options.working_prefix + ".softclips.fq.gz."	+ str_block_id;
			}
			if(!boost::filesystem::exists(file_name)) {
				break;
			}
			local_file_names.push_back(file_name);
			++n_blocks;
		}
		string file_name = options.prefix + ".softclips.fq.gz.bak";
		if(!options.working_dir.empty()) {
			file_name = options.working_prefix + ".softclips.fq.gz.bak";
		}
		if(boost::filesystem::exists(file_name)) {
			local_file_names.push_back(file_name);
		}
	});
	tasks.push_back([&]{
			auto& local_file_names = file_name_lists[1];
			int64_t n_blocks = 0;
			while(true) {
				string str_block_id = boost::lexical_cast<string>(n_blocks);
				string file_name = options.prefix + ".sr.1.fq.gz."	+ str_block_id;
				if(!options.working_dir.empty()) {
					file_name = options.working_prefix + ".sr.1.fq.gz."	+ str_block_id;
				}
				if(!boost::filesystem::exists(file_name)) {
					break;
				}
				local_file_names.push_back(file_name);
				++n_blocks;
			}
			string file_name = options.prefix + ".sr.1.fq.gz.bak";
			if(!options.working_dir.empty()) {
				file_name = options.working_prefix + ".sr.1.fq.gz.bak";
			}
			if(boost::filesystem::exists(file_name)) {
				local_file_names.push_back(file_name);
			}
		});
	tasks.push_back([&]{
				auto& local_file_names = file_name_lists[2];
				int64_t n_blocks = 0;
				while(true) {
					string str_block_id = boost::lexical_cast<string>(n_blocks);
					string file_name = options.prefix + ".sr.2.fq.gz."	+ str_block_id;
					if(!options.working_dir.empty()) {
						file_name = options.working_prefix + ".sr.2.fq.gz."	+ str_block_id;
					}
					if(!boost::filesystem::exists(file_name)) {
						break;
					}
					local_file_names.push_back(file_name);
					++n_blocks;
				}
				string file_name = options.prefix + ".sr.2.fq.gz.bak";
				if(!options.working_dir.empty()) {
					file_name = options.working_prefix + ".sr.2.fq.gz.bak";
				}
				if(boost::filesystem::exists(file_name)) {
					local_file_names.push_back(file_name);
				}
			});
	tasks.push_back([&]{
				auto& local_file_names = file_name_lists[3];
				int64_t n_blocks = 0;
				while(true) {
					string str_block_id = boost::lexical_cast<string>(n_blocks);
					string file_name = options.prefix + ".unmapped.fq.gz."	+ str_block_id;
					if(!options.working_dir.empty()) {
						file_name = options.working_prefix + ".unmapped.fq.gz."	+ str_block_id;
					}
					if(!boost::filesystem::exists(file_name)) {
						break;
					}
					local_file_names.push_back(file_name);
					++n_blocks;
				}
				string file_name = options.prefix + ".unmapped.fq.gz.bak";
				if(!options.working_dir.empty()) {
					file_name = options.working_prefix + ".unmapped.fq.gz.bak";
				}
				if(boost::filesystem::exists(file_name)) {
					local_file_names.push_back(file_name);
				}
			});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	vector<string> file_names;
	for(uint64_t l_id = 0; l_id < 4; ++l_id) {
		auto& local_file_names = file_name_lists[l_id];
		file_names.insert(file_names.end(), local_file_names.begin(), local_file_names.end());
	}
	castle::IOUtils::remove_files(file_names, n_cores);
}

} /* namespace meerkat */
