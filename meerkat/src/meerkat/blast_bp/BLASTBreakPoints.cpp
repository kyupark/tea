/*
 * BLASTBreakPoints.cpp
 *
 *  Created on: Jul 18, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#include "BLASTBreakPoints.hpp"

namespace meerkat {

BLASTBreakPoints::BLASTBreakPoints() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

BLASTBreakPoints::~BLASTBreakPoints() {
}

void BLASTBreakPoints::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}

//#	file format of prefix.sr.intra.refined
//#	del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, homology sizes
//#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion && insertion
//#	del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion && deletion (2 col), homology at break points
//#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, rchr (donor), ange of insertion (2 col), insert size, distance of deletion && insertion, homology at break points
//#	invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
//#	invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
//#	tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size, homology at break points
//#	file format of prefix.sr.inter.refine
//#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
//#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
//#	transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster, homology at break points

void BLASTBreakPoints::find_break_points() {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.find_break_points");
	checker.start();
	vector<function<void()> > tasks;
	string prefix = options.prefix;
	string db = options.reference_path;
	string unmappedfile = prefix + ".unmapped.fq.gz";
	string softclipfile = prefix + ".softclips.fq.gz";
	string blastdir = prefix + ".blast/";
	string srdir = prefix + ".sr/";

	string intra_filter = prefix + ".sr.intra.filtered";
	string intra_refine = prefix + ".intra.refined";
	string inter_filter = prefix + ".sr.inter.filtered";
	string inter_refine = prefix + ".inter.refined";
//	string intra_refine_coord = prefix + ".intra.refined.crd.sorted";
//	string inter_refine_coord = prefix + ".inter.refined.crd.sorted";
	string intra_refine_type = prefix + ".intra.refined.typ.sorted";
	string inter_refine_type = prefix + ".inter.refined.typ.sorted";

	if(!options.working_dir.empty()) {
		unmappedfile = options.working_prefix + ".unmapped.fq.gz";
		softclipfile = options.working_prefix + ".softclips.fq.gz";
		blastdir = options.working_prefix + ".blast/";
		srdir = options.working_prefix + ".sr/";
		intra_filter = options.working_prefix + ".sr.intra.filtered";
		intra_refine = options.working_prefix + ".intra.refined";
		inter_filter = options.working_prefix + ".sr.inter.filtered";
		inter_refine = options.working_prefix + ".inter.refined";
		intra_refine_type = options.working_prefix + ".intra.refined.typ.sorted";
		inter_refine_type = options.working_prefix + ".inter.refined.typ.sorted";
	}

	tasks.push_back([&] {
		if(boost::filesystem::exists(blastdir)) {
			boost::filesystem::remove_all(blastdir);
		}
		boost::filesystem::create_directories(blastdir);
	});

	tasks.push_back([&] {
		if(boost::filesystem::exists(srdir)) {
			boost::filesystem::remove_all(srdir);
		}
		boost::filesystem::create_directories(srdir);
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//#	generate query sequences

//#	%readsdetail: seq of each read. readsdetail{breakpoint name}{readname} = seq
	boost::unordered_map<string, string> readslist;
	boost::unordered_map<string, boost::unordered_map<string, string>> readsdetail;
	collect_cluster_name_read_name_map(readslist);
	collect_all_local_reads_detail(readsdetail, readslist);

	string n(100, 'N');
//#	generate subject sequences for intra chr events
	process_intra_chromosomal_events(readslist, readsdetail, n);
//	process_intra_chromosomal_events_alt(readslist, readsdetail, n);

//#	generate subject sequences for inter chr events
	process_inter_chromosomal_events(readslist, readsdetail, n);
//	process_inter_chromosomal_events_alt(readslist, readsdetail, n);

//#system "sort intra_refine -k 5,5 -k 6,6n > intra_refine_coord";
//#system "sort inter_refine -k 5,5 -k 6,6n > inter_refine_coord";

	tasks.push_back([&] {
		string sort_1_cmd = (boost::format("sort %s -k 1,1 -k 5,5 -k 6,6n > %s") % intra_refine % intra_refine_type).str();
		system(sort_1_cmd.c_str());
	});
	tasks.push_back([&] {
		string sort_2_cmd = (boost::format("sort %s -k 1,1 -k 5,5 -k 6,6n > %s") % inter_refine % inter_refine_type).str();
		system(sort_2_cmd.c_str());
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//#	system "rm intra_refine";
//#	system "rm inter_refine";
//	boost::filesystem::remove_all(blastdir);
	cout << checker;
}

void BLASTBreakPoints::collect_cluster_name_read_name_map(boost::unordered_map<string, string>& a_read_map) {
	string bp_readsfile = options.prefix + ".bp_reads";
	if(!options.working_dir.empty()) {
		bp_readsfile = options.working_prefix + ".bp_reads";
	}
	const char * delim_tab = "\t";
	vector<string> reads_name;
	string line;
	ifstream BPREADSRF(bp_readsfile, ios::binary);
	while (getline(BPREADSRF, line, '\n')) {
		string cluster_name = line;
		getline(BPREADSRF, line, '\n');
		castle::StringUtils::c_string_multi_split(line, delim_tab, reads_name);
//		cout << "[BLASTBreakPoints.collect_cluster_name_read_name_map] cluster: " << cluster_name << "\n";
//		cout << "[BLASTBreakPoints.collect_cluster_name_read_name_map] read names: " << line << "\n";
		for (auto a_read_name : reads_name) {
			a_read_map[a_read_name] = cluster_name;
		}
	}
}

void BLASTBreakPoints::collect_all_local_reads_detail(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.collect_all_local_reads_detail");
	checker.start();
//	string reads_detail_file = prefix + ".blast.target.reads";
//	if (boost::filesystem::exists(reads_detail_file)) {
//		const char* delim_tab = "\t";
//		string line;
//		vector<string> data;
//		ifstream in(reads_detail_file, ios::binary);
//		while (getline(in, line, '\n')) {
//			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
//			a_reads_detail[data[0]][data[1]] = data[2];
//		}
//		cout << checker;
//		return;
//	}

	vector<function<void()> > tasks;
	string unmappedfile = options.prefix + ".unmapped.fq.gz";
	string softclipfile = options.prefix + ".softclips.fq.gz";
	if(!options.working_dir.empty()) {
		unmappedfile = options.working_prefix + ".unmapped.fq.gz";
		softclipfile = options.working_prefix + ".softclips.fq.gz";
	}
	vector<boost::unordered_map<string, boost::unordered_map<string, string>>> local_readsdetail_list(2);
	collect_read_name_sequence_map(local_readsdetail_list[0], readslist, unmappedfile);
	collect_read_name_sequence_map(local_readsdetail_list[1], readslist, softclipfile);
	{
//		ofstream out(reads_detail_file, ios::binary);
		for (int64_t l_id = 0; l_id < 2; ++l_id) {
			auto& a_local_readsdetail = local_readsdetail_list[l_id];
			for (auto an_entry : a_local_readsdetail) {
				auto& the_first_key = an_entry.first;
				for (auto second_entry : an_entry.second) {
					auto& the_second_key = second_entry.first;
					auto& the_value = second_entry.second;
					a_reads_detail[the_first_key][the_second_key] = the_value;
//					out << the_first_key << "\t" << the_second_key << "\t" << the_value << "\n";
				}
			}
		}
	}
	cout << checker;
}

void BLASTBreakPoints::collect_read_name_sequence_map(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist, const string& a_gz_file_name) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.collect_read_name_sequence_map");
	checker.start();
	vector<function<void()> > tasks;
	BWACaller bc;
	bc.set_n_cores(n_cores);
	string input_name = a_gz_file_name + ".bak";
	int64_t file_size = castle::IOUtils::get_file_size(input_name);
	int64_t block_size = file_size / (double) n_cores;
	int64_t n_blocks = file_size / (double) block_size;
	vector<uint64_t> fq_blocks;
	castle::IOUtils::find_fastq_skip_points(fq_blocks, input_name, block_size, file_size, n_blocks, n_cores);
//	int64_t n_blocks = bc.split_FASTQ_single(input_name);
//	int64_t n_blocks = 0;
//	while (true) {
//		string str_block_id = boost::lexical_cast<string>(n_blocks);
//		string a_file_name = input_name + "." + str_block_id;
//		if (!boost::filesystem::exists(a_file_name)) {
//			break;
//		}
//		++n_blocks;
//	}
	cout << (boost::format("[BLASTBreakPoints.collect_read_name_sequence_map] # blocks: %d\n") % n_blocks).str();
	vector<boost::unordered_map<string, boost::unordered_map<string, string>>> reads_detail_list(n_blocks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			int64_t cur_pos = fq_blocks[block_id];
			int64_t next_pos = fq_blocks[block_id + 1];
			string str_block_id = boost::lexical_cast<string>(block_id);
			string line;
			auto& local_reads_detail = reads_detail_list[block_id];
			ifstream UMFQRF(input_name, ios::binary);
			if(0 != cur_pos) {
				UMFQRF.seekg(cur_pos, ios::beg);
			}
			while (getline(UMFQRF, line, '\n')) {
				cur_pos += line.size() + 1;
				string readname = line.substr(1);
				// get sequence
				getline(UMFQRF, line, '\n');
				int64_t length = line.size();
				cur_pos += line.size() + 1;
//				auto a_cluster_name_itr = readslist.find(readname);
//				if (readslist.end() != a_cluster_name_itr) {
//					local_reads_detail[a_cluster_name_itr->second][readname] = line;
//					//		#print "readslist{readname}\treadname\treadsdetail{readslist{readname}}{readname}\n";
//				} else {
				readname += "_" + boost::lexical_cast<string>(length);
				auto an_alt_cluster_name_itr = readslist.find(readname);
				if (readslist.end() != an_alt_cluster_name_itr) {
//					if(0 == block_id) {
//						cout << "[BLASTBreakPoints.collect_read_name_sequence_map] readname: " << readname << "\n";
//						cout << "[BLASTBreakPoints.collect_read_name_sequence_map] seq: " << line << "\n";
//					}
//					cout << readname << "\n";
					local_reads_detail[an_alt_cluster_name_itr->second][readname] = line;
				}
//				}

				// ignore '+' && quality values
				getline(UMFQRF, line, '\n');
				cur_pos += line.size() + 1;
				getline(UMFQRF, line, '\n');
				cur_pos += line.size() + 1;
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		auto& local_reads_detail = reads_detail_list[block_id];
		for (auto first_entry : local_reads_detail) {
			auto& first_key = first_entry.first;
			for (auto second_entry : first_entry.second) {
				auto& second_key = second_entry.first;
				auto& the_value = second_entry.second;
				a_reads_detail[first_key][second_key] = the_value;
			}
		}
	}
	cout << checker;
}

void BLASTBreakPoints::collect_read_name_sequence_map_alt(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist, const string& a_gz_file_name) {
	string line;
//	ifstream a_file(a_gz_file_name, ios::in | ios::binary);
//	boost::iostreams::filtering_istream UMFQRF;
//	UMFQRF.push(boost::iostreams::gzip_decompressor());
//	UMFQRF.push(a_file);
	igzstream UMFQRF(a_gz_file_name.c_str());

	while (getline(UMFQRF, line, '\n')) {
		string readname = line.substr(1);
		// get sequence
		getline(UMFQRF, line, '\n');
		int64_t length = line.size();
		readname += "_" + boost::lexical_cast<string>(length);
		auto a_cluster_name_itr = readslist.find(readname);
		if (readslist.end() != a_cluster_name_itr) {
			a_reads_detail[a_cluster_name_itr->second][readname] = line;
			//		#print "readslist{readname}\treadname\treadsdetail{readslist{readname}}{readname}\n";
		}
		// ignore '+' && quality values
		getline(UMFQRF, line, '\n');
		getline(UMFQRF, line, '\n');
	}
}

void BLASTBreakPoints::process_intra_chromosomal_events(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.process_intra_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	string intra_filter = options.prefix + ".sr.intra.filtered";
	if(!options.working_dir.empty()) {
		intra_filter = options.working_prefix + ".sr.intra.filtered";
	}
	vector<int64_t> block_positions;

	{
		block_positions.push_back(0);
		int64_t n_lines = castle::IOUtils::get_number_of_lines(intra_filter);
		int64_t BLOCK_SIZE = n_lines / (double) n_cores;
		string line;
		ifstream in(intra_filter);
		int64_t cur_pos = 0;
		int64_t n_reads = 0;
		int64_t n_total_reads = 0;
		while (getline(in, line, '\n')) {
			++n_reads;
			++n_total_reads;
			cur_pos += line.size() + 1;
			if (n_reads > BLOCK_SIZE) {
				block_positions.push_back(cur_pos);
				n_reads = 0;
			}
		}
		block_positions.push_back(numeric_limits<int64_t>::max());
	}
	int64_t n_blocks = block_positions.size() - 1;
	cout << (boost::format("[BLASTBreakPoints.process_intra_chromosomal_events] # blocks: %d\n") % n_blocks).str();
	vector<vector<EventEntry>> result_sr_list(n_blocks);
	vector<vector<EventEntry>> unsupport_del_list(n_blocks);
	vector<string> output_file_names(n_blocks);

	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string intra_refine = options.prefix + ".intra.refined." + str_block_id;
			string blast_dir = options.prefix + ".blast/" + str_block_id + "/";
			string sr_dir = options.prefix + ".sr/" + str_block_id + "/";
			if(!options.working_dir.empty()) {
				intra_refine = options.working_prefix + ".intra.refined." + str_block_id;
				blast_dir = options.working_prefix + ".blast/" + str_block_id + "/";
				sr_dir = options.working_prefix + ".sr/" + str_block_id + "/";
			}
			output_file_names[block_id] = intra_refine;
			int64_t cur_pos = block_positions[block_id];
			int64_t next_pos = block_positions[block_id + 1];
			string line;
			const char* delim_tab = "\t";
			const char* delim_slash = "/";
			vector<string> data;
			vector<string> mpd_id;

			IndexedFasta fai(options.reference_path.c_str());
			boost::filesystem::create_directories(blast_dir);
			boost::filesystem::create_directories(sr_dir);

			ofstream local_out(intra_refine, ios::binary);
			ifstream local_reader(intra_filter, ios::binary);
			local_reader.seekg(cur_pos, ios::beg);
			while(getline(local_reader, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[1], delim_slash, mpd_id);
				if ("del" == data[0]) {
					write_del(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if (string::npos != data[0].find("inssu")) {
					write_inssu(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if (string::npos != data[0].find("inssd")) {
					write_inssd(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if (string::npos != data[0].find("insou")) {
					write_insou(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if (string::npos != data[0].find("insod")) {
					write_insod(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("invers" == data[0]) {
					write_invers(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("del_invers" == data[0]) {
					write_del_invers(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("tandem_dup" == data[0]) {
					write_tandem_dup(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("invers_f" == data[0]) {
					write_invers_f(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("invers_r" == data[0]) {
					write_invers_r(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string intra_refine = options.prefix + ".intra.refined";
	if(!options.working_dir.empty()) {
		intra_refine = options.working_prefix + ".intra.refined";
	}
	ofstream INTRAREFINE(intra_refine, ios::out | ios::binary);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		INTRAREFINE << castle::IOUtils::read_fully(output_file_names[block_id]);
	}
	castle::IOUtils::remove_files(output_file_names, n_cores);
	cout << checker;
}

void BLASTBreakPoints::process_intra_chromosomal_events_alt(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.generate_intra_chromosomal_subject_sequences");
	checker.start();
	string intra_filter = options.prefix + ".sr.intra.filtered";
	string intra_refine = options.prefix + ".intra.refined";
	string blast_dir = options.prefix + ".blast/";
	if(!options.working_dir.empty()) {
		intra_filter = options.working_prefix + ".sr.intra.filtered";
		intra_refine = options.working_prefix + ".intra.refined";
		blast_dir = options.working_prefix + ".blast/";
	}
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	string line;
	vector<string> data;
	vector<string> mpd_id;
	int64_t i = 0;
//	int64_t selected_start = numeric_limits<int64_t>::max();
//	int64_t selected_end = selected_start;
//	int64_t selected_start = 0;
//	int64_t selected_end = numeric_limits<int64_t>::max();
//	int64_t selected_start = 0;
//	int64_t selected_end = selected_start + 100;
//	int64_t n_extra_selected_elements = 2;

	ofstream INTRAREFINE(intra_refine, ios::out | ios::binary);
	ifstream SRINTRAFTRF(intra_filter, ios::in | ios::binary);
	IndexedFasta fai(options.reference_path.c_str());
	while (getline(SRINTRAFTRF, line, '\n')) {
		++i;
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, mpd_id);
//		if ("del_insou" == data[0] && "20856_0/27077_0" == data[1]) {
//			selected_start = i;
//			selected_end = selected_start + n_extra_selected_elements;
//		} else {
//			if(0 == selected_start) {
//				continue;
//			}
//		}
//		if (i < selected_start) {
//			continue;
//		}
//
//		if (i > selected_end) {
//			break;
//		}
		cout << "intra: " << line << "\n";

		if ("del" == data[0]) {
			write_del(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if (string::npos != data[0].find("inssu")) {
			write_inssu(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if (string::npos != data[0].find("inssd")) {
			write_inssd(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if (string::npos != data[0].find("insou")) {
			write_insou(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if (string::npos != data[0].find("insod")) {
			write_insod(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("invers" == data[0]) {
			write_invers(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("del_invers" == data[0]) {
			write_del_invers(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("tandem_dup" == data[0]) {
			write_tandem_dup(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("invers_f" == data[0]) {
			write_invers_f(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("invers_r" == data[0]) {
			write_invers_r(blast_dir, INTRAREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		}
//		break;
	}
	cout << checker;
}

void BLASTBreakPoints::write_del(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;

	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_end + ref_is["rlu"]["selected"];

	const string& ref_id = data[4];
//	const bool debug = "chr21" == ref_id && ("24474417" == cur_start_str);
//	const bool debug = "chr21" == ref_id && ("24474417" == cur_start_str || "14051956" == cur_start_str);
	const bool debug = false;
	if(debug) {
		cout << "del here-0-a: positions: " << position1 << "/" << position2 << "/" << position3 << "/" << position4 << "/" << "\n";
	}
	string seq;
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_end_str;
	int64_t left_bound_blast = 0;
	int64_t right_bound_blast = 0;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		if(debug) {
			cout << "del here-0\n";
		}
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		if(debug) {
			cout << "del here-1\n";
		}
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n + fai.fetch_1_system(ref_id, position3, position4);
	}

	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	if(debug) {
		cout << "del here-2: " << bltoutfile << "\n";
	}
	left_bound_blast = position1 + ref_bps[0] - 1;
	if (is_covered) {
		if(debug) {
			cout << "del here-3\n";
		}
		right_bound_blast = position1 + ref_bps[1] - 1;
	} else {
		if(debug) {
			cout << "del here-4\n";
		}
		right_bound_blast = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}
	if (!ref_bps[0]) {
		if(debug) {
			cout << "del here-5\n";
		}
		left_bound_blast = cur_start;
	}
	if (!ref_bps[1]) {
		if(debug) {
			cout << "del here-6\n";
		}
		right_bound_blast = cur_end;
	}
	if (!ref_bps[0]) {
		if(debug) {
			cout << "del here-7\n";
		}
		ref_bps[2] = -1000;
	}
	int64_t eventsizedel = right_bound_blast - left_bound_blast - 1;
	if (eventsizedel > 0) {
		if (ref_bps[2] > 0) {
			if(debug) {
				cout << "del here-8: " << ref_bps[2] << "\n";
			}
			INTRAREFINE << (boost::format("del_ins\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t-\t-\t-\t%s\n")
			% data[1] % data[2] % data[3] % data[4] % left_bound_blast % right_bound_blast % eventsizedel % ref_bps[2]).str();
		}
		else {
			if(debug) {
				cout << "del here-9: " << ref_bps[2] << "\n";
			}
			int64_t homology = -ref_bps[2];
			if(!ref_bps[2]) {
				homology = 0;
			}
			INTRAREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
			% data[0] % data[1] % data[2] % data[3] % data[4] % left_bound_blast % right_bound_blast % eventsizedel % homology).str();
		}
	}
}

void BLASTBreakPoints::write_inssu(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	//cout << "inssu here-0\n";
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];
	string seq;
	//cout << "inssu here-1\n";
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		seq = fai.fetch_1_system(ref_id, position1, position4);
		//cout << "inssu here-2: seq_size: " << seq.size() << "\n";
	} else {

		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position3, position4);
		//cout << "inssu here-3: seq_size: " << seq.size() << "\n";
	}
	//cout << "inssu here-4\n";
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_mate_start_str;
	//cout << "inssu here-5\n";
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, false);
	//cout << "inssu here-6\n";
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
//		my (ref_bps, ref_bps2, left_bound_blast1, right_bound_blast1);
	//cout << "inssu here-7\n";
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "inssu here-8\n";
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "inssu here-9\n";
		right_bound_blast1 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "inssu here-10\n";
		right_bound_blast1 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}
	if (!ref_bps[0]) {
		//cout << "inssu here-11\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "inssu here-12\n";
		right_bound_blast1 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		//cout << "inssu here-13\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "inssu here-15\n";
		homology1 = 0;
	}
	ref_bps.clear();
	ref_bps2.clear();
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_end - ref_is["rlu"]["selected"];
	position4 = cur_mate_end + ref_is["rlu"]["selected"];
	is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "inssu here-16-a: positions: " << position1 << "," << position4 << "\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);

		//cout << "inssu here-16-b: seq_size: " << seq.size() << "\n";
	} else {
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position3, position4);
		//cout << "inssu here-17: seq_size: " << seq.size() << "\n";
	}
	name = ref_id + "__" + cur_end_str + "__" + ref_id + "__" + cur_mate_end_str;
	if (mpd_id.size() > 1) {
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, false);
	} else {
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, false);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);

	//cout << "inssu here-18\n";
	left_bound_blast2 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "inssu here-19\n";
		right_bound_blast2 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "inssu here-20\n";
		right_bound_blast2 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}
	if (!ref_bps[0]) {
		//cout << "inssu here-21\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "inssu here-22\n";
		right_bound_blast2 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "inssu here-23\n";
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "inssu here-25\n";
		homology2 = 0;
	}

	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast2 - right_bound_blast1 + 1;
	int64_t distance = left_bound_blast1 - right_bound_blast2;

	if (eventsizedel <= 10) {
		//cout << "inssu here-26\n";
		INTRAREFINE << (boost::format("inssu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[4] % right_bound_blast1 % right_bound_blast2 % eventsizeins % distance % homology1 % homology2).str();
	} else {
		//cout << "inssu here-27\n";
		INTRAREFINE << (boost::format("del_inssu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[4] % right_bound_blast1 % right_bound_blast2 % eventsizeins % distance % homology1 % homology2).str();
	}

}

void BLASTBreakPoints::write_inssd(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];

	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	string name;
	string bltoutfile;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
	string seq;

	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "inssd here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "inssd here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position3, position4);
	}
	name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_mate_start_str;
	if (mpd_id.size() > 1) {
		//cout << "inssd here-2\n";
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, true);
	} else {
		//cout << "inssd here-3\n";
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, true);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "inssd here-4\n";
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "inssd here-5\n";
		right_bound_blast1 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "inssd here-6\n";
		right_bound_blast1 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}

//		#print "ref_bps\n";
	if (!ref_bps[0]) {
		//cout << "inssd here-7\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "inssd here-8\n";
		right_bound_blast1 = cur_mate_start;
	}
	if ((ref_bps.size() > 0 && !ref_bps[0])) {
		//cout << "inssd here-9\n";
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];

	homology2 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "inssd here-11\n";
		homology2 = 0;
	}

	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_end - ref_is["rlu"]["selected"];
	position4 = cur_mate_end + ref_is["rlu"]["selected"];
	is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "inssd here-12\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "inssd here-13\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position3, position4);
	}
	name = ref_id + "__" + cur_end_str + "__" + ref_id + "__" + cur_mate_end_str;
	bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
//		my (ref_bps, ref_bps2, left_bound_blast2, right_bound_blast2);
	ref_bps.clear();
	ref_bps2.clear();
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "inssd here-14\n";
	left_bound_blast2 = position1 + ref_bps[0] - 1;
	if (is_covered) {

		//cout << "inssd here-15\n";
		right_bound_blast2 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "inssd here-16\n";
		right_bound_blast2 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}

//		#print "ref_bps\n";
	if (!ref_bps[0]) {
		//cout << "inssd here-17\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "inssd here-18\n";
		right_bound_blast2 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "inssd here-19\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "inssd here-21\n";
		homology1 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast2 - right_bound_blast1 + 1;
	int64_t distance = right_bound_blast1 - left_bound_blast2;
	if (eventsizedel <= 10) {
		//cout << "inssd here-22\n";
		INTRAREFINE << (boost::format("inssd\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[4] % right_bound_blast1 % right_bound_blast2 % eventsizeins % distance % homology1 % homology2).str();
	} else {
		//cout << "inssd here-23\n";
		INTRAREFINE
				<< (boost::format("del_inssd\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
						data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
						data[4] % right_bound_blast1 % right_bound_blast2 % eventsizeins % distance % homology1 % homology2).str();
	}

}

void BLASTBreakPoints::write_insou(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];
	//cout << "insou here-0a: positions: " << position1 << "/" << position2 << "/" << position3 << "/" << position4 << "\n";
	string seq;
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "insou here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "insou here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_mate_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, false);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "insou here-2\n";
		right_bound_blast1 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "insou here-3\n";
		right_bound_blast1 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		//cout << "insou here-4\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "insou here-5\n";
		right_bound_blast1 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "insou here-6\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];
	if (!ref_bps[2]) {
		homology1 = 0;
	}
	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_start - ref_is["rlu"]["selected"];
	position4 = cur_mate_start + ref_is["rlu"]["selected"];
	ref_bps.clear();
	ref_bps2.clear();
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "insou here-7\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "insou here-8\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	name = ref_id + "__" + cur_end_str + "__" + ref_id + "__" + cur_mate_start_str;
	if (mpd_id.size() > 1) {
		//cout << "insou here-9\n";
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, false);
	} else {
		//cout << "insou here-10\n";
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, false);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);

	left_bound_blast2 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "insou here-11\n";
		right_bound_blast2 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "insou here-12\n";
		right_bound_blast2 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	//cout << "insou here-12-a: positions: " << position1 << "/" << position2 << "/" << position3 << "/" << position4 << "\n";
	if (!ref_bps[0]) {
		//cout << "insou here-13\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "insou here-14\n";
		right_bound_blast2 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		//cout << "insou here-15\n";
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "insou here-16\n";
		homology2 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast1 - right_bound_blast2 + 1;
	int64_t distance = left_bound_blast1 - right_bound_blast1;

	if (eventsizedel <= 10) {
		if (eventsizedel <= 2 && abs(eventsizeins) <= 2) {
			if (right_bound_blast1 > left_bound_blast2) {
				//cout << "insou here-17\n";
				int64_t eventsizeinv = right_bound_blast1 - left_bound_blast2;
				INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
						data[2] % data[3] % data[4] % left_bound_blast2 % right_bound_blast1 % eventsizeinv % homology1 % homology2).str();
			} else {
				//cout << "insou here-18\n";
				int64_t eventsizeinv = left_bound_blast1 - right_bound_blast2;
				INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1]%
						data[2] % data[3] % data[4] % right_bound_blast2 % left_bound_blast1 % eventsizeinv % homology1 % homology2).str();
			}
		} else {
			//cout << "insou here-19\n";
			INTRAREFINE << (boost::format("insou\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
					data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
					data[4] % right_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
		}
	} else {
		//cout << "insou here-20\n";
		INTRAREFINE << (boost::format("del_insou\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[4] % right_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
	}

}

void BLASTBreakPoints::write_insod(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];

	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
	string seq;

	bool is_covered = covered(position1, position2, position3, position4);

	if (is_covered) {
		//cout << "insod here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "insod here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_mate_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "insod here-2\n";
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "insod here-3\n";
		right_bound_blast1 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "insod here-4\n";
		right_bound_blast1 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		//cout << "insod here-5\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "insod here-6\n";
		right_bound_blast1 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "insod here-7\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "insod here-8\n";
		homology1 = 0;
	}

	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_start - ref_is["rlu"]["selected"];
	position4 = cur_mate_start + ref_is["rlu"]["selected"];
	is_covered = covered(position1, position2, position3, position4);
	ref_bps.clear();
	ref_bps2.clear();
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	if (is_covered) {
		//cout << "insod here-9\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "insod here-10\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	name = ref_id + "__" + cur_end_str + "__" + ref_id + "__" + cur_mate_start_str;
	if (mpd_id.size() > 1) {
		//cout << "insod here-11\n";
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, true);
	} else {
		//cout << "insod here-12\n";
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, true);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);

	left_bound_blast2 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "insod here-13\n";
		right_bound_blast2 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "insod here-14\n";
		right_bound_blast2 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		//cout << "insod here-15\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "insod here-16\n";
		right_bound_blast2 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		//cout << "insod here-17\n";
		ref_bps[2] = -1000;
	}

	int64_t homology2 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "insod here-18\n";
		homology2 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast1 - right_bound_blast2 + 1;
	int64_t distance = right_bound_blast2 - left_bound_blast2;

	if (eventsizedel <= 10) {
		if (eventsizedel <= 2 && abs(eventsizeins) <= 2) {
			if (right_bound_blast1 > left_bound_blast2) {
				//cout << "insod here-19\n";
				int64_t eventsizeinv = right_bound_blast1 - left_bound_blast2;
				INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
						data[2] % data[3] % data[4] % left_bound_blast2 % right_bound_blast1 % eventsizeinv % homology1 % homology2).str();
			} else {
				//cout << "insod here-20\n";
				int64_t eventsizeinv = left_bound_blast1 - right_bound_blast2;
				INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
						data[2] % data[3] % data[4] % right_bound_blast2 % left_bound_blast1 % eventsizeinv % homology1 % homology2).str();
			}
		} else {
			//cout << "insod here-21\n";
			INTRAREFINE << (boost::format("insod\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
					data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
					data[4] % right_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
		}
	} else {
		//cout << "insod here-22\n";
		INTRAREFINE << (boost::format("del_insod\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[4] % right_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
	}

}

void BLASTBreakPoints::write_invers(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	auto& cut_sr = options.cut_sr;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
//			int64_t mate_start = boost::lexical_cast<int64_t>(data[10]);
//			int64_t mate_end = boost::lexical_cast<int64_t>(data[11]);
	int64_t position1 = cur_start - (ref_is["rlu"]["selected"] + cut_sr);
	int64_t position2 = cur_start + (ref_is["rlu"]["selected"] + cut_sr);
	int64_t position3 = cur_end - (ref_is["rlu"]["selected"] + cut_sr);
	int64_t position4 = cur_end + (ref_is["rlu"]["selected"] + cut_sr);
	const string& ref_id = data[4];
	string seq;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "invers here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "invers here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	if (ref_bps2.size() < 3) {
		//cout << "invers here-2\n";
		ref_bps2.resize(3);
	}
	if (is_covered) {
		//cout << "invers here-3\n";
		left_bound_blast1 = position1 + ref_bps[0] - 1;
		right_bound_blast1 = position1 + ref_bps[1] - 1;
		left_bound_blast2 = position1 + ref_bps2[0] - 1;
		right_bound_blast2 = position1 + ref_bps2[1] - 1;
	} else {
		//cout << "invers here-4\n";
		left_bound_blast1 = position1 + ref_bps[0] - 1;
		right_bound_blast1 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
		left_bound_blast2 = position1 + ref_bps2[0] - 1;
		right_bound_blast2 = position4 - (ref_bps2[1] - (position2 - position1 + 100)) + 2;
	}

//		#print "@{ref_bps}\n@{ref_bps2}\n";
	if (!ref_bps[0]) {
		//cout << "invers here-5\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "invers here-6\n";
		right_bound_blast1 = cur_end;
	}
	if (!ref_bps[0]) {
		//cout << "invers here-7\n";
		ref_bps[2] = -1000;
	}
	int64_t homology = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "invers here-8\n";
		homology = 0;
	}
	//cout << "invers here-8a: position: " << position1 << "/" << position2 << "/" << position3 << "/" << position4 << "\n";

	//cout << "invers here-8b: ref_bps: " << ref_bps[0] << "/" << ref_bps[1] << "/" << ref_bps[2] << "\n";
	//cout << "invers here-8c: ref_bps2: " << ref_bps2[0] << "/" << ref_bps2[1] << "/" << ref_bps2[2] << "\n";
	if (ref_bps2[1] && abs(ref_bps2[0] - ref_bps[0]) > 10 && abs(ref_bps2[1] - ref_bps[1]) > 10) {
		if (left_bound_blast1 < left_bound_blast2) {
			//cout << "invers here-9\n";
			int64_t eventsizedel = right_bound_blast2 - left_bound_blast1 - 1;
			int64_t eventsizeinv = right_bound_blast1 - left_bound_blast2;
			int64_t distance1 = left_bound_blast2 - left_bound_blast1;
			int64_t distance2 = right_bound_blast2 - right_bound_blast1;
			INTRAREFINE << (boost::format("del_invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[1] %
					data[2] % data[3] % data[4] % left_bound_blast1 % right_bound_blast2 % eventsizedel % //
					data[4] % left_bound_blast2 % right_bound_blast1 % eventsizeinv % distance1 % distance2 % homology).str();
		} else {
			//cout << "invers here-10\n";
			int64_t eventsizedel = right_bound_blast1 - left_bound_blast2 - 1;
			int64_t eventsizeinv = right_bound_blast2 - left_bound_blast1;
			int64_t distance1 = left_bound_blast1 - left_bound_blast2;
			int64_t distance2 = right_bound_blast1 - right_bound_blast2;
			INTRAREFINE << (boost::format("del_invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[1] %
					data[2] % data[3] % data[4] % left_bound_blast2 % right_bound_blast1 % eventsizedel % //
					data[4] % left_bound_blast1 % right_bound_blast2 % eventsizeinv % distance1 % distance2 % homology).str();
		}
	} else {
		//cout << "invers here-11\n";
		int64_t eventsizeinv = right_bound_blast1 - left_bound_blast1;
		INTRAREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[0] % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % right_bound_blast1 % eventsizeinv % homology).str();
	}

}

void BLASTBreakPoints::write_del_invers(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];
	string seq;
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast1 = 0;
	int64_t right_bound_blast1 = 0;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "del invers here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "del invers here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_mate_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "del invers here-2\n";
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "del invers here-3\n";
		right_bound_blast1 = position1 + ref_bps[1] - 1;
	} else {
		//cout << "del invers here-4\n";
		right_bound_blast1 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		//cout << "del invers here-5\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "del invers here-6\n";
		right_bound_blast1 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "del invers here-7\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "del invers here-8\n";
		homology1 = 0;
	}
	int64_t left_bound_blast2 = 0;
	int64_t right_bound_blast2 = 0;
	ref_bps.clear();
	ref_bps2.clear();
	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_start - ref_is["rlu"]["selected"];
	position4 = cur_mate_start + ref_is["rlu"]["selected"];
	is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "del invers here-9\n";
		seq = fai.fetch_1_system(ref_id, position3, position2);
	} else {
		//cout << "del invers here-10\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	name = ref_id + "__" + cur_end_str + "__" + ref_id + "__" + cur_mate_start_str;
	if (mpd_id.size() > 1) {
		//cout << "del invers here-11\n";
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, true);
	} else {
//		//cout << "del invers here-12\n";
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, true);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);

	if (is_covered) {
		//cout << "del invers here-13\n";
		left_bound_blast2 = position3 + ref_bps[1] - 1;
		right_bound_blast2 = position3 + ref_bps[0] - 1;
	} else {
		//cout << "del invers here-14\n";
		left_bound_blast2 = position1 + ref_bps[0] - 1;
		right_bound_blast2 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		//cout << "del invers here-15\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "del invers here-16\n";
		right_bound_blast2 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		//cout << "del invers here-17\n";
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "del invers here-18\n";
		homology2 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeinv = right_bound_blast1 - right_bound_blast2;
	int64_t distance1 = right_bound_blast2 - left_bound_blast1;
	int64_t distance2 = left_bound_blast2 - right_bound_blast1;

//		#print "left_bound_blast1\tleft_bound_blast2\tright_bound_blast1\tright_bound_blast2\neventsizedel\teventsizeinv\tdistance1\tdistance2\n";
	if (distance1 <= 10 && distance2 <= 10) {
		//cout << "del invers here-19\n";
		INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % homology1 % homology2).str();
	} else if (eventsizedel <= 10 && eventsizeinv <= 10) {
		//cout << "del invers here-20\n";
		int64_t eventsizeinv = right_bound_blast1 - left_bound_blast1;
		INTRAREFINE << (boost::format("invers\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % right_bound_blast1 % eventsizeinv % homology1 % homology2).str();
	} else if (distance1 < -10) {
		//cout << "del invers here-21\n";
		int64_t eventsizedel = left_bound_blast2 - right_bound_blast1 - 1;
		int64_t eventsizeins = left_bound_blast1 - right_bound_blast2 + 1;
		int64_t distance = right_bound_blast1 - left_bound_blast1;
		if (eventsizedel <= 10) {
			//cout << "del invers here-22\n";
			INTRAREFINE
					<< (boost::format("insou\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
							data[2] % data[3] % data[4] % right_bound_blast1 % left_bound_blast2 % eventsizedel % //
							data[4] % right_bound_blast2 % left_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
		} else {
			//cout << "del invers here-23\n";
			INTRAREFINE
					<< (boost::format("del_insou\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
							data[2] % data[3] % data[4] % right_bound_blast1 % left_bound_blast2 % eventsizedel % //
							data[4] % right_bound_blast2 % left_bound_blast1 % eventsizeins % distance % homology1
							% homology2).str();
		}
	} else if (distance2 < -10) {
		//cout << "del invers here-24\n";
		int64_t eventsizedel = right_bound_blast2 - left_bound_blast1 - 1;
		int64_t eventsizeins = right_bound_blast1 - left_bound_blast2 + 1;
		int64_t distance = left_bound_blast2 - right_bound_blast2;
		if (eventsizedel <= 10) {
			//cout << "del invers here-25\n";
			INTRAREFINE
					<< (boost::format("insod\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
							data[2] % data[3] % data[4] % left_bound_blast1 % right_bound_blast2 % eventsizedel %//
							data[4] % left_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1 % homology2).str();
		} else {
			//cout << "del invers here-26\n";
			INTRAREFINE
					<< (boost::format("del_insod\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
							data[2] % data[3] % data[4] % left_bound_blast1 % right_bound_blast2 % eventsizedel % //
							data[4] % left_bound_blast2 % right_bound_blast1 % eventsizeins % distance % homology1
							% homology2).str();
		}
	} else {
		//cout << "del invers here-27\n";
		INTRAREFINE
				<< (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[0] % data[1] %
						data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
						data[4] % right_bound_blast2 % right_bound_blast1 % eventsizeinv % distance1 % distance2
						% homology1 % homology2).str();
	}

}

void BLASTBreakPoints::write_tandem_dup(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
//			int64_t mate_start = boost::lexical_cast<int64_t>(data[10]);
//			int64_t mate_end = boost::lexical_cast<int64_t>(data[11]);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];
	string seq;
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast = 0;
	int64_t right_bound_blast = 0;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position3, position4);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
//		my (ref_bps, ref_bps2, left_bound_blast, right_bound_blast);

	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	left_bound_blast = position1 + ref_bps[0] - 1;
	if (is_covered) {
		right_bound_blast = position1 + ref_bps[1] - 1;
	} else {
		right_bound_blast = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	}
	if (!ref_bps[0]) {
		left_bound_blast = cur_start;
	}
	if (!ref_bps[1]) {
		right_bound_blast = cur_end;
	}
	if (!ref_bps[0]) {
		ref_bps[2] = -1000;
	}
	int64_t homology = -ref_bps[2];

	if (!ref_bps[2]) {
		homology = 0;
	}

	int64_t eventsizetransl = right_bound_blast - left_bound_blast;
	INTRAREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[0] %
			data[1] % data[2] % data[3] % data[4] % left_bound_blast % right_bound_blast % eventsizetransl % homology).str();
}

void BLASTBreakPoints::write_invers_f(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];

	string seq;
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast = 0;
	int64_t right_bound_blast = 0;

	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		//cout << "invers_f here-0\n";
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		//cout << "invers_f here-1\n";
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	//cout << "invers_f here-2\n";
	left_bound_blast = position1 + ref_bps[0] - 1;
	if (is_covered) {
		//cout << "invers_f here-3\n";
		right_bound_blast = position1 + ref_bps[1] - 1;
	} else {
		//cout << "invers_f here-4\n";
		right_bound_blast = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}

	if (!ref_bps[0]) {
		//cout << "invers_f here-5\n";
		left_bound_blast = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "invers_f here-6\n";
		right_bound_blast = cur_end;
	}
	if (!ref_bps[0]) {
		//cout << "invers_f here-7\n";
		ref_bps[2] = -1000;
	}
	//cout << "invers_f here-8\n";
	int64_t homology = -ref_bps[2];
	if (!ref_bps[2]) {
		//cout << "invers_f here-9\n";
		homology = 0;
	}
	int64_t eventsizeinv = right_bound_blast - left_bound_blast;
	INTRAREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[0] %
			data[1] % data[2] % data[3] % data[4] % left_bound_blast % right_bound_blast % eventsizeinv % homology).str();
}

void BLASTBreakPoints::write_invers_r(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_end + ref_is["rlu"]["selected"];
	const string& ref_id = data[4];
	string seq;
	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;
	int64_t left_bound_blast = 0;
	int64_t right_bound_blast = 0;
	bool is_covered = covered(position1, position2, position3, position4);
	if (is_covered) {
		seq = fai.fetch_1_system(ref_id, position1, position4);
	} else {
		seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(ref_id, position4, position3);
	}
	string name = ref_id + "__" + cur_start_str + "__" + ref_id + "__" + cur_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);

	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	left_bound_blast = position1 + ref_bps[0] - 1;
	if (is_covered) {
		right_bound_blast = position1 + ref_bps[1] - 1;
	} else {
		right_bound_blast = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	}
	if (!ref_bps[0]) {
		left_bound_blast = cur_start;
	}
	if (!ref_bps[1]) {
		right_bound_blast = cur_end;
	}
	if (!ref_bps[0]) {
		ref_bps[2] = -1000;
	}
	int64_t homology = -ref_bps[2];

	if (!ref_bps[2]) {
		homology = 0;
	}
	int64_t eventsizeinv = right_bound_blast - left_bound_blast;
	INTRAREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data[0] % data[1] %
			data[2] % data[3] % data[4] % left_bound_blast % right_bound_blast % eventsizeinv % homology).str();
}

void BLASTBreakPoints::process_inter_chromosomal_events(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.process_inter_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	string intra_filter = options.prefix + ".sr.inter.filtered";
	if(!options.working_dir.empty()) {
		intra_filter = options.working_prefix + ".sr.inter.filtered";
	}
	vector<int64_t> block_positions;

	{
		block_positions.push_back(0);
		int64_t n_lines = castle::IOUtils::get_number_of_lines(intra_filter);
		int64_t BLOCK_SIZE = n_lines / (double) n_cores;
		string line;
		ifstream in(intra_filter);
		int64_t cur_pos = 0;
		int64_t n_reads = 0;
		int64_t n_total_reads = 0;
		while (getline(in, line, '\n')) {
			++n_reads;
			++n_total_reads;
			cur_pos += line.size() + 1;
			if (n_reads > BLOCK_SIZE) {
				block_positions.push_back(cur_pos);
				n_reads = 0;
			}
		}
		block_positions.push_back(numeric_limits<int64_t>::max());
	}
	int64_t n_blocks = block_positions.size() - 1;
	cout << (boost::format("[BLASTBreakPoints.process_inter_chromosomal_events] # blocks: %d\n") % n_blocks).str();
	vector<vector<EventEntry>> result_sr_list(n_blocks);
	vector<vector<EventEntry>> unsupport_del_list(n_blocks);
	vector<string> output_file_names(n_blocks);

	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string inter_refine = options.prefix + ".inter.refined." + str_block_id;
			string blast_dir = options.prefix + ".blast/" + str_block_id + "/";
			string sr_dir = options.prefix + ".sr/" + str_block_id + "/";
			if(!options.working_dir.empty()) {
				inter_refine = options.working_prefix + ".inter.refined." + str_block_id;
				blast_dir = options.working_prefix + ".blast/" + str_block_id + "/";
				sr_dir = options.prefix + ".sr/" + str_block_id + "/";
			}

			int64_t cur_pos = block_positions[block_id];
			int64_t next_pos = block_positions[block_id + 1];
			string line;
			const char* delim_tab = "\t";
			const char* delim_slash = "/";
			vector<string> data;
			vector<string> mpd_id;

			IndexedFasta fai(options.reference_path.c_str());
			output_file_names[block_id] = inter_refine;

			boost::filesystem::create_directories(blast_dir);
			boost::filesystem::create_directories(sr_dir);

			ofstream local_out(inter_refine, ios::binary);
			ifstream local_reader(intra_filter, ios::binary);
			local_reader.seekg(cur_pos, ios::beg);
			while(getline(local_reader, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[1], delim_slash, mpd_id);
				if (string::npos != data[0].find("inss")) {
					write_inss(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if (string::npos != data[0].find("inso")) {
					write_inso(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				} else if ("transl_inter" == data[0]) {
					write_transl_inter(blast_dir, local_out, data, mpd_id, fai, readslist, readsdetail, n);
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string inter_refine = options.prefix + ".inter.refined";
	if(!options.working_dir.empty()) {
		inter_refine = options.working_prefix + ".inter.refined";
	}
	ofstream INTERREFINE(inter_refine, ios::out | ios::binary);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		INTERREFINE << castle::IOUtils::read_fully(output_file_names[block_id]);
	}
	castle::IOUtils::remove_files(output_file_names, n_cores);
	cout << checker;
}

void BLASTBreakPoints::process_inter_chromosomal_events_alt(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	castle::TimeChecker checker;
	checker.setTarget("BLASTBreakPoints.generate_inter_chromosomal_subject_sequences");
	checker.start();
	string inter_filter = options.prefix + ".sr.inter.filtered";
	string inter_refine = options.prefix + ".inter.refined";
	string blast_dir = options.prefix + ".blast/";
	if(!options.working_dir.empty()) {
		inter_filter = options.working_prefix + ".sr.inter.filtered";
		inter_refine = options.working_prefix + ".inter.refined";
		blast_dir = options.working_prefix + ".blast/";
	}
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	string line;
	vector<string> data;
	vector<string> mpd_id;
//	bool debug = false;
	IndexedFasta fai(options.reference_path.c_str());

	ifstream SRINTERFTRF(inter_filter, ios::in | ios::binary);
	ofstream INTERREFINE(inter_refine, ios::out | ios::binary);
	while (getline(SRINTERFTRF, line, '\n')) {
		cout << "inter: " << line << "\n";

		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, mpd_id);
//		if("transl_inter" == data[0] && "31245_0" == data[1]) {
//			debug = true;
//		}
//		if(!debug) {
//			continue;
//		}
		if (string::npos != data[0].find("inss")) {
			write_inss(blast_dir, INTERREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if (string::npos != data[0].find("inso")) {
			write_inso(blast_dir, INTERREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		} else if ("transl_inter" == data[0]) {
			write_transl_inter(blast_dir, INTERREFINE, data, mpd_id, fai, readslist, readsdetail, n);
		}
	}
	cout << checker;
}

void BLASTBreakPoints::write_inss(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;
	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];
	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);
	const string& ref_id = data[4];
	const string& mate_ref_id = data[8];

	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;

	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];
	//cout << "inss here-0\n";
	string seq = fai.fetch_1_system(ref_id, position1, position2);
	seq += n;
	seq += fai.fetch_1_system(mate_ref_id, position3, position4);
	string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	int64_t left_bound_blast1 = 0;
	//cout << "inss here-1\n";
	left_bound_blast1 = position1 + ref_bps[0] - 1;
	int64_t right_bound_blast1 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	//cout << "inss here-2\n";
	if (!ref_bps[0]) {
		//cout << "inss here-3\n";
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		//cout << "inss here-4\n";
		right_bound_blast1 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		//cout << "inss here-5\n";
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "inss here-7\n";
		homology1 = 0;
	}
	//cout << "inss here-7a\n";
	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_end - ref_is["rlu"]["selected"];
	position4 = cur_mate_end + ref_is["rlu"]["selected"];
	ref_bps.clear();
	ref_bps2.clear();
	seq = fai.fetch_1_system(ref_id, position1, position2);
	seq += n;
	seq += fai.fetch_1_system(mate_ref_id, position3, position4);
	name = ref_id + "__" + cur_end_str + "__" + mate_ref_id + "__" + cur_mate_end_str;
	if (mpd_id.size() > 1) {
		//cout << "inss here-9\n";
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, true);
	} else {
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, true);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	int64_t left_bound_blast2 = position1 + ref_bps[0] - 1;
	//cout << "inss here-11\n";
	int64_t right_bound_blast2 = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
	//cout << "inss here-12\n";
	if (!ref_bps[0]) {
		//cout << "inss here-13\n";
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		//cout << "inss here-14\n";
		right_bound_blast2 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		//cout << "inss here-15\n";
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];

	if (!ref_bps[2]) {
		//cout << "inss here-17\n";
		homology2 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast2 - right_bound_blast1 + 1;

	if (eventsizedel <= 10) {
		//cout << "inss here-18\n";
		INTERREFINE << (boost::format("inss\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[8] % right_bound_blast1 % right_bound_blast2 % eventsizeins % homology1 % homology2).str();
	} else {
		//cout << "inss here-19\n";
		INTERREFINE << (boost::format("del_inss\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel %//
				data[8] % right_bound_blast1 % right_bound_blast2 % eventsizeins % homology1 % homology2).str();
	}
}

void BLASTBreakPoints::write_inso(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;

	const string& cur_start_str = data[5];
	const string& cur_end_str = data[6];
	const string& cur_mate_start_str = data[9];
	const string& cur_mate_end_str = data[10];

	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_end = boost::lexical_cast<int64_t>(cur_end_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(cur_mate_end_str);

	int64_t position1 = cur_start - ref_is["rlu"]["selected"];
	int64_t position2 = cur_start + ref_is["rlu"]["selected"];
	int64_t position3 = cur_mate_end - ref_is["rlu"]["selected"];
	int64_t position4 = cur_mate_end + ref_is["rlu"]["selected"];

	const string& ref_id = data[4];
	const string& mate_ref_id = data[8];

	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;

	string seq = fai.fetch_1_system(ref_id, position1, position2);
	seq += n;
	seq += fai.fetch_1_system(mate_ref_id, position4, position3);
	string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_end_str;
	string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
	int64_t left_bound_blast1 = position1 + ref_bps[0] - 1;
	int64_t right_bound_blast1 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	if (!ref_bps[0]) {
		left_bound_blast1 = cur_start;
	}
	if (!ref_bps[1]) {
		right_bound_blast1 = cur_mate_end;
	}
	if (!ref_bps[0]) {
		ref_bps[2] = -1000;
	}
	int64_t homology1 = -ref_bps[2];

	if (!ref_bps[2]) {
		homology1 = 0;
	}
	position1 = cur_end - ref_is["rlu"]["selected"];
	position2 = cur_end + ref_is["rlu"]["selected"];
	position3 = cur_mate_start - ref_is["rlu"]["selected"];
	position4 = cur_mate_start + ref_is["rlu"]["selected"];
	ref_bps.clear();
	ref_bps2.clear();
	seq = fai.fetch_1_system(ref_id, position1, position2);
	seq += n;
	seq += fai.fetch_1_system(mate_ref_id, position4, position3);
	name = ref_id + "__" + cur_end_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
	if (mpd_id.size() > 1) {
		bltoutfile = run_blast(blast_dir, data, mpd_id[1], name, seq, readsdetail, true);
	} else {
		string a_tmp_mpd_id = "s" + mpd_id[0];
		bltoutfile = run_blast(blast_dir, data, a_tmp_mpd_id, name, seq, readsdetail, true);
	}
	parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);

	int64_t left_bound_blast2 = position1 + ref_bps[0] - 1;
	int64_t right_bound_blast2 = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
	if (!ref_bps[0]) {
		left_bound_blast2 = cur_end;
	}
	if (!ref_bps[1]) {
		right_bound_blast2 = cur_mate_start;
	}
	if (!ref_bps[0]) {
		ref_bps[2] = -1000;
	}
	int64_t homology2 = -ref_bps[2];

	if (!ref_bps[2]) {
		homology2 = 0;
	}
	int64_t eventsizedel = left_bound_blast2 - left_bound_blast1 - 1;
	int64_t eventsizeins = right_bound_blast1 - right_bound_blast2 + 1;

	if (eventsizedel <= 10) {
		INTERREFINE << (boost::format("inso\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel % //
				data[8] % right_bound_blast2 % right_bound_blast1 % eventsizeins % homology1 % homology2).str();
	} else {
		INTERREFINE << (boost::format("del_inso\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\n") % data[1] %
				data[2] % data[3] % data[4] % left_bound_blast1 % left_bound_blast2 % eventsizedel %//
				data[8] % right_bound_blast2 % right_bound_blast1 % eventsizeins % homology1 % homology2).str();
	}
}

void BLASTBreakPoints::write_transl_inter(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n) {
	auto& ref_is = options.is;

	const string& cur_start_str = data[5];
	const string& cur_mate_start_str = data[8];
	const string& cur_strand_str = data[6];
	const string& cur_mate_strand_str = data[9];

	int64_t cur_start = boost::lexical_cast<int64_t>(cur_start_str);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(cur_mate_start_str);
	int64_t strand = boost::lexical_cast<int64_t>(cur_strand_str);
	int64_t mate_strand = boost::lexical_cast<int64_t>(cur_mate_strand_str);

	const string& ref_id = data[4];
	const string& mate_ref_id = data[7];

	vector<int64_t> ref_bps;
	vector<int64_t> ref_bps2;

	if (strand == 1 && mate_strand == -1) {
		int64_t position1 = cur_start - ref_is["rlu"]["selected"];
		int64_t position2 = cur_start + ref_is["rlu"]["selected"];
		int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
		int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];

		string seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(mate_ref_id, position3, position4);
		string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
		string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
		parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
		int64_t left_bound_blast = position1 + ref_bps[0] - 1;
		int64_t right_bound_blast = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
		if (!ref_bps[0]) {
			left_bound_blast = cur_start;
		}
		if (!ref_bps[1]) {
			right_bound_blast = cur_mate_start;
		}
		if (!ref_bps[0]) {
			ref_bps[2] = -1000;
		}
		int64_t homology = -ref_bps[2];
		if (!ref_bps[2]) {
			homology = 0;
		}
		// type[0] % cluster_id[1] % mate_cluster_id[2] % n_supports[3] % n_mate_support[4]
//		% ref_id[5] % event_start[6] % event_end[7] % event_size_1[8] % mate_ref_id[9] % mate_event_start[10] % mate_event_end[11]
//		% event_size_2[12] % strand[13] % mate_strand[14]

		INTERREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% data[0] % data[1] % data[2] % data[3] % data[4] % left_bound_blast % data[6] % data[7] % right_bound_blast % data[9] % homology).str();
	} else if (strand == -1 && mate_strand == 1) {
		int64_t position1 = cur_start - ref_is["rlu"]["selected"];
		int64_t position2 = cur_start + ref_is["rlu"]["selected"];
		int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
		int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];

		string seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(mate_ref_id, position3, position4);
		string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
		string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
		parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
		int64_t left_bound_blast = position1 + ref_bps[0] - 1;
		int64_t right_bound_blast = position3 + (ref_bps[1] - (position2 - position1 + 100)) - 2;
		if (!ref_bps[0]) {
			left_bound_blast = cur_start;
		}
		if (!ref_bps[1]) {
			right_bound_blast = cur_mate_start;
		}
		if (!ref_bps[0]) {
			ref_bps[2] = -1000;
		}
		int64_t homology = -ref_bps[2];
		if (!ref_bps[2]) {
			homology = 0;
		}
		INTERREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% data[0] % data[1] % data[2] % data[3] % data[4] % left_bound_blast % data[6] % data[7] % right_bound_blast % data[9] % homology).str();
	} else if (strand == 1 && mate_strand == 1) {
		int64_t position1 = cur_start - ref_is["rlu"]["selected"];
		int64_t position2 = cur_start + ref_is["rlu"]["selected"];
		int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
		int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];

		string seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(mate_ref_id, position4, position3);
		string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
		string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
		parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
		int64_t left_bound_blast = position1 + ref_bps[0] - 1;
		int64_t right_bound_blast = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
		if (!ref_bps[0]) {
			left_bound_blast = cur_start;
		}
		if (!ref_bps[1]) {
			right_bound_blast = cur_mate_start;
		}
		if (!ref_bps[0]) {
			ref_bps[2] = -1000;
		}
		int64_t homology = -ref_bps[2];
		if (!ref_bps[2]) {
			homology = 0;
		}
		INTERREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% data[0] % data[1] % data[2] % data[3] % data[4] % left_bound_blast % data[6] % data[7] % right_bound_blast % data[9] % homology).str();
	} else if (strand == -1 && mate_strand == -1) {
		int64_t position1 = cur_start - ref_is["rlu"]["selected"];
		int64_t position2 = cur_start + ref_is["rlu"]["selected"];
		int64_t position3 = cur_mate_start - ref_is["rlu"]["selected"];
		int64_t position4 = cur_mate_start + ref_is["rlu"]["selected"];

		string seq = fai.fetch_1_system(ref_id, position1, position2);
		seq += n;
		seq += fai.fetch_1_system(mate_ref_id, position4, position3);
		string name = ref_id + "__" + cur_start_str + "__" + mate_ref_id + "__" + cur_mate_start_str;
		string bltoutfile = run_blast(blast_dir, data, mpd_id[0], name, seq, readsdetail, true);
		parse_blast(ref_bps, ref_bps2, bltoutfile, readslist, readsdetail);
		int64_t left_bound_blast = position1 + ref_bps[0] - 1;
		int64_t right_bound_blast = position4 - (ref_bps[1] - (position2 - position1 + 100)) + 2;
		if (!ref_bps[0]) {
			left_bound_blast = cur_start;
		}
		if (!ref_bps[1]) {
			right_bound_blast = cur_mate_start;
		}
		if (!ref_bps[0]) {
			ref_bps[2] = -1000;
		}
		int64_t homology = -ref_bps[2];
		if (!ref_bps[2]) {
			homology = 0;
		}
		INTERREFINE << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% data[0] % data[1] % data[2] % data[3] % data[4] % left_bound_blast % data[6] % data[7] % right_bound_blast % data[9] % homology).str();
	}
}

string BLASTBreakPoints::run_blast(const string blastdir, const vector<string>& ref_data, const string& mpd_id, const string& name, const string& seq, const boost::unordered_map<string, boost::unordered_map<string, string>>& ref_readsdetail, bool forward) {
	//cout << "run_blast here-0\n";
//	const string blastdir = options.prefix + ".blast/";
	//cout << "run_blast here-0-a:" << sbjctfile << "\n";
	const string delim_double_underscore = "__";
	const string sbjctfile = blastdir + mpd_id + ".sbj.fa";

	//cout << "run_blast here-1\n";
	ofstream SBJCTSEQ(sbjctfile, ios::binary);
	SBJCTSEQ << (boost::format(">%s\n%s\n") % name % seq).str();
	SBJCTSEQ.close();
	//cout << "run_blast here-2: " << name << "\n";
	string queryfile = blastdir + mpd_id + ".query.fa";
	string bltoutfile = blastdir + mpd_id + ".bltout";
	if(mpd_id.empty()) {
		cout << "[BLASTBreakPoints.run_blast] bltout file: " << bltoutfile << "\n";
		exit(0);
	}
	ofstream QRSEQ(queryfile, ios::binary);
	auto the_detail = ref_readsdetail.find(name);
	//	# flip the order of break points
	if (ref_readsdetail.end() == the_detail) {
		vector<castle::c_string_view> temp_bp;
		castle::StringUtils::tokenize(name, delim_double_underscore, temp_bp);
		swap(temp_bp[0], temp_bp[2]);
		swap(temp_bp[1], temp_bp[3]);
		string alt_name = castle::StringUtils::join(temp_bp, "__");
		the_detail = ref_readsdetail.find(alt_name);
	}
	if (ref_readsdetail.end() != the_detail) {
		//cout << "run_blast here-3\n";
//		auto& key = the_detail->first;
		auto& second_entry = the_detail->second;
		for (auto the_third_entry : second_entry) {
			auto& key = the_third_entry.first;
			auto& value = the_third_entry.second;
			QRSEQ << (boost::format(">%s\n%s\n") % key % value).str();
		}
	}
	QRSEQ.close();
	//cout << "run_blast here-6\n";
	// -i Input file for formatting -p F: nucleotide sequence;
	string formatdb_cmd = (boost::format("formatdb -i %s -p F") % sbjctfile).str();
	string blastall_cmd = (boost::format("blastall -p blastn -d %s -i %s -o %s -m 8 -F F") % sbjctfile % queryfile % bltoutfile).str();

	system(formatdb_cmd.c_str());
	//cout << "run_blast here-7\n";
	system(blastall_cmd.c_str());
	//cout << "run_blast here-8: " << castle::IOUtils::get_file_size(bltoutfile) << "\n";
	return bltoutfile;
}
//#	extract break points from blast output
void BLASTBreakPoints::parse_blast(vector<int64_t>& results, vector<int64_t>& results2, const string& bltoutfile, const boost::unordered_map<string, string>& ref_readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& ref_readsdetail) {
//	const bool debug = string::npos != bltoutfile.find("3065");
	const bool debug = false;
	if(debug) {
		cout << "parse_blast here-0\n";
	}
	auto& ref_is = options.is;
//	int64_t cut_sr = options.cut_sr;
	const int64_t support_reads = options.support_reads;
	const int64_t bp1 = ref_is["rlu"]["selected"];
	const int64_t bp2 = ref_is["rlu"]["selected"] * 3 + 100;
	string srfile = bltoutfile;
	boost::replace_all(srfile, "blast", "sr");
	boost::replace_all(srfile, "bltout", "srout");
	string sbjfafile = bltoutfile;
	boost::replace_all(sbjfafile, "bltout", "sbj.fa");
	const char* delim_amb_base = "N";
	const char* delim_tab_space = "\t ";
	string line;
	ifstream SBJFA(sbjfafile, ios::out | ios::binary);
	// discard read name
	getline(SBJFA, line, '\n');
	// read sequence
	getline(SBJFA, line, '\n');

	vector<string> ref_seq;
	castle::StringUtils::c_string_multi_split(line, delim_amb_base, ref_seq);
	if (1 == ref_seq.size()) {
		if(debug) {
			cout << "parse_blast here-1: " << sbjfafile << ":" << ref_seq.size() << "\n";
		}
	}

	ofstream SROUT(srfile, ios::out | ios::binary);

	boost::unordered_map<string, vector<BLASTEntry>> alg;
	map<int64_t, int64_t> freq_s1;
	map<int64_t, int64_t> freq_s2;
	map<int64_t, int64_t> freq_dist;
//	#	%freq_s1: frequency of break points on one region
//	#	%freq_s2: frequency of break points on the other region
//	#	%freq_dist: frequency of query break points distance, positive: unmatched bases, negative: homology
//	#	%alg: alignment details
//	#	alg{readname}{i}: reference to ith alignment detail


	{
		if(debug) {
			cout << "parse_blast here-3\n";
		}
		vector<string> data;
		ifstream BLTOUT(bltoutfile, ios::in | ios::binary);
		while (getline(BLTOUT, line, '\n')) {
			//		my (query, subject, identity, hitlength, mismatch, gap, qrstart, qrend, sbjstart, sbjend, evalue, score);
			if(debug) {
				cout << "parse_blast here-4: " << line << "\n";
			}
			castle::StringUtils::c_string_multi_split(line, delim_tab_space, data);
			BLASTEntry an_entry;
			an_entry.query_id = data[0];
			an_entry.subject_name = data[1];

			an_entry.percent_identity = boost::lexical_cast<double>(data[2]);
			an_entry.aln_length = boost::lexical_cast<int64_t>(data[3]);
			an_entry.n_mismatches = boost::lexical_cast<int64_t>(data[4]);
			an_entry.n_gap_opens = boost::lexical_cast<int64_t>(data[5]);
			an_entry.query_start = boost::lexical_cast<int64_t>(data[6]);
			an_entry.query_end = boost::lexical_cast<int64_t>(data[7]);
			an_entry.subject_start = boost::lexical_cast<int64_t>(data[8]);
			an_entry.subject_end = boost::lexical_cast<int64_t>(data[9]);
			an_entry.expect_value = boost::lexical_cast<double>(data[10]);
			an_entry.bit_score = boost::lexical_cast<double>(data[11]);
			alg[data[0]].push_back(an_entry);
			//		#		print "i\t@{alg{data[0]}{i}}\n";
		}
	}
	if(debug) {
		cout << "parse_blast here-5\n";
	}
	int64_t ref_seq_ori_det1 = 0;   // # determine the orientation of 1st ref seq, 0 not yet determined, 1 determined already
	int64_t ref_seq_ori_det2 = 0;
	int64_t ref_seq_ori1 = 1;    // # orientation of ref seq, 1 forward, 0 reverse
	int64_t ref_seq_ori2 = 1;
//	my @aligned_reads;
	vector<string> aligned_reads;
	for (auto& first_entry : alg) {
		auto& key = first_entry.first;
		auto& value = first_entry.second;
		if (value.size() < 2) {
			continue;
		}
		if(debug) {
			cout << "parse_blast key: " << key << "\n";
		}
//		# only 2 entries
		if (value.size() == 2) {
			if(debug) {
				cout << "parse_blast here-6\n";
			}
//			# not process reads not entirely used (such read should belong to other break point, caused by overlapping regions produced by small events)
			if ((value[1].query_start >= value[0].query_start && value[1].query_start <= value[0].query_end && value[1].query_end >= value[0].query_start && value[1].query_end <= value[0].query_end)
					|| (value[0].query_start >= value[1].query_start && value[0].query_start <= value[1].query_end && value[0].query_end >= value[1].query_start && value[0].query_end <= value[1].query_end)) {
				if(debug) {
					cout << "parse_blast here-7\n";
				}
				continue;
			}
//			# skip if not the entire read is used
//			#if (value[1].query_start > value[0].query_start)
//			#{
//			#	next if ((value[1].query_start-value[0].query_end)>5 && !(is_del));
//			#}
//			#else
//			#{
//			#	next if ((value[0].query_start-value[1].query_end)>5 && !(is_del));
//			#}
//			# 1st half read on top
			if (value[0].query_end > value[0].query_start && value[0].query_end < value[1].query_end) {
				if(debug) {
					cout << "parse_blast here-8: freq\n";
					cout << "parse_blast here-8-a: value0:" << value[0].subject_start << "/" << value[0].subject_end << "\n";
					cout << "parse_blast here-8-b: value1:" << value[1].subject_start << "/" << value[1].subject_end << "\n";
				}
				if (value[0].subject_end > value[1].subject_start) {
					if(debug) {
						cout << "parse_blast here-9: s1: " << value[1].subject_start << "/s2: " << value[0].subject_end << "\n";
					}
					++freq_s1[value[1].subject_start];
					++freq_s2[value[0].subject_end];
				} else {
					if(debug) {
						cout << "parse_blast here-10: s1: " << value[0].subject_end << "/s2: " << value[1].subject_start << "\n";
					}
					++freq_s1[value[0].subject_end];
					++freq_s2[value[1].subject_start];
				}

				int64_t dist = value[1].query_start - value[0].query_end - 1;
				if(debug) {
					cout << "parse_blast here-11: dist: " << value[1].query_start << "-" << value[0].query_end << " => " << dist << "\n";
				}
				++freq_dist[dist];
				if (!ref_seq_ori_det1) {
					if(debug) {
						cout << "parse_blast here-12\n";
					}
					if ((value[0].subject_start < value[1].subject_start && value[0].subject_start > value[0].subject_end) || (value[0].subject_start > value[1].subject_start && value[0].subject_start < value[0].subject_end)) {
						if(debug) {
							cout << "parse_blast here-13\n";
						}
						ref_seq[0] = castle::StringUtils::get_reverse_complement(ref_seq[0]);
						ref_seq_ori1 = 0;
					}
					ref_seq_ori_det1 = 1;
				}
				if (!ref_seq_ori_det2) {
					if(debug) {
						cout << "parse_blast here-14\n";
					}
					if ((value[0].subject_start < value[1].subject_start && value[1].subject_start > value[1].subject_end) || (value[0].subject_start > value[1].subject_start && value[1].subject_start < value[1].subject_end)) {
						if(debug) {
							cout << "parse_blast here-15\n";
						}
						if (ref_seq.size() > 1) {
							ref_seq[1] = castle::StringUtils::get_reverse_complement(ref_seq[1]);
						}
						ref_seq_ori2 = 0;
					}
					ref_seq_ori_det2 = 1;
				}
				int64_t num_space = 0;
				int64_t trim_size = 0;
//				my (space, num_space, trim_size);
				auto& a_query_id = value[0].query_id;
				auto a_read_itr = ref_readslist.find(a_query_id);

				string aligned_read;

				if (ref_readslist.end() != a_read_itr) {
					auto a_read_detail_itr = ref_readsdetail.find(a_read_itr->second);
					if (ref_readsdetail.end() != a_read_detail_itr) {
						auto the_read_itr = a_read_detail_itr->second.find(a_query_id);
						if (a_read_detail_itr->second.end() != the_read_itr) {
							aligned_read = the_read_itr->second;
						}
					}
				}
				if(debug) {
					cout << "parse_blast here-16: " << aligned_read << "\n";
				}
				if (value[0].subject_start < value[1].subject_start) {
					if(debug) {
						cout << "parse_blast here-17\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-18\n";
						}
						num_space = min(value[0].subject_start, value[0].subject_end) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-19\n";
						}
						num_space = ref_seq[0].size() - max(value[0].subject_start, value[0].subject_end);
					}
					trim_size = 1 - min(value[0].query_start, value[0].query_end);
				} else {
					if(debug) {
						cout << "parse_blast here-20\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-21\n";
						}
						num_space = min(value[1].subject_start, value[1].subject_end) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-22\n";
						}
						num_space = ref_seq[0].size() - max(value[1].subject_start, value[1].subject_end);
					}
					aligned_read = castle::StringUtils::get_reverse_complement(aligned_read);
					trim_size = 1 - min(value[0].query_start, value[0].query_end);
				}
				aligned_read = aligned_read.substr(-trim_size);
				if (0 != trim_size * num_space) {
					if(debug) {
						cout << "parse_blast here-23\n";
					}
					for (int64_t i = 0; i < -trim_size; ++i) {
						aligned_read = " " + aligned_read;
					}
				}
				if(debug) {
					cout << "parse_blast here-24\n";
				}
				for (int64_t i = 0; i < num_space; ++i) {
					aligned_read = " " + aligned_read;
				}
				if (!aligned_read.empty()) {
					if(debug) {
						cout << "parse_blast here-25\n";
					}
					aligned_reads.push_back(aligned_read);
				}
//				#print "value[0][0]\tnum_space\ttrim_size\tref_seq_ori1\tref_seq_ori2\n";
//				#print "value[1].query_start\tvalue[0].query_end\tdist\n";
			}

//			# 2nd half read on top
			if (value[1].query_end > value[1].query_start && value[1].query_end < value[0].query_end) {
				if(debug) {
					cout << "parse_blast here-26\n";
				}
				if (value[1].subject_end > value[0].subject_start) {
					if(debug) {
						cout << "parse_blast here-27: s1: " << value[0].subject_start << "/s2: " << value[1].subject_end << "\n";
					}
					++freq_s2[value[1].subject_end];
					++freq_s1[value[0].subject_start];
				} else {
					if(debug) {
						cout << "parse_blast here-28: s1: " << value[1].subject_end << "/s2: " << value[0].subject_start << "\n";
					}
					++freq_s1[value[1].subject_end];
					++freq_s2[value[0].subject_start];
				}
				int64_t dist = value[0].query_start - value[1].query_end - 1;
				if(debug) {
					cout << "parse_blast here-29: dist: " << value[0].query_start << "-" << value[1].query_end << "=>" << dist << "\n";
				}
				++freq_dist[dist];
				if (!ref_seq_ori_det1) {
					if(debug) {
						cout << "parse_blast here-30\n";
					}
					if ((value[1].subject_start < value[0].subject_start && value[1].subject_start > value[1].subject_end) || (value[1].subject_start > value[0].subject_start && value[1].subject_start < value[1].subject_end)) {
						if(debug) {
							cout << "parse_blast here-31\n";
						}
						if (ref_seq.size() > 0) {
							ref_seq[0] = castle::StringUtils::get_reverse_complement(ref_seq[0]);
						}
						ref_seq_ori1 = 0;
					}
					ref_seq_ori_det1 = 1;
				}
				if (!ref_seq_ori_det2) {
					if(debug) {
						cout << "parse_blast here-32\n";
					}
					if ((value[1].subject_start < value[0].subject_start && value[0].subject_start > value[0].subject_end) || (value[1].subject_start > value[0].subject_start && value[0].subject_start < value[0].subject_end)) {
						if(debug) {
							cout << "parse_blast here-33\n";
						}
						if (ref_seq.size() > 1) {
							ref_seq[1] = castle::StringUtils::get_reverse_complement(ref_seq[1]);
						}
						ref_seq_ori2 = 0;
					}
					ref_seq_ori_det2 = 1;
				}
//				my (space.substr(num_space, trim_size);
				int64_t num_space = 0;
				int64_t trim_size = 0;
//				string aligned_read = ref_readsdetail[ref_readslist[value[1].query_id]][value[1].query_id];
				auto a_query_id = value[1].query_id;
				auto a_read_itr = ref_readslist.find(a_query_id);
				string aligned_read;
				if (ref_readslist.end() != a_read_itr) {
					auto a_read_detail_itr = ref_readsdetail.find(a_read_itr->second);
					if (ref_readsdetail.end() != a_read_detail_itr) {
						auto the_read_itr = a_read_detail_itr->second.find(a_query_id);
						if (a_read_detail_itr->second.end() != the_read_itr) {
							aligned_read = the_read_itr->second;
						}
					}
				}
				if(debug) {
					cout << "parse_blast here-34: " << aligned_read << "\n";
				}
				if (value[1].subject_start < value[0].subject_start) {
					if(debug) {
						cout << "parse_blast here-35\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-36\n";
						}
						num_space = min(value[1].subject_start, value[1].subject_end) - 1;
						for (uint64_t base_id = 0; base_id < aligned_read.size(); ++base_id) {
							if ('N' == aligned_read[base_id]) {
								--num_space;
							} else {
								break;
							}
						}
					} else {
						if(debug) {
							cout << "parse_blast here-37\n";
						}
						num_space = ref_seq[0].size() - max(value[1].subject_start, value[1].subject_end);
					}
					trim_size = 1 - min(value[1].query_start, value[1].query_end);
				} else {
					if(debug) {
						cout << "parse_blast here-38\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-39\n";
						}
						num_space = min(value[0].subject_start, value[0].subject_end) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-40\n";
						}
						num_space = ref_seq[0].size() - max(value[0].subject_start, value[0].subject_end);
					}
					aligned_read = castle::StringUtils::get_reverse_complement(aligned_read);
					trim_size = 1 - min(value[1].query_start, value[1].query_end);
				}
				if (0 == trim_size * num_space) {
					if(debug) {
						cout << "parse_blast here-41\n";
					}
					aligned_read = aligned_read.substr(-trim_size);
				} else {
					if(debug) {
						cout << "parse_blast here-42\n";
					}
					aligned_read = aligned_read.substr(-trim_size);
					string an_empty_str(-trim_size, ' ');
					aligned_read = an_empty_str + aligned_read;
//					for (int64_t i = 0; i < -trim_size; ++i) {
//						aligned_read = " " + aligned_read;
//					}
				}
				if(debug) {
					cout << "parse_blast here-43\n";
				}
				if (num_space > 0) {
					string an_empty_str(num_space, ' ');
					aligned_read = an_empty_str + aligned_read;
				}
				if (!aligned_read.empty()) {
					if(debug) {
						cout << "parse_blast here-44\n";
					}
					aligned_reads.push_back(aligned_read);
				}

//				#print "value[0][0]\tnum_space\ttrim_size\tref_seq_ori1\tref_seq_ori2\n";
//				#print "value[0].query_start\tvalue[1].query_end\tdist\n";
			}
		}
//		# more than 2 entries
		else {
//			# get top 2 hits
			map<string, map<int32_t, BLASTEntry>> alg_top2;
//			my %alg_top2;
			int64_t maxk = value.size();
//			my maxk = keys %{value};
//			# initialize first entry
			alg_top2[key][0] = value[0];
//			@{alg_top2[key][0}} = @{value[0]};
//			value.erase(0);
//			delete value[0];
//			# get best hit for first entry
			for (int64_t k = 1; k < maxk; ++k) {
//				if(value.end() == value.find(k)) {
//					continue;
//				}
				if(debug) {
					cout << "parse_blast here-45\n";
				}
//			for(auto the_second_entry : value) {
				auto& a_BLAST_entry = value[k];
				if (alg_top2[key][0].query_start == a_BLAST_entry.query_start && alg_top2[key][0].query_end == a_BLAST_entry.query_end) {
					if(debug) {
						cout << "parse_blast here-46\n";
					}
//					# forward strand
					if (a_BLAST_entry.subject_end > a_BLAST_entry.subject_start) {
						if(debug) {
							cout << "parse_blast here-47\n";
						}
//						# first region
						if (abs(a_BLAST_entry.subject_end - bp1) < abs(a_BLAST_entry.subject_start - bp2)) {
							if(debug) {
								cout << "parse_blast here-48\n";
							}
							if (abs(a_BLAST_entry.subject_start - bp1) < abs(alg_top2[key][0].subject_start - bp1)) {
								if(debug) {
									cout << "parse_blast here-49\n";
								}
								alg_top2[key][0] = a_BLAST_entry;
							}
						}

//						# second region
						else {
							if(debug) {
								cout << "parse_blast here-50\n";
							}
							if (abs(a_BLAST_entry.subject_end - bp2) < abs(alg_top2[key][0].subject_end - bp2)) {
								if(debug) {
									cout << "parse_blast here-51\n";
								}
								alg_top2[key][0] = a_BLAST_entry;
							}
						}
					}
//					# reverse strand
					else {
						if(debug) {
							cout << "parse_blast here-52\n";
						}
//						# first region
						if (abs(a_BLAST_entry.subject_end - bp1) < abs(a_BLAST_entry.subject_start - bp2)) {
							if(debug) {
								cout << "parse_blast here-53: a blast: " << a_BLAST_entry.subject_start << "/" << a_BLAST_entry.subject_end << "\n";
								cout << "parse_blast here-53: top_0 " << alg_top2[key][0].subject_start << "/" << alg_top2[key][0].subject_end << "\n";
							}
							if (abs(a_BLAST_entry.subject_end - bp1) < abs(alg_top2[key][0].subject_end - bp1)) {
								if(debug) {
									cout << "parse_blast here-54\n";
								}
								alg_top2[key][0] = a_BLAST_entry;
							}
						}

//						# second region
						else {
							if(debug) {
								cout << "parse_blast here-55\n";
							}
							if (abs(a_BLAST_entry.subject_start - bp2) < abs(alg_top2[key][0].subject_start - bp2)) {
								if(debug) {
									cout << "parse_blast here-56\n";
								}
								alg_top2[key][0] = a_BLAST_entry;
							}
						}
					}
//					delete the_element_vec;
//					value.erase(k);
				} else {
					break;
				}
			}

//			# initialize second entry
			int64_t startk = 1;
			for (int64_t k = 1; k < maxk; ++k) {
//				if (value.end() == value.find(k)) {
//					continue;
//				}
				auto& a_BLAST_entry = value[k];
//			for(auto the_second_entry : value) {
//				int32_t k = the_second_entry.first;
//				auto& the_element_vec = the_second_entry.second;
				if (!a_BLAST_entry.query_id.empty()) {
					if(debug) {
						cout << "parse_blast here-57\n";
					}
					if ((a_BLAST_entry.query_start >= alg_top2[key][0].query_start && a_BLAST_entry.query_start <= alg_top2[key][0].query_end && a_BLAST_entry.query_end >= alg_top2[key][0].query_start && a_BLAST_entry.query_end <= alg_top2[key][0].query_end)
							|| (alg_top2[key][0].query_start >= a_BLAST_entry.query_start && alg_top2[key][0].query_start <= a_BLAST_entry.query_end && alg_top2[key][0].query_end >= a_BLAST_entry.query_start && alg_top2[key][0].query_end <= a_BLAST_entry.query_end)) {
						continue;
					}
					if(debug) {
						cout << "parse_blast here-58\n";
					}
					alg_top2[key][1] = a_BLAST_entry;
//					value.erase(k);
//					delete the_element_vec;
					startk = k + 1;
					break;
				}
			}

//			#print "@{alg_top2[key][0}}\n@{alg_top2[key][1}}\n";
//			# get best hit for second entry
			for (int64_t k = startk; k < maxk; ++k) {
//				if(value.end() == value.find(k)) {
//					continue;
//				}
				auto& a_BLAST_entry = value[k];
				if (alg_top2[key][1].query_start == a_BLAST_entry.query_start && alg_top2[key][1].query_end == a_BLAST_entry.query_end) {
					if(debug) {
						cout << "parse_blast here-59\n";
					}
//					# forward strand
					if (a_BLAST_entry.subject_end > a_BLAST_entry.subject_start) {
						if(debug) {
							cout << "parse_blast here-60\n";
						}
//						# first region
						if (abs(a_BLAST_entry.subject_end - bp1) < abs(a_BLAST_entry.subject_start - bp2)) {
							if(debug) {
								cout << "parse_blast here-61\n";
							}
							if (abs(a_BLAST_entry.subject_start - bp1) < abs(alg_top2[key][1].subject_start - bp1)) {
								if(debug) {
									cout << "parse_blast here-62\n";
								}
								alg_top2[key][1] = a_BLAST_entry;
							}
						}

//						# second region
						else {
							if(debug) {
								cout << "parse_blast here-63\n";
							}
							if (abs(a_BLAST_entry.subject_end - bp2) < abs(alg_top2[key][1].subject_end - bp2)) {
								if(debug) {
									cout << "parse_blast here-64\n";
								}
								alg_top2[key][1] = a_BLAST_entry;
							}
						}
					}

//					# reverse strand
					else {
						if(debug) {
							cout << "parse_blast here-65\n";
						}
//						# first region
						if (abs(a_BLAST_entry.subject_end - bp1) < abs(a_BLAST_entry.subject_start - bp2)) {
							if(debug) {
								cout << "parse_blast here-66\n";
							}
							if (abs(a_BLAST_entry.subject_end - bp1) < abs(alg_top2[key][1].subject_end - bp1)) {
								if(debug) {
									cout << "parse_blast here-67\n";
								}
								alg_top2[key][1] = a_BLAST_entry;
							}
						}
//						# second region
						else {
							if(debug) {
								cout << "parse_blast here-68\n";
							}
							if (abs(a_BLAST_entry.subject_start - bp2) < abs(alg_top2[key][1].subject_start - bp2)) {
								if(debug) {
									cout << "parse_blast here-69\n";
								}
								alg_top2[key][1] = a_BLAST_entry;
							}
						}
					}
//					delete the_element_vec;
//					value.erase(k);
				} else {
					if(debug) {
						cout << "parse_blast here-70\n";
					}
					break;
				}
			}

//			#print "@{alg_top2[key][0}}\n@{alg_top2[key][1}}\n";
//			# skip if not the entire read is used
//			#if (alg_top2[key][1].query_start > alg_top2[key][0].query_start)
//			#{
//			#	next if ((alg_top2[key][1].query_start-alg_top2[key][0].query_end)>5 && !(is_del));
//			#}
//			#else
//			#{
//			#	next if ((alg_top2[key][0].query_start-alg_top2[key][1].query_end)>5 && !(is_del));
//			#}
//			# 1st half read on top
			if (alg_top2[key][0].query_end > alg_top2[key][0].query_start && alg_top2[key][0].query_end < alg_top2[key][1].query_end) {
				if(debug) {
					cout << "parse_blast here-71-a: top2_key0: " << alg_top2[key][0].subject_start << "/" << alg_top2[key][0].subject_end << "\n";
					cout << "parse_blast here-71-b: top2_key1: " << alg_top2[key][1].subject_start << "/" << alg_top2[key][1].subject_end << "\n";
				}
				if (alg_top2[key][0].subject_end > alg_top2[key][1].subject_start) {
					if(debug) {
						cout << "parse_blast here-72: s1: " << alg_top2[key][1].subject_start << "/s2: " << alg_top2[key][0].subject_end << "\n";
					}
					++freq_s1[alg_top2[key][1].subject_start];
					++freq_s2[alg_top2[key][0].subject_end];
				} else {
					if(debug) {
						cout << "parse_blast here-73: s1: " << alg_top2[key][0].subject_end << "/s2: " << alg_top2[key][1].subject_start << "\n";
					}
					++freq_s1[alg_top2[key][0].subject_end];
					++freq_s2[alg_top2[key][1].subject_start];
				}
				int64_t dist = alg_top2[key][1].query_start - alg_top2[key][0].query_end - 1;
				if(debug) {
					cout << "parse_blast here-73-a:dist " << alg_top2[key][1].query_start << "-" << alg_top2[key][0].query_end << "=>" << dist << "\n";
				}
				++freq_dist[dist];
				if (!ref_seq_ori_det1) {
					if(debug) {
						cout << "parse_blast here-74\n";
					}
					if ((alg_top2[key][0].subject_start < alg_top2[key][1].subject_start && alg_top2[key][0].subject_start > alg_top2[key][0].subject_end) || (alg_top2[key][0].subject_start > alg_top2[key][1].subject_start && alg_top2[key][0].subject_start < alg_top2[key][0].subject_end)) {
						if(debug) {
							cout << "parse_blast here-75\n";
						}
						if (ref_seq.size() > 0) {
							ref_seq[0] = castle::StringUtils::get_reverse_complement(ref_seq[0]);
						}
						ref_seq_ori1 = 0;
					}
					ref_seq_ori_det1 = 1;
				}
				if (!ref_seq_ori_det2) {
					if(debug) {
						cout << "parse_blast here-76\n";
					}
					if ((alg_top2[key][0].subject_start < alg_top2[key][1].subject_start && alg_top2[key][1].subject_start > alg_top2[key][1].subject_end) || (alg_top2[key][0].subject_start > alg_top2[key][1].subject_start && alg_top2[key][1].subject_start < alg_top2[key][1].subject_end)) {
						if(debug) {
							cout << "parse_blast here-77\n";
						}
						if (ref_seq.size() > 1) {
							ref_seq[1] = castle::StringUtils::get_reverse_complement(ref_seq[1]);
						}
						ref_seq_ori2 = 0;
					}
					ref_seq_ori_det2 = 1;
				}
				if(debug) {
					cout << "parse_blast here-78\n";
				}
				int64_t num_space = 0;
				int64_t trim_size = 0;
//				my (space, num_space, trim_size);
//				string aligned_read = ref_readsdetail[ref_readslist[alg_top2[key][0].query_id]][alg_top2[key][0].query_id];
				auto a_query_id = alg_top2[key][0].query_id;
				auto a_read_itr = ref_readslist.find(a_query_id);
				string aligned_read;
				if (ref_readslist.end() != a_read_itr) {
					if(debug) {
						cout << "parse_blast here-79\n";
					}
					auto a_read_detail_itr = ref_readsdetail.find(a_read_itr->second);
					if (ref_readsdetail.end() != a_read_detail_itr) {
						auto the_read_itr = a_read_detail_itr->second.find(a_query_id);
						if (a_read_detail_itr->second.end() != the_read_itr) {
							aligned_read = the_read_itr->second;
						}
					}
				}
				if(debug) {
					cout << "parse_blast here-80: " << aligned_read << "\n";
				}

				if (alg_top2[key][0].subject_start < alg_top2[key][1].subject_start) {
					if(debug) {
						cout << "parse_blast here-81\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-82\n";
						}
						num_space = (min(alg_top2[key][0].subject_start, alg_top2[key][0].subject_end)) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-83\n";
						}
						num_space = ref_seq[0].size() - (max(alg_top2[key][0].subject_start, alg_top2[key][0].subject_end));
					}
					trim_size = 1 - (min(alg_top2[key][0].query_start, alg_top2[key][0].query_end));
				} else {
					if(debug) {
						cout << "parse_blast here-84\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-85\n";
						}
						num_space = (min(alg_top2[key][1].subject_start, alg_top2[key][1].subject_end)) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-86\n";
						}
						num_space = ref_seq[0].size() - (max(alg_top2[key][1].subject_start, alg_top2[key][1].subject_end));
					}
					aligned_read = castle::StringUtils::get_reverse_complement(aligned_read);
					trim_size = 1 - (min(alg_top2[key][0].query_start, alg_top2[key][0].query_end));
				}
				if (trim_size < 0 && static_cast<int64_t>(aligned_read.size()) > -trim_size) {
					aligned_read = aligned_read.substr(-trim_size);
				}
				if (0 != trim_size * num_space) {
					if(debug) {
						cout << "parse_blast here-88\n";
					}
					for (int64_t i = 0; i < -trim_size; ++i) {
						aligned_read = " " + aligned_read;
					}
				}
				if(debug) {
					cout << "parse_blast here-89\n";
				}
				for (int64_t i = 0; i < num_space; ++i) {
					aligned_read = " " + aligned_read;
				}
				if (!aligned_read.empty()) {
					if(debug) {
						cout << "parse_blast here-90\n";
					}
					aligned_reads.push_back(aligned_read);
				}

//				#print "alg_top2[key][0][0]\tnum_space\ttrim_size\tref_seq_ori1\tref_seq_ori2\n";
//				#print "alg_top2[key][1].query_start\talg_top2[key][0].query_end\tdist\n";
			}

//			# 2nd half read on top
			if (alg_top2[key][1].query_end > alg_top2[key][1].query_start && alg_top2[key][1].query_end < alg_top2[key][0].query_end) {
				if(debug) {
					cout << "parse_blast here-91-a: top2_key0: " << alg_top2[key][0].subject_start << "/" << alg_top2[key][0].subject_end << "\n";
					cout << "parse_blast here-91-b: top2_key1: " << alg_top2[key][1].subject_start << "/" << alg_top2[key][1].subject_end << "\n";
				}
				if (alg_top2[key][1].subject_end > alg_top2[key][0].subject_start) {
					if(debug) {
						cout << "parse_blast here-92: s1: " << alg_top2[key][0].subject_start << "/s2: " << alg_top2[key][1].subject_end << "\n";
					}
					++freq_s1[alg_top2[key][0].subject_start];
					++freq_s2[alg_top2[key][1].subject_end];
				} else {
					if(debug) {
						cout << "parse_blast here-93: s1: " << alg_top2[key][1].subject_end << "/s2: " << alg_top2[key][0].subject_start << "\n";
					}
					++freq_s1[alg_top2[key][1].subject_end];
					++freq_s2[alg_top2[key][0].subject_start];
				}
				int64_t dist = alg_top2[key][0].query_start - alg_top2[key][1].query_end - 1;
				if(debug) {
					cout << "parse_blast here-93-a:dist: " << alg_top2[key][0].query_start << "-" << alg_top2[key][1].query_end << "=>" << dist << "\n";
				}
				++freq_dist[dist];
				if (!ref_seq_ori_det1) {
					if(debug) {
						cout << "parse_blast here-94\n";
					}
					if ((alg_top2[key][1].subject_start < alg_top2[key][0].subject_start && alg_top2[key][1].subject_start > alg_top2[key][1].subject_end) || (alg_top2[key][1].subject_start > alg_top2[key][0].subject_start && alg_top2[key][1].subject_start < alg_top2[key][1].subject_end)) {
						if(debug) {
							cout << "parse_blast here-95\n";
						}
						if (ref_seq.size() > 0) {
							ref_seq[0] = castle::StringUtils::get_reverse_complement(ref_seq[0]);
						}
						ref_seq_ori1 = 0;
					}
					ref_seq_ori_det1 = 1;
				}
				if (!ref_seq_ori_det2) {
					if(debug) {
						cout << "parse_blast here-96\n";
					}
					if ((alg_top2[key][1].subject_start < alg_top2[key][0].subject_start && alg_top2[key][0].subject_start > alg_top2[key][0].subject_end) || (alg_top2[key][1].subject_start > alg_top2[key][0].subject_start && alg_top2[key][0].subject_start < alg_top2[key][0].subject_end)) {
						if(debug) {
							cout << "parse_blast here-97\n";
						}
						if (ref_seq.size() > 1) {
							ref_seq[1] = castle::StringUtils::get_reverse_complement(ref_seq[1]);
						}
						ref_seq_ori2 = 0;
					}
					ref_seq_ori_det2 = 1;
				}
//				my (space, num_space, trim_size);
				int64_t num_space = 0;
				int64_t trim_size = 0;
//				string aligned_read = ref_readsdetail[ref_readslist[alg_top2[key][1].query_id]][alg_top2[key][1].query_id];
				auto a_query_id = alg_top2[key][1].query_id;
				auto a_read_itr = ref_readslist.find(a_query_id);
				string aligned_read;
				if (ref_readslist.end() != a_read_itr) {
					if(debug) {
						cout << "parse_blast here-98\n";
					}
					auto a_read_detail_itr = ref_readsdetail.find(a_read_itr->second);
					if (ref_readsdetail.end() != a_read_detail_itr) {
						auto the_read_itr = a_read_detail_itr->second.find(a_query_id);
						if (a_read_detail_itr->second.end() != the_read_itr) {
							aligned_read = the_read_itr->second;
						}
					}
				}
				if (alg_top2[key][1].subject_start < alg_top2[key][0].subject_start) {
					if(debug) {
						cout << "parse_blast here-99\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-100\n";
						}
						num_space = (min(alg_top2[key][1].subject_start, alg_top2[key][1].subject_end)) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-101\n";
						}
						num_space = ref_seq[0].size() - (max(alg_top2[key][1].subject_start, alg_top2[key][1].subject_end));
					}
					trim_size = 1 - (min(alg_top2[key][1].query_start, alg_top2[key][1].query_end));
				} else {
					if(debug) {
						cout << "parse_blast here-102\n";
					}
					if (ref_seq_ori1) {
						if(debug) {
							cout << "parse_blast here-103\n";
						}
						num_space = (min(alg_top2[key][0].subject_start, alg_top2[key][0].subject_end)) - 1;
					} else {
						if(debug) {
							cout << "parse_blast here-104\n";
						}
						num_space = ref_seq[0].size() - (max(alg_top2[key][0].subject_start, alg_top2[key][0].subject_end));
					}
					aligned_read = castle::StringUtils::get_reverse_complement(aligned_read);
					trim_size = 1 - (min(alg_top2[key][1].query_start, alg_top2[key][1].query_end));
				}
				if (trim_size * num_space == 0) {
					if(debug) {
						cout << "parse_blast here-105\n";
					}
					aligned_read = aligned_read.substr(-trim_size);
				} else {
					if(debug) {
						cout << "parse_blast here-106\n";
					}
					aligned_read = aligned_read.substr(-trim_size);
					for (int64_t i = 0; i < -trim_size; ++i) {
						aligned_read = " " + aligned_read;
					}
				}
				if(debug) {
					cout << "parse_blast here-107\n";
				}
				for (int64_t i = 0; i < num_space; ++i) {
					aligned_read = " " + aligned_read;
				}
				if (!aligned_read.empty()) {
					if(debug) {
						cout << "parse_blast here-108\n";
					}
					aligned_reads.push_back(aligned_read);
				}

//				#print "alg_top2[key][0][0]\tnum_space\ttrim_size\tref_seq_ori1\tref_seq_ori2\n";
//				#print "alg_top2[key][0].query_start\talg_top2[key][1].query_end\tdist\n";
			}
		}
	}
	if(debug) {
		cout << "parse_blast here-109\n";
	}
//	int64_t_value_desc_sortedset sorted_freq_s1;
	vector<pair<int64_t, int64_t>> sorted_freq_s1;
	for (auto& an_entry : freq_s1) {
//		Int64Pair a_pair(an_entry.first, an_entry.second);
//		sorted_freq_s1.insert(a_pair);
		sorted_freq_s1.push_back(make_pair(an_entry.first, an_entry.second));
	}
	sort(sorted_freq_s1.begin(), sorted_freq_s1.end(), [&](const pair<int64_t, int64_t>& lhs, const pair<int64_t, int64_t>& rhs)->bool {
		if(lhs.second > rhs.second) {
			return true;
		} else if(lhs.second < rhs.second) {
			return false;
		}
		if(abs(lhs.first) < abs(rhs.first)) {
			return true;
		}
		return false;
	});

//	for (auto an_entry : freq_s1) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-109-a: freq_s1 " << i << ":" << an_entry.second << "\n";
//	}
//	for (auto an_entry : sorted_freq_s1) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-109-b: sorted_freq_s1 " << i << ":" << an_entry.second << "\n";
//	}

	int64_t fqi = 0;
//	foreach my i (sort {freq_s1[b} <=> freq_s1[a}} keys %freq_s1) {
	for (auto& an_entry : sorted_freq_s1) {
		int64_t i = an_entry.first;
		++fqi;
		if (fqi == 1) {
			results.push_back(i);
			continue;
		}
		if (abs(results[0] - i) <= 2) {
			continue;
		}
		if (an_entry.second >= support_reads) {
			results2.push_back(i);
		}

//		#print "s1 i\tfreq_s1[i}\n";
		break;
	}
	fqi = 0;
	if(debug) {
		cout << "parse_blast here-110\n";
	}

//	int64_t_value_desc_sortedset sorted_freq_s2;
	vector<pair<int64_t, int64_t>> sorted_freq_s2;
	for (auto an_entry : freq_s2) {
//		Int64Pair a_pair(an_entry.first, an_entry.second);
//		sorted_freq_s2.insert(a_pair);
		sorted_freq_s2.push_back(make_pair(an_entry.first, an_entry.second));
	}
	sort(sorted_freq_s2.begin(), sorted_freq_s2.end(), [&](const pair<int64_t, int64_t>& lhs, const pair<int64_t, int64_t>& rhs)->bool {
		if(lhs.second > rhs.second) {
			return true;
		} else if(lhs.second < rhs.second) {
			return false;
		}
		if(abs(lhs.first) > abs(rhs.first)) {
			return true;
		}
		return false;
	});

//	for (auto an_entry : freq_s2) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-110-a: freq_s2 " << i << ":" << an_entry.second << "\n";
//	}
//	for (auto an_entry : sorted_freq_s2) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-110-b: sorted_freq_s2 " << i << ":" << an_entry.second << "\n";
//	}
//	foreach my i (sort {freq_s2[b} <=> freq_s2[a}} keys %freq_s2) {
	for (auto& an_entry : sorted_freq_s2) {
		int64_t i = an_entry.first;
		++fqi;
		if (abs(results[1] - i) <= 2) {
			continue;
		}
		if (fqi == 1) {
			results.push_back(i);
			continue;
		}
		if (an_entry.second >= support_reads) {
			results2.push_back(i);
		}

//		#print "s2 i\tfreq_s2[i}\n";
		break;
	}
	if(debug) {
		cout << "parse_blast here-111\n";
	}
	vector<int64_t> raw_results = results;

//	#print "@results\n";
//	# adjust break points if query break points having overlap or gap
//	# adjust break point
	vector<pair<int64_t, int64_t>> sorted_freq_dist;
	for (auto& an_entry : freq_dist) {
		sorted_freq_dist.push_back(make_pair(an_entry.first, an_entry.second));
	}
	sort(sorted_freq_dist.begin(), sorted_freq_dist.end(), [&](const pair<int64_t, int64_t>& lhs, const pair<int64_t, int64_t>& rhs)->bool {
		if(lhs.second > rhs.second) {
			return true;
		} else if(lhs.second < rhs.second) {
			return false;
		}
		if(abs(lhs.first) < abs(rhs.first)) {
			return true;
		}
		return false;
	});
//	for (auto an_entry : freq_dist) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-111-a: freq_dist " << i << ":" << an_entry.second << "\n";
//	}
//	for (auto an_entry : sorted_freq_dist) {
//		int64_t i = an_entry.first;
//		cout << "parse_blast here-111-b: sorted_freq_dist " << i << ":" << an_entry.second << "\n";
//	}
	if(debug) {
		cout << "parse_blast here-112\n";
	}
//	foreach my i (sort {freq_dist{b} <=> freq_dist{a}} keys %freq_dist) {
	for (auto an_entry : sorted_freq_dist) {
		int64_t i = an_entry.first;
//		#print "i\t@results\n";
		if (i < 0) {
			if (i % 2) {
				results[0] -= -i / 2 + 0.5;
				results[1] += -i / 2 - 0.5;
			} else {
				results[0] -= -i / 2;
				results[1] += -i / 2;
			}
		}
		break;
	}

//	# get largest homology
	int64_t max_freq = 0;
	for (auto an_entry : sorted_freq_dist) {
		max_freq = an_entry.second;
		break;
	}
	if(debug) {
		cout << "parse_blast here-113 max_freq: " << max_freq << "\n";
	}
//	foreach my i (sort {freq_dist{b} <=> freq_dist{a}} keys %freq_dist) {
//		max_freq = freq_dist{i};
//		break;
//	}
	int64_t ins_size = 0;    //# positive number for insertion, negative number for homology
//	foreach my i (sort {freq_dist{b} <=> freq_dist{a}} keys %freq_dist) {
	for (auto an_entry : sorted_freq_dist) {
		int64_t i = an_entry.first;
//		#print "i\t@results\n";
		if(debug) {
			cout << "parse_blast here-113-a: " << an_entry.first << ":" << an_entry.second << "->" << ins_size << "\n";
		}
		if (an_entry.second < max_freq) {
			break;
		}
		if (!ins_size || i < ins_size) {
			ins_size = i;
		}
	}
	if (results.size() > 1) {
		sort(results.begin(), results.end());
	}
	if (results2.size() > 1) {
		sort(results2.begin(), results2.end());
	}
	if (0 != ins_size) {
		if(debug) {
			cout << "parse_blast here-115 ins_size: " << ins_size << "\n";
		}
		results.push_back(ins_size);
	}

	if (ref_seq.size() > 0) {
		SROUT << ref_seq[0] << "\n\n";
	}
	for (auto& an_entry : aligned_reads) {
		SROUT << an_entry << "\n";
	}

	int64_t pos1 = 0;
	if (raw_results.size() > 0) {
		if(debug) {
			cout << "parse_blast here-116-a: " << pos1 << "\n";
		}
		pos1 = raw_results[0] + 1;
		if(debug) {
			cout << "parse_blast here-116-b: " << pos1 << "\n";
		}
	}
	if (!ref_seq_ori1 && ref_seq.size() > 0) {
		if(debug) {
			cout << "parse_blast here-117-a: " << pos1 << "\n";
		}
		pos1 = ref_seq[0].size() - pos1 + 3;
		if(debug) {
			cout << "parse_blast here-117-b: " << pos1 << "\n";
		}
	}

	int64_t pos2 = 0;
	if (raw_results.size() > 1 && ref_seq.size() > 0) {
		if(debug) {
			cout << "parse_blast here-118-a: pos2: " << pos2 << "/" << raw_results[1] << "/" << ref_seq[0].size() << "\n";
		}
		pos2 = raw_results[1] - 100 - ref_seq[0].size();
		if(debug) {
			cout << "parse_blast here-118-b: pos2: " << pos2 << "\n";
		}
	}
	if (!ref_seq_ori2 && ref_seq.size() > 1) {
		if(debug) {
			cout << "parse_blast here-119-a: pos2: " << pos2 << "\n";
		}
		pos2 = ref_seq[1].size() - pos2 + 1;
		if(debug) {
			cout << "parse_blast here-119-b: pos2: " << pos2 << "\n";
		}
	}
	if (results.size() > 2) {
		pos2 -= results[2];
		if(debug) {
			cout << "parse_blast here-120: results[2]: " << results[2] << ", pos2: " << pos2 << "\n";
		}
	}
	// find the position of ref_seq[1].
	int64_t trim = pos2 - pos1;
//	if(trim < -100) {
//		trim -= 100;
//	}
	if (ref_seq.size() > 1) {
		if (trim > 0) {
//			string the_empty_space(100, ' ');
			if(debug) {
				cout << "parse_blast here-121: " << trim << "/" << ref_seq[1].size() << "\n";
			}

			if (trim < static_cast<int64_t>(ref_seq[1].size())) {
				ref_seq[1] = ref_seq[1].substr(trim);
			}
//			ref_seq[1] = the_empty_space + ref_seq[1];
		}
//		else if (trim < 0) {
//			string the_empty_space(-trim, ' ');
//			ref_seq[1] = the_empty_space + ref_seq[1];
//			if(debug) {
//				cout << "parse_blast here-122: " << pos1 << ", " << pos2 << "\n";
//				cout << "parse_blast here-123: " << trim << "\n";
//				cout << "parse_blast here-124: " << ref_seq[1] << "\n";
//			}
//
//		}
		if(debug) {
			cout << "parse_blast here-125\n";
		}
		SROUT << (boost::format("\n%s\n") % ref_seq[1]).str();
	}

	if (results.size() < 3) {
		results.resize(3);
	}
	if(debug) {
		cout << "parse_blast here-126: result " << results[0] << "/" << results[1] << "/" << results[2] << "\n";
		cout << "parse_blast here-126: raw_result " << raw_results[0] << "/" << raw_results[1] << "/" << raw_results[2] << "\n";
	}
//	#print "bltoutfile\n@results\n@results2\n";
//	#print "ref_seq_ori1\tref_seq_ori2\tpos1\tpos2\ttrim\n@aligned_reads\n";
//	return (\@results, \@results2);
}
} /* namespace meerkat */
