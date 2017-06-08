/*
 * SplitReadSVCaller.cpp
 *
 * Created on: Jul 8, 2016
 * Author: el174
 */
#include "SplitReadSVCaller.hpp"

namespace meerkat {

SplitReadSVCaller::SplitReadSVCaller() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

SplitReadSVCaller::~SplitReadSVCaller() {
}

void SplitReadSVCaller::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}
//#file format of prefix.sr.intra.out
//#del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size
//#del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion && insertion
//#del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion && deletion (2 col)
//#ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion && insertion
//#invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
//#invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
//#tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size

//#file format of prefix.sr.inter.out
//#del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size
//#ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size
//#transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster

//#left_bound_sr1: left bound in primary cluster
//#right_bound_sr1: right bound in primary cluster
//#left_bound_sr2: left bound in secondary cluster
//#right_bound_sr2: right bound in secondary cluster
//#primary cluster: cluster to initiate call in mpd
//#secondary cluster: cluster to accompany primary cluster in mpd

void SplitReadSVCaller::call_structural_variants() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.call_structural_variants");
	checker.start();
	BamTools::BamAlignment al;

	string bpinfofile = options.prefix + ".bp.info";
	string srintra_outfile = options.prefix + ".sr.intra.out";
	string srinter_outfile = options.prefix + ".sr.inter.out";
	string bp_readsfile = options.prefix + ".bp_reads";
	string sr_sortbam = options.prefix + ".sr.sorted.bam";

	if(!options.working_dir.empty()) {
		bpinfofile = options.working_prefix + ".bp.info";
		srintra_outfile = options.working_prefix + ".sr.intra.out";
		srinter_outfile = options.working_prefix + ".sr.inter.out";
		bp_readsfile = options.working_prefix + ".bp_reads";
		sr_sortbam = options.working_prefix + ".sr.sorted.bam";
	}
	if (boost::filesystem::exists(bp_readsfile)) {
		boost::filesystem::remove(bp_readsfile);
	}

	string a_path(sr_sortbam);
	string an_index_path(a_path);
	an_index_path += ".bai";
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}

	map<string, pair<int64_t, int64_t>> reverse_index_ref_id;
	const auto& ref_vec = reader.GetReferenceData();
	for (uint64_t ref_id = 0; ref_id < ref_vec.size(); ++ref_id) {
		reverse_index_ref_id[ref_vec[ref_id].RefName].first = ref_id;
		reverse_index_ref_id[ref_vec[ref_id].RefName].second = ref_vec[ref_id].RefLength;
	}

//# %cluster_region: index of corresponding region for each cluster
	map<string, CoordinateEntry> cluster_region;
	collect_bp_data(cluster_region);
	vector<EventEntry> result_sr;
	vector<EventEntry> unsupport_del;
	ofstream BPREAD(bp_readsfile, ios::binary);
//	detect_intra_chromosomal_events_alt(BPREAD, result_sr, unsupport_del, reader, reverse_index_ref_id, cluster_region);
	detect_intra_chromosomal_events(BPREAD, result_sr, unsupport_del, reader, reverse_index_ref_id, cluster_region);
//	detect_inter_chromosomal_events_alt(BPREAD, result_sr, unsupport_del, reader, reverse_index_ref_id, cluster_region);
	detect_inter_chromosomal_events(BPREAD, result_sr, unsupport_del, reader, reverse_index_ref_id, cluster_region);

	ofstream SRINTRAOUT(srintra_outfile, ios::binary);
	ofstream SRINTEROUT(srinter_outfile, ios::binary);
	for (auto& a_result : result_sr) {
		string toprint = a_result.sr_str();
		castle::StringUtils::trim(toprint);
		if (toprint.empty()) {
			continue;
		}
		if (string::npos != a_result.type.find("transl_inter") || "del_inss" == a_result.type || "del_inso" == a_result.type || "inss" == a_result.type || "inso" == a_result.type) {
			SRINTEROUT << toprint << "\n";
		} else {
			SRINTRAOUT << toprint << "\n";
		}
	}
}
void SplitReadSVCaller::collect_bp_data(map<string, CoordinateEntry>& cluster_region) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.collect_bp_data");
	checker.start();
	string bpinfofile = options.prefix + ".bp.info";
	if(!options.working_dir.empty()) {
		bpinfofile = options.working_prefix + ".bp.info";
	}
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	const char* delim_double_underscore = "__";
//	const char* delim_colon = ":";
	vector<string> data;
	vector<string> cl;
	vector<string> coords;
	string line;

	CoordinateEntry empty_entry;
	string empty_key;
	cluster_region[empty_key] = empty_entry;
	ifstream BPDTSRD(bpinfofile, ios::binary);
	while (getline(BPDTSRD, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[0], delim_slash, cl);
		castle::StringUtils::tokenize(data[1], delim_double_underscore, coords);
//		const bool debug = "1746_0/7577_0" == data[0];
//		if(debug) {
//			cout << "[SplitReadSVCaller.collect_bp_data] line: " << line << "\n";
		//cout << "[SplitReadSVCaller.collect_bp_data] data: " << data.size() << "/cl: " << cl.size() << "/coords:" << coords.size() << "\n";

//		}
		CoordinateEntry an_entry;
		an_entry.chr_bp_1 = coords[0];
		an_entry.start = boost::lexical_cast<int64_t>(coords[1]);
		an_entry.end = boost::lexical_cast<int64_t>(coords[2]);
		if (coords.size() > 3) {
			an_entry.chr_bp_2 = coords[3];
			an_entry.mate_start = boost::lexical_cast<int64_t>(coords[4]);
			an_entry.mate_end = boost::lexical_cast<int64_t>(coords[5]);
			an_entry.orientation = boost::lexical_cast<int32_t>(coords[6]);
		}
		an_entry.str = data[1];
		for (auto& a_cl : cl) {
//			if(debug) {
//				cout << "[SplitReadSVCaller.collect_bp_data] cl: " << a_cl << "\n";
//			}
			if(a_cl.empty()) {
				continue;
			} else {
				cluster_region[a_cl] = an_entry;
				if(string::npos == a_cl.rfind("_0")) {
//					cout << "[SplitReadSVCaller.collect_bp_data] cl: " << a_cl << "_0" << "\n";
					cluster_region[a_cl + "_0"] = an_entry;
				} else {
					auto the_suffix_pos = a_cl.rfind("_0");
					if(string::npos != the_suffix_pos) {
						auto temp_cl = a_cl.substr(0, the_suffix_pos);
						cluster_region[temp_cl] = an_entry;
					}
				}
			}
		}
	}
	cout << checker;
}

void SplitReadSVCaller::detect_intra_chromosomal_events(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.detect_intra_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}
	vector<int64_t> block_positions;

	{
		block_positions.push_back(0);
		int64_t n_lines = castle::IOUtils::get_number_of_lines(mpintra_outfile);
		int64_t BLOCK_SIZE = (n_lines + n_cores - 1) / (double) n_cores;
		string line;
		ifstream in(mpintra_outfile);
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
	cout << (boost::format("[SplitReadSVCaller.detect_intra_chromosomal_events] # blocks: %d\n") % n_blocks).str();
	vector<vector<EventEntry>> result_sr_list(n_blocks);
	vector<vector<EventEntry>> unsupport_del_list(n_blocks);
	vector<string> output_file_names(n_blocks);

	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string sr_sortbam = options.prefix + ".sr.sorted.bam";
			if(!options.working_dir.empty()) {
				sr_sortbam = options.working_prefix + ".sr.sorted.bam";
			}
			string a_path(sr_sortbam);
			string an_index_path(a_path);
			an_index_path += ".bai";
			BamTools::BamReader reader;
			if (!reader.Open(a_path, an_index_path)) {
				std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
				exit(1);
			}

			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = block_positions[block_id];
			int64_t next_pos = block_positions[block_id + 1];
			string line;
			const char* delim_tab = "\t";
			const char* delim_slash = "/";
			vector<string> data;
			vector<string> cl;
			vector<string> coords;
			vector<string> cluster_ids;
			vector<string> mpds;
			vector<string> temp_cols;
			string local_bpinfofile = options.prefix + ".bp_reads." + str_block_id;
			if(!options.working_dir.empty()) {
				local_bpinfofile = options.working_prefix + ".bp_reads." + str_block_id;
			}
			output_file_names[block_id] = local_bpinfofile;
			auto& local_result_sr = result_sr_list[block_id];
			auto& local_unsupport_del = unsupport_del_list[block_id];
			ofstream local_out(local_bpinfofile, ios::binary);
			ifstream local_reader(mpintra_outfile, ios::binary);
			local_reader.seekg(cur_pos, ios::beg);
			while(getline(local_reader, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
				if (cluster_region.end() == cluster_region.find(cl[0])) {
					continue;
				}
				map<int64_t, int32_t> bp_window;
				if (cl.size() > 1) {
					const auto& a_region_a = cluster_region.find(cl[0]);
					const auto& a_region_b = cluster_region.find(cl[1]);
					int64_t positiona1 = a_region_a->second.start;
					int64_t positionb1 = a_region_b->second.start;
					if (positiona1 > positionb1) {
						swap(cl[0], cl[1]);
					}
				}
				//# the window in candidate region where a break point is located
				//# bp_window{cur_start} = 1/2, 1 left window, 2 right window
				if (data[0] == "del") {
					detect_del(local_out, local_result_sr, local_unsupport_del, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds, temp_cols);
				} else if (string::npos != data[0].find("inssd")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_inssd(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (string::npos != data[0].find("inssu")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_inssu(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (string::npos != data[0].find("insod")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_insod(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (string::npos != data[0].find("insou")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_insou(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (data[0] == "invers" || data[0] == "del_invers") {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_invers(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (data[0] == "tandem_dup") {
					detect_tandem_dup(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (data[0] == "invers_f") {
					detect_invers_f(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (data[0] == "invers_r") {
					detect_invers_r(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		BPREAD << castle::IOUtils::read_fully(output_file_names[block_id]);
	}
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		auto& local_result_sr = result_sr_list[block_id];
		result_sr.insert(result_sr.end(), local_result_sr.begin(), local_result_sr.end());
		auto& local_unsupport_del = unsupport_del_list[block_id];
		unsupport_del.insert(unsupport_del.end(), local_unsupport_del.begin(), local_unsupport_del.end());
	}
	castle::IOUtils::remove_files(output_file_names, n_cores);
	cout << checker;
}
void SplitReadSVCaller::detect_intra_chromosomal_events_alt(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.detect_intra_chromosomal_events_alt");
	checker.start();
	bool debug = false;
	string line;
	int64_t i = 0;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	vector<string> data;
	vector<string> cl;
	vector<string> coords;
	vector<string> cluster_ids;
	vector<string> mpds;
	vector<string> temp_cols;
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}
	ifstream MPINTRASRD(mpintra_outfile, ios::binary);

	while (getline(MPINTRASRD, line, '\n')) {

		++i;
		//if (-1 != line_value && i < line_value) {
		//continue;
		//}
		//if (-1 != line_value && i > line_value) {
		//break;
		//}
		castle::StringUtils::tokenize(line, delim_tab, data);
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		if (cluster_region.end() == cluster_region.find(cl[0])) {
			continue;
		}
		cout << line << "\n";

	if(data[0] == "insou" && data[1] == "132200_0/277599_0") {
		debug = true;
	}
	if(!debug) {
		continue;
	}
		map<int64_t, int32_t> bp_window;
//		if (cl.size() > 1) {
//			const auto& a_region_a = cluster_region.find(cl[0]);
//			const auto& a_region_b = cluster_region.find(cl[1]);
//
//			int64_t positiona1 = a_region_a->second.start;
//			int64_t positionb1 = a_region_b->second.start;
////			if(debug) {
////				cout << "del_insou: " << cl[0] << "/" << cl[1] << "/\n";
////				cout << "del_insou: " << positiona1 << "/" << positionb1 << "\n";
////			}
//			if (positiona1 > positionb1) {
//				swap(cl[0], cl[1]);
//			}
//		}
		//# the window in candidate region where a break point is located
		//# bp_window{cur_start} = 1/2, 1 left window, 2 right window
		if (data[0] == "del") {
			detect_del(BPREAD, result_sr, unsupport_del, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds, temp_cols);
		} else if (string::npos != data[0].find("inssd")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_inssd(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (string::npos != data[0].find("inssu")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_inssu(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (string::npos != data[0].find("insod")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_insod(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (string::npos != data[0].find("insou")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_insou(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (data[0] == "invers" || data[0] == "del_invers") {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_invers(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (data[0] == "tandem_dup") {
			detect_tandem_dup(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (data[0] == "invers_f") {
			detect_invers_f(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (data[0] == "invers_r") {
			detect_invers_r(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		}
//		if (i > 1000) {
////		if (i > 1000 || ("inssu" == data[0] && "7891_0/23425_0" == data[1])) {
//			break;
//		}
	}
	cout << checker;
}
void SplitReadSVCaller::detect_del(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
		vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds, vector<string>& temp_cols) {
//	const bool debug = "3846" == data[1];
	const bool debug = false;
	auto& ref_is = options.is;
//	const char* delim_slash = "/";
	const char* delim_colon = ":";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	int64_t position1 = a_region_a->second.start;
	int64_t position2 = a_region_a->second.end;
	int64_t position3 = a_region_a->second.mate_start;
	int64_t position4 = a_region_a->second.mate_end;
	if (covered(cur_start, cur_start, position1, position2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, position3, position4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_end, cur_end, position1, position2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, position3, position4)) {
		bp_window[cur_end] = 2;
	}

	//#print "a_region_a->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\n";
	//my @alignments;
	//open( SAM,
	//"samtools_command view -X sr_sortbam a_region_a->second|"
	//);
	//while ( newline1 = <SAM> ) {
	//chomp newline1;
	//push @alignments, newline1;
	//}
	//close SAM;

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

	//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
	//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
	if(debug) {
		cout << "del here-0\n";
	}
	if (bp_window[cur_start] != bp_window[cur_end]) {
		discord_sr1.clear();
		if(debug) {
			cout << "del here-1\n";
		}
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			//int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			//int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand) {

				//#print "start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "del here-2\n";
			}
			vector<int64_t> ref_boundary1;
			vector<int64_t> ref_boundary2;
			vector<int64_t> ref_support_sr;
			map<int64_t, vector<string>> ref_bpread;
			sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

			//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
			//# large deletion with small close insertion
			if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
				if(debug) {
					cout << "del here-3\n";
				}
				//# inssu
				if (covered(ref_boundary1[0], ref_boundary1[1], ref_boundary1[2], ref_boundary1[3])
						&& ((ref_boundary2[0] >= (position2 - position1 + 101) && ref_boundary2[2] <= (position2 - position1 + 101) && ref_boundary1[0] <= (position2 - position1 + 101) && ref_boundary1[2] <= (position2 - position1 + 101))
								|| (ref_boundary2[2] >= (position2 - position1 + 101) && ref_boundary2[0] <= (position2 - position1 + 101) && ref_boundary1[0] <= (position2 - position1 + 101) && ref_boundary1[2] <= (position2 - position1 + 101)))) {
					if(debug) {
						cout << "del here-4\n";
					}
					if (ref_boundary2.size() > 2 && ref_boundary2[0] < ref_boundary2[2]) {
						if(debug) {
							cout << "del here-5\n";
						}
						int64_t left_bound_sr1 = 0;
						if(ref_boundary1.size() > 0) {
							left_bound_sr1 = position1 + ref_boundary1[0];
						}
						int64_t right_bound_sr1 = 0;
						if(ref_boundary2.size() > 1) {
							right_bound_sr1 = position1 + ref_boundary2[1] + options.cut_sr;
						}
						int64_t left_bound_sr2 = 0;
						if(ref_boundary1.size() > 3) {
							left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
						}
						int64_t right_bound_sr2 = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
						}
						int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
						int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "del here-5-a\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
							EventEntry an_entry;
							if (del_size > 10) {
								an_entry.type = "del_inssu";
							} else {
								an_entry.type = "inssu";
							}
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[1];
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = right_bound_sr1 + 1;
							an_entry.event_end = right_bound_sr2 + 1;
							an_entry.event_size_1 = del_size;
							an_entry.mate_ref_id = data[3];
							an_entry.mate_event_start = left_bound_sr1 + 1;
							an_entry.mate_event_end = left_bound_sr2 + 1;
							an_entry.event_size_2 = ins_size;
							result_sr.push_back(an_entry);

							//if (del_size > 10) {
							//result_sr[sri][0] = "del_inssu";
							//} else {
							//result_sr[sri][0] = "inssu";
							//}
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = right_bound_sr1;
							//result_sr[sri][6] = right_bound_sr2;
							//result_sr[sri][7] = del_size;
							//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
						}
					} else {
						if(debug) {
							cout << "del here-6\n";
						}
						int64_t left_bound_sr1 = 0;
						if (ref_boundary1.size() > 2) {
							left_bound_sr1 = position1 + ref_boundary1[2];
						}
						int64_t right_bound_sr1 = 0;
						if (ref_boundary2.size() > 3)
							right_bound_sr1 = position1 + ref_boundary2[3] + options.cut_sr;
						int64_t left_bound_sr2 = 0;
						if (ref_boundary1.size() > 1)
							left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
						int64_t right_bound_sr2 = 0;
						if (ref_boundary2.size() > 0) {
							right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
						}
						int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
						int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "del here-7\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();

							EventEntry an_entry;
							if (del_size > 10) {
								an_entry.type = "del_inssu";
							} else {
								an_entry.type = "inssu";
							}
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[1];
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = right_bound_sr1 + 1;
							an_entry.event_end = right_bound_sr2 + 1;
							an_entry.event_size_1 = del_size;
							an_entry.mate_ref_id = data[3];
							an_entry.mate_event_start = left_bound_sr1 + 1;
							an_entry.mate_event_end = left_bound_sr2 + 1;
							an_entry.event_size_2 = ins_size;
							result_sr.push_back(an_entry);

							//if (del_size > 10) {
							//result_sr[sri][0] = "del_inssu";
							//} else {
							//result_sr[sri][0] = "inssu";
							//}
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = right_bound_sr1;
							//result_sr[sri][6] = right_bound_sr2;
							//result_sr[sri][7] = del_size;
							//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
						}
					}
				}

				//# inssd
				else if (ref_boundary1.size() > 2 && ref_boundary2.size() > 3 && covered(ref_boundary2[0], ref_boundary2[1], ref_boundary2[2], ref_boundary2[3])
						&& ((ref_boundary1[0] <= (position2 - position1 + 101) && ref_boundary1[2] >= (position2 - position1 + 101) && ref_boundary2[0] >= (position2 - position1 + 101) && ref_boundary2[2] >= (position2 - position1 + 101))
								|| (ref_boundary1[2] <= (position2 - position1 + 101) && ref_boundary1[0] >= (position2 - position1 + 101) && ref_boundary2[0] >= (position2 - position1 + 101) && ref_boundary2[2] >= (position2 - position1 + 101)))) {
					if(debug) {
						cout << "del here-8\n";
					}
					if (ref_boundary1[0] < ref_boundary1[2]) {
						if(debug) {
							cout << "del here-9\n";
						}
						int64_t left_bound_sr1 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
						int64_t right_bound_sr1 = ref_boundary2[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
						int64_t left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
						int64_t right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
						int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
						int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "del here-9-a\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
							EventEntry an_entry;
							if (del_size > 10) {
								an_entry.type = "del_inssd";
							} else {
								an_entry.type = "inssd";
							}
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr2 + 1;
							an_entry.event_end = left_bound_sr1 + 1;
							an_entry.event_size_1 = del_size;
							an_entry.mate_ref_id = data[3];
							an_entry.mate_event_start = right_bound_sr2 + 1;
							an_entry.mate_event_end = right_bound_sr1 + 1;
							an_entry.event_size_2 = ins_size;
							result_sr.push_back(an_entry);
							//if (del_size > 10) {
							//result_sr[sri][0] = "del_inssd";
							//} else {
							//result_sr[sri][0] = "inssd";
							//}
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = left_bound_sr2;
							//result_sr[sri][6] = left_bound_sr1;
							//result_sr[sri][7] = del_size;
							//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
						}
					} else {
						if(debug) {
							cout << "del here-10\n";
						}
						int64_t left_bound_sr1 = 0;
						if(ref_boundary1.size() > 0) {
							left_bound_sr1 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
						}
						int64_t right_bound_sr1 = 0;
						if(ref_boundary2.size() > 1) {
							right_bound_sr1 = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
						}
						int64_t left_bound_sr2 = 0;
						if(ref_boundary1.size() > 3) {
							left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
						}
						int64_t right_bound_sr2 = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
						}
						int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
						int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "del here-11\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
							EventEntry an_entry;
							if (del_size > 10) {
								an_entry.type = "del_inssd";
							} else {
								an_entry.type = "inssd";
							}
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr2 + 1;
							an_entry.event_end = left_bound_sr1 + 1;
							an_entry.event_size_1 = del_size;
							an_entry.mate_ref_id = data[3];
							an_entry.mate_event_start = right_bound_sr2 + 1;
							an_entry.mate_event_end = right_bound_sr1 + 1;
							an_entry.event_size_2 = ins_size;
							result_sr.push_back(an_entry);
							//if (del_size > 10) {
							//result_sr[sri][0] = "del_inssd";
							//} else {
							//result_sr[sri][0] = "inssd";
							//}
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = left_bound_sr2;
							//result_sr[sri][6] = left_bound_sr1;
							//result_sr[sri][7] = del_size;
							//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							//
						}
					}
				}

				//# del
				else {
					if(debug) {
						cout << "del here-12\n";
					}
					//# re-process deletion to require direction
					discord_sr1.clear();
					reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
					BamTools::BamAlignment al;
					//int64_t prev_bam_pos = backend_bgzf.Tell();
					while (reader.GetNextAlignmentBasic(al)) {
						int64_t start = al.Position;
						int32_t strand = 1;
						if (al.IsReverseStrand()) {
							strand = -1;
						}
//	string mseqid = ref_vec[al.MateRefID].RefName;
						int64_t mstart = al.MatePosition;
						int32_t mstrand = 1;
						if (al.IsMateReverseStrand()) {
							mstrand = -1;
						}
						int64_t isize = al.InsertSize;

						if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (position2 - position1) && mstart >= (position2 - position1 + 100)) {
							if ((al.IsFirstMate() && strand == 1) || (al.IsSecondMate() && strand == -1)) {

								//#print "start mstart\tisize\tstrand mstrand\n";
								discord_sr1.push_back(al);
							}
						}
					}
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;
					sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

					//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
					if(debug) {
						cout << "del here-12-b: " << ref_boundary1.size() << "/" << ref_boundary2.size() << "/" << ref_support_sr.size() << "/" << ref_bpread.size() << "/" << discord_sr1.size() << "\n";
					}
					int64_t left_bound_sr = 0;
					int64_t right_bound_sr = 0;
					if (ref_boundary1.size() > 1 && ref_boundary2.size() > 0 && ref_boundary1[1] < (position2 - position1 + 101) && ref_boundary2[0] >= (position2 - position1 + 101)) {
						left_bound_sr = position1 + ref_boundary1[1] + options.cut_sr;
						right_bound_sr = ref_boundary2[0] - (position2 - position1 + 101) + position3;
					}
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 2 && ref_boundary1[3] < (position2 - position1 + 101) && ref_boundary2[2] >= (position2 - position1 + 101)) {
						left_bound_sr = position1 + ref_boundary1[3] + options.cut_sr;
						right_bound_sr = ref_boundary2[2] - (position2 - position1 + 101) + position3;
					}
					int64_t del_size = right_bound_sr - left_bound_sr - 1;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr && right_bound_sr) {
						if(debug) {
							cout << "del here-13\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = del_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = ref_support_sr[0];
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr;
						//result_sr[sri][6] = right_bound_sr;
						//result_sr[sri][7] = "del_size\n";
					} else {
						if(debug) {
							cout << "del here-14\n";
						}
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.ref_id = data[3];
						an_entry.event_start = cur_start + 1;
						an_entry.event_end = cur_end + 1;
						an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
						castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
						for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
							an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
						}
						unsupport_del.push_back(an_entry);
						//unsupport_del[usdi][0] = data[0];
						//unsupport_del[usdi][1] = data[1];
						//unsupport_del[usdi][2] = data[2];
						//unsupport_del[usdi][3] = data[3];
						//unsupport_del[usdi][4] = cur_start;
						//unsupport_del[usdi][5] = cur_end;
						//unsupport_del[usdi][6] = data[6];
						//unsupport_del[usdi][7] = "data[7]\n";
					}
				}
			}

			//# del
			else {
				if(debug) {
					cout << "del here-15\n";
				}
				//# re-process deletion to require direction
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (position2 - position1) && mstart >= (position2 - position1 + 100)) {
						if ((al.IsFirstMate() && strand == 1) || (al.IsSecondMate() && strand == -1)) {

							//#print "start mstart\tisize\tstrand mstrand\n";
							discord_sr1.push_back(al);
						}
					}
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if(debug) {
					cout << "del here-16\n";
				}
				int64_t left_bound_sr = 0;
				if (ref_boundary1.size() > 1 && ref_boundary1[1] < (position2 - position1 + 101)) {
					left_bound_sr = position1 + ref_boundary1[1] + options.cut_sr;
				}

				int64_t right_bound_sr = 0;
				if (ref_boundary2.size() > 0 && ref_boundary2[0] >= (position2 - position1 + 101)) {
					right_bound_sr = ref_boundary2[0] - (position2 - position1 + 101) + position3;
				}
				if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
					left_bound_sr = cur_start;
				}
				if (!ref_boundary2.empty() && !ref_boundary2[0]) {
					right_bound_sr = cur_end;
				}
				int64_t del_size = right_bound_sr - left_bound_sr - 1;
				if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr && right_bound_sr) {
					if(debug) {
						cout << "del here-17\n";
					}
					string reads = castle::StringUtils::join(ref_bpread[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[0];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr + 1;
					an_entry.event_end = right_bound_sr + 1;
					an_entry.event_size_1 = del_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = ref_support_sr[0];
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = left_bound_sr;
					//result_sr[sri][6] = right_bound_sr;
					//result_sr[sri][7] = "del_size\n";
				} else {
					if(debug) {
						cout << "del here-18\n";
					}
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.ref_id = data[3];
					an_entry.event_start = cur_start;
					an_entry.event_end = cur_end;
					an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
					castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
					for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
						an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
					}
					unsupport_del.push_back(an_entry);
					//unsupport_del[usdi][0] = data[0];
					//unsupport_del[usdi][1] = data[1];
					//unsupport_del[usdi][2] = data[2];
					//unsupport_del[usdi][3] = data[3];
					//unsupport_del[usdi][4] = cur_start;
					//unsupport_del[usdi][5] = cur_end;
					//unsupport_del[usdi][6] = data[6];
					//unsupport_del[usdi][7] = "data[7]\n";
				}
			}
		} else {
			if(debug) {
				cout << "del here-19\n";
			}
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.ref_id = data[3];
			an_entry.event_start = cur_start;
			an_entry.event_end = cur_end;
			an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
			castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
			for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
				an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
			}
			unsupport_del.push_back(an_entry);
			//unsupport_del[usdi][0] = data[0];
			//unsupport_del[usdi][1] = data[1];
			//unsupport_del[usdi][2] = data[2];
			//unsupport_del[usdi][3] = data[3];
			//unsupport_del[usdi][4] = cur_start;
			//unsupport_del[usdi][5] = cur_end;
			//unsupport_del[usdi][6] = data[6];
			//unsupport_del[usdi][7] = "data[7]\n";
		}
	} else {
		if(debug) {
			cout << "del here-20\n";
		}
		if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1) {
			if(debug) {
				cout << "del here-21\n";
			}
			discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (position2 - position1) && mstart < (position2 - position1)) {

					//#print "start mstart\tisize\tstrand mstrand flags\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "del here-22\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "del here-23\n";
					}
					//# inssu
					if (ref_boundary1.size() > 3 && covered(ref_boundary1[0], ref_boundary1[1], ref_boundary1[2], ref_boundary1[3])) {
						if (ref_boundary2.size() > 2 && ref_boundary2[0] < ref_boundary2[2]) {
							int64_t left_bound_sr1 = position1 + ref_boundary1[0];
							int64_t right_bound_sr1 = position1 + ref_boundary2[1] + options.cut_sr;
							int64_t left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
							int64_t right_bound_sr2 = position1 + ref_boundary2[2];
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-24\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssu";
								//} else {
								//result_sr[sri][0] = "inssu";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = right_bound_sr1;
								//result_sr[sri][6] = right_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
							}
						} else {
							if(debug) {
								cout << "del here-25\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = position1 + ref_boundary1[2];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = position1 + ref_boundary2[3] + options.cut_sr;
							}

							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = position1 + ref_boundary2[0];
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-25-a\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssu";
								//} else {
								//result_sr[sri][0] = "inssu";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = right_bound_sr1;
								//result_sr[sri][6] = right_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
							}
						}
					}

					//# inssd
					else if (ref_boundary2.size() > 3 && covered(ref_boundary2[0], ref_boundary2[1], ref_boundary2[2], ref_boundary2[3])) {
						if(debug) {
							cout << "del here-26\n";
						}
						if (ref_boundary1.size() > 2 && ref_boundary1[0] < ref_boundary1[2]) {
							int64_t left_bound_sr1 = position1 + ref_boundary1[2];
							int64_t right_bound_sr1 = position1 + ref_boundary2[3] + options.cut_sr;
							int64_t left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
							int64_t right_bound_sr2 = position1 + ref_boundary2[0];
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-26-a\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssd";
								//} else {
								//result_sr[sri][0] = "inssd";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr2;
								//result_sr[sri][6] = left_bound_sr1;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							}
						} else {
							if(debug) {
								cout << "del here-27\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = position1 + ref_boundary1[0];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = position1 + ref_boundary2[1] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = position1 + ref_boundary2[2];
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-27-a\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssd";
								//} else {
								//result_sr[sri][0] = "inssd";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr2;
								//result_sr[sri][6] = left_bound_sr1;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							}
						}
					} else {
						if(debug) {
							cout << "del here-28\n";
						}
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 1) {
							left_bound_sr = position1 + ref_boundary1[1] + options.cut_sr;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 0) {
							right_bound_sr = position1 + ref_boundary2[0];
						}
						if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
							left_bound_sr = cur_start;
						}
						if (!ref_boundary2.empty() && !ref_boundary2[0]) {
							right_bound_sr = cur_end;
						}
						int64_t del_size = right_bound_sr - left_bound_sr - 1;
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "del here-28-a\n";
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = data[0];
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = del_size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = data[0];
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = ref_support_sr[0];
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = left_bound_sr;
							//result_sr[sri][6] = right_bound_sr;
							//result_sr[sri][7] = "del_size\n";
						} else {
							if(debug) {
								cout << "del here-29\n";
							}
							EventEntry an_entry;
							an_entry.type = data[0];
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.ref_id = data[3];
							an_entry.event_start = cur_start;
							an_entry.event_end = cur_end;
							an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
							castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
							for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
								an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
							}
							unsupport_del.push_back(an_entry);
							//unsupport_del[usdi][0] = data[0];
							//unsupport_del[usdi][1] = data[1];
							//unsupport_del[usdi][2] = data[2];
							//unsupport_del[usdi][3] = data[3];
							//unsupport_del[usdi][4] = cur_start;
							//unsupport_del[usdi][5] = cur_end;
							//unsupport_del[usdi][6] = data[6];
							//unsupport_del[usdi][7] = "data[7]\n";
						}
					}
				} else {
					if(debug) {
						cout << "del here-30\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 1) {
						left_bound_sr = position1 + ref_boundary1[1] + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 0) {
						right_bound_sr = position1 + ref_boundary2[0];
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						left_bound_sr = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr = cur_end;
					}
					int64_t del_size = right_bound_sr - left_bound_sr - 1;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "del here-30-a\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] %
								(right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = del_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = ref_support_sr[0];
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr;
						//result_sr[sri][6] = right_bound_sr;
						//result_sr[sri][7] = "del_size\n";
					} else {
						if(debug) {
							cout << "del here-31\n";
						}
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.ref_id = data[3];
						an_entry.event_start = cur_start;
						an_entry.event_end = cur_end;
						an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
						castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
						for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
							an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
						}
						unsupport_del.push_back(an_entry);
						//unsupport_del[usdi][0] = data[0];
						//unsupport_del[usdi][1] = data[1];
						//unsupport_del[usdi][2] = data[2];
						//unsupport_del[usdi][3] = data[3];
						//unsupport_del[usdi][4] = cur_start;
						//unsupport_del[usdi][5] = cur_end;
						//unsupport_del[usdi][6] = data[6];
						//unsupport_del[usdi][7] = "data[7]\n";
					}
				}
			} else {
				if(debug) {
					cout << "del here-32\n";
				}
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.ref_id = data[3];
				an_entry.event_start = cur_start;
				an_entry.event_end = cur_end;
				an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
				for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
					an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
				}
				unsupport_del.push_back(an_entry);
				//unsupport_del[usdi][0] = data[0];
				//unsupport_del[usdi][1] = data[1];
				//unsupport_del[usdi][2] = data[2];
				//unsupport_del[usdi][3] = data[3];
				//unsupport_del[usdi][4] = cur_start;
				//unsupport_del[usdi][5] = cur_end;
				//unsupport_del[usdi][6] = data[6];
				//unsupport_del[usdi][7] = "data[7]\n";
			}
		}

		//# both break points in window 2
		else {
			if(debug) {
				cout << "del here-32\n";
			}
			discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start >= (position2 - position1 + 100) && mstart >= (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "del here-33\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "del here-34\n";
					}
					//# inssu
					if (ref_boundary1.size() > 3 && covered(ref_boundary1[0], ref_boundary1[1], ref_boundary1[2], ref_boundary1[3])) {
						if(debug) {
							cout << "del here-35\n";
						}
						if (ref_boundary2.size() > 2 && ref_boundary2[0] < ref_boundary2[2]) {
							if(debug) {
								cout << "del here-36\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-36-a\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssu";
								//} else {
								//result_sr[sri][0] = "inssu";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = right_bound_sr1;
								//result_sr[sri][6] = right_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
							}
						} else {
							if(debug) {
								cout << "del here-37\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-38\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssu";
								//} else {
								//result_sr[sri][0] = "inssu";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = right_bound_sr1;
								//result_sr[sri][6] = right_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
							}
						}
					}
					//# inssd
					else if (ref_boundary2.size() > 3 && covered(ref_boundary2[0], ref_boundary2[1], ref_boundary2[2], ref_boundary2[3])) {
						if(debug) {
							cout << "del here-39\n";
						}
						if (ref_boundary1.size() > 3 && ref_boundary1[0] < ref_boundary1[2]) {
							if(debug) {
								cout << "del here-40\n";
							}
							int64_t left_bound_sr1 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							int64_t right_bound_sr1 = ref_boundary2[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							int64_t left_bound_sr2 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							int64_t right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-41\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssd";
								//} else {
								//result_sr[sri][0] = "inssd";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr2;
								//result_sr[sri][6] = left_bound_sr1;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							}
						} else {
							if(debug) {
								cout << "del here-42\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "del here-42-a\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//if (del_size > 10) {
								//result_sr[sri][0] = "del_inssd";
								//} else {
								//result_sr[sri][0] = "inssd";
								//}
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr2;
								//result_sr[sri][6] = left_bound_sr1;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
								//
							}
						}
					} else {
						if(debug) {
							cout << "del here-43\n";
						}
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 1) {
							left_bound_sr = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 0) {
							right_bound_sr = ref_boundary2[0] - (position2 - position1 + 101) + position3;
						}
						if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
							left_bound_sr = cur_start;
						}
						if (!ref_boundary2.empty() && !ref_boundary2[0]) {
							right_bound_sr = cur_end;
						}
						int64_t del_size = right_bound_sr - left_bound_sr - 1;
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "del here-43-a\n";
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = data[0];
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = del_size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = data[0];
							//result_sr[sri][1] = data[1];
							//result_sr[sri][2] = data[2];
							//result_sr[sri][3] = ref_support_sr[0];
							//result_sr[sri][4] = data[3];
							//result_sr[sri][5] = left_bound_sr;
							//result_sr[sri][6] = right_bound_sr;
							//result_sr[sri][7] = "del_size\n";
						} else {
							if(debug) {
								cout << "del here-44\n";
							}
							EventEntry an_entry;
							an_entry.type = data[0];
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.ref_id = data[3];
							an_entry.event_start = cur_start;
							an_entry.event_end = cur_end;
							an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
							castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
							for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
								an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
							}
							unsupport_del.push_back(an_entry);
							//unsupport_del[usdi][0] = data[0];
							//unsupport_del[usdi][1] = data[1];
							//unsupport_del[usdi][2] = data[2];
							//unsupport_del[usdi][3] = data[3];
							//unsupport_del[usdi][4] = cur_start;
							//unsupport_del[usdi][5] = cur_end;
							//unsupport_del[usdi][6] = data[6];
							//unsupport_del[usdi][7] = "data[7]\n";
						}
					}
				} else {
					if(debug) {
						cout << "del here-45\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 1) {
						left_bound_sr = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 0) {
						right_bound_sr = ref_boundary2[0] - (position2 - position1 + 101) + position3;
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						left_bound_sr = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr = cur_end;
					}
					int64_t del_size = right_bound_sr - left_bound_sr - 1;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "del here-46\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = del_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = ref_support_sr[0];
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr;
						//result_sr[sri][6] = right_bound_sr;
						//result_sr[sri][7] = "del_size\n";
					} else {
						if(debug) {
							cout << "del here-47\n";
						}
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.ref_id = data[3];
						an_entry.event_start = cur_start;
						an_entry.event_end = cur_end;
						an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
						castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
						for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
							an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
						}
						unsupport_del.push_back(an_entry);
						//unsupport_del[usdi][0] = data[0];
						//unsupport_del[usdi][1] = data[1];
						//unsupport_del[usdi][2] = data[2];
						//unsupport_del[usdi][3] = data[3];
						//unsupport_del[usdi][4] = cur_start;
						//unsupport_del[usdi][5] = cur_end;
						//unsupport_del[usdi][6] = data[6];
						//unsupport_del[usdi][7] = "data[7]\n";
					}
				}
			} else {
				if(debug) {
					cout << "del here-48\n";
				}
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.ref_id = data[3];
				an_entry.event_start = cur_start;
				an_entry.event_end = cur_end;
				an_entry.event_size_1 = boost::lexical_cast<int64_t>(data[6]);
				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, temp_cols);
				for (uint64_t c_id = 0; c_id < temp_cols.size(); ++c_id) {
					an_entry.cbps.push_back(boost::lexical_cast<int64_t>(temp_cols[c_id]));
				}
				unsupport_del.push_back(an_entry);
				//unsupport_del[usdi][0] = data[0];
				//unsupport_del[usdi][1] = data[1];
				//unsupport_del[usdi][2] = data[2];
				//unsupport_del[usdi][3] = data[3];
				//unsupport_del[usdi][4] = cur_start;
				//unsupport_del[usdi][5] = cur_end;
				//unsupport_del[usdi][6] = data[6];
				//unsupport_del[usdi][7] = "data[7]\n";
			}
		}
	}

}
void SplitReadSVCaller::detect_inssd(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "15646_0/16210_0" == data[1];
	const bool debug = false;
	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);

	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	//int32_t orientationa = a_region_a->second.orientation;

	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;

	//int32_t orientationb = a_region_b->second.orientation;
	if(debug) {
		cout << "inssd positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
		cout << "inssd positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
	}

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;

	if (covered(cur_end, cur_end, positiona1, positiona2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, positiona3, positiona4)) {
		bp_window[cur_end] = 2;
	}

	if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
		bp_window[cur_mate_end] = 1;
	}

	if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
		bp_window[cur_mate_end] = 2;
	}

	if (covered(cur_start, cur_start, positionb1, positionb2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, positionb3, positionb4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
		bp_window[cur_mate_start] = 1;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
		bp_window[cur_mate_start] = 2;
	}
	int64_t n_bp_zeros = 0;
	if (0 == bp_window[cur_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_end]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_end]) {
		++n_bp_zeros;
	}
	if (n_bp_zeros >= 2) {
		swap(cl[0], cl[1]);

		const auto& a_region_a = cluster_region.find(cl[0]);
		positiona1 = a_region_a->second.start;
		positiona2 = a_region_a->second.end;
		positiona3 = a_region_a->second.mate_start;
		positiona4 = a_region_a->second.mate_end;
		//orientationa = a_region_a->second.orientation;

		const auto& a_region_b = cluster_region.find(cl[1]);
		positionb1 = a_region_b->second.start;
		positionb2 = a_region_b->second.end;
		positionb3 = a_region_b->second.mate_start;
		positionb4 = a_region_b->second.mate_end;

		//orientationb = a_region_b->second.orientation;
		if(debug) {
			cout << "inssd positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
			cout << "inssd positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
		}

		const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
		const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

		local_ref_id = rev_index->second.first;
		local_ref_end = rev_index->second.second;
		local_ref_id_second = rev_index_second->second.first;
		local_ref_end_second = rev_index_second->second.second;

		if (covered(cur_end, cur_end, positiona1, positiona2)) {
			bp_window[cur_end] = 1;
		}
		if (covered(cur_end, cur_end, positiona3, positiona4)) {
			bp_window[cur_end] = 2;
		}

		if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
			bp_window[cur_mate_end] = 1;
		}

		if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
			bp_window[cur_mate_end] = 2;
		}

		if (covered(cur_start, cur_start, positionb1, positionb2)) {
			bp_window[cur_start] = 1;
		}
		if (covered(cur_start, cur_start, positionb3, positionb4)) {
			bp_window[cur_start] = 2;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
			bp_window[cur_mate_start] = 1;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
			bp_window[cur_mate_start] = 2;
		}
	}

	if(debug) {
		cout << "inssd bp_window cur start/end, cur mate_start/end: " << bp_window[cur_start] << "/" << bp_window[cur_end] << "/" << bp_window[cur_mate_start] << "/" << bp_window[cur_mate_end] << "\n";
	}
	//#print "a_region_a->second\na_region_b->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\ncur_mate_start\tbp_window{cur_mate_start}\ncur_mate_end\tbp_window{cur_mate_end}\n";
	if (a_region_a->second == a_region_b->second) {
		if(debug) {
			cout << "inssd here-0\n";
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
		//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
//	discord_sr1.clear();
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			//int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			//int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand) {

				//#print "start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "inssd here-1\n";
			}
			vector<int64_t> ref_boundary1;
			vector<int64_t> ref_boundary2;
			vector<int64_t> ref_support_sr;
			map<int64_t, vector<string>> ref_bpread;
			sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 1);

			//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
			if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
				if(debug) {
					cout << "inssd here-2\n";
				}
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;

				if (bp_window[cur_end] == 1 && ref_boundary1[2] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-0\n";
					}
					left_bound_sr1 = positiona1 + ref_boundary1[2];
				}

				if (bp_window[cur_end] == 2 && ref_boundary1[2] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-1\n";
					}
					left_bound_sr1 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (bp_window[cur_mate_end] == 1 && ref_boundary2[3] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-2\n";
					}
					right_bound_sr1 = positiona1 + ref_boundary2[3] + options.cut_sr;
				}

				if (ref_boundary2.size() > 3 && bp_window[cur_mate_end] == 2 && ref_boundary2[3] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-3\n";
					}
					right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (bp_window[cur_start] == 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-4\n";
					}
					left_bound_sr2 = positiona1 + ref_boundary1[1] + options.cut_sr;
				}

				if (bp_window[cur_start] == 2 && ref_boundary1[1] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-5\n";
					}
					left_bound_sr2 = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary2.size() > 0 && bp_window[cur_mate_start] == 1 && ref_boundary2[0] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-6\n";
					}
					right_bound_sr2 = positiona1 + ref_boundary2[0];
				}

				if (ref_boundary2.size() > 0 && bp_window[cur_mate_start] == 2 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd-d here-7\n";
					}
					right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
					if(debug) {
						cout << "inssd-d here-8\n";
					}
					left_bound_sr1 = cur_end;
				}
				if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
					if(debug) {
						cout << "inssd-d here-9\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
				if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
					if(debug) {
						cout << "inssd-d here-10\n";
					}
					left_bound_sr2 = cur_start;
				}
				if (!ref_boundary2.empty() && !ref_boundary2[0]) {
					if(debug) {
						cout << "inssd-d here-11\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
				int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
				int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
					if(debug) {
						cout << "inssd here-3\n";
					}
					string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
					string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[0];
					an_entry.n_mate_support = ref_support_sr[1];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr2 + 1;
					an_entry.event_end = left_bound_sr1 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = right_bound_sr2 + 1;
					an_entry.mate_event_end = right_bound_sr1 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = left_bound_sr2;
					//result_sr[sri][6] = left_bound_sr1;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					//
				} else {
					if(debug) {
						cout << "inssd-d here-12\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "inssd-d here-13\n";
						}
						int64_t left_bound_sr = left_bound_sr1;
						int64_t right_bound_sr = right_bound_sr1;
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (ref_bpread.end() != ref_bpread.find(1)) {
							if(debug) {
								cout << "inssd-d here-14\n";
							}
							reads = castle::StringUtils::join(ref_bpread[1], "\t");
						} else {
							if(debug) {
								cout << "inssd-d here-15\n";
							}
							left_bound_sr = positiona1 + ref_boundary1[0];
							right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "inssd here-1\n";
							}
							int64_t size = right_bound_sr - left_bound_sr;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "tandem_dup";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";

						}
					}
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "inssd here-4\n";
						}
						int64_t left_bound_sr = left_bound_sr2;
						int64_t right_bound_sr = right_bound_sr2;
						int64_t size = right_bound_sr - left_bound_sr - 1;
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "del";
						an_entry.cluster_id = cluster_id2;
						an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
						an_entry.n_supports = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "del\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";

					}
				}
			} else {
				if(debug) {
					cout << "inssd here-5\n";
				}
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;

				if (ref_boundary1.size() > 0 && bp_window[cur_end] == 1 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-a\n";
					}
					left_bound_sr1 = positiona1 + ref_boundary1[0];
				}

				if (ref_boundary1.size() > 0 && bp_window[cur_end] == 2 && ref_boundary1[0] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-b\n";
					}
					left_bound_sr1 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary2.size() > 1 && bp_window[cur_mate_end] == 1 && ref_boundary2[1] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-c\n";
					}
					right_bound_sr1 = positiona1 + ref_boundary2[1] + options.cut_sr;
				}

				if (ref_boundary2.size() > 1 && bp_window[cur_mate_end] == 2 && ref_boundary2[1] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-d\n";
					}
					right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary1.size() > 3 && bp_window[cur_start] == 1 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-e\n";
					}
					left_bound_sr2 = positiona1 + ref_boundary1[3] + options.cut_sr;
				}

				if (ref_boundary1.size() > 3 && bp_window[cur_start] == 2 && ref_boundary1[3] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-f\n";
					}
					left_bound_sr2 = ref_boundary1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary2.size() > 2 && bp_window[cur_mate_start] == 1 && ref_boundary2[2] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-g\n";
					}
					right_bound_sr2 = positiona1 + ref_boundary2[2];
				}

				if (ref_boundary2.size() > 2 && bp_window[cur_mate_start] == 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssd here-5-h\n";
					}
					right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
					if(debug) {
						cout << "inssd here-5-i\n";
					}
					left_bound_sr1 = cur_end;
				}
				if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
					if(debug) {
						cout << "inssd here-5-j\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
				if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
					if(debug) {
						cout << "inssd here-5-k\n";
					}
					left_bound_sr2 = cur_start;
				}
				if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
					if(debug) {
						cout << "inssd here-5-l\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
				int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
				int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
					if(debug) {
						cout << "inssd here-6\n";
					}
					string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
					string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[1];
					an_entry.n_mate_support = ref_support_sr[0];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr2 + 1;
					an_entry.event_end = left_bound_sr1 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = right_bound_sr2 + 1;
					an_entry.mate_event_end = right_bound_sr1 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = left_bound_sr2;
					//result_sr[sri][6] = left_bound_sr1;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					//
				} else {
					if(debug) {
						cout << "inssd here-7\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						if(debug) {
							cout << "inssd here-8\n";
						}
						int64_t left_bound_sr = left_bound_sr1;
						int64_t right_bound_sr = right_bound_sr1;
						int64_t size = right_bound_sr - left_bound_sr;
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "tandem_dup";
						an_entry.cluster_id = cluster_id1;
						an_entry.n_supports = mpd1;
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
					}
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "inssd here-9\n";
						}
						int64_t left_bound_sr = left_bound_sr2;
						int64_t right_bound_sr = right_bound_sr2;
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (ref_bpread.end() != ref_bpread.find(1)) {
							if(debug) {
								cout << "inssd-d here-16\n";
							}
							reads = castle::StringUtils::join(ref_bpread[1], "\t");
						} else {
							if(debug) {
								cout << "inssd-d here-17\n";
							}
							int64_t tmp_left_bound_sr = 0;
							if(ref_boundary1.size() > 1) {
								tmp_left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t tmp_right_bound_sr = 0;
							if(ref_boundary2.size() > 0) {
								tmp_right_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if (tmp_left_bound_sr > 0) {
								left_bound_sr = tmp_left_bound_sr;
							}
							if (tmp_right_bound_sr > 0) {
								right_bound_sr = tmp_right_bound_sr;
							}
						}
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "inssd here-10\n";
							}
							int64_t size = right_bound_sr - left_bound_sr - 1;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "del";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "del\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				}
			}
		}
	}
	// a_region_a->second != a_region_b->second
	else {
		if(debug) {
			cout << "inssd here-11\n";
		}
		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (bp_window[cur_end] == 1 && bp_window[cur_mate_end] == 2) {
			if(debug) {
				cout << "inssd-d here-18\n";
			}
			discord_sr1.clear();
//	reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			const int64_t first_end = positiona2 - positiona1;
			const int64_t min_mate_start = first_end + 100;

			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
//			if (debug) {
//				cout << "debug position: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "/" << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
//				cout << "debug position: a2-a1: " << (positiona2 - positiona1) << "/a2-a1 + 100:" << (positiona2 - positiona1 + 100) << "\n";
//			}
//			int64_t debug_i = 0;
			while (reader.GetNextAlignmentBasic(al)) {
//				++debug_i;
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;
//				if (debug) {
//					cout << "debug aln-" << debug_i << ":" << strand << "/" << mstrand << "/" << start << "/" << mstart << "/" << isize << "\n";
//				}
				if (isize > 100 && strand == mstrand && start < first_end && mstart > min_mate_start) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if(debug) {
				cout << "inssd here-12:" << local_ref_id << ":" << local_ref_end << "/" << discord_sr1.size() << "\n";
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssd here-13\n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
//				for(auto& e : src_return1) {
//					cout << "inssd here-13 src_return1: " << e << "\n";
//				}
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "inssd-d here-17\n";
					}
					left_bound_sr1 = positiona1 + src_return1[0];
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "inssd-d here-18\n";
					}
					right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssd-d here-19\n";
					}
					left_bound_sr1 = cur_end;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "inssd-d here-20\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
				if(debug) {
					cout << "inssd here-13 left_bound_sr1: " << left_bound_sr1 << "/right_bound_sr1: " << right_bound_sr1 << "\n";
				}
			}
		} else if (bp_window[cur_end] == 1 && bp_window[cur_mate_end] == 1) {
			if(debug) {
				cout << "inssd here-14\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssd here-15\n";
				}
//	vector<int64_t> src_return1;
//	map<int64_t, vector<string>> bpread;
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				if(src_return1.size() > 0) {
					left_bound_sr1 = positiona1 + src_return1[0];
				}
				if(src_return1.size() > 3) {
					right_bound_sr1 = positiona1 + src_return1[3] + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssd-d here-21\n";
					}
					left_bound_sr1 = cur_end;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "inssd-d here-22\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
			}
		} else if (bp_window[cur_end] == 2 && bp_window[cur_mate_end] == 2) {
			if(debug) {
				cout << "inssd here-16\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssd here-17\n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "inssd-d here-18\n";
					}
					left_bound_sr1 = src_return1[0] - (positiona2 - positiona1 + 101) + positiona3;
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "inssd-d here-19\n";
					}
					right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssd-d here-20\n";
					}
					left_bound_sr1 = cur_end;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "inssd-d here-21\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
			}
		}

		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_b->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		if (bp_window[cur_start] == 1 && bp_window[cur_mate_start] == 2) {
			if(debug) {
				cout << "inssd here-18\n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				if(debug) {
					cout << "inssd here-19:\n";
				}
//				for(auto& e : src_return2) {
//					cout << "inssd here-19 src_return2: " << e << "\n";
//				}
				if (src_return2.size() > 1) {
					left_bound_sr2 = positionb1 + src_return2[1] + options.cut_sr;
				}
				if (src_return2.size() > 2) {
					right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssd-d here-22\n";
					}
					left_bound_sr2 = cur_start;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssd-d here-23\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
				if(debug) {
					cout << "inssd here-19 left_bound_sr2: " << left_bound_sr2 << "/right_bound_sr1: " << right_bound_sr2 << "\n";
				}
			}
		} else if (bp_window[cur_start] == 1 && bp_window[cur_mate_start] == 1) {
			if(debug) {
				cout << "inssd here-20\n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (positionb2 - positionb1) && mstart < (positionb2 - positionb1)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssd here-21\n";
				}
//	vector<int64_t> src_return2;
//	map<int64_t, vector<string>> bpread;
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				if(src_return2.size() > 1) {
					left_bound_sr2 = positionb1 + src_return2[1] + options.cut_sr;
				}
				if (src_return2.size() > 2) {
					if(debug) {
						cout << "inssd-d here-24\n";
					}
					right_bound_sr2 = positionb1 + src_return2[2];
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssd-d here-25\n";
					}
					left_bound_sr2 = cur_start;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssd-d here-26\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
			}
		} else if (bp_window[cur_start] == 2 && bp_window[cur_mate_start] == 2) {
			if(debug) {
				cout << "inssd here-22\n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start > (positionb2 - positionb1 + 100) && mstart > (positionb2 - positionb1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssd here-23\n";
				}
//	vector<int64_t> src_return2;
//	map<int64_t, vector<string>> bpread;
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				if(src_return2.size() > 1) {
					left_bound_sr2 = src_return2[1] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
				}
				if(src_return2.size() > 2) {
					right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssd-d here-27\n";
					}
					left_bound_sr2 = cur_start;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssd-d here-28\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
			}
		}

		int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
		int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
			if(debug) {
				cout << "inssd here-24\n";
			}
			string reads0 = castle::StringUtils::join(bpread_2[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_1[0], "\t");
			BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = data[3];
			an_entry.event_start = left_bound_sr2 + 1;
			an_entry.event_end = left_bound_sr1 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = data[3];
			an_entry.mate_event_start = right_bound_sr2 + 1;
			an_entry.mate_event_end = right_bound_sr1 + 1;
			an_entry.event_size_2 = ins_size;
			result_sr.push_back(an_entry);
			//result_sr[sri][0] = data[0];
			//result_sr[sri][1] = data[1];
			//result_sr[sri][2] = data[2];
			//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
			//result_sr[sri][4] = data[3];
			//result_sr[sri][5] = left_bound_sr2;
			//result_sr[sri][6] = left_bound_sr1;
			//result_sr[sri][7] = del_size;
			//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
		} else {
			if(debug) {
				cout << "inssd here-25\n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
				if(debug) {
					cout << "inssd here-26\n";
				}
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "tandem_dup";
				an_entry.cluster_id = cluster_id1;
				an_entry.n_supports = mpd1;
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
				if(debug) {
					cout << "inssd here-27\n";
				}
				int64_t left_bound_sr = left_bound_sr2;
				int64_t right_bound_sr = right_bound_sr2;
				int64_t size = right_bound_sr - left_bound_sr - 1;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "del";
				an_entry.cluster_id = cluster_id2;
				an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
				an_entry.n_supports = src_return2[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "del\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
			}
		}
	}

}

void SplitReadSVCaller::detect_inssu(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "16110_0/6390_0" == data[1];
	const bool debug = false;
	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	//int32_t orientationa = a_region_a->second.orientation;

	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;

	if(debug) {
		cout << "inssu positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
		cout << "inssu positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
	}

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;

	//int32_t orientationb = a_region_b->second.orientation;
	if (covered(cur_mate_start, cur_mate_start, positiona1, positiona2)) {
		bp_window[cur_mate_start] = 1;
	}
	if (covered(cur_mate_start, cur_mate_start, positiona3, positiona4)) {
		bp_window[cur_mate_start] = 2;
	}
	if (covered(cur_start, cur_start, positiona3, positiona4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_start, cur_start, positiona1, positiona2)) {
		bp_window[cur_start] = 1;
	}

	if (covered(cur_mate_end, cur_mate_end, positionb1, positionb2)) {
		bp_window[cur_mate_end] = 1;
	}
	if (covered(cur_mate_end, cur_mate_end, positionb3, positionb4)) {
		bp_window[cur_mate_end] = 2;
	}
	if (covered(cur_end, cur_end, positionb1, positionb2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, positionb3, positionb4)) {
		bp_window[cur_end] = 2;
	}

	int64_t n_bp_zeros = 0;
	if (0 == bp_window[cur_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_end]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_end]) {
		++n_bp_zeros;
	}
	if (n_bp_zeros >= 2) {
		swap(cl[0], cl[1]);

		const auto& a_region_a = cluster_region.find(cl[0]);
		positiona1 = a_region_a->second.start;
		positiona2 = a_region_a->second.end;
		positiona3 = a_region_a->second.mate_start;
		positiona4 = a_region_a->second.mate_end;
		//orientationa = a_region_a->second.orientation;

		const auto& a_region_b = cluster_region.find(cl[1]);
		positionb1 = a_region_b->second.start;
		positionb2 = a_region_b->second.end;
		positionb3 = a_region_b->second.mate_start;
		positionb4 = a_region_b->second.mate_end;

		//orientationb = a_region_b->second.orientation;
		if(debug) {
			cout << "inssu positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
			cout << "inssu positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
		}

		const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
		const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

		local_ref_id = rev_index->second.first;
		local_ref_end = rev_index->second.second;
		local_ref_id_second = rev_index_second->second.first;
		local_ref_end_second = rev_index_second->second.second;

		if (covered(cur_mate_start, cur_mate_start, positiona1, positiona2)) {
			bp_window[cur_mate_start] = 1;
		}
		if (covered(cur_mate_start, cur_mate_start, positiona3, positiona4)) {
			bp_window[cur_mate_start] = 2;
		}
		if (covered(cur_start, cur_start, positiona3, positiona4)) {
			bp_window[cur_start] = 2;
		}
		if (covered(cur_start, cur_start, positiona1, positiona2)) {
			bp_window[cur_start] = 1;
		}

		if (covered(cur_mate_end, cur_mate_end, positionb1, positionb2)) {
			bp_window[cur_mate_end] = 1;
		}
		if (covered(cur_mate_end, cur_mate_end, positionb3, positionb4)) {
			bp_window[cur_mate_end] = 2;
		}
		if (covered(cur_end, cur_end, positionb1, positionb2)) {
			bp_window[cur_end] = 1;
		}
		if (covered(cur_end, cur_end, positionb3, positionb4)) {
			bp_window[cur_end] = 2;
		}
	}

	//#print "a_region_a->second\na_region_b->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\ncur_mate_start\tbp_window{cur_mate_start}\ncur_mate_end\tbp_window{cur_mate_end}\n";
	if (a_region_a->second == a_region_b->second) {
		const int64_t min_isize = (ref_is["rlu"]["selected"] - options.cut_sr + 5);
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if(debug) {
			cout << "inssu here-0: \n";
		}
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			//int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			//int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > min_isize && strand == mstrand) {
				//#print "readname\tstart mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}

		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "inssu here-1: \n";
			}
			vector<int64_t> ref_boundary1;
			vector<int64_t> ref_boundary2;
			vector<int64_t> ref_support_sr;
//	map<int64_t, vector<string>> ref_bpread;
			sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, bpread_1, discord_sr1, 1, 2);
//			cout << "inssu: ref_boundary1\n";
//			for(auto& e: ref_boundary1) {
//				cout << e << "\n";
//			}
//			cout << "inssu: ref_boundary2\n";
//			for(auto& e: ref_boundary2) {
//				cout << e << "\n";
//			}
//			cout << "inssu: ref_support_sr\n";
//			for(auto& e: ref_support_sr) {
//				cout << e << "\n";
//			}
//			cout << "inssu: bpread_1\n";
//			for(auto& e: bpread_1) {
//				cout << e.first << "\n";
//				for(auto& s: e.second) {
//					cout << s << "\n";
//				}
//			}
			//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n@{ref_bpread{0}}\n@{ref_bpread{1}}\n";
			if (ref_boundary2.size() > 2 && ref_boundary2[2] > ref_boundary2[0]) {
				if(debug) {
					cout << "inssu here-2: \n";
				}
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;
				//cout << "here-5: \n";
				if (ref_boundary1.size() > 0 && bp_window[cur_mate_start] == 1 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-a: \n";
					}
					left_bound_sr1 = positiona1 + ref_boundary1[0];
				}

				if (ref_boundary1.size() > 0 && bp_window[cur_mate_start] == 2 && ref_boundary1[0] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-b: \n";
					}
					left_bound_sr1 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary2.size() > 1 && bp_window[cur_start] == 1 && ref_boundary2[1] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-c: \n";
					}
					right_bound_sr1 = positiona1 + ref_boundary2[1] + options.cut_sr;
				}

				if (ref_boundary2.size() > 1 && bp_window[cur_start] == 2 && ref_boundary2[1] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-d: \n";
					}
					right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary1.size() > 3 && bp_window[cur_mate_end] == 1 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-e: \n";
					}
					left_bound_sr2 = positiona1 + ref_boundary1[3] + options.cut_sr;
				}

				if (ref_boundary1.size() > 3 && bp_window[cur_mate_end] == 2 && ref_boundary1[3] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-f: \n";
					}
					left_bound_sr2 = ref_boundary1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary2.size() > 2 && bp_window[cur_end] == 1 && ref_boundary2[2] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-g: \n";
					}
					right_bound_sr2 = positiona1 + ref_boundary2[2];
				}

				if (ref_boundary2.size() > 2 && bp_window[cur_end] == 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-2-h: \n";
					}
					right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
					if(debug) {
						cout << "inssu here-2-i: \n";
					}
					left_bound_sr1 = cur_mate_start;
				}
				if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
					if(debug) {
						cout << "inssu here-2-j: \n";
					}
					right_bound_sr1 = cur_start;
				}
				if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
					if(debug) {
						cout << "inssu here-2-k: \n";
					}
					left_bound_sr2 = cur_mate_end;
				}
				if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
					if(debug) {
						cout << "inssu here-2-l: \n";
					}
					right_bound_sr2 = cur_end;
				}
				int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
				int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
					string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
					string reads1 = castle::StringUtils::join(bpread_1[1], "\t");
					string a_line = (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
					if(debug) {
						cout << "inssu here-3: " << a_line << "\n";
					}

					BPREAD << a_line;
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[0];
					an_entry.n_mate_support = ref_support_sr[1];
					an_entry.ref_id = data[3];
					an_entry.event_start = right_bound_sr1 + 1;
					an_entry.event_end = right_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = left_bound_sr1 + 1;
					an_entry.mate_event_end = left_bound_sr2 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = right_bound_sr1;
					//result_sr[sri][6] = right_bound_sr2;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
				} else {
					if(debug) {
						cout << "inssu here-4: \n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						if(debug) {
							cout << "inssu here-5: \n";
						}
						int64_t left_bound_sr = left_bound_sr1;
						int64_t right_bound_sr = right_bound_sr1;
						int64_t size = right_bound_sr - left_bound_sr;
						string reads = castle::StringUtils::join(bpread_1[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr - 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "tandem_dup";
						an_entry.cluster_id = cluster_id1;
						an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr - 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
					}
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "inssu here-6: " << left_bound_sr1 << ":" << left_bound_sr2 << "/" << right_bound_sr1 << ":" << right_bound_sr2 << "\n";
						}
						int64_t left_bound_sr = left_bound_sr2;
						int64_t right_bound_sr = right_bound_sr2;
						int64_t size = right_bound_sr - left_bound_sr - 1;
						string reads = castle::StringUtils::join(bpread_1[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr - 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "del";
						an_entry.cluster_id = cluster_id2;
						an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
						an_entry.n_supports = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr - 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//"del\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
					}
				}
			} else {
				if(debug) {
					cout << "inssu here-7: \n";
				}
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;

				if(ref_boundary1.size() > 2 && bp_window[cur_mate_start] == 1 && ref_boundary1[2] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-a: \n";
					}
					left_bound_sr1 = positiona1 + ref_boundary1[2];
				}

				if (ref_boundary1.size() > 2 && bp_window[cur_mate_start] == 2 && ref_boundary1[2] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-b: \n";
					}
					left_bound_sr1 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary2.size() > 3 && bp_window[cur_start] == 1 && ref_boundary2[3] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-c: \n";
					}
					right_bound_sr1 = positiona1 + ref_boundary2[3] + options.cut_sr;
				}

				if (ref_boundary2.size() > 3 && bp_window[cur_start] == 2 && ref_boundary2[3] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-d: \n";
					}
					right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary1.size() > 1 && bp_window[cur_mate_end] == 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-e: \n";
					}
					left_bound_sr2 = positiona1 + ref_boundary1[1] + options.cut_sr;
				}

				if (ref_boundary1.size() > 1 && bp_window[cur_mate_end] == 2 && ref_boundary1[1] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-f: \n";
					}
					left_bound_sr2 = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}

				if (ref_boundary2.size() > 0 && bp_window[cur_end] == 1 && ref_boundary2[0] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-g: \n";
					}
					right_bound_sr2 = positiona1 + ref_boundary2[0];
				}

				if (ref_boundary2.size() > 0 && bp_window[cur_end] == 2 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "inssu here-7-h: \n";
					}
					right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
				}

				if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
					if(debug) {
						cout << "inssu here-7-i: \n";
					}
					left_bound_sr1 = cur_mate_start;
				}

				if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
					if(debug) {
						cout << "inssu here-7-j: \n";
					}
					right_bound_sr1 = cur_start;
				}

				if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
					if(debug) {
						cout << "inssu here-7-k: \n";
					}
					left_bound_sr2 = cur_mate_end;
				}

				if (!ref_boundary2.empty() && !ref_boundary2[0]) {
					if(debug) {
						cout << "inssu here-7-l: \n";
					}
					right_bound_sr2 = cur_end;
				}

				int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
				int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;

				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
					if(debug) {
						cout << "inssu here-8: \n";
					}
					string reads0 = castle::StringUtils::join(bpread_1[1], "\t");
					string reads1 = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[1];
					an_entry.n_mate_support = ref_support_sr[0];
					an_entry.ref_id = data[3];
					an_entry.event_start = right_bound_sr1 + 1;
					an_entry.event_end = right_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = left_bound_sr1 + 1;
					an_entry.mate_event_end = left_bound_sr2 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = right_bound_sr1;
					//result_sr[sri][6] = right_bound_sr2;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";
				} else {
					if(debug) {
						cout << "inssu here-9: \n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						if(debug) {
							cout << "inssu here-10: \n";
						}
						int64_t left_bound_sr = left_bound_sr1;
						int64_t right_bound_sr = right_bound_sr1;
						int64_t size = right_bound_sr - left_bound_sr;
						string reads = castle::StringUtils::join(bpread_1[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "tandem_dup";
						an_entry.cluster_id = cluster_id1;
						an_entry.n_supports = mpd1;
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
					}
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "inssu here-11: \n";
						}
						int64_t left_bound_sr = left_bound_sr2;
						int64_t right_bound_sr = right_bound_sr2;
						int64_t size = right_bound_sr - left_bound_sr - 1;
						string reads = castle::StringUtils::join(bpread_1[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "del";
						an_entry.cluster_id = cluster_id2;
						an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "del\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
					}
				}
			}
		}
	}
	// a_region_a->second != a_region_b->second
	else {
		if(debug) {
			cout << "inssu here-12: \n";
		}

		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (bp_window[cur_mate_start] == 1 && bp_window[cur_start] == 2) {
			if(debug) {
				cout << "inssu here-13: \n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-14: \n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "inssu here-14-a: \n";
					}
					left_bound_sr1 = positiona1 + src_return1[0];
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "inssu here-14-b: \n";
					}
					right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssu here-14-c: \n";
					}
					left_bound_sr1 = cur_mate_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "inssu here-14-d: \n";
					}
					right_bound_sr1 = cur_start;
				}
			}
		} else if (bp_window[cur_mate_start] == 1 && bp_window[cur_start] == 1) {
			if(debug) {
				cout << "inssu here-15: \n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-16: \n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "inssu here-16-a: \n";
					}
					left_bound_sr1 = positiona1 + src_return1[0];
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "inssu here-16-b: \n";
					}
					right_bound_sr1 = positiona1 + src_return1[3] + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssu here-16-c: \n";
					}
					left_bound_sr1 = cur_mate_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					//cout << "here-16-d: \n";
					right_bound_sr1 = cur_start;
				}
			}
		} else if (bp_window[cur_mate_start] == 2 && bp_window[cur_start] == 2) {
			if(debug) {
				cout << "inssu here-17: \n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-18: \n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "inssu here-18-a: \n";
					}
					left_bound_sr1 = src_return1[0] - (positiona2 - positiona1 + 101) + positiona3;
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "inssu here-18-b: \n";
					}
					right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "inssu here-18-c: \n";
					}
					left_bound_sr1 = cur_mate_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "inssu here-18-d: \n";
					}
					right_bound_sr1 = cur_start;
				}
			}
		}

		// the cl[1] part
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		if (bp_window[cur_mate_end] == 1 && bp_window[cur_end] == 2) {
			if(debug) {
				cout << "inssu here-19: \n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {

				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;
				if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {
					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-20: \n";
				}
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);

				if (src_return2.size() > 1) {
					if(debug) {
						cout << "inssu here-20-a: \n";
					}
					left_bound_sr2 = positionb1 + src_return2[1] + options.cut_sr;
				}
				if (src_return2.size() > 2) {
					if(debug) {
						cout << "inssu here-20-b: \n";
					}
					right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssu here-20-c: \n";
					}
					left_bound_sr2 = cur_mate_end;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssu here-20-d: \n";
					}
					right_bound_sr2 = cur_end;
				}
			}
		} else if (bp_window[cur_mate_end] == 1 && bp_window[cur_end] == 1) {
			if(debug) {
				cout << "inssu here-21: \n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (positionb2 - positionb1) && mstart < (positionb2 - positionb1)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-22: \n";
				}
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				left_bound_sr2 = positionb1 + src_return2[1] + options.cut_sr;
				if (src_return2.size() > 2) {
					if(debug) {
						cout << "inssu here-22-a: \n";
					}
					right_bound_sr2 = positionb1 + src_return2[2];
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssu here-22-b: \n";
					}
					left_bound_sr2 = cur_mate_end;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssu here-22-c: \n";
					}
					right_bound_sr2 = cur_end;
				}
			}
		} else if (bp_window[cur_mate_end] == 2 && bp_window[cur_end] == 2) {
			if(debug) {
				cout << "inssu here-23: \n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start > (positionb2 - positionb1 + 100) && mstart > (positionb2 - positionb1 + 100)) {
					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inssu here-24: \n";
				}
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				if (src_return2.size() > 1) {
					if(debug) {
						cout << "inssu here-24-a: \n";
					}
					left_bound_sr2 = src_return2[1] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
				}
				if (src_return2.size() > 2) {
					if(debug) {
						cout << "inssu here-24-b: \n";
					}
					right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
						cout << "inssu here-24-c: \n";
					}
					left_bound_sr2 = cur_mate_end;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
						cout << "inssu here-24-d: \n";
					}
					right_bound_sr2 = cur_end;
				}
			}
		}
		//cout << "here-53-1: " << src_return1.size() << "/" << src_return2.size() << "\n";
		int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
		int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
			if(debug) {
				cout << "inssu here-25: \n";
			}
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = data[3];
			an_entry.event_start = right_bound_sr1 + 1;
			an_entry.event_end = right_bound_sr2 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = data[3];
			an_entry.mate_event_start = left_bound_sr1 + 1;
			an_entry.mate_event_end = left_bound_sr2 + 1;
			an_entry.event_size_2 = ins_size;
			result_sr.push_back(an_entry);
			//result_sr[sri][0] = data[0];
			//result_sr[sri][1] = data[1];
			//result_sr[sri][2] = data[2];
			//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
			//result_sr[sri][4] = data[3];
			//result_sr[sri][5] = right_bound_sr1;
			//result_sr[sri][6] = right_bound_sr2;
			//result_sr[sri][7] = del_size;
			//result_sr[sri][8] = "data[3]\tleft_bound_sr1\tleft_bound_sr2\tins_size\n";

		} else {
			if(debug) {
				cout << "inssu here-26: \n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
				if(debug) {
					cout << "inssu here-27: \n";
				}
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr - 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "tandem_dup";
				an_entry.cluster_id = cluster_id1;
				an_entry.n_supports = mpd1;
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr - 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "tandem_dup\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
				//
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
				if(debug) {
					cout << "inssu here-28: \n";
				}
				int64_t left_bound_sr = left_bound_sr2;
				int64_t right_bound_sr = right_bound_sr2;
				int64_t size = right_bound_sr - left_bound_sr - 1;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr - 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "del";
				an_entry.cluster_id = cluster_id2;
				an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
				an_entry.n_supports = src_return2[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr - 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "del\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
			}
		}
	}
}

void SplitReadSVCaller::detect_insod(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "132539_0/198512_0" == data[1];
	const bool debug = false;
//	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	int32_t orientationa = a_region_a->second.orientation;

	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;
	int32_t orientationb = a_region_b->second.orientation;
	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;
	if(debug) {
		cout << "insod positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
		cout << "insod positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
	}

	if (covered(cur_start, cur_start, positiona1, positiona2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, positiona3, positiona4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
		bp_window[cur_mate_end] = 1;
	}
	if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
		bp_window[cur_mate_end] = 2;
	}
	if (covered(cur_end, cur_end, positionb1, positionb2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, positionb3, positionb4)) {
		bp_window[cur_end] = 2;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
		bp_window[cur_mate_start] = 1;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
		bp_window[cur_mate_start] = 2;
	}

	int64_t n_bp_zeros = 0;
	if (0 == bp_window[cur_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_end]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_end]) {
		++n_bp_zeros;
	}
	if (n_bp_zeros >= 2) {
		swap(cl[0], cl[1]);

		const auto& a_region_a = cluster_region.find(cl[0]);
		positiona1 = a_region_a->second.start;
		positiona2 = a_region_a->second.end;
		positiona3 = a_region_a->second.mate_start;
		positiona4 = a_region_a->second.mate_end;
		//orientationa = a_region_a->second.orientation;

		const auto& a_region_b = cluster_region.find(cl[1]);
		positionb1 = a_region_b->second.start;
		positionb2 = a_region_b->second.end;
		positionb3 = a_region_b->second.mate_start;
		positionb4 = a_region_b->second.mate_end;

		//orientationb = a_region_b->second.orientation;
		if(debug) {
			cout << "insod positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
			cout << "insod positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
		}

		const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
		const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

		local_ref_id = rev_index->second.first;
		local_ref_end = rev_index->second.second;
		local_ref_id_second = rev_index_second->second.first;
		local_ref_end_second = rev_index_second->second.second;

		if (covered(cur_start, cur_start, positiona1, positiona2)) {
			bp_window[cur_start] = 1;
		}
		if (covered(cur_start, cur_start, positiona3, positiona4)) {
			bp_window[cur_start] = 2;
		}
		if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
			bp_window[cur_mate_end] = 1;
		}
		if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
			bp_window[cur_mate_end] = 2;
		}
		if (covered(cur_end, cur_end, positionb1, positionb2)) {
			bp_window[cur_end] = 1;
		}
		if (covered(cur_end, cur_end, positionb3, positionb4)) {
			bp_window[cur_end] = 2;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
			bp_window[cur_mate_start] = 1;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
			bp_window[cur_mate_start] = 2;
		}
	}

	//#print "a_region_a->second\na_region_b->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\ncur_mate_start\tbp_window{cur_mate_start}\ncur_mate_end\tbp_window{cur_mate_end}\n";
	if (a_region_a->second == a_region_b->second) {
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if(debug) {
			cout << "insod here-0\n";
		}
		if (bp_window[cur_start] != bp_window[cur_mate_end] || bp_window[cur_end] != bp_window[cur_mate_start]) {
			if(debug) {
				cout << "insod here-1\n";
			}
			if (orientationa == 1) {
				if(debug) {
					cout << "insod here-2:" << local_ref_start << "/" << local_ref_end << "\n";
				}
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "1 start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
			}
//	discord_sr2.clear();
			if(debug) {
				cout << "insod here-3:" << local_ref_start << "/" << local_ref_end << "\n";
			}
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				//int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				//int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand) {

					//#print "2 start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}

			//# small deletion, large insertion, close, 3 break points in one window
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads) && discord_sr2.size() >= static_cast<uint64_t>(options.support_reads) && bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1) {
				if(debug) {
					cout << "insod here-4:\n";
				}
//	vector<int64_t> src_return1;
//	vector<int64_t> src_return2;
//	map<int64_t, vector<string>> bpread_1;
//	map<int64_t, vector<string>> bpread_2;
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

				//#print "src_return1\nsrc_return2\n";
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;

				if (src_return1.size() > 1 && src_return1[1] < (positiona2 - positiona1 + 101)) {
					left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
				}

				if (src_return1.size() > 2 && src_return1[2] >= (positiona2 - positiona1 + 101)) {
					right_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
				}

				if (src_return2.size() > 0 && src_return2[0] < (positiona2 - positiona1 + 101)) {
					left_bound_sr2 = positiona1 + src_return2[0];
				}

				if (src_return2.size() > 2 && src_return2[2] < (positiona2 - positiona1 + 101)) {
					right_bound_sr2 = positiona1 + src_return2[2];
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					left_bound_sr1 = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					right_bound_sr1 = cur_mate_end;
				}
				if (!src_return2.empty() && !src_return2[0]) {
					left_bound_sr2 = cur_end;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					right_bound_sr2 = cur_mate_start;
				}
				int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
				int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

				if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1
						&& right_bound_sr2) {
					if(debug) {
						cout << "insod here-5:\n";
					}
					string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
					string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.n_mate_support = src_return2[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr1 + 1;
					an_entry.event_end = left_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = right_bound_sr2 + 1;
					an_entry.mate_event_end = right_bound_sr1 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = left_bound_sr1;
					//result_sr[sri][6] = left_bound_sr2;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
				} else {
					if(debug) {
						cout << "insod here-6:\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					//string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					//string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						if(debug) {
							cout << "insod here-7:\n";
						}
						int64_t left_bound_sr = left_bound_sr1;
						int64_t right_bound_sr = right_bound_sr1;
						int64_t size = right_bound_sr - left_bound_sr;
						string reads = castle::StringUtils::join(bpread_1[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "invers_f";
						an_entry.cluster_id = cluster_id1;
						an_entry.n_supports = mpd1;
						an_entry.n_mate_support = src_return1[4];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
//	if (src_return2.size() > 4 && src_return2.size() > 4 && src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
					//#int64_t left_bound_sr = left_bound_sr2;//#int64_t right_bound_sr = right_bound_sr2;//#int64_t size = right_bound_sr - left_bound_sr;//#my reads = join("\t", @{bpread_2{0}});//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#
//	}
				}
			}
			//# small deletion, small insertion, far, 2 break points in one window, the other 2 in the other window
			else if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "insod here-8:\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
					if(debug) {
						cout << "insod here-9:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}

					if (ref_boundary2.size() > 0 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr1 = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
					}

					if (ref_boundary1.size() > 2 && ref_boundary1[2] < (positiona2 - positiona1 + 101)) {
						left_bound_sr2 = positiona1 + ref_boundary1[2];
					}

					if (ref_boundary2.size() > 3 && ref_boundary2[3] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr2 = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && left_bound_sr2 && right_bound_sr1 && right_bound_sr2) {
						if(debug) {
							cout << "insod here-10:\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insod here-11:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
							if(debug) {
								cout << "insod here-12:\n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
							if(debug) {
								cout << "insod here-13:\n";
							}
							int64_t left_bound_sr = left_bound_sr2;
							int64_t right_bound_sr = right_bound_sr2;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_r";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				} else {
					if(debug) {
						cout << "insod here-14:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 3 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}

					if (ref_boundary2.size() > 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr1 = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
					}

					if (ref_boundary1.size() > 0 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
						left_bound_sr2 = positiona1 + ref_boundary1[0];
					}

					if (ref_boundary2.size() > 1 && ref_boundary2[1] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr2 = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insod here-15:\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insod here-16:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
							if(debug) {
								cout << "insod here-17:\n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
							if(debug) {
								cout << "insod here-18:\n";
							}
							int64_t left_bound_sr = left_bound_sr2 + 1;
							int64_t right_bound_sr = right_bound_sr2 - 1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % left_bound_sr % data[3] % right_bound_sr % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_r";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr;
							an_entry.event_end = right_bound_sr;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				}
			}
			//# large deletion, small insertion, close, bp_info 0, || small deletion, small insertion, far, bp_info 0
			else if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "insod here-19:\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr2, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
					if(debug) {
						cout << "insod here-20:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101))
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;

					if (ref_boundary2.size() > 1 && ref_boundary2[1] >= (positiona2 - positiona1 + 101))
						right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;

					if (ref_boundary1.size() > 2 && ref_boundary1[2] >= (positiona2 - positiona1 + 101))
						left_bound_sr2 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;

					if (ref_boundary2.size() > 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101))
						right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insod here-21:\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insod here-22:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
							if(debug) {
								cout << "insod here-22-a: \n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
						//if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						//#int64_t left_bound_sr = left_bound_sr2;//#int64_t right_bound_sr = right_bound_sr2;//#int64_t size = right_bound_sr - left_bound_sr;//#my reads = join("\t", @{ref_bpread{1}});//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#
						//}
					}
				} else {
					if(debug) {
						cout << "insod here-23:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 3 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}

					if (ref_boundary2.size() > 3 && ref_boundary2[3] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}

					if (ref_boundary1.size() > 0 && ref_boundary1[0] >= (positiona2 - positiona1 + 101)) {
						left_bound_sr2 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
					}

					if (ref_boundary2.size() > 0 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
						right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						left_bound_sr2 = cur_end;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insod here-24:\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insod here-25:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
							if(debug) {
								cout << "insod here-25-a:\n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
						//if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						//#int64_t left_bound_sr = left_bound_sr2;//#int64_t right_bound_sr = right_bound_sr2;//#int64_t size = right_bound_sr - left_bound_sr;//#my reads = join("\t", @{ref_bpread{0}});//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#
						//}
					}
				}
			}
		}

		//# small deletion, small insertion, close, 4 break points in one window
		else {
			if(debug) {
				cout << "insod here-26:\n";
			}
			if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1 && bp_window[cur_mate_end] == 1) {
				if(debug) {
					cout << "insod here-27:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1 + 101) && mstart < (positiona2 - positiona1 + 101)) {

						//#print "i start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-28:\n";
					}
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);
					if(debug) {
						cout << "insod here-28: ref_boundary1\n";
					}
//					for(auto& e: ref_boundary1) {
//						cout << "insod here-28: " << e << "\n";
//					}
//					cout << "insod here-28: ref_boundary2\n";
//					for(auto& e: ref_boundary2) {
//						cout << "insod here-28: " << e << "\n";
//					}
//					cout << "insod here-28: ref_support_sr\n";
//					for(auto& e: ref_support_sr) {
//						cout << "insod here-28: " << e << "\n";
//					}
//					cout << "insod here-28: bpread_1\n";
//					for(auto& e: bpread_1) {
//						cout << "insod here-28: " << e.first << "\n";
//						for(auto& s: e.second) {
//							cout << "insod here-28: " << s << "\n";
//						}
//					}

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n";
					//# 2 cluster event
					if(debug) {
						cout << "insod here-29:\n";
					}
					if (ref_support_sr.size() > 1) {
						if(debug) {
							cout << "insod here-30:\n";
						}
						if (ref_boundary1.size() > 3 && ref_boundary1[1] < ref_boundary1[3]) {
							if(debug) {
								cout << "insod here-31:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = positiona1 + ref_boundary2[1] + options.cut_sr;
							}
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = positiona1 + ref_boundary1[2];
							}
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = positiona1 + ref_boundary2[2];
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								left_bound_sr1 = cur_start;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								right_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								left_bound_sr2 = cur_end;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								right_bound_sr2 = cur_mate_start;
							}
							int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insod here-32:\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = left_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = data[0];
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr1;
								//result_sr[sri][6] = left_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							} else {
								if(debug) {
									cout << "insod here-33:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insod here-33-a:\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insod here-34:\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
							}
						} else {
							if(debug) {
								cout << "insod here-35:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;

							if (ref_boundary1.size() > 3 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
								if(debug) {
									cout << "insod-d here-35-a:\n";
								}
								left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
							}

							if (ref_boundary2.size() > 3 && ref_boundary2[3] < (positiona2 - positiona1 + 101)) {
								if(debug) {
									cout << "insod-d here-35-b:\n";
								}
								right_bound_sr1 = positiona1 + ref_boundary2[3] + options.cut_sr;
							}

							if (ref_boundary1.size() > 0 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
								if(debug) {
									cout << "insod-d here-35-c:\n";
								}
								left_bound_sr2 = positiona1 + ref_boundary1[0];
							}

							if (ref_boundary2.size() > 0 && ref_boundary2[0] < (positiona2 - positiona1 + 101)) {
								if(debug) {
									cout << "insod-d here-35-d:\n";
								}
								right_bound_sr2 = positiona1 + ref_boundary2[0];
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
									cout << "insod-d here-35-e:\n";
								}
								left_bound_sr1 = cur_start;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
									cout << "insod-d here-35-f:\n";
								}
								right_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
									cout << "insod-d here-35-g:\n";
								}
								left_bound_sr2 = cur_end;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
									cout << "insod-d here-35-h:\n";
								}
								right_bound_sr2 = cur_mate_start;
							}
							int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insod here-35-a:\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = left_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = data[0];
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr1;
								//result_sr[sri][6] = left_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							} else {
								if(debug) {
									cout << "insod here-36:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insod here-36-a:\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insod here-37:\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
							}
						}
					} else if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "insod here-38:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						int64_t left_bound_sr = 0;
						if (ref_boundary1.size() > 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
							left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
						}
						int64_t right_bound_sr = 0;
						if (ref_boundary2.size() > 1 && ref_boundary2[1] < (positiona2 - positiona1 + 101)) {
							right_bound_sr = positiona1 + ref_boundary2[1] + options.cut_sr;
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "insod here-38-a:\n";
							}
							int64_t size = right_bound_sr - left_bound_sr;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				}
			}
			//# 4 break points in window 2
			else {
				if(debug) {
					cout << "insod here-39:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start >= (positiona2 - positiona1 + 101) && mstart >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insod here-40:\n";
						}
						//#print "i start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-41:\n";
					}
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n";
					//# 2 cluster event
					if (ref_support_sr.size() > 1) {
						if(debug) {
							cout << "insod here-42:\n";
						}
						if (ref_boundary1.size() > 1 && ref_boundary1[1] < ref_boundary1[3]) {
							if(debug) {
								cout << "insod here-43:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								left_bound_sr1 = cur_start;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								right_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								left_bound_sr2 = cur_end;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								right_bound_sr2 = cur_mate_start;
							}
							int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insod here-44:\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = left_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = data[0];
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr1;
								//result_sr[sri][6] = left_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							} else {
								if(debug) {
									cout << "insod here-45:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insod here-45-a:\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insod here-46:\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
							}
						} else {
							if(debug) {
								cout << "insod here-47:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = ref_boundary1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}

							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								left_bound_sr1 = cur_start;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								right_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								left_bound_sr2 = cur_end;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								right_bound_sr2 = cur_mate_start;
							}
							int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insod here-48:\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = left_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = data[0];
								//result_sr[sri][1] = data[1];
								//result_sr[sri][2] = data[2];
								//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
								//result_sr[sri][4] = data[3];
								//result_sr[sri][5] = left_bound_sr1;
								//result_sr[sri][6] = left_bound_sr2;
								//result_sr[sri][7] = del_size;
								//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
							} else {
								if(debug) {
									cout << "insod here-49:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insod here-50:\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insod here-51:\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
								}
							}
						}
					} else if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "insod here-52:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						int64_t left_bound_sr = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						int64_t right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "insod here-52-a:\n";
							}
							int64_t size = right_bound_sr - left_bound_sr;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				}
			}
		}
	} else {
		if(debug) {
			cout << "insod here-53:\n";
		}
		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		if (bp_window[cur_start] != bp_window[cur_mate_end]) {
			if(debug) {
				cout << "insod here-54:\n";
			}
			if (orientationa == 1) {
				if(debug) {
					cout << "insod here-55:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-56:\n";
					}
					sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if(src_return1.size() > 2) {
						right_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (src_return1.size() > 2 && !src_return1[2]) {
						right_bound_sr1 = cur_mate_end;
					}
				}
			} else {
				if(debug) {
					cout << "insod here-57:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-58:\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						right_bound_sr1 = cur_mate_end;
					}
				}
			}
		} else {
			if(debug) {
				cout << "insod here-59:\n";
			}

			if (bp_window[cur_start] == 1 && bp_window[cur_mate_end] == 1) {
				if(debug) {
					cout << "insod here-60:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-61:\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = positiona1 + src_return1[3] + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						right_bound_sr1 = cur_mate_end;
					}
				}
			}
			if (bp_window[cur_start] == 2 && bp_window[cur_mate_end] == 2) {
				if(debug) {
					cout << "insod here-62:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-63:\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = src_return1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_start;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						right_bound_sr1 = cur_mate_end;
					}
				}
			}
		}

		// cl[1] parts
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		if (bp_window[cur_end] != bp_window[cur_mate_start]) {
			if(debug) {
				cout << "insod here-64:\n";
			}
			if (orientationb == 1) {
				if(debug) {
					cout << "insod here-65:\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {
						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-66:\n";
					}
					sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if(src_return2.size() > 3) {
						right_bound_sr2 = positionb4 - (src_return2[3] - (positionb2 - positionb1 + 101) + options.cut_sr);
					}
					if (src_return1.size() > 0 && !src_return1[0]) {
						left_bound_sr2 = cur_end;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						right_bound_sr2 = cur_mate_start;
					}
				}
			} else {
				if(debug) {
					cout << "insod here-67:\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-68:\n";
					}
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if(src_return2.size() > 2) {
						right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (!src_return2.empty() && !src_return2[0]) {
						left_bound_sr2 = cur_end;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						right_bound_sr2 = cur_mate_start;
					}
				}
			}
		} else {
			if(debug) {
				cout << "insod here-69:\n";
			}
			if (bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1) {
				if(debug) {
					cout << "insod here-70:\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart < (positionb2 - positionb1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);
					if(debug) {
						cout << "insod here-71:" << src_return2.size() << "/" << bpread_2.size() << "/" << discord_sr1.size() << "\n";
					}

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 2) {
						right_bound_sr2 = positionb1 + src_return2[2];
					}
					if (!src_return2.empty() && !src_return2[0]) {
						left_bound_sr2 = cur_end;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						right_bound_sr2 = cur_mate_start;
					}
				}
				if(debug) {
					cout << "insod here-71-b:\n";
				}
			}
			if (bp_window[cur_end] == 2 && bp_window[cur_mate_start] == 2) {
				if(debug) {
					cout << "insod here-72-a:\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positionb2 - positionb1 + 100) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insod here-72:\n";
					}
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if(src_return2.size() > 0) {
						left_bound_sr2 = src_return2[0] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if(src_return2.size() > 2) {
						right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (!src_return2.empty() && !src_return2[0]) {
						left_bound_sr2 = cur_end;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						right_bound_sr2 = cur_mate_start;
					}
				}
			}
		}

		int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
		int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			string a_line = (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads1).str();
			if(debug) {
				cout << "insod here-73:" << a_line;
			}
			BPREAD << a_line;
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = data[3];
			an_entry.event_start = left_bound_sr1 + 1;
			an_entry.event_end = left_bound_sr2 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = data[3];
			an_entry.mate_event_start = right_bound_sr2 + 1;
			an_entry.mate_event_end = right_bound_sr1 + 1;
			an_entry.event_size_2 = ins_size;
			result_sr.push_back(an_entry);
			//result_sr[sri][0] = data[0];
			//result_sr[sri][1] = data[1];
			//result_sr[sri][2] = data[2];
			//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
			//result_sr[sri][4] = data[3];
			//result_sr[sri][5] = left_bound_sr1;
			//result_sr[sri][6] = left_bound_sr2;
			//result_sr[sri][7] = del_size;
			//result_sr[sri][8] = "data[3]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
		} else {
			if(debug) {
				cout << "insod here-74:\n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1.size() > 4 && src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
				if(debug) {
					cout << "insod here-74-a:\n";
				}
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "invers_f";
				an_entry.cluster_id = cluster_id1;
				an_entry.n_supports = mpd1;
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
				if(debug) {
					cout << "insod here-75:\n";
				}
				int64_t left_bound_sr = left_bound_sr2;
				int64_t right_bound_sr = right_bound_sr2;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "invers_r";
				an_entry.cluster_id = cluster_id2;
				an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
				an_entry.n_supports = src_return2[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
			}
		}
	}

}

void SplitReadSVCaller::detect_insou(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "965_0/1127_0" == data[1];
	const bool debug = false;

	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);
	if(debug) {
		cout << "insou-d pre here-0:\n";
	}
	const auto& a_region_a = cluster_region.find(cl[0]);
	if(debug) {
		cout << "insou-d pre here-1:\n";
	}
	const auto& a_region_b = cluster_region.find(cl[1]);
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	int32_t orientationa = a_region_a->second.orientation;

	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;
	int32_t orientationb = a_region_b->second.orientation;
	if(debug) {
		cout << "insou-d pre here-2: " << cl[0] << "\n";
	}
	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	if(debug) {
		cout << "insou-d pre here-3:\n";
	}
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;
	if(debug) {
		cout << "insou-d pre here-4:\n";
	}
	if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
		bp_window[cur_mate_end] = 1;
	}
	if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
		bp_window[cur_mate_end] = 2;
	}
	if (covered(cur_start, cur_start, positiona1, positiona2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, positiona3, positiona4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
		bp_window[cur_mate_start] = 1;
	}
	if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
		bp_window[cur_mate_start] = 2;
	}
	if (covered(cur_end, cur_end, positionb1, positionb2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, positionb3, positionb4)) {
		bp_window[cur_end] = 2;
	}

	if(debug) {
		cout << "insou-d pre here-5:\n";
	}
	int64_t n_bp_zeros = 0;
	if (0 == bp_window[cur_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_end]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_end]) {
		++n_bp_zeros;
	}
	if(debug) {
		cout << "insou-d pre here-6:\n";
	}
	if (n_bp_zeros >= 2) {
		if(debug) {
			cout << "insou-d pre here-7:\n";
		}
		swap(cl[0], cl[1]);

		const auto& a_region_a = cluster_region.find(cl[0]);
		positiona1 = a_region_a->second.start;
		positiona2 = a_region_a->second.end;
		positiona3 = a_region_a->second.mate_start;
		positiona4 = a_region_a->second.mate_end;
		//orientationa = a_region_a->second.orientation;
		if(debug) {
			cout << "insou-d pre here-8:\n";
		}
		const auto& a_region_b = cluster_region.find(cl[1]);
		positionb1 = a_region_b->second.start;
		positionb2 = a_region_b->second.end;
		positionb3 = a_region_b->second.mate_start;
		positionb4 = a_region_b->second.mate_end;

		//orientationb = a_region_b->second.orientation;
		if(debug) {
		cout << "insou positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
		cout << "insou positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
		cout << "insou-d pre here-9:\n";
		}
		const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
		if(debug) {
			cout << "insou-d pre here-10:\n";
		}
		const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

		local_ref_id = rev_index->second.first;
		local_ref_end = rev_index->second.second;
		local_ref_id_second = rev_index_second->second.first;
		local_ref_end_second = rev_index_second->second.second;
		if(debug) {
			cout << "insou-d pre here-11:\n";
		}
		if (covered(cur_mate_end, cur_mate_end, positiona1, positiona2)) {
			bp_window[cur_mate_end] = 1;
		}
		if (covered(cur_mate_end, cur_mate_end, positiona3, positiona4)) {
			bp_window[cur_mate_end] = 2;
		}
		if (covered(cur_start, cur_start, positiona1, positiona2)) {
			bp_window[cur_start] = 1;
		}
		if (covered(cur_start, cur_start, positiona3, positiona4)) {
			bp_window[cur_start] = 2;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb1, positionb2)) {
			bp_window[cur_mate_start] = 1;
		}
		if (covered(cur_mate_start, cur_mate_start, positionb3, positionb4)) {
			bp_window[cur_mate_start] = 2;
		}
		if (covered(cur_end, cur_end, positionb1, positionb2)) {
			bp_window[cur_end] = 1;
		}
		if (covered(cur_end, cur_end, positionb3, positionb4)) {
			bp_window[cur_end] = 2;
		}
	}

	//#print "a_region_a->second\na_region_b->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\ncur_mate_start\tbp_window{cur_mate_start}\ncur_mate_end\tbp_window{cur_mate_end}\n";
	if (a_region_a->second == a_region_b->second) {
		if(debug) {
			cout << "insou-d here-0:\n";
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);

		if (bp_window[cur_start] != bp_window[cur_mate_end] || bp_window[cur_end] != bp_window[cur_mate_start]) {
			if(debug) {
				cout << "insou-d here-1:\n";
			}
			discord_sr1.clear();
			if (orientationa == 1) {
				if(debug) {
					cout << "insou-d here-2:\n";
				}
				discord_sr1.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {
						//#print "1 start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
			}
			discord_sr2.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				//int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				//int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand) {

					//#print "2 start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}

			//# small deletion, large insertion, close, 3 break points in one window
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads) && discord_sr2.size() >= static_cast<uint64_t>(options.support_reads) && bp_window[cur_start] == 2 && bp_window[cur_mate_end] == 2) {

				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);

				//#print "src_return1\nsrc_return2\n"
				if(debug) {
					cout << "insou-d here-3: src_return1: " << src_return1[0] << ", " << src_return1[1] << ", " << src_return1[2]
							<< ", src_return2: " << src_return2[0] << ", " << src_return2[1] << ", " << src_return2[2] << "\n";
				}
				int64_t left_bound_sr1 = 0;
				int64_t right_bound_sr1 = 0;
				int64_t left_bound_sr2 = 0;
				int64_t right_bound_sr2 = 0;

				if (src_return1.size() > 2 && src_return1[2] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "insou-d here-4: " << (positiona2 - positiona1 + 101) << "\n";
					}
					left_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
				}

				if (src_return1.size() > 0 && src_return1[0] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "insou-d here-5:\n";
					}
					right_bound_sr1 = positiona4 - (src_return1[0] - (positiona2 - positiona1 + 101));
				}

				if (src_return2.size() > 0 && src_return2[0] < (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "insou-d here-6:\n";
					}
					left_bound_sr2 = positiona1 + src_return2[0];
				}

				if (src_return2.size() > 3 && src_return2[3] >= (positiona2 - positiona1 + 101)) {
					if(debug) {
						cout << "insou-d here-7:\n";
					}
					right_bound_sr2 = positiona4 - (src_return2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "insou-d here-8:\n";
					}
					left_bound_sr1 = cur_mate_end;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "insou-d here-9:\n";
					}
					right_bound_sr1 = cur_start;
				}
				if (!src_return2.empty() && !src_return2[0]) {
					if(debug) {
						cout << "insou-d here-10:\n";
					}
					left_bound_sr2 = cur_mate_start;
				}
				if (src_return2.size() > 3 && !src_return2[3]) {
					if(debug) {
						cout << "insou-d here-11:\n";
					}
					right_bound_sr2 = cur_end;
				}
				int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
				int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

				if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2
						&& right_bound_sr2) {
					if(debug) {
						cout << "insou here-0\n";
					}
					string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
					string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.n_mate_support = src_return2[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = right_bound_sr1 + 1;
					an_entry.event_end = right_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = data[3];
					an_entry.mate_event_start = left_bound_sr2 + 1;
					an_entry.mate_event_end = left_bound_sr1 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = data[0];
					//result_sr[sri][1] = data[1];
					//result_sr[sri][2] = data[2];
					//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
					//result_sr[sri][4] = data[3];
					//result_sr[sri][5] = right_bound_sr1;
					//result_sr[sri][6] = right_bound_sr2;
					//result_sr[sri][7] = del_size;
					//result_sr[sri][8] = "data[3]\tleft_bound_sr2\tleft_bound_sr1\tins_size\n";
				} else {
					if(debug) {
						cout << "insou-d here-12:\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					//string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					//string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					//int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
//					if (src_return1.size() > 4 && src_return1.size() > 4 && src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
					//#int64_t left_bound_sr = left_bound_sr1;//#int64_t right_bound_sr = right_bound_sr1;//#int64_t size = right_bound_sr - left_bound_sr;//#my reads = join("\t", @{bpread_1{0}});//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#
//					}
					if (src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insou here-1\n";
						}
						int64_t left_bound_sr = left_bound_sr2;
						int64_t right_bound_sr = right_bound_sr2;
						int64_t size = right_bound_sr - left_bound_sr;
						string reads = castle::StringUtils::join(bpread_2[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = "invers_r";
						an_entry.cluster_id = cluster_id2;
						an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
						an_entry.n_supports = src_return2[4];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
					}
				}
			}
			//# small deletion, small insertion, far, 2 break points in one window, the other 2 in the other window
			else if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "insou-d here-13:\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 2);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary2.size() > 2 && ref_boundary2[2] > ref_boundary2[0]) {
					if(debug) {
						cout << "insou-d here-14:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 3 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-15:\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}

					if (ref_boundary2.size() > 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-16:\n";
						}
						right_bound_sr1 = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
					}

					if (ref_boundary1.size() > 0 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-17:\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[0];
					}

					if (ref_boundary2.size() > 1 && ref_boundary2[1] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-18:\n";
						}
						right_bound_sr2 = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}

					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
							cout << "insou-d here-19:\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
							cout << "insou-d here-20:\n";
						}
						right_bound_sr1 = cur_start;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
							cout << "insou-d here-21:\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "insou-d here-22:\n";
						}
						right_bound_sr2 = cur_end;
					}
					int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
					int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insou here-2\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 - 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = right_bound_sr1 + 1;
						an_entry.event_end = right_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = left_bound_sr2 + 1;
						an_entry.mate_event_end = left_bound_sr1 - 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = right_bound_sr1;
						//result_sr[sri][6] = right_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tleft_bound_sr2\tleft_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insou-d here-23:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);

						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "insou-d here-24:\n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_bpread.end() != ref_bpread.find(1)) {
								reads = castle::StringUtils::join(ref_bpread[1], "\t");
							} else {
								if(ref_boundary1.size() > 1) {
									left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
								}
								if(ref_boundary2.size() > 0) {
									right_bound_sr = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
								}
							}
							if (left_bound_sr && right_bound_sr) {
								if(debug) {
									cout << "insou here-3\n";
								}
								int64_t size = right_bound_sr - left_bound_sr;
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_f";
								an_entry.cluster_id = cluster_id1;
								an_entry.n_supports = mpd1;
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
							}
						}
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
							if(debug) {
								cout << "insou here-4\n";
							}
							int64_t left_bound_sr = left_bound_sr2;
							int64_t right_bound_sr = right_bound_sr2;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_r";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				} else {
					if(debug) {
						cout << "insou-d here-25:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-26:\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}

					if (ref_boundary2.size() > 0 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-27:\n";
						}
						right_bound_sr1 = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
					}

					if (ref_boundary1.size() > 2 && ref_boundary1[2] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-28:\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[2];
					}

					if (ref_boundary2.size() > 3 && ref_boundary2[3] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-29:\n";
						}
						right_bound_sr2 = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}

					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "insou-d here-30:\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
							cout << "insou-d here-31:\n";
						}
						right_bound_sr1 = cur_start;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
							cout << "insou-d here-32:\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
							cout << "insou-d here-33:\n";
						}

						right_bound_sr2 = cur_end;
					}
					int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
					int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insou here-5\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = right_bound_sr1 + 1;
						an_entry.event_end = right_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = left_bound_sr2 + 1;
						an_entry.mate_event_end = left_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = right_bound_sr1;
						//result_sr[sri][6] = right_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tleft_bound_sr2\tleft_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insou-d here-34:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
							if(debug) {
								cout << "insou here-6\n";
							}
							int64_t left_bound_sr = left_bound_sr1;
							int64_t right_bound_sr = right_bound_sr1;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "insou-d here-35:\n";
							}
							int64_t left_bound_sr = left_bound_sr2;
							int64_t right_bound_sr = right_bound_sr2;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_bpread.end() != ref_bpread.find(1)) {
								reads = castle::StringUtils::join(ref_bpread[1], "\t");
							} else {
								if(ref_boundary1.size() > 0) {
									left_bound_sr = positiona1 + ref_boundary1[0];
								}
								if(ref_boundary2.size() > 1) {
									right_bound_sr = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
								}
							}
							if (left_bound_sr && right_bound_sr) {
								if(debug) {
									cout << "insou here-7\n";
								}
								int64_t size = right_bound_sr - left_bound_sr;
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
								//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
							}
						}
					}
				}
			}
			//# large deletion, small insertion, close, 3 break points in one window
			else if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "insou-d here-36:\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr2, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary2.size() > 2 && ref_boundary2[2] > ref_boundary2[0]) {
					if(debug) {
						cout << "insou-d here-37:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 1 && ref_boundary1[1] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-38:\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}

					if (ref_boundary2.size() > 1 && ref_boundary2[1] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-39:\n";
						}
						right_bound_sr1 = positiona1 + ref_boundary2[1] + options.cut_sr;
					}

					if (ref_boundary1.size() > 2 && ref_boundary1[2] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-40:\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[2];
					}

					if (ref_boundary2.size() > 2 && ref_boundary2[2] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-41:\n";
						}
						right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "insou-d here-42:\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "insou-d here-43:\n";
						}
						right_bound_sr1 = cur_start;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
							cout << "insou-d here-44:\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
							cout << "insou-d here-45:\n";
						}
						right_bound_sr2 = cur_end;
					}
					int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
					int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insou here-8\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = data[3];
						an_entry.event_start = right_bound_sr1 + 1;
						an_entry.event_end = right_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = left_bound_sr2 + 1;
						an_entry.mate_event_end = left_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = right_bound_sr1;
						//result_sr[sri][6] = right_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[3]\tleft_bound_sr2\tleft_bound_sr1\tins_size\n";
					} else {
						if(debug) {
							cout << "insou-d here-46:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						//string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						//string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						//int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
//						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						//#int64_t left_bound_sr = left_bound_sr1;//#int64_t right_bound_sr = right_bound_sr1;//#my reads = join("\t", @{ref_bpread{0}});//#if (ref_bpread{1})//#{//#reads = join("\t", @{ref_bpread{1}});//#}//#else//#{//#left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;//#right_bound_sr = positiona4 - (ref_boundary2[0] - (positiona2-positiona1+101));//#}//#if (left_bound_sr && right_bound_sr)//#{//#int64_t size = right_bound_sr - left_bound_sr;//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#//#}
//						}
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
							if(debug) {
								cout << "insou here-9\n";
							}
							int64_t left_bound_sr = left_bound_sr2;
							int64_t right_bound_sr = right_bound_sr2;
							int64_t size = right_bound_sr - left_bound_sr;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_r";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "invers_r\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";
						}
					}
				} else {
					if(debug) {
						cout << "insou-d here-47:\n";
					}
					int64_t left_bound_sr1 = 0;
					int64_t right_bound_sr1 = 0;
					int64_t left_bound_sr2 = 0;
					int64_t right_bound_sr2 = 0;

					if (ref_boundary1.size() > 3 && ref_boundary1[3] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-48:\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}

					if (ref_boundary2.size() > 3 && ref_boundary2[3] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-49:\n";
						}
						right_bound_sr1 = positiona1 + ref_boundary2[3] + options.cut_sr;
					}

					if (ref_boundary1.size() > 0 && ref_boundary1[0] < (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-50:\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[0];
					}

					if (ref_boundary2.size() > 0 && ref_boundary2[0] >= (positiona2 - positiona1 + 101)) {
						if(debug) {
							cout << "insou-d here-51:\n";
						}
						right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
							cout << "insou-d here-52:\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
							cout << "insou-d here-53:\n";
						}
						right_bound_sr1 = cur_start;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
							cout << "insou-d here-54:\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
							cout << "insou-d here-55:\n";
						}
						right_bound_sr2 = cur_end;
					}
					int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
					int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
						if(debug) {
							cout << "insou here-10\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = right_bound_sr1 + 1;
						an_entry.event_end = right_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = data[3];
						an_entry.mate_event_start = left_bound_sr2 + 1;
						an_entry.mate_event_end = left_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
					} else {
						if(debug) {
							cout << "insou-d here-56:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						//string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						//string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						//int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
//						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						//#int64_t left_bound_sr = left_bound_sr1;//#int64_t right_bound_sr = right_bound_sr1;//#int64_t size = right_bound_sr - left_bound_sr;//#my reads = join("\t", @{ref_bpread{0}});//#BPREAD <<  "data[3]__left_bound_sr__data[3]__right_bound_sr\nreads\n";//#result_sr[sri][0] = "invers_f\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\tright_bound_sr\tsize\n";//#
//						}
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "insou-d here-57:\n";
							}
							int64_t left_bound_sr = left_bound_sr2;
							int64_t right_bound_sr = right_bound_sr2;
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_bpread.end() != ref_bpread.find(1)) {
								reads = castle::StringUtils::join(ref_bpread[1], "\t");
							} else {
								if(ref_boundary1.size() > 0) {
									left_bound_sr = positiona1 + ref_boundary1[0];
								}
								if(ref_boundary2.size() > 1) {
									right_bound_sr = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
								}
							}
							if (left_bound_sr && right_bound_sr) {
								if(debug) {
									cout << "insou here-11\n";
								}
								int64_t size = right_bound_sr - left_bound_sr;
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
						}
					}
				}
			}
		}

		//# small deletion, small insertion, close, 4 break points in one window
		else {
			if(debug) {
				cout << "insou-d here-58:\n";
			}

			if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1 && bp_window[cur_mate_end] == 1) {
				if(debug) {
					cout << "insou-d here-59:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-60:\n";
					}
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n@{ref_support_sr}\n";
					//# 2 cluster event
					if (ref_support_sr.size() > 1) {
						if(debug) {
							cout << "insou-d here-61:\n";
						}
						if (ref_boundary2.size() > 3 && ref_boundary2[1] < ref_boundary2[3]) {
							if(debug) {
								cout << "insou-d here-62:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = positiona1 + ref_boundary2[1] + options.cut_sr;
							}
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = positiona1 + ref_boundary1[2];
							}
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = positiona1 + ref_boundary2[2];
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
									cout << "insou-d here-63:\n";
								}
								left_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
									cout << "insou-d here-64:\n";
								}
								right_bound_sr1 = cur_start;
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
									cout << "insou-d here-65:\n";
								}
								left_bound_sr2 = cur_mate_start;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
									cout << "insou-d here-66:\n";
								}
								right_bound_sr2 = cur_end;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insou here-12\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = left_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							} else {
								if(debug) {
									cout << "insou-d here-67:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insou here-13\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insou here-14\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							if(debug) {
								cout << "insou-d here-68:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = positiona1 + ref_boundary2[3] + options.cut_sr;
							}
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = positiona1 + ref_boundary1[0];
							}
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = positiona1 + ref_boundary2[0];
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
									cout << "insou-d here-69:\n";
								}
								left_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
									cout << "insou-d here-70:\n";
								}
								right_bound_sr1 = cur_start;
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
									cout << "insou-d here-71:\n";
								}
								left_bound_sr2 = cur_mate_start;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
									cout << "insou-d here-72:\n";
								}
								right_bound_sr2 = cur_end;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insou here-15\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = left_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							} else {
								if(debug) {
									cout << "insou-d here-73:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insou here-16\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insou here-17\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					} else if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "insou-d here-74:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 1) {
							left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 1) {
							right_bound_sr = positiona1 + ref_boundary2[1] + options.cut_sr;
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "insou here-18\n";
							}
							int64_t size = right_bound_sr - left_bound_sr;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
						}
					}
				}
			}

			//# 4 break points in window 2
			else {
				if(debug) {
					cout << "insou-d here-75:\n";
				}
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start >= (positiona2 - positiona1 + 100) && mstart >= (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-76:\n";
					}
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n@{ref_support_sr}\n";
					//# 2 cluster event
					if (ref_support_sr.size() > 1) {
						if(debug) {
							cout << "insou-d here-77:\n";
						}
						if (ref_boundary2.size() > 3 && ref_boundary2[1] < ref_boundary2[3]) {
							if(debug) {
								cout << "insou-d here-78:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
									cout << "insou-d here-79:\n";
								}
								left_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
									cout << "insou-d here-80:\n";
								}
								right_bound_sr1 = cur_start;
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
									cout << "insou-d here-81:\n";
								}
								left_bound_sr2 = cur_mate_start;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
									cout << "insou-d here-82:\n";
								}
								right_bound_sr2 = cur_end;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insou here-19\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = left_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							} else {
								if(debug) {
									cout << "insou-d here-83:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insou here-20\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insou here-21\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							if(debug) {
								cout << "insou-d here-84:\n";
							}
							int64_t left_bound_sr1 = 0;
							int64_t right_bound_sr1 = 0;
							int64_t left_bound_sr2 = 0;
							int64_t right_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = ref_boundary1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
									cout << "insou-d here-85:\n";
								}
								left_bound_sr1 = cur_mate_end;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
									cout << "insou-d here-86:\n";
								}
								right_bound_sr1 = cur_start;
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
									cout << "insou-d here-87:\n";
								}
								left_bound_sr2 = cur_mate_start;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
									cout << "insou-d here-88:\n";
								}
								right_bound_sr2 = cur_end;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
								if(debug) {
									cout << "insou here-22\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = data[0];
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = left_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							} else {
								if(debug) {
									cout << "insou-d here-89:\n";
								}
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
									if(debug) {
										cout << "insou here-23\n";
									}
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
									if(debug) {
										cout << "insou here-24\n";
									}
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					} else if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "insou-d here-90:\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						//string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						//string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						//int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 1) {
							left_bound_sr = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 1) {
							right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						if (left_bound_sr && right_bound_sr) {
							if(debug) {
								cout << "insou here-25\n";
							}
							int64_t size = right_bound_sr - left_bound_sr;
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "invers_f";
							an_entry.cluster_id = cluster_id1;
							an_entry.n_supports = mpd1;
							an_entry.n_mate_support = ref_support_sr[0];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = size;
							result_sr.push_back(an_entry);
						}
					}
				}
			}
		}
	} else {
		if(debug) {
			cout << "insou-d here-91\n";
		}
		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (bp_window[cur_start] != bp_window[cur_mate_end]) {
			if(debug) {
				cout << "insou-d here-92\n";
			}
			if (orientationa == 1) {
				if(debug) {
					cout << "insou-d here-93\n";
				}
				discord_sr1.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-94\n";
					}
					sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1];
					}
					if(src_return1.size() > 2) {
						right_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_mate_end;
					}
					if (src_return1.size() > 2 && !src_return1[2]) {
						right_bound_sr1 = cur_start;
					}
				}
			} else {
				if(debug) {
					cout << "insou-d here-95\n";
				}
				discord_sr1.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-96\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						if(debug) {
							cout << "insou-d here-97\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						if(debug) {
							cout << "insou-d here-98\n";
						}
						right_bound_sr1 = cur_start;
					}
				}
			}
		} else {
			if (bp_window[cur_start] == 1 && bp_window[cur_mate_end] == 1) {
				if(debug) {
					cout << "insou-d here-99\n";
				}
				discord_sr1.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-100\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = positiona1 + src_return1[3] + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						if(debug) {
							cout << "insou-d here-101\n";
						}
						left_bound_sr1 = cur_mate_end;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						if(debug) {
							cout << "insou-d here-102\n";
						}
						right_bound_sr1 = cur_start;
					}
				}
			} else if (bp_window[cur_start] == 2 && bp_window[cur_mate_end] == 2) {
				if(debug) {
					cout << "insou-d here-103\n";
				}
				discord_sr1.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-104\n";
					}
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if(src_return1.size() > 1) {
						left_bound_sr1 = src_return1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if(src_return1.size() > 3) {
						right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						left_bound_sr1 = cur_mate_end;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						right_bound_sr1 = cur_start;
					}
				}
			}
		}

		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_b->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		//int64_t local_ref_id = reverse_index_ref_id[a_region_a->second.str].first;
		//int64_t local_ref_start = 0;
//	int64_t local_ref_end = reverse_index_ref_id[a_region_a->second.str].second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

		//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
		//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
		if (bp_window[cur_end] != bp_window[cur_mate_start]) {
			if(debug) {
				cout << "insou-d here-105\n";
			}
			if (orientationb == 1) {
				if(debug) {
					cout << "insou-d here-106\n";
				}
				discord_sr2.clear();

				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-107\n";
					}
					sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
					if(src_return2.size() > 0) {
						left_bound_sr2 = positionb1 + src_return2[0] + options.cut_sr;
					}
					if(src_return2.size() > 3) {
						right_bound_sr2 = positionb4 - (src_return2[3] - (positionb2 - positionb1 + 101));
					}
					if (!src_return2.empty() && !src_return2[0]) {
						if(debug) {
							cout << "insou-d here-108\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (src_return2.size() > 3 && !src_return2[3]) {
						if(debug) {
							cout << "insou-d here-109\n";
						}
						right_bound_sr2 = cur_end;
					}
				}
			} else {
				if(debug) {
					cout << "insou-d here-110\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-111\n";
					}
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						if(debug) {
							cout << "insou-d here-112\n";
						}
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 2) {
						if(debug) {
							cout << "insou-d here-113\n";
						}
						right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (!src_return2.empty() && !src_return2[0]) {
						if(debug) {
							cout << "insou-d here-114\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						if(debug) {
							cout << "insou-d here-116\n";
						}
						right_bound_sr2 = cur_end;
					}
				}
			}
		} else {
			if (bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1) {
				if(debug) {
					cout << "insou-d here-117\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart < (positionb2 - positionb1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-118\n";
					}
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						if(debug) {
							cout << "insou-d here-119\n";
						}
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 2) {
						if(debug) {
							cout << "insou-d here-120\n";
						}
						right_bound_sr2 = positionb1 + src_return2[2];
					}
					if (!src_return2.empty() && !src_return2[0]) {
						if(debug) {
							cout << "insou-d here-121\n";
						}
						left_bound_sr2 = cur_mate_start;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						if(debug) {
							cout << "insou-d here-122\n";
						}
						right_bound_sr2 = cur_end;
					}
				}
			} else if (bp_window[cur_end] == 2 && bp_window[cur_mate_start] == 2) {
				if(debug) {
					cout << "insou-d here-123\n";
				}
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positionb2 - positionb1 + 100) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					if(debug) {
						cout << "insou-d here-124\n";
					}
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if(src_return2.size() > 0) {
						left_bound_sr2 = src_return2[0] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if(src_return2.size() > 2) {
						right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (!src_return2.empty() && !src_return2[0]) {
						cout << "insou-d here-125\n";
						left_bound_sr2 = cur_mate_start;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						cout << "insou-d here-126\n";
						right_bound_sr2 = cur_end;
					}
				}
			}
		}

		int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
		int64_t ins_size = left_bound_sr1 - left_bound_sr2 + 1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1 && left_bound_sr2 && right_bound_sr2) {
			if(debug) {
				cout << "insou here-26\n";
			}
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = data[3];
			an_entry.event_start = right_bound_sr1 + 1;
			an_entry.event_end = right_bound_sr2 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = data[3];
			an_entry.mate_event_start = left_bound_sr2 + 1;
			an_entry.mate_event_end = left_bound_sr1 + 1;
			an_entry.event_size_2 = ins_size;
			result_sr.push_back(an_entry);
		} else {
			if(debug) {
				cout << "insou-d here-127\n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				string a_line = (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				if(debug) {
					cout << "insou here-27: " << a_line << "\n";
				}
				BPREAD << a_line;
				EventEntry an_entry;
				an_entry.type = "invers_f";
				an_entry.cluster_id = cluster_id1;
				an_entry.n_supports = mpd1;
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads && left_bound_sr2 && right_bound_sr2) {
				if(debug) {
					cout << "insou here-28\n";
				}
				int64_t left_bound_sr = left_bound_sr2 + 1;
				int64_t right_bound_sr = right_bound_sr2 - 1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % left_bound_sr % data[3] % right_bound_sr % reads).str();
				EventEntry an_entry;
				an_entry.type = "invers_r";
				an_entry.cluster_id = cluster_id2;
				an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
				an_entry.n_supports = src_return2[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr;
				an_entry.event_end = right_bound_sr;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
			}
		}
	}

}

void SplitReadSVCaller::detect_invers(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;
	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	int32_t orientationa = a_region_a->second.orientation;

	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;

	//int32_t orientationb = a_region_b->second.orientation;
	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;

	if (covered(cur_start, cur_start, positiona1, positiona2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, positiona3, positiona4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_mate_start, cur_mate_start, positiona1, positiona2)) {
		bp_window[cur_mate_start] = 1;
	}
	if (covered(cur_mate_start, cur_mate_start, positiona3, positiona4)) {
		bp_window[cur_mate_start] = 2;
	}
	if (covered(cur_end, cur_end, positionb1, positionb2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, positionb3, positionb4)) {
		bp_window[cur_end] = 2;
	}
	if (covered(cur_mate_end, cur_mate_end, positionb1, positionb2)) {
		bp_window[cur_mate_end] = 1;
	}
	if (covered(cur_mate_end, cur_mate_end, positionb3, positionb4)) {
		bp_window[cur_mate_end] = 2;
	}

	int64_t n_bp_zeros = 0;
	if (0 == bp_window[cur_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_end]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_start]) {
		++n_bp_zeros;
	}
	if (0 == bp_window[cur_mate_end]) {
		++n_bp_zeros;
	}
	if (n_bp_zeros >= 2) {
		swap(cl[0], cl[1]);

		const auto& a_region_a = cluster_region.find(cl[0]);
		positiona1 = a_region_a->second.start;
		positiona2 = a_region_a->second.end;
		positiona3 = a_region_a->second.mate_start;
		positiona4 = a_region_a->second.mate_end;
		//orientationa = a_region_a->second.orientation;

		const auto& a_region_b = cluster_region.find(cl[1]);
		positionb1 = a_region_b->second.start;
		positionb2 = a_region_b->second.end;
		positionb3 = a_region_b->second.mate_start;
		positionb4 = a_region_b->second.mate_end;

		//orientationb = a_region_b->second.orientation;

//		if(debug) {
			//cout << "invers positiona: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 << "\n";
			//cout << "invers positionb: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 << "\n";
//		}

		const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
		const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);
		local_ref_id = rev_index->second.first;
		local_ref_start = 0;
		local_ref_end = rev_index->second.second;
		local_ref_id_second = rev_index_second->second.first;
		local_ref_start_second = 0;
		local_ref_end_second = rev_index_second->second.second;

		if (covered(cur_start, cur_start, positiona1, positiona2)) {
			bp_window[cur_start] = 1;
		}
		if (covered(cur_start, cur_start, positiona3, positiona4)) {
			bp_window[cur_start] = 2;
		}
		if (covered(cur_mate_start, cur_mate_start, positiona1, positiona2)) {
			bp_window[cur_mate_start] = 1;
		}
		if (covered(cur_mate_start, cur_mate_start, positiona3, positiona4)) {
			bp_window[cur_mate_start] = 2;
		}
		if (covered(cur_end, cur_end, positionb1, positionb2)) {
			bp_window[cur_end] = 1;
		}
		if (covered(cur_end, cur_end, positionb3, positionb4)) {
			bp_window[cur_end] = 2;
		}
		if (covered(cur_mate_end, cur_mate_end, positionb1, positionb2)) {
			bp_window[cur_mate_end] = 1;
		}
		if (covered(cur_mate_end, cur_mate_end, positionb3, positionb4)) {
			bp_window[cur_mate_end] = 2;
		}
	}

	//#print "a_region_a->second\na_region_b->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\ncur_mate_start\tbp_window{cur_mate_start}\ncur_mate_end\tbp_window{cur_mate_end}\n";
	if (a_region_a->second == a_region_b->second) {
		if (bp_window[cur_start] != bp_window[cur_mate_start] || bp_window[cur_end] != bp_window[cur_mate_end]) {
			//my @alignments;
			//open( SAM,
			//"samtools_command view -X sr_sortbam a_region_a->second|"
			//);
			//while ( newline1 = <SAM> ) {
			//chomp newline1;
			//push @alignments, newline1;
			//}
			//close SAM;

//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

			//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
			//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
			if (orientationa == 1) {
//	discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;
					sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

					//#print "@{ref_boundary1}\t\n@{ref_boundary2}\n@{ref_support_sr}\n";

					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						int64_t left_bound_sr = positiona1 + int((ref_boundary1[1] + ref_boundary1[2] + options.cut_sr) / 2);
						int64_t right_bound_sr = positiona4 - (int((ref_boundary2[1] + ref_boundary2[2] + options.cut_sr) / 2) - (positiona2 - positiona1 + 101));
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers here-0\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);
						} else {
							castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
							string& cluster_id1 = cluster_ids[0];
							string& cluster_id2 = cluster_ids[1];
							castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
							string& mpd_str_1 = mpds[0];
							string& mpd_str_2 = mpds[1];
							int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
							int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
								//cout << "invers here-1\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_f";
								an_entry.cluster_id = cluster_id1;
								an_entry.n_supports = mpd1;
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
							if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-2\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 1) {
								left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 3) {
								right_bound_sr1 = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101));
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 2) {
								left_bound_sr2 = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 0) {
								right_bound_sr2 = positiona4 - ((ref_boundary2[0] + options.cut_sr) - (positiona2 - positiona1 + 101));
							}
							int64_t del_size = right_bound_sr1 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr2 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr1 - right_bound_sr2;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-3\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr2 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-4\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-5\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 3) {
								left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 1) {
								right_bound_sr1 = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101));
							}

							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 0) {
								left_bound_sr2 = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 2) {
								right_bound_sr2 = positiona4 - ((ref_boundary2[2] + options.cut_sr) - (positiona2 - positiona1 + 101));
							}
							int64_t del_size = right_bound_sr1 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr2 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr1 - right_bound_sr2;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-6\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr2 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-7\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr2;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									if (ref_bpread.end() != ref_bpread.find(1)) {
										reads = castle::StringUtils::join(ref_bpread[1], "\t");
									} else {
										if (ref_boundary1.size() > 1) {
											left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
										}
									}
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-8\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					}
				}
			} else {
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\t\n@{ref_boundary2}\n@{ref_support_sr}\n";
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						int64_t left_bound_sr = 0;
						if (ref_boundary1.size() > 2) {
							left_bound_sr = positiona1 + int((ref_boundary1[1] + ref_boundary1[2] + options.cut_sr) / 2);
						}
						int64_t right_bound_sr = 0;
						if (ref_boundary2.size() > 2) {
							right_bound_sr = int((ref_boundary2[1] + ref_boundary2[2] + options.cut_sr) / 2 - (positiona2 - positiona1 + 101)) + positiona3;
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers here-9\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);

						} else {
							castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
							string& cluster_id1 = cluster_ids[0];
							string& cluster_id2 = cluster_ids[1];
							castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
							string& mpd_str_1 = mpds[0];
							string& mpd_str_2 = mpds[1];
							int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
							int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
								//cout << "invers here-10\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_f";
								an_entry.cluster_id = cluster_id1;
								an_entry.n_supports = mpd1;
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
							if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-11\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);

							}
						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 1) {
								left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 2) {
								left_bound_sr2 = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] + options.cut_sr - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t del_size = right_bound_sr1 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr2 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr1 - right_bound_sr2;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-12\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr2 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-13\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-14\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 3) {
								left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 0) {
								left_bound_sr2 = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] + options.cut_sr - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t del_size = right_bound_sr1 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr2 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr1 - right_bound_sr2;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-15\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr2 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-16\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									if (ref_bpread.end() != ref_bpread.find(1)) {
										reads = castle::StringUtils::join(ref_bpread[1], "\t");
									} else {
										if (ref_boundary2.size() > 1) {
											left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
										}
									}
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-17\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					}
				}
			}
		}

		//# 4 break points in one window
		// !(bp_window[cur_start] != bp_window[cur_mate_start] || bp_window[cur_end] != bp_window[cur_mate_end])
		else {
			//my @alignments;
			//open( SAM,
			//"samtools_command view -X sr_sortbam a_region_a->second"
			//);
			//while ( newline1 = <SAM> ) {
			//chomp newline1;
			//push @alignments, newline1;
			//}
			//close SAM;

			if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1 && bp_window[cur_mate_start] == 1 && bp_window[cur_mate_end] == 1) {
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n@{ref_support_sr}\n";
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						int64_t left_bound_sr = 0;
						if (ref_boundary1.size() > 2) {
							left_bound_sr = positiona1 + ref_boundary1[2];
						}
						int64_t right_bound_sr = 0;
						if (ref_boundary2.size() > 2) {
							right_bound_sr = positiona1 + ref_boundary2[2];
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers here-18\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);

						} else {
							castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
							string& cluster_id1 = cluster_ids[0];
							string& cluster_id2 = cluster_ids[1];
							castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
							string& mpd_str_1 = mpds[0];
							string& mpd_str_2 = mpds[1];
							int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
							int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
								//cout << "invers here-19\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_f";
								an_entry.cluster_id = cluster_id1;
								an_entry.n_supports = mpd1;
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
							if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-20\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);

							}
						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 1) {
								left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 0) {
								right_bound_sr1 = positiona1 + ref_boundary2[0] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 2) {
								left_bound_sr2 = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 3) {
								right_bound_sr2 = positiona1 + ref_boundary2[3];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-21\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-22\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-23\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 3) {
								left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 2) {
								right_bound_sr1 = positiona1 + ref_boundary2[2] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 0) {
								left_bound_sr2 = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 1) {
								right_bound_sr2 = positiona1 + ref_boundary2[1];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-24\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-25\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									if (ref_bpread.end() != ref_bpread.find(1)) {
										reads = castle::StringUtils::join(ref_bpread[1], "\t");
									} else {
										if (ref_boundary1.size() > 1) {
											left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
										}
										if (ref_boundary2.size() > 0) {
											right_bound_sr = positiona1 + ref_boundary2[0] + options.cut_sr;
										}
									}
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-26\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					}
				}
			}
			if (bp_window[cur_start] == 2 && bp_window[cur_end] == 2 && bp_window[cur_mate_start] == 2 && bp_window[cur_mate_end] == 2) {
				//cout << "invers here-27\n";
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					//cout << "invers here-28\n";
					vector<int64_t> ref_boundary1;
					vector<int64_t> ref_boundary2;
					vector<int64_t> ref_support_sr;
					map<int64_t, vector<string>> ref_bpread;

					sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

					//#print "@{ref_boundary1}\n@{ref_boundary2}\n@{ref_support_sr}\n";
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						//cout << "invers here-29\n";
						int64_t left_bound_sr = 0;
						if (ref_boundary1.size() > 2) {
							left_bound_sr = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
						}
						int64_t right_bound_sr = 0;
						if (ref_boundary2.size() > 2) {
							right_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers here-30\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);

						} else {
							//cout << "invers here-31\n";
							castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
							string& cluster_id1 = cluster_ids[0];
							string& cluster_id2 = cluster_ids[1];
							castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
							string& mpd_str_1 = mpds[0];
							string& mpd_str_2 = mpds[1];
							int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
							int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
								//cout << "invers here-32\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_f";
								an_entry.cluster_id = cluster_id1;
								an_entry.n_supports = mpd1;
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);
							}
							if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-33\n";
								int64_t size = right_bound_sr - left_bound_sr;
								string reads = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
								EventEntry an_entry;
								an_entry.type = "invers_r";
								an_entry.cluster_id = cluster_id2;
								an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
								an_entry.n_supports = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr + 1;
								an_entry.event_end = right_bound_sr + 1;
								an_entry.event_size_1 = size;
								result_sr.push_back(an_entry);

							}
						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							//cout << "invers here-34\n";
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 1) {
								left_bound_sr1 = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 0) {
								right_bound_sr1 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 2) {
								left_bound_sr2 = ref_boundary1[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 3) {
								right_bound_sr2 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-35\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							} else {
								//cout << "invers here-36\n";
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-37\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-38\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[1], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						} else {
							int64_t left_bound_sr1 = 0;
							if (ref_boundary1.size() > 3) {
								left_bound_sr1 = ref_boundary1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if (ref_boundary2.size() > 2) {
								right_bound_sr1 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if (ref_boundary1.size() > 0) {
								left_bound_sr2 = ref_boundary1[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t right_bound_sr2 = 0;
							if (ref_boundary2.size() > 1) {
								right_bound_sr2 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;
							//cout << "invers here-39\n";
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers here-40\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							} else {
								castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
								string& cluster_id1 = cluster_ids[0];
								string& cluster_id2 = cluster_ids[1];
								castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
								string& mpd_str_1 = mpds[0];
								string& mpd_str_2 = mpds[1];
								int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
								int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
								//cout << "invers here-41\n";
								if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
									//cout << "invers here-42\n";
									int64_t left_bound_sr = left_bound_sr1;
									int64_t right_bound_sr = right_bound_sr1;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									if (ref_bpread.end() != ref_bpread.find(1)) {
										reads = castle::StringUtils::join(ref_bpread[1], "\t");
									} else {
										if (ref_boundary1.size() > 1) {
											left_bound_sr = ref_boundary1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
										}
										if (ref_boundary2.size() > 0) {
											right_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
										}
									}
									int64_t size = right_bound_sr - left_bound_sr;
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_f";
									an_entry.cluster_id = cluster_id1;
									an_entry.n_supports = mpd1;
									an_entry.n_mate_support = ref_support_sr[0];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);
								}
								if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
									//cout << "invers here-43\n";
									int64_t left_bound_sr = left_bound_sr2;
									int64_t right_bound_sr = right_bound_sr2;
									int64_t size = right_bound_sr - left_bound_sr;
									string reads = castle::StringUtils::join(ref_bpread[0], "\t");
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "invers_r";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = data[3];
									an_entry.event_start = left_bound_sr + 1;
									an_entry.event_end = right_bound_sr + 1;
									an_entry.event_size_1 = size;
									result_sr.push_back(an_entry);

								}
							}
						}
					}
				}
			}
		}
	} else {
		//cout << "invers here-44\n";
		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

//		int64_t local_ref_id = reverse_index_ref_id[a_region_a->second.str].first;
//		int64_t local_ref_start = 0;
//		int64_t local_ref_end = reverse_index_ref_id[a_region_a->second.str].second;
//		int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
//		int64_t local_ref_start_second = 0;
//		int64_t local_ref_end_second = reverse_index_ref_id[a_region_b->second.str].second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

		//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
		//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
		if (bp_window[cur_start] != bp_window[cur_mate_start]) {
			//cout << "invers here-45\n";
			if (orientationa == 1) {
				//cout << "invers here-46\n";
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
//	vector<int64_t> src_return1;
//	map<int64_t, vector<string>> bpread;
					sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
					if (src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if (src_return1.size() > 2) {
						right_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
					}
				}
			} else {
				//cout << "invers here-47\n";
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
					if (src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if (src_return1.size() > 2) {
						right_bound_sr1 = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
				}
			}
		} else {
			//cout << "invers here-48\n";
			if (bp_window[cur_start] == 1 && bp_window[cur_mate_start] == 1) {
				//cout << "invers here-49\n";
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart < (positiona2 - positiona1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if (src_return1.size() > 1) {
						left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
					if (src_return1.size() > 3) {
						right_bound_sr1 = positiona1 + src_return1[3] + options.cut_sr;
					}
				}
			}
			if (bp_window[cur_start] == 2 && bp_window[cur_mate_start] == 2) {
				//cout << "invers here-50\n";
				discord_sr1.clear();
				reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positiona2 - positiona1 + 100) && mstart > (positiona2 - positiona1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr1.push_back(al);
					}
				}
				if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

					//#print "src_return1\n";
					if (src_return1.size() > 1) {
						left_bound_sr1 = src_return1[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 3) {
						right_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
				}
			}
		}

		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_b->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		//int64_t local_ref_id = reverse_index_ref_id[a_region_a->second.str].first;
		//int64_t local_ref_start = 0;
//	int64_t local_ref_end = reverse_index_ref_id[a_region_a->second.str].second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

		//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
		//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
		if (bp_window[cur_end] != bp_window[cur_mate_end]) {
			if (orientationa == 1) {
				//cout << "invers here-51\n";
				discord_sr2.clear();

				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 3) {
						right_bound_sr2 = positionb4 - (src_return2[3] - (positionb2 - positionb1 + 101) + options.cut_sr);
					}
				}
			} else {
				//cout << "invers here-52\n";
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 3) {
						right_bound_sr2 = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
					}
				}
			}
		} else {
			if (bp_window[cur_end] == 1 && bp_window[cur_mate_end] == 1) {
				//cout << "invers here-53\n";
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart < (positionb2 - positionb1)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (!src_return2.empty()) {
						left_bound_sr2 = positionb1 + src_return2[0];
					}
					if (src_return2.size() > 2) {
						right_bound_sr2 = positionb1 + src_return2[2];
					}
				}
			} else if (bp_window[cur_end] == 2 && bp_window[cur_mate_end] == 2) {
				//cout << "invers here-54\n";
				discord_sr2.clear();
				BamTools::BamAlignment al;
				//int64_t prev_bam_pos = backend_bgzf.Tell();
				while (reader.GetNextAlignmentBasic(al)) {
					int64_t start = al.Position;
					int32_t strand = 1;
					if (al.IsReverseStrand()) {
						strand = -1;
					}
//	string mseqid = ref_vec[al.MateRefID].RefName;
					int64_t mstart = al.MatePosition;
					int32_t mstrand = 1;
					if (al.IsMateReverseStrand()) {
						mstrand = -1;
					}
					int64_t isize = al.InsertSize;

					if (isize > 0 && strand != mstrand && start > (positionb2 - positionb1 + 100) && mstart > (positionb2 - positionb1 + 100)) {

						//#print "start mstart\tisize\tstrand mstrand\n";
						discord_sr2.push_back(al);
					}
				}
				if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
					sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);

					//#print "src_return2\n";
					if (src_return2.size() > 0) {
						left_bound_sr2 = src_return2[0] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (src_return2.size() > 2) {
						right_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
				}
			}
		}

		int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
		int64_t inv_size = right_bound_sr1 - left_bound_sr2;
		int64_t distance1 = left_bound_sr2 - left_bound_sr1;
		int64_t distance2 = right_bound_sr2 - right_bound_sr1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads) {
			//cout << "invers here-55\n";
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = data[3];
			an_entry.event_start = left_bound_sr1 + 1;
			an_entry.event_end = right_bound_sr2 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = data[3];
			an_entry.mate_event_start = left_bound_sr2 + 1;
			an_entry.mate_event_end = right_bound_sr1 + 1;
			an_entry.event_size_2 = inv_size;
			an_entry.distance_1 = distance1;
			an_entry.distance_2 = distance2;
			result_sr.push_back(an_entry);
		} else {
			//cout << "invers here-56\n";
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				//cout << "invers here-57\n";
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "invers_f";
				an_entry.cluster_id = cluster_id1;
				an_entry.n_supports = mpd1;
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads) {
				//cout << "invers here-58\n";
				int64_t left_bound_sr = left_bound_sr2;
				int64_t right_bound_sr = right_bound_sr2;
				int64_t size = right_bound_sr - left_bound_sr;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = "invers_r";
				an_entry.cluster_id = cluster_id2;
				an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
				an_entry.n_supports = src_return2[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
			}
		}
	}

}

void SplitReadSVCaller::detect_tandem_dup(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "4774_0" == data[1];
	const bool debug = false;
	auto& ref_is = options.is;
//	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	int64_t position1 = a_region_a->second.start;
	int64_t position2 = a_region_a->second.end;
	int64_t position3 = a_region_a->second.mate_start;
	int64_t position4 = a_region_a->second.mate_end;
	//int32_t orientation = a_region_a->second.orientation;
	if (covered(cur_start, cur_start, position1, position2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, position3, position4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_end, cur_end, position1, position2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, position3, position4)) {
		bp_window[cur_end] = 2;
	}
	//#print "a_region_a->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\n";
	//my @alignments;
	//open( SAM,
	//"samtools_command view -X sr_sortbam a_region_a->second|"
	//);
	//while ( newline1 = <SAM> ) {
	//chomp newline1;
	//push @alignments, newline1;
	//}
	//close SAM;

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

	//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
	//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
	reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
	if (bp_window[cur_start] != bp_window[cur_end]) {
		if(debug) {
			cout << "tandem_dup here-0\n";
		}

		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

				//#print "start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "tandem_dup here-1\n";
			}
			sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
			int64_t left_bound_sr = 0;
			if (src_return1.size() > 0) {
				if(debug) {
					cout << "tandem_dup here-2\n";
				}
				left_bound_sr = position1 + src_return1[0];
			}
			int64_t right_bound_sr = 0;
			if (src_return1.size() > 3) {
				if(debug) {
					cout << "tandem_dup here-3\n";
				}
				right_bound_sr = src_return1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
			}
			if (src_return1.size() > 0 && !src_return1[0]) {
				if(debug) {
					cout << "tandem_dup here-4\n";
				}
				left_bound_sr = cur_start;
			}
			if (src_return1.size() > 3 && !src_return1[3]) {
				if(debug) {
					cout << "tandem_dup here-5\n";
				}
				right_bound_sr = cur_end;
			}
			int64_t size = right_bound_sr - left_bound_sr;
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				if(debug) {
					cout << "tandem_dup here-6\n";
				}
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = data[3];
				an_entry.event_start = left_bound_sr + 1;
				an_entry.event_end = right_bound_sr + 1;
				an_entry.event_size_1 = size;
				result_sr.push_back(an_entry);
			}
		}
	} else {
		if(debug) {
			cout << "tandem_dup here-7\n";
		}
		if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1) {
			if(debug) {
				cout << "tandem_dup here-8\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start < (position2 - position1) && mstart < (position2 - position1)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "tandem_dup here-9\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "tandem_dup here-10\n";
					}
					//# inssu
					if (ref_boundary1.size() > 3 && covered(ref_boundary1[0], ref_boundary1[1], ref_boundary1[2], ref_boundary1[3])) {
						if(debug) {
							cout << "tandem_dup here-11\n";
						}
						if (ref_boundary2.size() > 0 && ref_boundary2[0] < ref_boundary2[2]) {
							if(debug) {
								cout << "tandem_dup here-12\n";
							}

							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = position1 + ref_boundary1[0];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = position1 + ref_boundary2[1] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = position1 + ref_boundary2[2];
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-13\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						} else {
							if(debug) {
								cout << "tandem_dup here-14\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = position1 + ref_boundary1[2];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = position1 + ref_boundary2[3] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = position1 + ref_boundary2[0];
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-15\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						}
					}

					//# inssd
					if (ref_boundary2.size() > 3 && covered(ref_boundary2[0], ref_boundary2[1], ref_boundary2[2], ref_boundary2[3])) {
						if(debug) {
							cout << "tandem_dup here-16\n";
						}
						if (ref_boundary1.size() > 2 && ref_boundary1[0] < ref_boundary1[2]) {
							if(debug) {
								cout << "tandem_dup here-17\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = position1 + ref_boundary1[2];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = position1 + ref_boundary2[3] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = position1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = position1 + ref_boundary2[0];
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-18\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						} else {
							if(debug) {
								cout << "tandem_dup here-19\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = position1 + ref_boundary1[0];
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = position1 + ref_boundary2[1] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = position1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = position1 + ref_boundary2[2];
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-20\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						}
					}
				}

				//# single cluster
				else {
					if(debug) {
						cout << "tandem_dup here-21\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 0) {
						left_bound_sr = position1 + ref_boundary1[0];
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 1) {
						right_bound_sr = position1 + ref_boundary2[1] + options.cut_sr;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
							cout << "tandem_dup here-22\n";
						}
						left_bound_sr = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "tandem_dup here-23\n";
						}
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "tandem_dup here-24\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		} else if (bp_window[cur_start] == 2 && bp_window[cur_end] == 2) {
			if(debug) {
				cout << "tandem_dup here-25\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > (ref_is["rlu"]["selected"] - options.cut_sr + 5) && strand == mstrand && start > (position2 - position1 + 100) && mstart > (position2 - position1 + 100)) {
					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "tandem_dup here-26\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 0);
				if(debug) {
					cout << "tandem_dup here-26: ref_boundary1\n";
				}
//				for (auto& e : ref_boundary1) {
//					cout << "tandem_dup here-26:" << e << "\n";
//				}
//				cout << "tandem_dup here-26: ref_boundary2\n";
//				for (auto& e : ref_boundary2) {
//					cout << "tandem_dup here-26:" << e << "\n";
//				}
//				cout << "tandem_dup here-26: ref_support_sr\n";
//				for (auto& e : ref_support_sr) {
//					cout << "tandem_dup here-26:" << e << "\n";
//				}
//				cout << "tandem_dup here-26: bpread_1\n";
//				for (auto& e : bpread_1) {
//					cout << "tandem_dup here-26:" << e.first << "\n";
//					for (auto& s : e.second) {
//						cout << "tandem_dup here-26:" << s << "\n";
//					}
//				}
				//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "tandem_dup here-27\n";
					}
					//# inssu
					if (ref_boundary1.size() > 3 && covered(ref_boundary1[0], ref_boundary1[1], ref_boundary1[2], ref_boundary1[3])) {
						if(debug) {
							cout << "tandem_dup here-28\n";
						}
						if (ref_boundary2.size() > 0 && ref_boundary2[0] < ref_boundary2[2]) {
							if(debug) {
								cout << "tandem_dup here-29\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-30\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						} else {
							if(debug) {
								cout << "tandem_dup here-31\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - right_bound_sr1 - 1;
							int64_t ins_size = left_bound_sr2 - left_bound_sr1 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-32\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (right_bound_sr1 + 1) % data[3] % (left_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssu";
								} else {
									an_entry.type = "inssu";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[1];
								an_entry.n_mate_support = ref_support_sr[0];
								an_entry.ref_id = data[3];
								an_entry.event_start = right_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr1 + 1;
								an_entry.mate_event_end = left_bound_sr2 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						}
					}

					//# inssd
					if (ref_boundary2.size() > 3 && covered(ref_boundary2[0], ref_boundary2[1], ref_boundary2[2], ref_boundary2[3])) {
						if(debug) {
							cout << "tandem_dup here-33\n";
						}
						if (ref_boundary1.size() > 2 && ref_boundary1[0] < ref_boundary1[2]) {
							if(debug) {
								cout << "tandem_dup here-34\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr1 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr1 = ref_boundary2[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr2 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr2 = ref_boundary2[0] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-35\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						} else {
							if(debug) {
								cout << "tandem_dup here-36\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr1 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr1 = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr2 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr2 = ref_boundary2[2] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = left_bound_sr1 - left_bound_sr2 - 1;
							int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "tandem_dup here-37\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr2 + 1) % data[3] % (right_bound_sr2 + 1) % reads0 % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads1).str();
								EventEntry an_entry;
								if (del_size > 10) {
									an_entry.type = "del_inssd";
								} else {
									an_entry.type = "inssd";
								}
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr2 + 1;
								an_entry.event_end = left_bound_sr1 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = right_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = ins_size;
								result_sr.push_back(an_entry);
							}
						}
					}
				}

				//# single cluster
				else {
					if(debug) {
						cout << "tandem_dup here-38\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 0) {
						left_bound_sr = ref_boundary1[0] - (position2 - position1 + 101) + position3;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 1) {
						right_bound_sr = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
							cout << "tandem_dup here-39\n";
						}
						left_bound_sr = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "tandem_dup here-40\n";
						}
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "tandem_dup here-41\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		}
	}

}

void SplitReadSVCaller::detect_invers_f(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "211466_0" == data[1];
	const bool debug = false;
	auto& ref_is = options.is;
//	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);

	const auto& a_region_a = cluster_region.find(cl[0]);

	int64_t position1 = a_region_a->second.start;
	int64_t position2 = a_region_a->second.end;
	int64_t position3 = a_region_a->second.mate_start;
	int64_t position4 = a_region_a->second.mate_end;
	int32_t orientation = a_region_a->second.orientation;
	if (covered(cur_start, cur_start, position1, position2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, position3, position4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_end, cur_end, position1, position2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, position3, position4)) {
		bp_window[cur_end] = 2;
	}
	//#print "a_region_a->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\n";
	//my @alignments;
	//open( SAM,
	//"samtools_command view -X sr_sortbam a_region_a->second|"
	//);
	//while ( newline1 = <SAM> ) {
	//chomp newline1;
	//push @alignments, newline1;
	//}
	//close SAM;

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
//	int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
	if(debug) {
		cout << "invers_f-d first here-0: " << cl[0] << ": " << local_ref_end << "\n";
		cout << "invers_f-d first here-1: " << position1 << "/" << position2 << "/" << position3 << "/" << position4 << "\n";
	}

	reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
	if (bp_window[cur_start] != bp_window[cur_end]) {
		if(debug) {
			cout << "invers_f-d here-0\n";
		}
		if (orientation == 1) {
			if(debug) {
				cout << "invers_f-d here-1\n";
			}
//	discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {
					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "invers_f-d here-2\n";
					}
					left_bound_sr = position1 + src_return1[1] + options.cut_sr;
				}
				int64_t right_bound_sr = 0;
				if(src_return1.size() > 2) {
					right_bound_sr = position4 - (src_return1[2] - (position2 - position1 + 101));
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
						cout << "invers_f-d here-3\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					right_bound_sr = cur_end;
				}
				int64_t size = right_bound_sr - left_bound_sr;
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "invers_f here-0\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
					an_entry.n_mate_support = src_return1[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr + 1;
					an_entry.event_end = right_bound_sr + 1;
					an_entry.event_size_1 = size;
					result_sr.push_back(an_entry);
				}
			}
		} else {
			if(debug) {
				cout << "invers_f-d here-5\n";
			}
//	discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "invers_f-d here-6\n";
				}
				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "invers_f-d here-7\n";
					}
					left_bound_sr = position1 + src_return1[1] + options.cut_sr;
				}

				int64_t right_bound_sr = 0;
				if (src_return1.size() > 2) {
					if(debug) {
						cout << "invers_f-d here-8\n";
					}
					right_bound_sr = src_return1[2] - (position2 - position1 + 101) + position3;
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
						cout << "invers_f-d here-9\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "invers_f-d here-10\n";
					}
					right_bound_sr = cur_end;
				}
				int64_t size = right_bound_sr - left_bound_sr;
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "invers_f here-1\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
					an_entry.n_mate_support = src_return1[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr + 1;
					an_entry.event_end = right_bound_sr + 1;
					an_entry.event_size_1 = size;
					result_sr.push_back(an_entry);
				}
			}
		}
	} else {
		if(debug) {
			cout << "invers_f-d here-11\n";
		}
		if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1) {
			if(debug) {
				cout << "invers_f-d here-12\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart < (position2 - position1)) {
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "invers_f-d here-13\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "invers_f-d here-14\n";
					}
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						if(debug) {
							cout << "invers_f-d here-15\n";
						}
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 2) {
							left_bound_sr = position1 + ref_boundary1[2];
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr = position1 + ref_boundary2[2];
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "invers_f here-2\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);
						}
					} else {
						if(debug) {
							cout << "invers_f-d here-16\n";
						}
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							if(debug) {
								cout << "invers_f-d here-17\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = position1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr1 = position1 + ref_boundary2[0] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = position1 + ref_boundary1[2];
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr2 = position1 + ref_boundary2[3];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "invers_f here-3\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						} else {
							if(debug) {
								cout << "invers_f-d here-18\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = position1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr1 = position1 + ref_boundary2[2] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = position1 + ref_boundary1[0];
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr2 = position1 + ref_boundary2[1];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "invers_f here-4\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						}
					}
				}
				//# single cluster
				else {
					if(debug) {
						cout << "invers_f-d here-19\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 1) {
						left_bound_sr = position1 + ref_boundary1[1] + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 1) {
						right_bound_sr = position1 + ref_boundary2[1] + options.cut_sr;
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "invers_f-d here-20\n";
						}
						left_bound_sr = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "invers_f-d here-21\n";
						}
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "invers_f here-5: " << data[1] << "/" << data[2] << "\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		}
		if (bp_window[cur_start] == 2 && bp_window[cur_end] == 2) {
			if(debug) {
				cout << "invers_f-d here-22\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start > (position2 - position1 + 100) && mstart > (position2 - position1 + 100)) {
					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "invers_f-d here-23\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "invers_f-d here-24\n";
					}
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						if(debug) {
							cout << "invers_f-d here-25\n";
						}
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 2) {
							left_bound_sr = ref_boundary1[2] - (position2 - position1 + 101) + position3;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr = ref_boundary2[2] - (position2 - position1 + 101) + position3;
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "invers_f here-6\n";
							}
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);
						}
					} else {
						if(debug) {
							cout << "invers_f-d here-26\n";
						}
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							if(debug) {
								cout << "invers_f-d here-27\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr1 = ref_boundary2[0] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr2 = ref_boundary2[3] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "invers_f here-7\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);
							}
						} else {
							if(debug) {
								cout << "invers_f-d here-28\n";
							}
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr1 = ref_boundary2[2] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr2 = ref_boundary2[1] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								if(debug) {
									cout << "invers_f here-8\n";
								}
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						}
					}
				}
				//# single cluster
				else {
					if(debug) {
						cout << "invers_f-d here-29\n";
					}
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 1) {
						left_bound_sr = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 1) {
						right_bound_sr = ref_boundary2[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "invers_f-d here-30\n";
						}
						left_bound_sr = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
							cout << "invers_f-d here-31\n";
						}
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
							cout << "invers_f here-9\n";
						}
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.n_supports = boost::lexical_cast<int64_t>(data[2]);
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		}
	}
}

void SplitReadSVCaller::detect_invers_r(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {

	auto& ref_is = options.is;
//	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;
//cout << "here-1\n";
	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	const auto& a_region_a = cluster_region.find(cl[0]);
	int64_t position1 = a_region_a->second.start;
	int64_t position2 = a_region_a->second.end;
	int64_t position3 = a_region_a->second.mate_start;
	int64_t position4 = a_region_a->second.mate_end;
	int32_t orientation = a_region_a->second.orientation;
	if (covered(cur_start, cur_start, position1, position2)) {
		bp_window[cur_start] = 1;
	}
	if (covered(cur_start, cur_start, position3, position4)) {
		bp_window[cur_start] = 2;
	}
	if (covered(cur_end, cur_end, position1, position2)) {
		bp_window[cur_end] = 1;
	}
	if (covered(cur_end, cur_end, position3, position4)) {
		bp_window[cur_end] = 2;
	}
//cout << "here-2\n";
//#print "a_region_a->second\ncur_start\tbp_window{cur_start}\ncur_end\tbp_window{cur_end}\n";
//my @alignments;
//open( SAM,
//"samtools_command view -X sr_sortbam a_region_a->second|"
//);
//while ( newline1 = <SAM> ) {
//chomp newline1;
//push @alignments, newline1;
//}
//close SAM;

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
//	int64_t local_ref_start_second = a_region_b->second.start;
//	int64_t local_ref_end_second = a_region_b->second.end;
//	int64_t local_ref_size = reverse_index_ref_id[a_region_a->second.str].second;

//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
	if (bp_window[cur_start] != bp_window[cur_end]) {
		if (orientation == 1) {
			//cout << "here-3\n";
//discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
//cout << "here-4\n";
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				//cout << "here-5\n";
//vector<int64_t> src_return;
//map<int64_t, vector<string>> bpread;
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 0) {
					left_bound_sr = position1 + src_return1[0];
				}
				int64_t right_bound_sr = 0;
				if(src_return1.size() > 3) {
					right_bound_sr = position4 - (src_return1[3] - (position2 - position1 + 101) + options.cut_sr);
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					right_bound_sr = cur_end;
				}
				int64_t size = right_bound_sr - left_bound_sr;
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					//cout << "invers_r here-6\n";
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr + 1;
					an_entry.event_end = right_bound_sr + 1;
					an_entry.event_size_1 = size;
					result_sr.push_back(an_entry);
				}
			}
		} else {
			//cout << "here-7\n";
			discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
//cout << "here-8\n";
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				//cout << "here-9\n";
//vector<int64_t> src_return;
//map<int64_t, vector<string>> bpread;
				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 0) {
					left_bound_sr = position1 + src_return1[0];
				}
				int64_t right_bound_sr = 0;
				if (src_return1.size() > 3) {
					right_bound_sr = src_return1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					right_bound_sr = cur_end;
				}
//cout << "here-10\n";
				int64_t size = right_bound_sr - left_bound_sr;
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					//cout << "invers_r here-10-a\n";
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = data[3];
					an_entry.event_start = left_bound_sr + 1;
					an_entry.event_end = right_bound_sr + 1;
					an_entry.event_size_1 = size;
					result_sr.push_back(an_entry);
				}
			}
		}

	} else {
		if (bp_window[cur_start] == 1 && bp_window[cur_end] == 1) {
			discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart < (position2 - position1)) {

//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 2) {
							left_bound_sr = position1 + ref_boundary1[2];
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr = position1 + ref_boundary2[2];
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers_r here-10-a\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);

						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = position1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr1 = position1 + ref_boundary2[0] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = position1 + ref_boundary1[2];
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr2 = position1 + ref_boundary2[3];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers_r here-11-a\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						} else {
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = position1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr1 = position1 + ref_boundary2[2] + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = position1 + ref_boundary1[0];
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr2 = position1 + ref_boundary2[1];
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers_r here-12-a\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						}
					}
				}

//# single cluster
				else {
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 0) {
						left_bound_sr = position1 + ref_boundary1[0];
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 0) {
						right_bound_sr = position1 + ref_boundary2[0];
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						left_bound_sr = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						//cout << "invers_r here-13-a\n";
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		}
		if (bp_window[cur_start] == 2 && bp_window[cur_end] == 2) {
			discord_sr1.clear();
			reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
			BamTools::BamAlignment al;
//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start > (position2 - position1 + 100) && mstart > (position2 - position1 + 100)) {

//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
//# 2 cluster event
				if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
					if (ref_boundary1.size() > 3 && ref_boundary2.size() > 3 &&
							((abs(ref_boundary1[2] - ref_boundary1[1]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[2] - ref_boundary2[1]) < (ref_is["rlu"]["selected"] - options.cut_sr))
							|| (abs(ref_boundary1[0] - ref_boundary1[3]) < (ref_is["rlu"]["selected"] - options.cut_sr) && abs(ref_boundary2[0] - ref_boundary2[3]) < (ref_is["rlu"]["selected"] - options.cut_sr)))) {
						int64_t left_bound_sr = 0;
						if(ref_boundary1.size() > 2) {
							left_bound_sr = ref_boundary1[2] - (position2 - position1 + 101) + position3;
						}
						int64_t right_bound_sr = 0;
						if(ref_boundary2.size() > 2) {
							right_bound_sr = ref_boundary2[2] - (position2 - position1 + 101) + position3;
						}
						int64_t inv_size = right_bound_sr - left_bound_sr;
						if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
							//cout << "invers_r here-14-a\n";
							string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
							string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\t%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads0 % reads1).str();
							EventEntry an_entry;
							an_entry.type = "invers";
							an_entry.cluster_id = data[1];
							an_entry.mate_cluster_id = data[2];
							an_entry.n_supports = ref_support_sr[0];
							an_entry.n_mate_support = ref_support_sr[1];
							an_entry.ref_id = data[3];
							an_entry.event_start = left_bound_sr + 1;
							an_entry.event_end = right_bound_sr + 1;
							an_entry.event_size_1 = inv_size;
							result_sr.push_back(an_entry);

						}
					} else {
						if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr1 = ref_boundary1[1] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr1 = ref_boundary2[0] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 2) {
								left_bound_sr2 = ref_boundary1[2] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 3) {
								right_bound_sr2 = ref_boundary2[3] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers_r here-15-a\n";
								string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						} else {
							int64_t left_bound_sr1 = 0;
							if(ref_boundary1.size() > 3) {
								left_bound_sr1 = ref_boundary1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t right_bound_sr1 = 0;
							if(ref_boundary2.size() > 2) {
								right_bound_sr1 = ref_boundary2[2] - (position2 - position1 + 101) + position3 + options.cut_sr;
							}
							int64_t left_bound_sr2 = 0;
							if(ref_boundary1.size() > 0) {
								left_bound_sr2 = ref_boundary1[0] - (position2 - position1 + 101) + position3;
							}
							int64_t right_bound_sr2 = 0;
							if(ref_boundary2.size() > 1) {
								right_bound_sr2 = ref_boundary2[1] - (position2 - position1 + 101) + position3;
							}
							int64_t del_size = right_bound_sr2 - left_bound_sr1 - 1;
							int64_t inv_size = right_bound_sr1 - left_bound_sr2;
							int64_t distance1 = left_bound_sr2 - left_bound_sr1;
							int64_t distance2 = right_bound_sr2 - right_bound_sr1;

							if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
								//cout << "invers_r here-16-a\n";
								string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
								string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
								BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr1 + 1) % data[3] % (right_bound_sr1 + 1) % reads0 % data[3] % (right_bound_sr2 + 1) % data[3] % (left_bound_sr2 + 1) % reads1).str();
								EventEntry an_entry;
								an_entry.type = "del_invers";
								an_entry.cluster_id = data[1];
								an_entry.mate_cluster_id = data[2];
								an_entry.n_supports = ref_support_sr[0];
								an_entry.n_mate_support = ref_support_sr[1];
								an_entry.ref_id = data[3];
								an_entry.event_start = left_bound_sr1 + 1;
								an_entry.event_end = right_bound_sr2 + 1;
								an_entry.event_size_1 = del_size;
								an_entry.mate_ref_id = data[3];
								an_entry.mate_event_start = left_bound_sr2 + 1;
								an_entry.mate_event_end = right_bound_sr1 + 1;
								an_entry.event_size_2 = inv_size;
								an_entry.distance_1 = distance1;
								an_entry.distance_2 = distance2;
								result_sr.push_back(an_entry);

							}
						}
					}
				}

//# single cluster
				else {
					int64_t left_bound_sr = 0;
					if(ref_boundary1.size() > 0) {
						left_bound_sr = ref_boundary1[0] - (position2 - position1 + 101) + position3;
					}
					int64_t right_bound_sr = 0;
					if(ref_boundary2.size() > 0) {
						right_bound_sr = ref_boundary2[0] - (position2 - position1 + 101) + position3;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						left_bound_sr = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						right_bound_sr = cur_end;
					}
					int64_t size = right_bound_sr - left_bound_sr;
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						//cout << "invers_r here-17-a\n";
						string reads = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % data[3] % (left_bound_sr + 1) % data[3] % (right_bound_sr + 1) % reads).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.ref_id = data[3];
						an_entry.event_start = left_bound_sr + 1;
						an_entry.event_end = right_bound_sr + 1;
						an_entry.event_size_1 = size;
						result_sr.push_back(an_entry);
					}
				}
			}
		}
	}

}

void SplitReadSVCaller::detect_inter_chromosomal_events(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.detect_inter_chromosomal_events");
	checker.start();
	vector<function<void()> > tasks;
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
	vector<int64_t> block_positions;

	{
		block_positions.push_back(0);
		int64_t n_lines = castle::IOUtils::get_number_of_lines(mpinter_outfile);
		int64_t BLOCK_SIZE = (n_lines + n_cores - 1) / (double) n_cores;
		string line;
		ifstream in(mpinter_outfile);
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
	cout << (boost::format("[SplitReadSVCaller.detect_inter_chromosomal_events] # blocks: %d\n") % n_blocks).str();
	vector<vector<EventEntry>> result_sr_list(n_blocks);
	vector<string> output_file_names(n_blocks);

	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string sr_sortbam = options.prefix + ".sr.sorted.bam";
			if(!options.working_dir.empty()) {
				sr_sortbam = options.working_prefix + ".sr.sorted.bam";
			}
			string a_path(sr_sortbam);
			string an_index_path(a_path);
			an_index_path += ".bai";
			BamTools::BamReader reader;
			if (!reader.Open(a_path, an_index_path)) {
				std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
				exit(1);
			}

			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = block_positions[block_id];
			int64_t next_pos = block_positions[block_id + 1];
			string line;
			const char* delim_tab = "\t";
			const char* delim_slash = "/";
			vector<string> data;
			vector<string> cl;
			vector<string> coords;
			vector<string> cluster_ids;
			vector<string> mpds;
			string local_bpinfofile = options.prefix + ".bp_reads." + str_block_id;
			if(!options.working_dir.empty()) {
				local_bpinfofile = options.working_prefix + ".bp_reads." + str_block_id;
			}
			output_file_names[block_id] = local_bpinfofile;
			auto& local_result_sr = result_sr_list[block_id];
			ofstream local_out(local_bpinfofile, ios::binary);
			ifstream local_reader(mpinter_outfile, ios::binary);
			local_reader.seekg(cur_pos, ios::beg);
			while(getline(local_reader, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
				map<int64_t, int32_t> bp_window;
				if (cluster_region.end() == cluster_region.find(cl[0])) {
					continue;
				}
				if (string::npos != data[0].find("inss")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_inss(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (string::npos != data[0].find("inso")) {
					if (cluster_region.end() == cluster_region.find(cl[1])) {
						continue;
					}
					detect_inso(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				} else if (data[0] == "transl_inter") {
					detect_transl_inter(local_out, local_result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
				}
//				if("inss" == data[0] && "86860_0/26254_0" == data[1]) {
//					cout << local_result_sr.back().sr_str() << "\n";
//				}
				if(cur_pos >= next_pos) {
					break;
				}
			}

		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		BPREAD << castle::IOUtils::read_fully(output_file_names[block_id]);
	}
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		auto& local_result_sr = result_sr_list[block_id];
		result_sr.insert(result_sr.end(), local_result_sr.begin(), local_result_sr.end());
	}
	castle::IOUtils::remove_files(output_file_names, n_cores);
	//# add deletion clusters without split read support into other events if supporting complex deletion events
	auto& ref_is = options.is;
	for (auto& del : unsupport_del) {
		//#string toprint = join ("\t", @del); print "toprint";
		for (auto& a_result : result_sr) {
			if (string::npos != a_result.type.find("del")) {
				if (del.ref_id == a_result.ref_id && abs(del.event_start - a_result.event_start) < ref_is["rlu"]["selected"] && abs(del.event_end - a_result.event_end) < ref_is["rlu"]["selected"]) {
					a_result.cluster_id += "/" + del.cluster_id;
					a_result.mate_cluster_id += "/" + del.mate_cluster_id;
				}
			}
		}
	}
	cout << checker;
}
void SplitReadSVCaller::detect_inter_chromosomal_events_alt(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadSVCaller.detect_inter_chromosomal_events_alt");
	checker.start();
	auto& ref_is = options.is;
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
//	int64_t n_elems = 100;
	string line;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	vector<string> data;
	vector<string> cl;
	vector<string> coords;
	vector<string> cluster_ids;
	vector<string> mpds;

	int64_t i = 0;
	ifstream MPINTERSRD(mpinter_outfile, ios::binary);
	while (getline(MPINTERSRD, line, '\n')) {
		++i;
		//if (-1 != line_value && i < line_value) {
		//continue;
		//}
		//if (-1 != line_value && i > line_value) {
		//break;
		//}

		castle::StringUtils::tokenize(line, delim_tab, data);
//		if("del_inso" != data[0] || "91157_0/3464_0" != data[1]) {
//			continue;
//		}
		cout << line << "\n";
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		map<int64_t, int32_t> bp_window;
		if (cluster_region.end() == cluster_region.find(cl[0])) {
			continue;
		}
		if (string::npos != data[0].find("inss")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_inss(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (string::npos != data[0].find("inso")) {
			if (cluster_region.end() == cluster_region.find(cl[1])) {
				continue;
			}
			detect_inso(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		} else if (data[0] == "transl_inter") {
			detect_transl_inter(BPREAD, result_sr, reader, reverse_index_ref_id, cluster_region, bp_window, data, cl, cluster_ids, mpds);
		}
		cout << result_sr.back().sr_str() << "\n";
//		if(i > n_elems) {
//			break;
//		}
	}

	//# add deletion clusters without split read support into other events if supporting complex deletion events
	for (auto& del : unsupport_del) {
		//#string toprint = join ("\t", @del); print "toprint";
		for (auto& a_result : result_sr) {
			if (string::npos != a_result.type.find("del")) {
				if (del.ref_id == a_result.ref_id && abs(del.event_start - a_result.event_start) < ref_is["rlu"]["selected"] && abs(del.event_end - a_result.event_end) < ref_is["rlu"]["selected"]) {
					a_result.cluster_id += "/" + del.cluster_id;
					a_result.mate_cluster_id += "/" + del.mate_cluster_id;
				}
			}
		}
	}
	cout << checker;
}
void SplitReadSVCaller::detect_inss(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "28982_0/28902_0" == data[1];
	const bool debug = false;
//	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);
	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);
	string chr_bp_a1 = a_region_a->second.chr_bp_1;
//	int64_t chr_bp_a1_id = reverse_index_ref_id[chr_bp_a1].first;
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	string chr_bp_a2 = a_region_a->second.chr_bp_2;
//	int64_t chr_bp_a2_id = reverse_index_ref_id[chr_bp_a2].first;
	int64_t positiona3 = a_region_a->second.mate_start;
	//int64_t positiona4 = a_region_a->second.mate_end;
	//int32_t orientationa = a_region_a->second.orientation;

	string chr_bp_b1 = a_region_b->second.chr_bp_1;
//	int64_t chr_bp_b1_id = reverse_index_ref_id[chr_bp_b1].first;
	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	string chr_bp_b2 = a_region_b->second.chr_bp_2;
//	int64_t chr_bp_b2_id = reverse_index_ref_id[chr_bp_b2].first;
	int64_t positionb3 = a_region_b->second.mate_start;
	//int64_t positionb4 = a_region_b->second.mate_end;
	//int32_t orientationb = a_region_b->second.orientation;

	string& ref_id = data[3];
	string& mate_ref_id = data[7];


	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;

	if(debug) {
		cout << "inss: " << cur_start << "/" << cur_end << "/" << cur_mate_start << "/" << cur_mate_end <<
				", chrbp: " << chr_bp_a1 << "/" << chr_bp_a2 << "/" << chr_bp_b1 << "/" << chr_bp_b2 <<
				", position_a: " << positiona1 << "/" << positiona2 << "/" << positiona3 <<
				", position_b: " << positionb1 << "/" << positionb2 << "/" << positionb3 <<
				", local_ref:" << local_ref_id << "/" << local_ref_start << "/" << local_ref_end <<
				", local_ref_second:" << local_ref_id_second << "/" << local_ref_start_second << "/" << local_ref_end_second << "\n" ;
	}

	string& local_chr_bp_a1 = mate_ref_id;
	string& local_chr_bp_a2 = ref_id;

	if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
		if (a_region_a->second == a_region_b->second) {
			cerr << "inss i\n";
		}
		//#(positiona1, positiona2, positiona3, positiona4) = (positiona3, positiona4, positiona1, positiona2);
	}
	//#string& local_chr_bp_1 = data[7];
	//string& local_chr_bp_2 = data[3];

	//if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2){
	//	#(positionb1, positionb2, positionb3, positionb4) = (positionb3, positionb4, positionb1, positionb2);}//#print "a_region_a->second\na_region_b->second\npositiona1, positiona2, positiona3, positiona4\npositionb1, positionb2, positionb3, positionb4\n";

	if (a_region_a->second == a_region_b->second) {
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);

		if(debug) {
		cout << "inss here-0\n";
		}
		discord_sr1.clear();
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {
				discord_sr1.push_back(al);
			}
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
			cout << "inss here-1\n";
			}
			vector<int64_t> ref_boundary1;
			vector<int64_t> ref_boundary2;
			vector<int64_t> ref_support_sr;
			map<int64_t, vector<string>> ref_bpread;
			sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 1);

			if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
				if(debug) {
				cout << "inss here-2\n";
				}

				int64_t left_bound_sr1 = 0;
				if (ref_boundary1.size() > 1) {
					left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
				}
				int64_t right_bound_sr1 = 0;
				if (ref_boundary2.size() > 0) {
					right_bound_sr1 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
				}
				int64_t left_bound_sr2 = 0;
				if (ref_boundary1.size() > 2) {
					left_bound_sr2 = positiona1 + ref_boundary1[2];
				}
				int64_t right_bound_sr2 = 0;
				if (ref_boundary2.size() > 3) {
					right_bound_sr2 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				string& local_chr_bp_a1 = mate_ref_id;
				string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {

					if (ref_boundary2.size() > 1) {
						if(debug) {
						cout << "inss here-4\n";
						}
						left_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (ref_boundary1.size() > 0) {
						if(debug) {
						cout << "inss here-5\n";
						}
						right_bound_sr1 = positiona1 + ref_boundary1[0];
					}
					if (ref_boundary2.size() > 2) {
						if(debug) {
						cout << "inss here-6\n";
						}
						left_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (ref_boundary1.size() > 3) {
						if(debug) {
							cout << "inss here-7\n";
						}
						right_bound_sr2 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}
				}
				if (ref_id < mate_ref_id) {

					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "inss here-8\n";
						}

						left_bound_sr1 = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
													cout << "inss here-9\n";
												}
						right_bound_sr1 = cur_mate_start;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
													cout << "inss here-10\n";
												}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
							cout << "inss here-11\n";
						}
						cout << "inss here-7\n";
						right_bound_sr2 = cur_mate_end;
					}
				} else {
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "inss here-12\n";
						}
						left_bound_sr1 = cur_end;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
							cout << "inss here-13\n";
						}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
							cout << "inss here-14\n";
						}
						left_bound_sr2 = cur_start;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
							cout << "inss here-15\n";
						}
						right_bound_sr2 = cur_mate_start;
					}
				}
				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
					if(debug) {
						cout << "inss here-16\n";
					}
					string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
					string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr2 - right_bound_sr1 + 1;
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[0];
					an_entry.n_mate_support = ref_support_sr[1];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr1 + 1;
					an_entry.event_end = left_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr1 + 1;
					an_entry.mate_event_end = right_bound_sr2 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
				} else {
					if(debug) {
						cout << "inss here-17\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if (ref_id < mate_ref_id) {
							if(debug) {
								cout << "inss here-18\n";
							}
							int64_t left_bound_sr = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = mate_ref_id;
							string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 1) {
									if(debug) {
										cout << "inss here-19\n";
									}
									left_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 0) {
									if(debug) {
										cout << "inss here-20\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[0];
								}

							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
									cout << "inss here-21\n";
								}
								left_bound_sr = cur_start;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
									cout << "inss here-22\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();

							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id1;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = ref_id;
							an_entry.event_start = left_bound_sr + 1;
							an_entry.strand = 1;
							an_entry.mate_ref_id = mate_ref_id;
							an_entry.mate_event_start = right_bound_sr + 1;
							an_entry.mate_strand = -1;
							result_sr.push_back(an_entry);
						} else {
							if(debug) {
								cout << "inss here-23\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 1) {
								if(debug) {
															cout << "inss here-24\n";
														}
								cout << "inss here-22\n";
								left_bound_sr = positiona1 + ref_boundary1[1];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 0) {
								if(debug) {
									cout << "inss here-25\n";
								}
								right_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 1) {
									if(debug) {
																		cout << "inss here-26\n";
																	}
									left_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 0) {
									if(debug) {
																		cout << "inss here-27\n";
																	}
									right_bound_sr = positiona1 + ref_boundary1[0] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
																	cout << "inss here-28\n";
																}
								cout << "inss here-26\n";
								left_bound_sr = cur_end;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
																	cout << "inss here-29\n";
																}
								cout << "inss here-27\n";
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id1;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = mate_ref_id;
							an_entry.event_start = right_bound_sr + 1;
							an_entry.strand = 1;
							an_entry.mate_ref_id = ref_id;
							an_entry.mate_event_start = left_bound_sr + 1;
							an_entry.mate_strand = -1;
							result_sr.push_back(an_entry);
						}
					}
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {

						if (ref_id < mate_ref_id) {
							if(debug) {
																cout << "inss here-30\n";
															}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 2) {
								if(debug) {
																	cout << "inss here-31\n";
																}
								left_bound_sr = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 3) {
								if(debug) {
									cout << "inss here-32\n";
								}
								right_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 2) {
									if(debug) {
																		cout << "inss here-33\n";
																	}
									left_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 3) {
									if(debug) {
										cout << "inss here-34\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
																	cout << "inss here-35\n";
																}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
																	cout << "inss here-36\n";
																}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = ref_id;
							an_entry.event_start = left_bound_sr + 1;
							an_entry.strand = -1;
							an_entry.mate_ref_id = mate_ref_id;
							an_entry.mate_event_start = right_bound_sr + 1;
							an_entry.mate_strand = 1;
							result_sr.push_back(an_entry);
						} else {
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 2) {
								if(debug) {
																	cout << "inss here-38\n";
																}
								left_bound_sr = positiona1 + ref_boundary1[2] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 3) {
								if(debug) {
																	cout << "inss here-39\n";
																}
								right_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = mate_ref_id;
							string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 2) {
									if(debug) {
																		cout << "inss here-40\n";
																	}
									left_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 3) {
									if(debug) {
																		cout << "inss here-41\n";
																	}
									right_bound_sr = positiona1 + ref_boundary1[3];
								}
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
																	cout << "inss here-42\n";
																}
								left_bound_sr = cur_start;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
																	cout << "inss here-43\n";
																}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = mate_ref_id;
							an_entry.event_start = right_bound_sr + 1;
							an_entry.strand = -1;
							an_entry.mate_ref_id = ref_id;
							an_entry.mate_event_start = left_bound_sr + 1;
							an_entry.mate_strand = 1;
							result_sr.push_back(an_entry);
						}
					}
				}
			} else {
				int64_t left_bound_sr1 = 0;
				if (ref_boundary1.size() > 3) {
					if(debug) {
						cout << "inss here-44\n";
					}
					left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
				}
				int64_t right_bound_sr1 = 0;
				if (ref_boundary2.size() > 2) {
					if(debug) {
						cout << "inss here-45\n";
					}
					right_bound_sr1 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
				}
				int64_t left_bound_sr2 = 0;
				if (ref_boundary1.size() > 0) {
					if(debug) {
											cout << "inss here-46\n";
										}
					left_bound_sr2 = positiona1 + ref_boundary1[0];
				}
				int64_t right_bound_sr2 = 0;
				if (ref_boundary2.size() > 1) {
					if(debug) {
											cout << "inss here-47\n";
										}
					right_bound_sr2 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				string& local_chr_bp_a1 = mate_ref_id;
				string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if (ref_boundary2.size() > 3) {
						if(debug) {
												cout << "inss here-48\n";
											}
						left_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (ref_boundary1.size() > 2) {
						if(debug) {
												cout << "inss here-49\n";
											}
						right_bound_sr1 = positiona1 + ref_boundary1[2];
					}
					if (ref_boundary2.size() > 0) {
						if(debug) {
												cout << "inss here-50\n";
											}
						left_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (ref_boundary1.size() > 1) {
						if(debug) {
												cout << "inss here-51\n";
											}
						right_bound_sr2 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}
				}
				if (ref_id < mate_ref_id) {
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
							cout << "inss here-52\n";
						}
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
													cout << "inss here-53\n";
												}
						right_bound_sr1 = cur_mate_start;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
													cout << "inss here-54\n";
												}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
													cout << "inss here-55\n";
												}
						right_bound_sr2 = cur_mate_end;
					}
				} else {
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
													cout << "inss here-56\n";
												}
						left_bound_sr1 = cur_end;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
													cout << "inss here-57\n";
												}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
													cout << "inss here-58\n";
												}
						left_bound_sr2 = cur_start;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
													cout << "inss here-59\n";
												}
						right_bound_sr2 = cur_mate_start;
					}
				}
				if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
					string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
					string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr2 - right_bound_sr1 + 1;
					string a_line = (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
					if(debug) {
						cout << "inss here-60\n";
					}
					BPREAD << a_line;
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = ref_support_sr[1];
					an_entry.n_mate_support = ref_support_sr[0];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr1 + 1;
					an_entry.event_end = left_bound_sr2 + 1;
					an_entry.event_size_1 = del_size;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr1 + 1;
					an_entry.mate_event_end = right_bound_sr2 + 1;
					an_entry.event_size_2 = ins_size;
					result_sr.push_back(an_entry);
				} else {
					if(debug) {
						cout << "inss here-61\n";
					}
					castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
					string& cluster_id1 = cluster_ids[0];
					string& cluster_id2 = cluster_ids[1];
					castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
					string& mpd_str_1 = mpds[0];
					string& mpd_str_2 = mpds[1];
					int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
					int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
					if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
						if (data[3] < data[7]) {
							if(debug) {
													cout << "inss here-62\n";
												}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 3) {
								if(debug) {
														cout << "inss here-63\n";
													}
								left_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 2) {
								if(debug) {
														cout << "inss here-64\n";
													}
								right_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 3) {
									if(debug) {
															cout << "inss here-65\n";
														}
									left_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 2) {
									if(debug) {
															cout << "inss here-66\n";
														}
									right_bound_sr = positiona1 + ref_boundary1[2];
								}
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
														cout << "inss here-67\n";
													}
								left_bound_sr = cur_start;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
														cout << "inss here-68\n";
													}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id1;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = ref_id;
							an_entry.event_start = left_bound_sr + 1;
							an_entry.strand = 1;
							an_entry.mate_ref_id = mate_ref_id;
							an_entry.mate_event_start = right_bound_sr + 1;
							an_entry.mate_strand = -1;
							result_sr.push_back(an_entry);
						} else {
							if(debug) {
													cout << "inss here-69\n";
												}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 3) {
								if(debug) {
									cout << "inss here-70\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[3];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 2) {
								if(debug) {
									cout << "inss here-71\n";
								}
								right_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							string& local_chr_bp_a1 = mate_ref_id;
							string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 3) {
									if(debug) {
																		cout << "inss here-72\n";
																	}
									left_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 2) {
									if(debug) {
																		cout << "inss here-73\n";
																	}
									right_bound_sr = positiona1 + ref_boundary1[2] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
																	cout << "inss here-74\n";
																}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
																	cout << "inss here-75\n";
																}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();

							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id1;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
							an_entry.n_supports = ref_support_sr[1];
							an_entry.ref_id = mate_ref_id;
							an_entry.event_start = right_bound_sr + 1;
							an_entry.strand = 1;
							an_entry.mate_ref_id = ref_id;
							an_entry.mate_event_start = left_bound_sr + 1;
							an_entry.mate_strand = -1;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[7]\tright_bound_sr\t1\tdata[3]\tleft_bound_sr\t-1\n";
						}
					}
					if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
						if(debug) {
															cout << "inss here-76\n";
														}
						if (ref_id < mate_ref_id) {
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 0) {
								if(debug) {
																	cout << "inss here-77\n";
																}
								left_bound_sr = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 1) {
								if(debug) {
									cout << "inss here-78\n";
								}
								right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}

							string& local_chr_bp_a1 = mate_ref_id;
							string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 0) {
									if(debug) {
																		cout << "inss here-79\n";
																	}
									left_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 1) {
									if(debug) {
																		cout << "inss here-80\n";
																	}
									right_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
																	cout << "inss here-81\n";
																}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
																	cout << "inss here-82\n";
																}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = ref_id;
							an_entry.event_start = left_bound_sr + 1;
							an_entry.strand = -1;
							an_entry.mate_ref_id = mate_ref_id;
							an_entry.mate_event_start = right_bound_sr + 1;
							an_entry.mate_strand = 1;
							result_sr.push_back(an_entry);
						} else {
							if(debug) {
								cout << "inss here-83\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 0) {
								if(debug) {
																cout << "inss here-84\n";
															}
								left_bound_sr = positiona1 + ref_boundary1[0] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 1) {
								if(debug) {
																cout << "inss here-85\n";
															}
								right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = mate_ref_id;
							string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 0) {
									if(debug) {
										cout << "inss here-86\n";
									}
									left_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 1) {
									if(debug) {
										cout << "inss here-87\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[1];
								}
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
																cout << "inss here-88\n";
															}
								left_bound_sr = cur_start;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
									cout << "inss here-89\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
							EventEntry an_entry;
							an_entry.type = "transl_inter";
							an_entry.cluster_id = cluster_id2;
							an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
							an_entry.n_supports = ref_support_sr[0];
							an_entry.ref_id = mate_ref_id;
							an_entry.event_start = right_bound_sr + 1;
							an_entry.strand = -1;
							an_entry.mate_ref_id = ref_id;
							an_entry.mate_event_start = left_bound_sr + 1;
							an_entry.mate_strand = 1;
							result_sr.push_back(an_entry);
							//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[7]\tright_bound_sr\t-1\tdata[3]\tleft_bound_sr\t1\n";

						}
					}
				}
			}
		}
	} else {

		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		//int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
		//int64_t local_ref_size_second = reverse_index_ref_id[a_region_b->second.str].second;
		discord_sr1.clear();
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

				//#print "1 start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if(debug) {
			cout << "inss here-90: " << discord_sr1.size() << "\n";
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "inss here-90\n";
			}
			sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
			if (ref_id < mate_ref_id) {
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "inss here-91\n";
					}
					left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
				}
				if (src_return1.size() > 2) {
					if(debug) {
						cout << "inss here-92\n";
					}
					right_bound_sr1 = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3;
				}
				string& local_chr_bp_a1 = mate_ref_id;
				string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if (src_return1.size() > 3) {
						if(debug) {
								cout << "inss here-93\n";
							}
						left_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 0) {
						if(debug) {
									cout << "inss here-94\n";
								}
						right_bound_sr1 = positiona1 + src_return1[0];
					}
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
						cout << "inss here-95\n";
					}
					left_bound_sr1 = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "inss here-96\n";
					}
					right_bound_sr1 = cur_mate_start;
				}
			} else {
				if (src_return1.size() > 0) {
					if(debug) {
								cout << "inss here-97\n";
							}
					left_bound_sr2 = positiona1 + src_return1[0];
				}
				if (src_return1.size() > 3) {
					if(debug) {
								cout << "inss here-98\n";
							}
					right_bound_sr2 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				string& local_chr_bp_a1 = mate_ref_id;
				string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if (src_return1.size() > 2) {
						cout << "inss here-103\n";
						left_bound_sr2 = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
					if (src_return1.size() > 1) {
						if(debug) {
							cout << "inss here-99\n";
						}
						right_bound_sr2 = positiona1 + src_return1[1] + options.cut_sr;
					}
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
								cout << "inss here-100\n";
							}
					left_bound_sr2 = cur_end;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
								cout << "inss here-101\n";
							}
					right_bound_sr2 = cur_mate_end;
				}
			}
		}

		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_b->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		discord_sr2.clear();
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		//BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {
				discord_sr2.push_back(al);
			}
		}
		if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {

			sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
			if (ref_id < mate_ref_id) {
				if(debug) {
					cout << "inss here-102\n";
				}
				if (!src_return2.empty()) {
					if(debug) {
						cout << "inss here-103\n";
					}
					left_bound_sr2 = positionb1 + src_return2[0];
				}
				if(debug) {
					cout << "inss here-104\n";
				}
				if(src_return2.size() > 3) {
					right_bound_sr2 = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
				}
				string& local_chr_bp_a1 = mate_ref_id;
				string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if (src_return2.size() > 2) {
						if(debug) {
							cout << "inss here-105\n";
						}
						left_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (src_return2.size() > 1) {
						if(debug) {
							cout << "inss here-106\n";
						}
						right_bound_sr2 = positionb1 + src_return2[1] + options.cut_sr;
					}
				}
				if (!src_return2.empty() && !src_return2[0]) {
					if(debug) {
						cout << "inss here-107\n";
					}
					left_bound_sr2 = cur_end;
				}
				if (src_return2.size() > 3 && !src_return2[3]) {
					if(debug) {
						cout << "inss here-108\n";
					}
					right_bound_sr2 = cur_mate_end;
				}
			} else {
				if(debug) {
					cout << "inss here-109\n";
				}
				if(src_return2.size() > 1) {
					left_bound_sr1 = positionb1 + src_return2[1] + options.cut_sr;
				}
				if(src_return2.size() > 2) {
					right_bound_sr1 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
				}
				string& local_chr_bp_1 = mate_ref_id;
				string& local_chr_bp_2 = ref_id;

				if (chr_bp_b1 == local_chr_bp_1 && chr_bp_b2 == local_chr_bp_2) {
					if (src_return2.size() > 3) {
						if(debug) {
							cout << "inss here-110\n";
						}
						left_bound_sr1 = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
					}
					if (src_return2.size() > 0) {
						if(debug) {
													cout << "inss here-111\n";
												}
						right_bound_sr1 = positionb1 + src_return2[0];
					}
				}
				if (src_return2.size() > 1 && !src_return2[1]) {
					if(debug) {
												cout << "inss here-112\n";
											}
					left_bound_sr1 = cur_start;
				}
				if (src_return2.size() > 2 && !src_return2[2]) {
					if(debug) {
												cout << "inss here-113\n";
											}
					right_bound_sr1 = cur_mate_start;
				}
			}
		}

		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads) {
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			if (ref_id < mate_ref_id) {
				if(debug) {
					cout << "inss here-114\n";
				}
				int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
				int64_t ins_size = right_bound_sr2 - right_bound_sr1 + 1;
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.n_supports = src_return1[4];
				an_entry.n_mate_support = src_return2[4];
				an_entry.ref_id = ref_id;
				an_entry.event_start = left_bound_sr1 + 1;
				an_entry.event_end = left_bound_sr2 + 1;
				an_entry.event_size_1 = del_size;
				an_entry.mate_ref_id = mate_ref_id;
				an_entry.mate_event_start = right_bound_sr1 + 1;
				an_entry.mate_event_end = right_bound_sr2 + 1;
				an_entry.event_size_2 = ins_size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = data[0];
				//result_sr[sri][1] = data[1];
				//result_sr[sri][2] = data[2];
				//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
				//result_sr[sri][4] = data[3];
				//result_sr[sri][5] = left_bound_sr1;
				//result_sr[sri][6] = left_bound_sr2;
				//result_sr[sri][7] = del_size;
				//result_sr[sri][8] = "data[7]\tright_bound_sr1\tright_bound_sr2\tins_size\n";

			} else {
				if(debug) {
					cout << "inss here-115\n";
				}
				swap(reads0, reads1);
				int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
				int64_t ins_size = right_bound_sr2 - right_bound_sr1 + 1;
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.n_supports = src_return2[4];
				an_entry.n_mate_support = src_return1[4];
				an_entry.ref_id = ref_id;
				an_entry.event_start = left_bound_sr1 + 1;
				an_entry.event_end = left_bound_sr2 + 1;
				an_entry.event_size_1 = del_size;
				an_entry.mate_ref_id = mate_ref_id;
				an_entry.mate_event_start = right_bound_sr1 + 1;
				an_entry.mate_event_end = right_bound_sr2 + 1;
				an_entry.event_size_2 = ins_size;
				result_sr.push_back(an_entry);
				//result_sr[sri][0] = data[0];
				//result_sr[sri][1] = data[1];
				//result_sr[sri][2] = data[2];
				//result_sr[sri][3] = "src_return2[4]/src_return1[4]";
				//result_sr[sri][4] = data[3];
				//result_sr[sri][5] = left_bound_sr1;
				//result_sr[sri][6] = left_bound_sr2;
				//result_sr[sri][7] = del_size;
				//result_sr[sri][8] = "data[7]\tright_bound_sr1\tright_bound_sr2\tins_size\n";

			}
		} else {
			if(debug) {
				cout << "inss here-116\n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				if (ref_id < mate_ref_id) {
					int64_t left_bound_sr = 0;
					if (src_return1.size() > 1) {
						if(debug) {
							cout << "inss here-117\n";
						}
						left_bound_sr = positiona1 + src_return1[1] + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if (src_return1.size() > 2) {
						if(debug) {
							cout << "inss here-118\n";
						}
						right_bound_sr = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
					string& local_chr_bp_a1 = mate_ref_id;
					string& local_chr_bp_a2 = ref_id;

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if (src_return1.size() > 3) {
							if(debug) {
								cout << "inss here-119\n";
							}
							left_bound_sr = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						if (src_return1.size() > 0) {
							if(debug) {
								cout << "inss here-120\n";
							}
							right_bound_sr = positiona1 + src_return1[0];
						}
					}
					if (src_return1.size() > 1 && !src_return1[1]) {
						if(debug) {
							cout << "inss here-121\n";
						}
						left_bound_sr = cur_start;
					}
					if (src_return1.size() > 2 && !src_return1[2]) {
						if(debug) {
														cout << "inss here-122\n";
													}
						right_bound_sr = cur_mate_start;
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id1;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = 1;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr + 1;
					an_entry.mate_strand = -1;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[3]\tleft_bound_sr\t1\tdata[7]\tright_bound_sr\t-1\n";

				} else {
					int64_t left_bound_sr = 0;
					if (src_return1.size() > 0) {
						if(debug) {
							cout << "inss here-123\n";
						}
						left_bound_sr = positiona1 + src_return1[0];
					}
					int64_t right_bound_sr = 0;
					if (src_return1.size() > 3) {
						if(debug) {
							cout << "inss here-124\n";
						}
						right_bound_sr = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					string& local_chr_bp_a1 = mate_ref_id;
					string& local_chr_bp_a2 = ref_id;

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if (src_return1.size() > 2) {
							if(debug) {
								cout << "inss here-125\n";
							}
							left_bound_sr = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3;
						}
						if (src_return1.size() > 1) {
							if(debug) {
								cout << "inss here-126\n";
							}
							right_bound_sr = positiona1 + src_return1[1] + options.cut_sr;
						}
					}
					if (src_return1.size() > 0 && !src_return1[0]) {
						if(debug) {
							cout << "inss here-127\n";
						}
						left_bound_sr = cur_end;
					}
					if (src_return1.size() > 3 && !src_return1[3]) {
						if(debug) {
							cout << "inss here-128\n";
						}
						right_bound_sr = cur_mate_end;
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id1;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = mate_ref_id;
					an_entry.event_start = right_bound_sr + 1;
					an_entry.strand = 1;
					an_entry.mate_ref_id = ref_id;
					an_entry.mate_event_start = left_bound_sr + 1;
					an_entry.mate_strand = -1;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tsrc_return1[4]\tdata[7]\tright_bound_sr\t1\tdata[3]\tleft_bound_sr\t-1\n";

				}
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads) {
				if (ref_id < mate_ref_id) {
					if(debug) {
						cout << "inss here-129\n";
					}

					int64_t left_bound_sr = 0;
					if (src_return2.size() > 0) {
						if(debug) {
							cout << "inss here-130\n";
						}
						left_bound_sr = positionb1 + src_return2[0];
					}
					int64_t right_bound_sr = 0;
					if (src_return2.size() > 3) {
						if(debug) {
							cout << "inss here-131\n";
						}
						right_bound_sr = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
					}
					string& local_chr_bp_1 = mate_ref_id;
					string& local_chr_bp_2 = ref_id;

					if (chr_bp_b1 == local_chr_bp_1 && chr_bp_b2 == local_chr_bp_2) {
						if (src_return2.size() > 2) {
							if(debug) {
								cout << "inss here-132\n";
							}
							left_bound_sr = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
						}
						if (src_return2.size() > 1) {
							if(debug) {
								cout << "inss here-133\n";
							}
							right_bound_sr = positionb1 + src_return2[1] + options.cut_sr;
						}
					}
					if (!src_return2.empty() && !src_return2[0]) {
						if(debug) {
							cout << "inss here-134\n";
						}
						left_bound_sr = cur_end;
					}
					if (src_return2.size() > 3 && !src_return2[3]) {
						if(debug) {
							cout << "inss here-135\n";
						}
						right_bound_sr = cur_mate_end;
					}
					string reads = castle::StringUtils::join(bpread_2[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id2;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
					an_entry.n_supports = src_return2[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = -1;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr + 1;
					an_entry.mate_strand = 1;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[3]\tleft_bound_sr\t-1\tdata[7]\tright_bound_sr\t1\n";

				} else {
					if(debug) {
						cout << "inss here-136\n";
					}
					int64_t left_bound_sr = 0;
					if (src_return2.size() > 1) {
						if(debug) {
							cout << "inss here-137\n";
						}
						left_bound_sr = positionb1 + src_return2[1] + options.cut_sr;
					}
					int64_t right_bound_sr = 0;
					if (src_return2.size() > 2) {
						if(debug) {
							cout << "inss here-138\n";
						}
						right_bound_sr = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					string& local_chr_bp_1 = mate_ref_id;
					string& local_chr_bp_2 = ref_id;

					if (chr_bp_b1 == local_chr_bp_1 && chr_bp_b2 == local_chr_bp_2) {
						if (src_return2.size() > 3) {
							if(debug) {
								cout << "inss here-139\n";
							}
							left_bound_sr = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3 + options.cut_sr;
						}
						if (src_return2.size() > 0) {
							if(debug) {
								cout << "inss here-140\n";
							}
							right_bound_sr = positionb1 + src_return2[0];
						}
					}
					if (src_return2.size() > 1 && !src_return2[1]) {
						if(debug) {
														cout << "inss here-141\n";
													}
						left_bound_sr = cur_start;
					}
					if (src_return2.size() > 2 && !src_return2[2]) {
						if(debug) {
							cout << "inss here-142\n";
						}
						right_bound_sr = cur_mate_start;
					}
					string reads = castle::StringUtils::join(bpread_2[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id2;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
					an_entry.n_supports = src_return2[4];
					an_entry.ref_id = mate_ref_id;
					an_entry.event_start = right_bound_sr + 1;
					an_entry.strand = -1;
					an_entry.mate_ref_id = ref_id;
					an_entry.mate_event_start = left_bound_sr + 1;
					an_entry.mate_strand = 1;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tsrc_return2[4]\tdata[7]\tright_bound_sr\t-1\tdata[3]\tleft_bound_sr\t1\n";
				}
			}
		}
	}
}

void SplitReadSVCaller::detect_inso(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {
//	const bool debug = "17814_0/5021_0" == data[1];
	const bool debug = false;
//	auto& ref_is = options.is;
	const char* delim_slash = "/";
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;

	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_end = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_start = boost::lexical_cast<int64_t>(data[8]);
	int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& a_region_b = cluster_region.find(cl[1]);

	string chr_bp_a1 = a_region_a->second.chr_bp_1;
	int64_t positiona1 = a_region_a->second.start;
	int64_t positiona2 = a_region_a->second.end;
	string chr_bp_a2 = a_region_a->second.chr_bp_2;
	int64_t positiona3 = a_region_a->second.mate_start;
	int64_t positiona4 = a_region_a->second.mate_end;
	int32_t orientationa = a_region_a->second.orientation;

	string chr_bp_b1 = a_region_b->second.chr_bp_1;
	int64_t positionb1 = a_region_b->second.start;
	int64_t positionb2 = a_region_b->second.end;
	string chr_bp_b2 = a_region_b->second.chr_bp_2;
	int64_t positionb3 = a_region_b->second.mate_start;
	int64_t positionb4 = a_region_b->second.mate_end;
	int32_t orientationb = a_region_b->second.orientation;

	const string& ref_id = data[3];
	const string& mate_ref_id = data[7];

	const string& local_chr_bp_a1 = mate_ref_id;
	const string& local_chr_bp_a2 = ref_id;

	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);
	const auto& rev_index_second = reverse_index_ref_id.find(a_region_b->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
	int64_t local_ref_id_second = rev_index_second->second.first;
	int64_t local_ref_start_second = 0;
	int64_t local_ref_end_second = rev_index_second->second.second;

	if(debug) {
		cout << "inso: cl: " << cl[0] << "/" << cl[1] <<
		   ", cur_pos: " << cur_start << "/" << cur_end << "/" << cur_mate_start << "/" << cur_mate_end <<
				", chrbp: " << chr_bp_a1 << "/" << chr_bp_a2 << "/" << chr_bp_b1 << "/" << chr_bp_b2 <<
				", position_a: " << positiona1 << "/" << positiona2 << "/" << positiona3 << "/" << positiona4 <<
				", position_b: " << positionb1 << "/" << positionb2 << "/" << positionb3 << "/" << positionb4 <<
				", local_ref:" << local_ref_id << "/" << local_ref_start << "/" << local_ref_end <<
				", local_ref_second:" << local_ref_id_second << "/" << local_ref_start_second << "/" << local_ref_end_second << "\n" ;
	}

	if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
		if (a_region_a->second == a_region_b->second) {
			cerr << "inso i\n";
		}
		//#(positiona1, positiona2, positiona3, positiona4) = (positiona3, positiona4, positiona1, positiona2);
	}
	//#string& local_chr_bp_1 = data[7];
	//string& local_chr_bp_2 = data[3];
	//if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2){
	//#(positionb1, positionb2, positionb3, positionb4) = (positionb3, positionb4, positionb1, positionb2);}//#print "a_region_a->second\na_region_b->second\npositiona1, positiona2, positiona3, positiona4\npositionb1, positionb2, positionb3, positionb4\n";

	if (a_region_a->second == a_region_b->second) {
		if(debug) {
			cout << "inso here-1\n";
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);

		if (orientationa == 1) {
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;
				if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {
					//#print "readname\tstart mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if(debug) {
				cout << "inso here-2: " << discord_sr1.size() << "\n";
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inso here-3\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;
				sr_cluster(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
					int64_t left_bound_sr1 = 0;
					if (ref_boundary1.size() > 1) {
						if(debug) {
							cout << "inso here-4\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}

					int64_t right_bound_sr1 = 0;
					if (ref_boundary2.size() > 0) {
						if(debug) {
							cout << "inso here-5\n";
						}
						right_bound_sr1 = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
					}
					int64_t left_bound_sr2 = 0;
					if (ref_boundary1.size() > 2) {
						if(debug) {
							cout << "inso here-6\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[2];
					}
					int64_t right_bound_sr2 = 0;
					if (ref_boundary2.size() > 3) {
						if(debug) {
							cout << "inso here-7\n";
						}
						right_bound_sr2 = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}
					string& local_chr_bp_a1 = data[7];
					string& local_chr_bp_a2 = data[3];

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if (ref_boundary2.size() > 0) {
							if(debug) {
								cout << "inso here-8\n";
							}
							left_bound_sr1 = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
						}
						if (ref_boundary1.size() > 1) {
							if(debug) {
								cout << "inso here-9\n";
							}
							right_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
						}
						if (ref_boundary2.size() > 3) {
							if(debug) {
								cout << "inso here-10\n";
							}
							left_bound_sr2 = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
						}
						if (ref_boundary1.size() > 2) {
							if(debug) {
								cout << "inso here-11\n";
							}
							right_bound_sr2 = positiona1 + ref_boundary1[2];
						}
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
							cout << "inso here-12\n";
						}
						left_bound_sr1 = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
							cout << "inso here-13\n";
						}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
							cout << "inso here-14\n";
						}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
							cout << "inso here-15\n";
						}
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
						if(debug) {
							cout << "inso here-16\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = ref_id;
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = mate_ref_id;
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[7]\tright_bound_sr2\tright_bound_sr1\tins_size\n";

					} else {
						if(debug) {
							cout << "inso here-17\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
								cout << "inso here-18\n";
							}
							int64_t left_bound_sr = 0;
							if(ref_boundary1.size() > 1) {
								left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if(ref_boundary2.size() > 0) {
								right_bound_sr = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 0) {
									if(debug) {
										cout << "inso here-19\n";
									}
									left_bound_sr = positiona4 - (ref_boundary2[0] - (positiona2 - positiona1 + 101));
								}
								if (ref_boundary1.size() > 1) {
									if(debug) {
										cout << "inso here-20\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
									cout << "inso here-21\n";
								}
								left_bound_sr = cur_start;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
									cout << "inso here-22\n";
								}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
									cout << "inso here-23\n";
								}
								if(ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);

								//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\t1\tdata[7]\tright_bound_sr\t1\n";
								}
							} else {
								if(debug) {
									cout << "inso here-24\n";
								}
								if(ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = 1;
								}
								//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[7]\tright_bound_sr\t1\tdata[3]\tleft_bound_sr\t1\n";
							}
						}
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
								cout << "inso here-25\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 2) {
								if(debug) {
									cout << "inso here-26\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 3) {
								if(debug) {
									cout << "inso here-27\n";
								}
								right_bound_sr = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 3) {
									if(debug) {
										cout << "inso here-28\n";
									}
									left_bound_sr = positiona4 - (ref_boundary2[3] - (positiona2 - positiona1 + 101) + options.cut_sr);
								}
								if (ref_boundary1.size() > 2) {
									if(debug) {
										cout << "inso here-29\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[2];
								}
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
									cout << "inso here-30\n";
								}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
									cout << "inso here-31\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
									cout << "inso here-32\n";
								}
								if(ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[3]\tleft_bound_sr\t-1\tdata[7]\tright_bound_sr\t-1\n";
								}
							} else {
								if(debug) {
									cout << "inso here-33\n";
								}
								if(ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tref_support_sr[1]\tdata[7]\tright_bound_sr\t-1\tdata[3]\tleft_bound_sr\t-1\n";
								}
							}
						}
					}
				} else {
					if(debug) {
					cout << "inso here-34\n";
					}
					int64_t left_bound_sr1 = 0;
					if (ref_boundary1.size() > 3) {
						if(debug) {
						cout << "inso here-35\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}
					int64_t right_bound_sr1 = 0;
					if (ref_boundary2.size() > 2) {
						if(debug) {
						cout << "inso here-36\n";
						}
						right_bound_sr1 = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
					}
					int64_t left_bound_sr2 = 0;
					if (ref_boundary1.size() > 0) {
						if(debug) {
						cout << "inso here-37\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[0];
					}
					int64_t right_bound_sr2 = 0;
					if (ref_boundary2.size() > 1) {
						if(debug) {
						cout << "inso here-38\n";
						}
						right_bound_sr2 = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
					}
					const string& local_chr_bp_a1 = mate_ref_id;
					const string& local_chr_bp_a2 = ref_id;

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if (ref_boundary2.size() > 2) {
							if(debug) {
							cout << "inso here-39\n";
							}
							left_bound_sr1 = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
						}
						if (ref_boundary1.size() > 3) {
							if(debug) {
							cout << "inso here-40\n";
							}
							right_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
						}
						if (ref_boundary2.size() > 1) {
							if(debug) {
							cout << "inso here-41\n";
							}
							left_bound_sr2 = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
						}
						if (ref_boundary1.size() > 0) {
							if(debug) {
							cout << "inso here-42\n";
							}
							right_bound_sr2 = positiona1 + ref_boundary1[0];
						}
					}
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
						cout << "inso here-43\n";
						}
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
						cout << "inso here-44\n";
						}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
						cout << "inso here-45\n";
						}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
						cout << "inso here-46\n";
						}
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads && left_bound_sr1 && right_bound_sr1) {
						if(debug) {
						cout << "inso here-47\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = ref_id;
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = mate_ref_id;
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
					} else {
						if(debug) {
						cout << "inso here-48\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
							cout << "inso here-49\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 3) {
								if(debug) {
								cout << "inso here-50\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 2) {
								if(debug) {
								cout << "inso here-51\n";
								}
								right_bound_sr = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if(debug) {
								cout << "inso here-52\n";
								}
								if (ref_boundary2.size() > 2) {
									if(debug) {
									cout << "inso here-53\n";
									}
									left_bound_sr = positiona4 - (ref_boundary2[2] - (positiona2 - positiona1 + 101));
								}
								if (ref_boundary1.size() > 3) {
									if(debug) {
									cout << "inso here-54\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
								cout << "inso here-55\n";
								}
								left_bound_sr = cur_start;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
								cout << "inso here-56\n";
								}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
								cout << "inso here-57\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[3]\tleft_bound_sr\t1\tdata[7]\tright_bound_sr\t1\n";
								}
							} else {
								if(debug) {
								cout << "inso here-58\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[1]\tdata[7]\tright_bound_sr\t1\tdata[3]\tleft_bound_sr\t1\n";
								}
							}
						}
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
							cout << "inso here-59\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 0) {
								if(debug) {
								cout << "inso here-60\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 1) {
								if(debug) {
								cout << "inso here-61\n";
								}
								right_bound_sr = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
							}
							const string& local_chr_bp_a1 = mate_ref_id;
							const string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if(debug) {
								cout << "inso here-62\n";
								}
								if (ref_boundary2.size() > 1) {
									if(debug) {
									cout << "inso here-63\n";
									}
									left_bound_sr = positiona4 - (ref_boundary2[1] - (positiona2 - positiona1 + 101) + options.cut_sr);
								}
								if (ref_boundary1.size() > 0) {
									if(debug) {
									cout << "inso here-64\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[0];
								}
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
								cout << "inso here-65\n";
								}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
								cout << "inso here-66\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
								cout << "inso here-67\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr - 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[3]\tleft_bound_sr\t-1\tdata[7]\tright_bound_sr\t-1\n";
								}
							} else {
								if(debug) {
								cout << "inso here-68\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr - 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id2\tmpd2\tref_support_sr[0]\tdata[7]\tright_bound_sr\t-1\tdata[3]\tleft_bound_sr\t-1\n";
								}

							}
						}
					}
				}
			}
		} else {
			if(debug) {
			cout << "inso here-69\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "readname\tstart mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
				cout << "inso here-70\n";
				}
				vector<int64_t> ref_boundary1;
				vector<int64_t> ref_boundary2;
				vector<int64_t> ref_support_sr;
				map<int64_t, vector<string>> ref_bpread;

				sr_cluster_1(ref_boundary1, ref_boundary2, ref_support_sr, ref_bpread, discord_sr1, 1);

				//#print "@ref_boundary1\n@ref_boundary2\n@ref_support_sr\n";
				if (ref_boundary1.size() > 2 && ref_boundary1[2] > ref_boundary1[0]) {
					if(debug) {
					cout << "inso here-71\n";
					}
					int64_t left_bound_sr1 = 0;
					if (ref_boundary1.size() > 1) {
						if(debug) {
						cout << "inso here-72\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[1] + options.cut_sr;
					}
					int64_t right_bound_sr1 = 0;
					if (ref_boundary2.size() > 0) {
						if(debug) {
						cout << "inso here-73\n";
						}
						right_bound_sr1 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
					}
					int64_t left_bound_sr2 = 0;
					if (ref_boundary1.size() > 2) {
						if(debug) {
						cout << "inso here-74\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[2];
					}
					int64_t right_bound_sr2 = 0;
					if (ref_boundary2.size() > 3) {
						if(debug) {
						cout << "inso here-75\n";
						}
						right_bound_sr2 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					string& local_chr_bp_a1 = data[7];
					string& local_chr_bp_a2 = data[3];

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if(debug) {
						cout << "inso here-76\n";
						}
						if (ref_boundary2.size() > 1) {
							if(debug) {
							cout << "inso here-77\n";
							}
							left_bound_sr1 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						if (ref_boundary1.size() > 0) {
							if(debug) {
							cout << "inso here-78\n";
							}
							right_bound_sr1 = positiona1 + ref_boundary1[0];
						}
						if (ref_boundary2.size() > 2) {
							if(debug) {
							cout << "inso here-79\n";
							}
							left_bound_sr2 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
						}
						if (ref_boundary1.size() > 3) {
							if(debug) {
							cout << "inso here-80\n";
							}
							right_bound_sr2 = positiona1 + ref_boundary1[3] + options.cut_sr;
						}
					}
					if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
						if(debug) {
						cout << "inso here-81\n";
						}
						left_bound_sr1 = cur_start;
					}
					if (!ref_boundary2.empty() && !ref_boundary2[0]) {
						if(debug) {
						cout << "inso here-82\n";
						}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
						if(debug) {
						cout << "inso here-83\n";
						}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
						if(debug) {
						cout << "inso here-84\n";
						}
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
						if(debug) {
						cout << "inso here-85\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[0], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[1], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[0];
						an_entry.n_mate_support = ref_support_sr[1];
						an_entry.ref_id = ref_id;
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = mate_ref_id;
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[0]/ref_support_sr[1]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[7]\tright_bound_sr2\tright_bound_sr1\tins_size\n";

					} else {
						if(debug) {
						cout << "inso here-86\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
							cout << "inso here-87\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 1) {
								if(debug) {
								cout << "inso here-88\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 0) {
								if(debug) {
								cout << "inso here-89\n";
								}
								right_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if(debug) {
								cout << "inso here-90\n";
								}
								if (ref_boundary2.size() > 1) {
									if(debug) {
									cout << "inso here-91\n";
									}
									left_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 0) {
									if(debug) {
									cout << "inso here-92\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[0];
								}
							}
							if (ref_boundary1.size() > 1 && !ref_boundary1[1]) {
								if(debug) {
								cout << "inso here-93\n";
								}
								left_bound_sr = cur_start;
							}
							if (!ref_boundary2.empty() && !ref_boundary2[0]) {
								if(debug) {
								cout << "inso here-94\n";
								}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (data[3] < data[7]) {
								if(debug) {
								cout << "inso here-95\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[3]\tleft_bound_sr\t1\tdata[7]\tright_bound_sr\t1\n";
								}
							} else {
								if(debug) {
								cout << "inso here-96\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
									//result_sr[sri][0] = "transl_inter\tcluster_id1\tmpd1\tref_support_sr[0]\tdata[7]\tright_bound_sr\t1\tdata[3]\tleft_bound_sr\t1\n";
								}
							}
						}
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
							cout << "inso here-97\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 2) {
								if(debug) {
								cout << "inso here-98\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[2];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 3) {
								if(debug) {
								cout << "inso here-99\n";
								}
								right_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 2) {
									if(debug) {
									cout << "inso here-100\n";
									}
									left_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 3) {
									if(debug) {
									cout << "inso here-101\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 2 && !ref_boundary1[2]) {
								if(debug) {
								cout << "inso here-102\n";
								}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 3 && !ref_boundary2[3]) {
								if(debug) {
								cout << "inso here-103\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
								cout << "inso here-104\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
								}
							} else {
								if(debug) {
								cout << "inso here-105\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
								}
							}
						}
					}
				} else {
					if(debug) {
					cout << "inso here-106\n";
					}
					int64_t left_bound_sr1 = 0;
					if (ref_boundary1.size() > 3) {
						if(debug) {
						cout << "inso here-107\n";
						}
						left_bound_sr1 = positiona1 + ref_boundary1[3] + options.cut_sr;
					}
					int64_t right_bound_sr1 = 0;
					if (ref_boundary2.size() > 2) {
						if(debug) {
						cout << "inso here-108\n";
						}
						right_bound_sr1 = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
					}
					int64_t left_bound_sr2 = 0;
					if (ref_boundary1.size() > 0) {
						if(debug) {
						cout << "inso here-109\n";
						}
						left_bound_sr2 = positiona1 + ref_boundary1[0];
					}
					int64_t right_bound_sr2 = 0;
					if (ref_boundary2.size() > 1) {
						if(debug) {
						cout << "inso here-110\n";
						}
						right_bound_sr2 = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					string& local_chr_bp_a1 = data[7];
					string& local_chr_bp_a2 = data[3];

					if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
						if (ref_boundary2.size() > 3) {
							if(debug) {
							cout << "inso here-111\n";
							}
							left_bound_sr1 = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
						}
						if (ref_boundary1.size() > 2) {
							if(debug) {
							cout << "inso here-112\n";
							}
							right_bound_sr1 = positiona1 + ref_boundary1[2];
						}
						if (ref_boundary2.size() > 0) {
							if(debug) {
							cout << "inso here-113\n";
							}
							left_bound_sr2 = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
						}
						if (ref_boundary1.size() > 1) {
							if(debug) {
							cout << "inso here-114\n";
							}
							right_bound_sr2 = positiona1 + ref_boundary1[1] + options.cut_sr;
						}
					}
					if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
						if(debug) {
						cout << "inso here-115\n";
						}
						left_bound_sr1 = cur_start;
					}
					if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
						if(debug) {
						cout << "inso here-116\n";
						}
						right_bound_sr1 = cur_mate_end;
					}
					if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
						if(debug) {
						cout << "inso here-117\n";
						}
						left_bound_sr2 = cur_end;
					}
					if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
						if(debug) {
						cout << "inso here-118\n";
						}
						right_bound_sr2 = cur_mate_start;
					}
					int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
					int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
					if (ref_support_sr.size() > 1 && ref_support_sr[0] >= options.support_reads && ref_support_sr[1] >= options.support_reads) {
						if(debug) {
						cout << "inso here-119\n";
						}
						string reads0 = castle::StringUtils::join(ref_bpread[1], "\t");
						string reads1 = castle::StringUtils::join(ref_bpread[0], "\t");
						BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
						EventEntry an_entry;
						an_entry.type = data[0];
						an_entry.cluster_id = data[1];
						an_entry.mate_cluster_id = data[2];
						an_entry.n_supports = ref_support_sr[1];
						an_entry.n_mate_support = ref_support_sr[0];
						an_entry.ref_id = ref_id;
						an_entry.event_start = left_bound_sr1 + 1;
						an_entry.event_end = left_bound_sr2 + 1;
						an_entry.event_size_1 = del_size;
						an_entry.mate_ref_id = mate_ref_id;
						an_entry.mate_event_start = right_bound_sr2 + 1;
						an_entry.mate_event_end = right_bound_sr1 + 1;
						an_entry.event_size_2 = ins_size;
						result_sr.push_back(an_entry);
						//result_sr[sri][0] = data[0];
						//result_sr[sri][1] = data[1];
						//result_sr[sri][2] = data[2];
						//result_sr[sri][3] = "ref_support_sr[1]/ref_support_sr[0]";
						//result_sr[sri][4] = data[3];
						//result_sr[sri][5] = left_bound_sr1;
						//result_sr[sri][6] = left_bound_sr2;
						//result_sr[sri][7] = del_size;
						//result_sr[sri][8] = "data[7]\tright_bound_sr2\tright_bound_sr1\tins_size\n";

					} else {
						if(debug) {
						cout << "inso here-120\n";
						}
						castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
						string& cluster_id1 = cluster_ids[0];
						string& cluster_id2 = cluster_ids[1];
						castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
						string& mpd_str_1 = mpds[0];
						string& mpd_str_2 = mpds[1];
						int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
						int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
						if (ref_support_sr.size() > 1 && ref_support_sr[1] >= options.support_reads) {
							if(debug) {
							cout << "inso here-121\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 3) {
								if(debug) {
								cout << "inso here-122\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[3] + options.cut_sr;
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 2) {
								if(debug) {
								cout << "inso here-123\n";
								}
								right_bound_sr = ref_boundary2[2] - (positiona2 - positiona1 + 101) + positiona3;
							}
							string& local_chr_bp_a1 = data[7];
							string& local_chr_bp_a2 = data[3];

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if(debug) {
								cout << "inso here-124\n";
								}
								if (ref_boundary2.size() > 3) {
									if(debug) {
									cout << "inso here-125\n";
									}

									left_bound_sr = ref_boundary2[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
								}
								if (ref_boundary1.size() > 2) {
									if(debug) {
									cout << "inso here-126\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[2];
								}
							}
							if (ref_boundary1.size() > 3 && !ref_boundary1[3]) {
								if(debug) {
								cout << "inso here-127\n";
								}
								left_bound_sr = cur_start;
							}
							if (ref_boundary2.size() > 2 && !ref_boundary2[2]) {
								if(debug) {
								cout << "inso here-128\n";
								}
								right_bound_sr = cur_mate_end;
							}
							string reads = castle::StringUtils::join(ref_bpread[1], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
								cout << "inso here-129\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
								}

							} else {
								if(debug) {
								cout << "inso here-130\n";
								}
								if (ref_support_sr.size() > 1) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id1;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
									an_entry.n_supports = ref_support_sr[1];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = 1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = 1;
									result_sr.push_back(an_entry);
								}

							}
						}
						if (ref_support_sr.size() > 0 && ref_support_sr[0] >= options.support_reads) {
							if(debug) {
							cout << "inso here-131\n";
							}
							int64_t left_bound_sr = 0;
							if (ref_boundary1.size() > 0) {
								if(debug) {
								cout << "inso here-132\n";
								}
								left_bound_sr = positiona1 + ref_boundary1[0];
							}
							int64_t right_bound_sr = 0;
							if (ref_boundary2.size() > 1) {
								if(debug) {
								cout << "inso here-133\n";
								}
								right_bound_sr = ref_boundary2[1] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
							}
							const string& local_chr_bp_a1 = mate_ref_id;
							const string& local_chr_bp_a2 = ref_id;

							if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
								if (ref_boundary2.size() > 0) {
									if(debug) {
									cout << "inso here-134\n";
									}
									left_bound_sr = ref_boundary2[0] - (positiona2 - positiona1 + 101) + positiona3;
								}
								if (ref_boundary1.size() > 1) {
									if(debug) {
									cout << "inso here-135\n";
									}
									right_bound_sr = positiona1 + ref_boundary1[1] + options.cut_sr;
								}
							}
							if (ref_boundary1.size() > 0 && !ref_boundary1[0]) {
								if(debug) {
								cout << "inso here-136\n";
								}
								left_bound_sr = cur_end;
							}
							if (ref_boundary2.size() > 1 && !ref_boundary2[1]) {
								if(debug) {
								cout << "inso here-137\n";
								}
								right_bound_sr = cur_mate_start;
							}
							string reads = castle::StringUtils::join(ref_bpread[0], "\t");
							if (ref_id < mate_ref_id) {
								if(debug) {
								cout << "inso here-138\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = ref_id;
									an_entry.event_start = left_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = mate_ref_id;
									an_entry.mate_event_start = right_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
								}
							} else {
								if(debug) {
								cout << "inso here-139\n";
								}
								if (ref_support_sr.size() > 0) {
									BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr + 1) % ref_id % (left_bound_sr + 1) % reads).str();
									EventEntry an_entry;
									an_entry.type = "transl_inter";
									an_entry.cluster_id = cluster_id2;
									an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
									an_entry.n_supports = ref_support_sr[0];
									an_entry.ref_id = mate_ref_id;
									an_entry.event_start = right_bound_sr + 1;
									an_entry.strand = -1;
									an_entry.mate_ref_id = ref_id;
									an_entry.mate_event_start = left_bound_sr + 1;
									an_entry.mate_strand = -1;
									result_sr.push_back(an_entry);
								}
							}
						}
					}
				}
			}
		}
	} else {
		if(debug) {
		cout << "inso here-140\n";
		}
		int64_t left_bound_sr1 = 0;
		int64_t right_bound_sr1 = 0;
		int64_t left_bound_sr2 = 0;
		int64_t right_bound_sr2 = 0;
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (orientationa == 1) {
			if(debug) {
			cout << "inso here-141\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
				cout << "inso here-142\n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);

				//#print "src_return1\n";
				if (src_return1.size() > 1) {
					if(debug) {
					cout << "inso here-143\n";
					}
					left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
				}
				if (src_return1.size() > 2) {
					if(debug) {
					cout << "inso here-144\n";
					}
					right_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
				}
				const string& local_chr_bp_a1 = mate_ref_id;
				const string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if(debug) {
					cout << "inso here-145\n";
					}
					if (src_return1.size() > 2) {
						if(debug) {
						cout << "inso here-146\n";
						}
						left_bound_sr1 = positiona4 - (src_return1[2] - (positiona2 - positiona1 + 101));
					}
					if (src_return1.size() > 1) {
						if(debug) {
						cout << "inso here-147\n";
						}
						right_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
					}
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
					cout << "inso here-148\n";
					}
					left_bound_sr1 = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
					cout << "inso here-149\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
			}
		} else {
			if(debug) {
			cout << "inso here-150\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (positiona2 - positiona1) && mstart > (positiona2 - positiona1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
				cout << "inso here-151\n";
				}
				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);

				//#print "src_return1\n";
				if (src_return1.size() > 1) {
					if(debug) {
					cout << "inso here-152\n";
					}
					left_bound_sr1 = positiona1 + src_return1[1] + options.cut_sr;
				}
				if (src_return1.size() > 2) {
					if(debug) {
					cout << "inso here-153\n";
					}
					right_bound_sr1 = src_return1[2] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
				}
				const string& local_chr_bp_a1 = mate_ref_id;
				const string& local_chr_bp_a2 = ref_id;

				if (chr_bp_a1 == local_chr_bp_a1 && chr_bp_a2 == local_chr_bp_a2) {
					if(debug) {
					cout << "inso here-154\n";
					}
					if (src_return1.size() > 3) {
						if(debug) {
						cout << "inso here-155\n";
						}
						left_bound_sr1 = src_return1[3] - (positiona2 - positiona1 + 101) + positiona3 + options.cut_sr;
					}
					if (src_return1.size() > 0) {
						if(debug) {
						cout << "inso here-156\n";
						}
						right_bound_sr1 = positiona1 + src_return1[0] + options.cut_sr;
					}
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
					cout << "inso here-157\n";
					}
					left_bound_sr1 = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "inso here-158\n";
					}
					right_bound_sr1 = cur_mate_end;
				}
			}
		}

		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_b->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		// cl[1] parts
		reader.SetRegion(local_ref_id_second, local_ref_start_second, local_ref_id_second, local_ref_end_second);
		if (orientationb == 1) {
			if(debug) {
				cout << "inso here-159\n";
			}
			discord_sr2.clear();

			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inso here-160\n";
				}
				sr_cluster_no_cl(src_return2, bpread_2, discord_sr2);
				if (!src_return2.empty()) {
					if(debug) {
						cout << "inso here-161\n";
					}
					left_bound_sr2 = positionb1 + src_return2[0];
				}
				if (src_return2.size() > 3) {
					if(debug) {
						cout << "inso here-162\n";
					}
					right_bound_sr2 = positionb4 - (src_return2[3] - (positionb2 - positionb1 + 101) + options.cut_sr);
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_b1 == local_chr_bp_1 && chr_bp_b2 == local_chr_bp_2) {
					if(debug) {
						cout << "inso here-163\n";
					}
					if (src_return2.size() > 3) {
						if(debug) {
							cout << "inso here-164\n";
						}
						left_bound_sr2 = positionb4 - (src_return2[3] - (positionb2 - positionb1 + 101) + options.cut_sr);
					}
					if (src_return2.size() > 0) {
						if(debug) {
							cout << "inso here-165\n";
						}
						right_bound_sr2 = positionb1 + src_return2[0];
					}
				}
				if (!src_return2.empty() && !src_return2[0]) {
					if(debug) {
						cout << "inso here-166\n";
					}
					left_bound_sr2 = cur_end;
				}
				if (src_return2.size() > 3 && !src_return2[3]) {
					if(debug) {
						cout << "inso here-167\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
			}
		} else {
			if(debug) {
				cout << "inso here-168\n";
			}
			discord_sr2.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (positionb2 - positionb1) && mstart > (positionb2 - positionb1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr2.push_back(al);
				}
			}
			if (discord_sr2.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "inso here-169\n";
				}
				sr_cluster_1_no_cl(src_return2, bpread_2, discord_sr2);
				if (!src_return2.empty()) {
					if(debug) {
						cout << "inso here-170\n";
					}
					left_bound_sr2 = positionb1 + src_return2[0];
				}
				if (src_return2.size() > 3) {
					if(debug) {
						cout << "inso here-171\n";
					}
					right_bound_sr2 = src_return2[3] - (positionb2 - positionb1 + 101) + positionb3;
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_b1 == local_chr_bp_1 && chr_bp_b2 == local_chr_bp_2) {
					cout << "inso here-172\n";
					if (src_return2.size() > 2) {
						if(debug) {
							cout << "inso here-173\n";
						}
						left_bound_sr2 = src_return2[2] - (positionb2 - positionb1 + 101) + positionb3;
					}
					if (src_return2.size() > 1) {
						if(debug) {
							cout << "inso here-174\n";
						}
						right_bound_sr2 = positionb1 + src_return2[1];
					}
				}
				if (!src_return2.empty() && !src_return2[0]) {
					if(debug) {
						cout << "inso here-175\n";
					}
					left_bound_sr2 = cur_end;
				}
				if (src_return2.size() > 3 && !src_return2[3]) {
					if(debug) {
						cout << "inso here-176\n";
					}
					right_bound_sr2 = cur_mate_start;
				}
			}
		}

		int64_t del_size = left_bound_sr2 - left_bound_sr1 - 1;
		int64_t ins_size = right_bound_sr1 - right_bound_sr2 + 1;
		if (src_return1.size() > 4 && src_return2.size() > 4 && src_return1[4] >= options.support_reads && src_return2[4] >= options.support_reads) {
			if(debug) {
				cout << "inso here-177\n";
			}
			string reads0 = castle::StringUtils::join(bpread_1[0], "\t");
			string reads1 = castle::StringUtils::join(bpread_2[0], "\t");
			BPREAD << (boost::format("%s__%s__%s__%s\n%s\n%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr1 + 1) % mate_ref_id % (right_bound_sr1 + 1) % reads0 % ref_id % (left_bound_sr2 + 1) % mate_ref_id % (right_bound_sr2 + 1) % reads1).str();
			EventEntry an_entry;
			an_entry.type = data[0];
			an_entry.cluster_id = data[1];
			an_entry.mate_cluster_id = data[2];
			an_entry.n_supports = src_return1[4];
			an_entry.n_mate_support = src_return2[4];
			an_entry.ref_id = ref_id;
			an_entry.event_start = left_bound_sr1 + 1;
			an_entry.event_end = left_bound_sr2 + 1;
			an_entry.event_size_1 = del_size;
			an_entry.mate_ref_id = mate_ref_id;
			an_entry.mate_event_start = right_bound_sr2 + 1;
			an_entry.mate_event_end = right_bound_sr1 + 1;
			an_entry.event_size_2 = ins_size;
			result_sr.push_back(an_entry);
			//result_sr[sri][0] = data[0];
			//result_sr[sri][1] = data[1];
			//result_sr[sri][2] = data[2];
			//result_sr[sri][3] = "src_return1[4]/src_return2[4]";
			//result_sr[sri][4] = data[3];
			//result_sr[sri][5] = left_bound_sr1;
			//result_sr[sri][6] = left_bound_sr2;
			//result_sr[sri][7] = del_size;
			//result_sr[sri][8] = "data[7]\tright_bound_sr2\tright_bound_sr1\tins_size\n";
		} else {
			if(debug) {
				cout << "inso here-178\n";
			}
			castle::StringUtils::c_string_multi_split(data[1], delim_slash, cluster_ids);
			string& cluster_id1 = cluster_ids[0];
			string& cluster_id2 = cluster_ids[1];
			castle::StringUtils::c_string_multi_split(data[2], delim_slash, mpds);
			string& mpd_str_1 = mpds[0];
			string& mpd_str_2 = mpds[1];
			int64_t mpd1 = boost::lexical_cast<int64_t>(mpd_str_1);
			int64_t mpd2 = boost::lexical_cast<int64_t>(mpd_str_2);
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				if(debug) {
					cout << "inso here-179\n";
				}
				int64_t left_bound_sr = left_bound_sr1;
				int64_t right_bound_sr = right_bound_sr1;
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				if (ref_id < mate_ref_id) {
					if(debug) {
						cout << "inso here-189\n";
					}
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr - 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id1;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = 1;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr - 1;
					an_entry.mate_strand = 1;
					result_sr.push_back(an_entry);

				} else {
					if(debug) {
						cout << "inso here-190\n";
					}
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr - 1) % ref_id % (left_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id1;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd1);
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = mate_ref_id;
					an_entry.event_start = right_bound_sr - 1;
					an_entry.strand = 1;
					an_entry.mate_ref_id = ref_id;
					an_entry.mate_event_start = left_bound_sr + 1;
					an_entry.mate_strand = 1;
					result_sr.push_back(an_entry);

				}
			}
			if (src_return2.size() > 4 && src_return2[4] >= options.support_reads) {
				if(debug) {
					cout << "inso here-191\n";
				}
				int64_t left_bound_sr = left_bound_sr2;
				int64_t right_bound_sr = right_bound_sr2;
				string reads = castle::StringUtils::join(bpread_2[0], "\t");
				if (ref_id < mate_ref_id) {
					if(debug) {
						cout << "inso here-192\n";
					}
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr - 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id2;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
					an_entry.n_supports = src_return2[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = -1;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr - 1;
					an_entry.mate_strand = -1;
					result_sr.push_back(an_entry);
				} else {
					if(debug) {
						cout << "inso here-193\n";
					}
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % mate_ref_id % (right_bound_sr - 1) % ref_id % (left_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = "transl_inter";
					an_entry.cluster_id = cluster_id2;
					an_entry.mate_cluster_id = boost::lexical_cast<string>(mpd2);
					an_entry.n_supports = src_return2[4];
					an_entry.ref_id = mate_ref_id;
					an_entry.event_start = right_bound_sr - 1;
					an_entry.strand = -1;
					an_entry.mate_ref_id = ref_id;
					an_entry.mate_event_start = left_bound_sr + 1;
					an_entry.mate_strand = -1;
					result_sr.push_back(an_entry);
				}
			}
		}
	}

}

void SplitReadSVCaller::detect_transl_inter(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window, vector<string>& data,
		vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds) {

//	auto& ref_is = options.is;
	vector<BamTools::BamAlignment> discord_sr1;
	vector<BamTools::BamAlignment> discord_sr2;
	vector<int64_t> src_return1;
	vector<int64_t> src_return2;
	map<int64_t, vector<string>> bpread_1;
	map<int64_t, vector<string>> bpread_2;
	int64_t cur_start = boost::lexical_cast<int64_t>(data[4]);
	int64_t cur_strand = boost::lexical_cast<int64_t>(data[5]);
	int64_t cur_mate_strand = boost::lexical_cast<int64_t>(data[8]);
	//int64_t cur_mate_end = boost::lexical_cast<int64_t>(data[9]);
	const string& ref_id = data[3];
	const string& mate_ref_id = data[6];

	const auto& a_region_a = cluster_region.find(cl[0]);
	const auto& rev_index = reverse_index_ref_id.find(a_region_a->second.str);

	int64_t local_ref_id = rev_index->second.first;
	int64_t local_ref_start = 0;
	int64_t local_ref_end = rev_index->second.second;
//	int64_t local_ref_end =
//	int64_t local_ref_id_second = reverse_index_ref_id[a_region_b->second.str].first;
//	int64_t local_ref_start_second = 0;
//	int64_t local_ref_end_second = reverse_index_ref_id[a_region_b->second.str].second;

//	const bool debug = "18846_0" == data[1];
	const bool debug = false;
	if(debug) {
		cout << "transl_inter: cl: " << cl[0] << ", cur_start: " << cur_start << ", strand: " << cur_strand << "/" << cur_mate_strand <<
				", local_ref:" << local_ref_id << "/" << local_ref_start << "/" << local_ref_end << "\n";
	}


	if (cur_strand == 1 && cur_mate_strand == -1) {
		if(debug) {
			cout << "transl_inter here-0\n";
		}
		string chr_bp_1 = a_region_a->second.chr_bp_1;
		int64_t position1 = a_region_a->second.start;
		int64_t position2 = a_region_a->second.end;
		string chr_bp_2 = a_region_a->second.chr_bp_2;
		int64_t position3 = a_region_a->second.mate_start;
		int64_t position4 = a_region_a->second.mate_end;
		const string& local_chr_bp_1 = mate_ref_id;
		const string& local_chr_bp_2 = ref_id;

		if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
			swap(position1, position3);
			swap(position2, position4);
			//(position1, position2, position3, position4) = (position3, position4, position1, position2);
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		discord_sr1.clear();

//		int64_t tmp_local_ref_end = position2 - position1 + ref_is["rlu"]["selected"];
//		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, tmp_local_ref_end);
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);

		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		int64_t n_aln = 0;
		while (reader.GetNextAlignmentBasic(al)) {
			++n_aln;
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

				//#print "start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if(debug) {
			cout << "transl_inter here-1: " << discord_sr1.size() << "\n";
		}

		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
			int64_t left_bound_sr = 0;
			if (src_return1.size() > 1) {
				if(debug) {
					cout << "transl_inter here-2\n";
				}
				left_bound_sr = position1 + src_return1[1] + options.cut_sr;
			}
			int64_t right_bound_sr = 0;
			if (src_return1.size() > 2) {
				if(debug) {
					cout << "transl_inter here-3\n";
				}
				right_bound_sr = src_return1[2] - (position2 - position1 + 101) + position3;
			}
			const string& local_chr_bp_1 = mate_ref_id;
			const string& local_chr_bp_2 = ref_id;

			if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
				if(debug) {
					cout << "transl_inter here-4\n";
				}
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "transl_inter here-5\n";
					}
					left_bound_sr = position1 + src_return1[1];
				}
				if (src_return1.size() > 2) {
					if(debug) {
						cout << "transl_inter here-6\n";
					}
					right_bound_sr = src_return1[2] - (position2 - position1 + 101) + position3 + options.cut_sr;
				}
			}
			if (src_return1.size() > 1 && !src_return1[1]) {
				if(debug) {
					cout << "transl_inter here-7\n";
				}
				left_bound_sr = cur_start;
			}
			if (src_return1.size() > 2 && !src_return1[2]) {
				if(debug) {
					cout << "transl_inter here-8\n";
				}
				right_bound_sr = boost::lexical_cast<int64_t>(data[8]);
			}
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				if(debug) {
					cout << "transl_inter here-9\n";
				}
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.n_supports = src_return1[4];
				an_entry.ref_id = ref_id;
				an_entry.event_start = left_bound_sr + 1;
				an_entry.strand = cur_strand;
				an_entry.mate_ref_id = mate_ref_id;
				an_entry.mate_event_start = right_bound_sr + 1;
				an_entry.mate_strand = cur_mate_strand;
				result_sr.push_back(an_entry);
				////EventEntry an_entry;
				//an_entry.type = data[0];
				//an_entry.the_id_i = data[1];
				//an_entry.the_id_k = data[2];
				//an_entry.support_i = src_return[4];
				//an_entry.seqid_i = data[3];
				//an_entry.starts_map_start = left_bound_sr;
				//an_entry.starts_map_end = cur_end;
				//an_entry.mseqid_i = data[6];
				//an_entry.mstarts_map_start = right_bound_sr;
				//an_entry.mstarts_map_end = cur_mate_start;
				//result_sr.push_back(an_entry);
				//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";

			}
		}
	} else if (cur_strand == -1 && cur_mate_strand == 1) {
		if(debug) {
			cout << "transl_inter here-10\n";
		}
		string chr_bp_1 = a_region_a->second.chr_bp_1;
		int64_t position1 = a_region_a->second.start;
		int64_t position2 = a_region_a->second.end;
		string chr_bp_2 = a_region_a->second.chr_bp_2;
		int64_t position3 = a_region_a->second.mate_start;
		int64_t position4 = a_region_a->second.mate_end;
		//int32_t orientation = a_region_a->second.orientation;
		const string& local_chr_bp_1 = mate_ref_id;
		const string& local_chr_bp_2 = ref_id;

		if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
			swap(position1, position3);
			swap(position2, position4);
			//(position1, position2, position3, position4) = (position3, position4, position1, position2);
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;

		discord_sr1.clear();
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		BamTools::BamAlignment al;
		//int64_t prev_bam_pos = backend_bgzf.Tell();
		while (reader.GetNextAlignmentBasic(al)) {
			int64_t start = al.Position;
			int32_t strand = 1;
			if (al.IsReverseStrand()) {
				strand = -1;
			}
//	string mseqid = ref_vec[al.MateRefID].RefName;
			int64_t mstart = al.MatePosition;
			int32_t mstrand = 1;
			if (al.IsMateReverseStrand()) {
				mstrand = -1;
			}
			int64_t isize = al.InsertSize;

			if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

				//#print "start mstart\tisize\tstrand mstrand\n";
				discord_sr1.push_back(al);
			}
		}
		if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
			if(debug) {
				cout << "transl_inter here-11\n";
			}
//	vector<int64_t> src_return;
//	map<int64_t, vector<string>> bpread;
			sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
			int64_t left_bound_sr = 0;
			if (src_return1.size() > 0) {
				if(debug) {
					cout << "transl_inter here-12\n";
				}
				left_bound_sr = position1 + src_return1[0];
			}
			int64_t right_bound_sr = 0;
			if (src_return1.size() > 3) {
				if(debug) {
					cout << "transl_inter here-13\n";
				}
				right_bound_sr = src_return1[3] - (position2 - position1 + 101) + options.cut_sr + position3;
			}
			const string& local_chr_bp_1 = mate_ref_id;
			const string& local_chr_bp_2 = ref_id;

			if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
				if(debug) {
					cout << "transl_inter here-14\n";
				}
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "transl_inter here-15\n";
					}
					left_bound_sr = position1 + src_return1[0] + options.cut_sr;
				}
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "transl_inter here-16\n";
					}
					right_bound_sr = src_return1[3] - (position2 - position1 + 101) + position3;
				}
			}
			if (src_return1.size() > 0 && !src_return1[0]) {
				if(debug) {
					cout << "transl_inter here-17\n";
				}
				left_bound_sr = cur_start;
			}
			if (src_return1.size() > 3 && !src_return1[3]) {
				if(debug) {
					cout << "transl_inter here-18\n";
				}
				right_bound_sr = boost::lexical_cast<int64_t>(data[7]);
			}
			if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
				if(debug) {
					cout << "transl_inter here-19\n";
				}
				string reads = castle::StringUtils::join(bpread_1[0], "\t");
				BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
				EventEntry an_entry;
				an_entry.type = data[0];
				an_entry.cluster_id = data[1];
				an_entry.mate_cluster_id = data[2];
				an_entry.n_supports = src_return1[4];
				an_entry.ref_id = ref_id;
				an_entry.event_start = left_bound_sr + 1;
				an_entry.strand = cur_strand;
				an_entry.mate_ref_id = mate_ref_id;
				an_entry.mate_event_start = right_bound_sr + 1;
				an_entry.mate_strand = cur_mate_strand;
				result_sr.push_back(an_entry);
				//EventEntry an_entry;
				//an_entry.type = data[0];
				//an_entry.the_id_i = data[1];
				//an_entry.the_id_k = data[2];
				//an_entry.support_i = src_return[4];
				//an_entry.seqid_i = data[3];
				//an_entry.starts_map_start = left_bound_sr;
				//an_entry.starts_map_end = cur_end;
				//an_entry.mseqid_i = data[6];
				//an_entry.mstarts_map_start = right_bound_sr;
				//an_entry.mstarts_map_end = cur_mate_start;
				//result_sr.push_back(an_entry);
				//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";

			}
		}
	} else if (cur_strand == 1 && cur_mate_strand == 1) {
		if(debug) {
			cout << "transl_inter here-20\n";
		}
		string chr_bp_1 = a_region_a->second.chr_bp_1;
		int64_t position1 = a_region_a->second.start;
		int64_t position2 = a_region_a->second.end;
		string chr_bp_2 = a_region_a->second.chr_bp_2;
		int64_t position3 = a_region_a->second.mate_start;
		int64_t position4 = a_region_a->second.mate_end;
		int32_t orientation = a_region_a->second.orientation;
		const string& local_chr_bp_1 = mate_ref_id;
		const string& local_chr_bp_2 = ref_id;
		if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
			swap(position1, position3);
			swap(position2, position4);
			//(position1, position2, position3, position4) = (position3, position4, position1, position2);
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (orientation == 1) {
			if(debug) {
				cout << "transl_inter here-21\n";
			}
			discord_sr1.clear();

			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "transl_inter here-22\n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "transl_inter here-23\n";
					}
					left_bound_sr = position1 + src_return1[1] + options.cut_sr;
				}
				int64_t right_bound_sr = 0;
				if (src_return1.size() > 2) {
					if(debug) {
						cout << "transl_inter here-24\n";
					}
					right_bound_sr = position4 - (src_return1[2] - (position2 - position1 + 101));
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
					if(debug) {
						cout << "transl_inter here-25\n";
					}
					if (src_return1.size() > 1) {
						if(debug) {
							cout << "transl_inter here-26\n";
						}
						left_bound_sr = position1 + src_return1[1];
					}
					if (src_return1.size() > 2) {
						if(debug) {
							cout << "transl_inter here-27\n";
						}
						right_bound_sr = position4 - (src_return1[2] - (position2 - position1 + 101) + options.cut_sr);
					}
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
						cout << "transl_inter here-28\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "transl_inter here-29\n";
					}
					right_bound_sr = boost::lexical_cast<int64_t>(data[8]);
				}
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "transl_inter here-30\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr - 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = cur_strand;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr - 1;
					an_entry.mate_strand = cur_mate_strand;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";
				}
			}
		} else {
			if(debug) {
				cout << "transl_inter here-31\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "transl_inter here-32\n";
				}
//	vector<int64_t> src_return;
//	map<int64_t, vector<string>> bpread;
				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 1) {
					if(debug) {
						cout << "transl_inter here-33\n";
					}
					left_bound_sr = position1 + src_return1[1] + options.cut_sr;
				}
				int64_t right_bound_sr = 0;
				if (src_return1.size() > 2) {
					if(debug) {
						cout << "transl_inter here-34\n";
					}
					right_bound_sr = src_return1[2] - (position2 - position1 + 101) + position3;
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
					if(debug) {
						cout << "transl_inter here-35\n";
					}
					if (src_return1.size() > 1) {
						if(debug) {
							cout << "transl_inter here-36\n";
						}
						left_bound_sr = position1 + src_return1[1];
					}
					if (src_return1.size() > 2) {
						if(debug) {
							cout << "transl_inter here-37\n";
						}
						right_bound_sr = src_return1[2] - (position2 - position1 + 101) + position3 + options.cut_sr;
					}
				}
				if (src_return1.size() > 1 && !src_return1[1]) {
					if(debug) {
						cout << "transl_inter here-38\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 2 && !src_return1[2]) {
					if(debug) {
						cout << "transl_inter here-39\n";
					}
					right_bound_sr = boost::lexical_cast<int64_t>(data[8]);
				}
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "transl_inter here-40\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = cur_strand;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr + 1;
					an_entry.mate_strand = cur_mate_strand;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";
				}
			}
		}

	} else if (cur_strand == -1 && cur_mate_strand == -1) {
		if(debug) {
			cout << "transl_inter here-41\n";
		}
		string chr_bp_1 = a_region_a->second.chr_bp_1;
		int64_t position1 = a_region_a->second.start;
		int64_t position2 = a_region_a->second.end;
		string chr_bp_2 = a_region_a->second.chr_bp_2;
		int64_t position3 = a_region_a->second.mate_start;
		int64_t position4 = a_region_a->second.mate_end;
		int32_t orientation = a_region_a->second.orientation;
		const string& local_chr_bp_1 = mate_ref_id;
		const string& local_chr_bp_2 = ref_id;
		if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
			swap(position1, position3);
			swap(position2, position4);
			//(position1, position2, position3, position4) = (position3, position4, position1, position2);
		}
		//my @alignments;
		//open( SAM,
		//"samtools_command view -X sr_sortbam a_region_a->second|"
		//);
		//while ( newline1 = <SAM> ) {
		//chomp newline1;
		//push @alignments, newline1;
		//}
		//close SAM;
		reader.SetRegion(local_ref_id, local_ref_start, local_ref_id, local_ref_end);
		if (orientation == 1) {
			if(debug) {
				cout << "transl_inter here-42\n";
			}
			discord_sr1.clear();

			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 100 && strand == mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "transl_inter here-43\n";
				}
				sr_cluster_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "transl_inter here-44\n";
					}
					left_bound_sr = position1 + src_return1[0];
				}
				int64_t right_bound_sr = 0;
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "transl_inter here-45\n";
					}
					right_bound_sr = position4 - (src_return1[3] - (position2 - position1 + 101) + options.cut_sr);
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
					if(debug) {
						cout << "transl_inter here-46\n";
					}
					if (src_return1.size() > 0) {
						if(debug) {
							cout << "transl_inter here-47\n";
						}
						left_bound_sr = position1 + src_return1[0] + options.cut_sr;
					}
					if (src_return1.size() > 3) {
						if(debug) {
							cout << "transl_inter here-48\n";
						}
						right_bound_sr = position4 - (src_return1[3] - (position2 - position1 + 101));
					}
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "transl_inter here-49\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "transl_inter here-50\n";
					}
					right_bound_sr = boost::lexical_cast<int64_t>(data[7]);
				}

				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "transl_inter here-51\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = cur_strand;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr + 1;
					an_entry.mate_strand = cur_mate_strand;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";

				}
			}
		} else {
			if(debug) {
				cout << "transl_inter here-52\n";
			}
			discord_sr1.clear();
			BamTools::BamAlignment al;
			//int64_t prev_bam_pos = backend_bgzf.Tell();
			while (reader.GetNextAlignmentBasic(al)) {
				int64_t start = al.Position;
				int32_t strand = 1;
				if (al.IsReverseStrand()) {
					strand = -1;
				}
//	string mseqid = ref_vec[al.MateRefID].RefName;
				int64_t mstart = al.MatePosition;
				int32_t mstrand = 1;
				if (al.IsMateReverseStrand()) {
					mstrand = -1;
				}
				int64_t isize = al.InsertSize;

				if (isize > 0 && strand != mstrand && start < (position2 - position1) && mstart > (position2 - position1 + 100)) {

					//#print "start mstart\tisize\tstrand mstrand\n";
					discord_sr1.push_back(al);
				}
			}
			if (discord_sr1.size() >= static_cast<uint64_t>(options.support_reads)) {
				if(debug) {
					cout << "transl_inter here-53\n";
				}
//	vector<int64_t> src_return;
//	map<int64_t, vector<string>> bpread;
				sr_cluster_1_no_cl(src_return1, bpread_1, discord_sr1);
				int64_t left_bound_sr = 0;
				if (src_return1.size() > 0) {
					if(debug) {
						cout << "transl_inter here-54\n";
					}
					left_bound_sr = position1 + src_return1[0];
				}
				int64_t right_bound_sr = 0;
				if (src_return1.size() > 3) {
					if(debug) {
						cout << "transl_inter here-55\n";
					}
					right_bound_sr = src_return1[3] - (position2 - position1 + 101) + position3 + options.cut_sr;
				}
				const string& local_chr_bp_1 = mate_ref_id;
				const string& local_chr_bp_2 = ref_id;

				if (chr_bp_1 == local_chr_bp_1 && chr_bp_2 == local_chr_bp_2) {
					if(debug) {
						cout << "transl_inter here-56\n";
					}
					if (src_return1.size() > 0) {
						if(debug) {
							cout << "transl_inter here-57\n";
						}
						left_bound_sr = position1 + src_return1[0] + options.cut_sr;
					}
					if (src_return1.size() > 3) {
						if(debug) {
							cout << "transl_inter here-58\n";
						}
						right_bound_sr = src_return1[3] - (position2 - position1 + 101) + position3;
					}
				}
				if (src_return1.size() > 0 && !src_return1[0]) {
					if(debug) {
						cout << "transl_inter here-59\n";
					}
					left_bound_sr = cur_start;
				}
				if (src_return1.size() > 3 && !src_return1[3]) {
					if(debug) {
						cout << "transl_inter here-60\n";
					}
					right_bound_sr = boost::lexical_cast<int64_t>(data[7]);
				}
				if (src_return1.size() > 4 && src_return1[4] >= options.support_reads) {
					if(debug) {
						cout << "transl_inter here-61\n";
					}
					string reads = castle::StringUtils::join(bpread_1[0], "\t");
					BPREAD << (boost::format("%s__%s__%s__%s\n%s\n") % ref_id % (left_bound_sr + 1) % mate_ref_id % (right_bound_sr + 1) % reads).str();
					EventEntry an_entry;
					an_entry.type = data[0];
					an_entry.cluster_id = data[1];
					an_entry.mate_cluster_id = data[2];
					an_entry.n_supports = src_return1[4];
					an_entry.ref_id = ref_id;
					an_entry.event_start = left_bound_sr + 1;
					an_entry.strand = cur_strand;
					an_entry.mate_ref_id = mate_ref_id;
					an_entry.mate_event_start = right_bound_sr + 1;
					an_entry.mate_strand = cur_mate_strand;
					result_sr.push_back(an_entry);
					//result_sr[sri][0] = "data[0]\tdata[1]\tdata[2]\tsrc_return[4]\tdata[3]\tleft_bound_sr\tcur_end\tdata[6]\tright_bound_sr\tcur_mate_start\n";

				}
			}
		}

	}

}
// cluster discordant split reads
//cl2: # there are 2 clusters in this region
//# 0 no break, 1 break on left side, 2 break on right side
void SplitReadSVCaller::sr_cluster(vector<int64_t>& boundary1, vector<int64_t>& boundary2, vector<int64_t>& support_sr, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& ref_discord_sr, int64_t cl2, int32_t internal_break) {
	auto& ref_is = options.is;
	vector<string> temp_cols;
	map<int64_t, vector<BamTools::BamAlignment>> cluster;
	int64_t count = 0;
	const char* delim_underscore = "_";

	for (auto& al : ref_discord_sr) {
		string& readname = al.Name;
		int64_t seqid_num = al.RefID;
		int64_t start = al.Position;
		int64_t mseqid_num = al.MateRefID;
		int64_t mstart = al.MatePosition;
		int64_t isize = al.InsertSize;
		if (seqid_num != mseqid_num) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(readname, delim_underscore, temp_cols);
		int64_t readlength = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
		int64_t topush = -1;
//		//cout << "sr_cluster: read " << readname << ": " << start << "/" << mstart << "/" << isize << "\n";
		for (int64_t i = count - 5; i <= count; ++i) {
			bool terminate_SR = false;
			if (cluster.end() == cluster.find(i)) {
				continue;
			}

			for (auto& m_al : cluster[i]) {
				string& readnamet = m_al.Name;
				int64_t seqid_numt = m_al.RefID;
				int64_t startt = m_al.Position;
				int64_t mseqid_numt = m_al.MateRefID;
				int64_t mstartt = m_al.MatePosition;
				int64_t isizet = m_al.InsertSize;
				if (seqid_numt != mseqid_numt) {
					continue;
				}
				castle::StringUtils::c_string_multi_split(readnamet, delim_underscore, temp_cols);
				int64_t readlengtht = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
//				//cout << "sr_cluster: pair " << readnamet << ": " << startt << "/" << mstartt << "/" << isizet << "\n";
//				if(abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"]
//																											&& abs((isize - isizet) - (readlength - readlengtht)) <= 2) {
//					topush = i;
//					terminate_SR = true;
//					break;
//					if(abs((isize - isizet) - (readlength - readlengtht)) != abs((isizet - isize) - (readlengtht - readlength))) {
//						cout << "sr_cluster: condition: " << abs((isize - isizet) - (readlength - readlengtht)) <<
//							abs((isizet - isize) - (readlengtht - readlength)) << "\n";
//					}
				if (abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"] && readlength >= readlengtht && abs((isize - isizet) - (readlength - readlengtht)) <= 2) {
					topush = i;
					terminate_SR = true;
					break;
				}
				if (abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"] && readlength < readlengtht && abs((isizet - isize) - (readlengtht - readlength)) <= 2) {
					topush = i;
					terminate_SR = true;
					break;
				}
//				}

//				}
			}
			if (terminate_SR) {
				break;
			}
		}
		if (-1 == topush) {
			++count;
			cluster[count].push_back(al);
		} else {
			cluster[topush].push_back(al);
		}
	}

	int64_t_value_desc_sortedset clustersize;
	for (int64_t i = 1; i <= count; ++i) {
		Int64Pair a_pair(i, cluster[i].size());
		clustersize.insert(a_pair);
	}
	//cout << "sr_cluster: clustersize.size(): " << clustersize.size() << "/count: " << count << "\n";
//	# split 1 cluster into 2, split at largest distance read
	if (0 != internal_break) {
		bool tosplit = (1 == clustersize.size());
		if (tosplit) {
			int64_t clusteri = clustersize.begin()->first;
			int64_t laststart = -1;
			int64_t lastmstart = -1;
			int64_t max_dis = 0;
			string max_name;
			int64_t max_pos = 0;
			map<int64_t, vector<BamTools::BamAlignment>> cluster_temp;
			// the break point resides in the first read
			if (internal_break == 1) {
//				# find largest distance
				for (auto& m_alt : cluster[clusteri]) {
					string& readname = m_alt.Name;
					int64_t start = m_alt.Position;

					if (-1 != laststart) {
						int64_t distance = abs(start - laststart);
						if (distance > max_dis && max_dis > ref_is["rlu"]["selected"]) {
							max_dis = distance;
							max_name = readname;
							max_pos = start;
						}
					}
					laststart = start;
				}

//				# split one cluster into two
				bool cluster2 = false;
				for (auto& m_alt : cluster[clusteri]) {
					string& readname = m_alt.Name;
					int64_t start = m_alt.Position;

					if (readname == max_name && start == max_pos) {
						cluster2 = true;
					}
					if (cluster2) {
						cluster_temp[2].push_back(m_alt);
					} else {
						cluster_temp[1].push_back(m_alt);
					}
				}
				cluster.swap(cluster_temp);
			}
			// the break point resides in the second read
			else if (internal_break == 2) {
				map<string, BamTools::BamAlignment> sorted_mstart_ob;
				// sort by mstart by ascending order
				string_int64_t_value_sortedset sorted_mstart_t;
				for (auto& m_alt : cluster[clusteri]) {
					string& readname = m_alt.Name;
					int64_t mstart = m_alt.MatePosition;
					sorted_mstart_ob[readname] = m_alt;
					StringInt64Pair a_pair(readname, mstart);
					sorted_mstart_t.insert(a_pair);
				}

//				# find largest distance
				for (auto& m_readname_entry : sorted_mstart_t) {
					auto& m_alt = sorted_mstart_ob[m_readname_entry.first];
					auto& readname = m_alt.Name;
					auto& mstart = m_alt.MatePosition;
					if (-1 != lastmstart) {
						int64_t distance = abs(mstart - lastmstart);
						if (distance > max_dis && max_dis > ref_is["rlu"]["selected"]) {
							max_dis = distance;
							max_name = readname;
							max_pos = mstart;
						}
					}
					lastmstart = mstart;
				}
				//cout << "sr_cluster: max_name: " << max_name << "/max_distance:"  << max_dis << "\n";

//				# split one cluster into two
				bool cluster2 = false;
				for (auto& m_readname_entry : sorted_mstart_t) {
					auto& m_alt = sorted_mstart_ob[m_readname_entry.first];
					auto& readname = m_alt.Name;
					auto& mstart = m_alt.MatePosition;
					if (readname == max_name && mstart == max_pos) {
						cluster2 = true;
					}
					if (cluster2) {
						cluster_temp[2].push_back(m_alt);
					} else {
						cluster_temp[1].push_back(m_alt);
					}
				}
				cluster.swap(cluster_temp);
			}

//			# update clustersize
			clustersize.clear();
			for (int64_t i = 1; i <= 2; ++i) {
				if (cluster.end() == cluster.find(i)) {
					continue;
				}
				Int64Pair a_pair(i, cluster[i].size());
				clustersize.insert(a_pair);
//				$clustersize{i} = @{ $cluster{i} } if ( $cluster{i} );
			}
		}
	}

//	my ( @boundary1, @boundary2, @support_sr, %bpread );
	int64_t k = 0;
	//cout << "sr_cluster: cluster data start\n";
//	for (auto& cluster_entry : clustersize) {
//		int64_t i = cluster_entry.first;
//		//cout << "sr_cluster: " << cluster_entry.first << ": " << cluster_entry.second << "\n";
//		for (auto& m_alt : cluster[i]) {
//			string& readname = m_alt.Name;
//			int64_t start = m_alt.Position;
//			int64_t mstart = m_alt.MatePosition;
//			//cout << "sr_cluster: " << readname << ": " << start << "/" << mstart << "\n";
//		}
//	}
	//cout << "sr_cluster: cluster data end\n";
	int64_t selected_i = -1;
	for (auto& cluster_entry : clustersize) {
//	  my i ( sort { $clustersize{$b} <=> $clustersize{$a} } keys %clustersize )
//	{
		vector<int64_t> start_list;
		vector<int64_t> mstart_list;
//		my ( @start, @mstart );
		int64_t i = cluster_entry.first;
		if (-1 != selected_i) {
			int64_t next_i = selected_i + 1;
			//cout << "sr_cluster: selected_i: " << selected_i << "\n";
			if (clustersize.size() > 1) {
				auto the_next = clustersize.begin();
				++the_next;
				if (the_next->second > 1) {
					next_i = the_next->first;
				} else {
					for (; next_i <= count; ++next_i) {
						if (cluster.end() == cluster.find(next_i)) {
							continue;
						}
						bool found_discrepancy = false;
						for (auto& m_alt : cluster[next_i]) {
							int64_t start = m_alt.Position;
							int64_t mstart = m_alt.MatePosition;
							//cout << "sr_cluster: boundary 1,2 : " << boundary1.back() << "," << boundary2.back() << "/" << start << "/" << mstart << "\n";
							if (boundary1.back() > start || boundary2.back() > mstart) {
								found_discrepancy = true;
								break;
							}
						}
						if (!found_discrepancy) {
							i = next_i;
							break;
						}
					}
				}
			}

//			i = count;
		}
		for (auto& m_alt : cluster[i]) {
//			my @data1    = split( /\t/, $_ );
//			my $readname = $data1[0];
//			my start    = $data1[3];
//			my mstart   = $data1[7];
			string& readname = m_alt.Name;
			int64_t start = m_alt.Position;
			int64_t mstart = m_alt.MatePosition;
//			#print "i start mstart\n";
			start_list.push_back(start);
			mstart_list.push_back(mstart);
			bpread[k].push_back(readname);
			//cout << "sr_cluster: insert before sorting: " << readname << "/" << start << "/" << mstart << "\n";
		}
		sort(start_list.begin(), start_list.end());
		sort(mstart_list.begin(), mstart_list.end());
		if ((start_list.back() - start_list[0]) < ref_is["rlu"]["selected"] && (mstart_list.back() - mstart_list[0]) < ref_is["rlu"]["selected"]) {
			int64_t local_support_sr = cluster[i].size();
			if (1 == cl2) {
				boundary1.push_back(start_list[0]);
				boundary1.push_back(start_list.back());
				boundary2.push_back(mstart_list[0]);
				boundary2.push_back(mstart_list.back());
				support_sr.push_back(local_support_sr);
				selected_i = i;
			}
		} else {
			bpread[k].clear();
			continue;
		}
		if (0 != k) {
			break;
		}
		++k;
	}
	if(boundary1.size() < 2) {
		boundary1.resize(2);
	}
	if(boundary2.size() < 2) {
		boundary2.resize(2);
	}
	if(support_sr.size() < 1) {
		support_sr.resize(1);
	}
}
// cl2 == 0, internal_break = 0
void SplitReadSVCaller::sr_cluster_no_cl(vector<int64_t>& src_return, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& discord_sr) {
	auto& ref_is = options.is;
//	int64_t cut_sr = options.cut_sr;
//	vector<string> data1;
	vector<string> temp_cols;
	map<int64_t, vector<BamTools::BamAlignment>> cluster;
	int64_t count = 0;
	const char* delim_underscore = "_";
//	const char* delim_tab = "\t";

	for (auto& al : discord_sr) {
		string& readname = al.Name;
		int64_t seqid_num = al.RefID;
		int64_t start = al.Position;
		int64_t mseqid_num = al.MateRefID;
		int64_t mstart = al.MatePosition;
		int64_t isize = al.InsertSize;
		if (seqid_num != mseqid_num) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(readname, delim_underscore, temp_cols);
		int64_t readlength = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
//		bool incluster = false;
		int64_t topush = -1;
		for (int64_t i = count - 5; i <= count; ++i) {
			bool terminate_SR = false;
			for (auto& m_al : cluster[i]) {
				string& readnamet = m_al.Name;
				int64_t seqid_numt = m_al.RefID;
				int64_t startt = m_al.Position;
				int64_t mseqid_numt = m_al.MateRefID;
				int64_t mstartt = m_al.MatePosition;
				int64_t isizet = m_al.InsertSize;
				if (seqid_numt != mseqid_numt) {
					continue;
				}
				castle::StringUtils::c_string_multi_split(readnamet, delim_underscore, temp_cols);
				int64_t readlengtht = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);

				if (readlength >= readlengtht && abs((isize - isizet) - (readlength - readlengtht)) <= 2 && abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"]) {
					topush = i;
					terminate_SR = true;
					break;
				}
				if (readlength < readlengtht && abs((isizet - isize) - (readlengtht - readlength)) <= 2 && abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"]) {
					topush = i;
					terminate_SR = true;
					break;
				}
			}
			if (terminate_SR) {
				break;
			}
		}
		if (-1 == topush) {
			++count;
			cluster[count].push_back(al);
		} else {
			cluster[topush].push_back(al);
		}
	}

	int64_t_value_desc_sortedset clustersize;
	for (int64_t i = 1; i <= count; ++i) {
		Int64Pair a_pair(i, cluster[i].size());
		clustersize.insert(a_pair);
	}
	//cout << "sr_cluster_no_cl: cluster data start\n";
//	for (auto& cluster_entry : clustersize) {
//		int64_t i = cluster_entry.first;
//		//cout << "sr_cluster_no_cl: " << cluster_entry.first << ": " << cluster_entry.second << "\n";
//		for (auto& m_alt : cluster[i]) {
//			string& readname = m_alt.Name;
//			int64_t start = m_alt.Position;
//			int64_t mstart = m_alt.MatePosition;
//			//cout << "sr_cluster_no_cl: " << readname << ": " << start << "/" << mstart << "\n";
//		}
//	}
	//cout << "sr_cluster_no_cl: cluster data end\n";
	int64_t k = 0;
	for (auto& cluster_entry : clustersize) {
		vector<int64_t> start_list;
		vector<int64_t> mstart_list;
		int64_t i = cluster_entry.first;
		for (auto& m_alt : cluster[i]) {
			string& readname = m_alt.Name;
			int64_t start = m_alt.Position;
			int64_t mstart = m_alt.MatePosition;
			start_list.push_back(start);
			mstart_list.push_back(mstart);
			bpread[k].push_back(readname);
		}
		sort(start_list.begin(), start_list.end());
		sort(mstart_list.begin(), mstart_list.end());
		if ((start_list.back() - start_list[0]) < ref_is["rlu"]["selected"] && (mstart_list.back() - mstart_list[0]) < ref_is["rlu"]["selected"]) {
			int64_t local_support_sr = cluster[i].size();
			src_return.push_back(start_list[0]);
			src_return.push_back(start_list.back());
			src_return.push_back(mstart_list[0]);
			src_return.push_back(mstart_list.back());
			src_return.push_back(local_support_sr);
			return;
		} else {
			bpread[k].clear();
			continue;
		}
		if (0 != k) {
			break;
		}
		++k;
	}
//	if(src_return.size() < 4) {
//		src_return.resize(4);
//	}
}
// cluster discordant split reads in one region (for inversion only)
void SplitReadSVCaller::sr_cluster_1(vector<int64_t>& boundary1, vector<int64_t>& boundary2, vector<int64_t>& support_sr, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& discord_sr, int64_t cl2) {
	auto& ref_is = options.is;
//	int64_t cut_sr         = options.cut_sr;
//	int64_t sr_dis_cutoff  = ref_is["rlu"]["selected"] - cut_sr;
	vector<string> temp_cols;
	map<int64_t, vector<BamTools::BamAlignment>> cluster;
//	my %cluster;
	int64_t count = 0;
	const char* delim_underscore = "_";

	for (auto& al : discord_sr) {
		string& readname = al.Name;
		int64_t seqid_num = al.RefID;
		int64_t start = al.Position;
		int64_t mseqid_num = al.MateRefID;
		int64_t mstart = al.MatePosition;
		if (seqid_num != mseqid_num) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(readname, delim_underscore, temp_cols);
		int64_t readlength = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
		int64_t topush = -1;
		for (int64_t i = count - 5; i <= count; ++i) {
			bool terminate_SR1 = false;
			if (cluster.end() == cluster.find(i) || cluster[i].empty()) {
				continue;
			}

			auto& m_al = cluster[i].back();
			string& readnamet = m_al.Name;
			int64_t startt = m_al.Position;
			int64_t mstartt = m_al.MatePosition;
			castle::StringUtils::c_string_multi_split(readnamet, delim_underscore, temp_cols);
			int64_t readlengtht = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
			if (abs(abs(start - startt + mstart - mstartt) - abs(readlength - readlengtht)) <= 2 && abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"]) {
				topush = i;
				terminate_SR1 = true;
				break;
			}
			if (terminate_SR1) {
				break;
			}
		}
		if (-1 == topush) {
			++count;
			cluster[count].push_back(al);
		} else {
			cluster[topush].push_back(al);
		}
	}

	int64_t_value_desc_sortedset clustersize;
	for (int64_t i = 1; i <= count; ++i) {
		Int64Pair a_pair(i, cluster[i].size());
		clustersize.insert(a_pair);
	}

	//cout << "sr_cluster_1: cluster data start\n";
//	for (auto& cluster_entry : clustersize) {
//		int64_t i = cluster_entry.first;
//		//cout << "sr_cluster_1: " << cluster_entry.first << ": " << cluster_entry.second << "\n";
//		for (auto& m_alt : cluster[i]) {
//			string& readname = m_alt.Name;
//			int64_t start = m_alt.Position;
//			int64_t mstart = m_alt.MatePosition;
//			//cout << "sr_cluster_1: " << readname << ": " << start << "/" << mstart << "\n";
//		}
//	}
	//cout << "sr_cluster_1: cluster data end\n";

	int64_t k = 0;
	int64_t selected_i = -1;
	for (auto& cluster_entry : clustersize) {
		vector<int64_t> start_list;
		vector<int64_t> mstart_list;
		vector<int64_t> isize_list;
		int64_t i = cluster_entry.first;
		if (cluster.end() == cluster.find(i)) {
			continue;
		}
		if (-1 != selected_i) {
			int64_t next_i = selected_i + 1;
//				//cout << "sr_cluster_1: selected_i: " << selected_i << "\n";

			if (clustersize.size() > 1) {
				auto the_next = clustersize.begin();
				++the_next;
				if (the_next->second > 1) {
					next_i = the_next->first;
				} else {
					for (; next_i <= count; ++next_i) {
						if (cluster.end() == cluster.find(next_i)) {
							continue;
						}
						bool found_discrepancy = false;
						for (auto& m_alt : cluster[next_i]) {
							int64_t start = m_alt.Position;
							int64_t mstart = m_alt.MatePosition;
							//cout << "sr_cluster: boundary 1,2 : " << boundary1.back() << "," << boundary2.back() << "/" << start << "/" << mstart << "\n";
							if (boundary1.back() > start || boundary2.back() > mstart) {
								found_discrepancy = true;
								break;
							}
						}
						if (!found_discrepancy) {
							i = next_i;
							break;
						}
					}
				}
			}

			//			i = count;
		}
		for (auto& m_alt : cluster[i]) {
			string& readname = m_alt.Name;
			int64_t start = m_alt.Position;
			int64_t mstart = m_alt.MatePosition;
			start_list.push_back(start);
			mstart_list.push_back(mstart);
			bpread[k].push_back(readname);
		}
		sort(start_list.begin(), start_list.end());
		sort(mstart_list.begin(), mstart_list.end());
		if ((start_list.back() - start_list[0]) < ref_is["rlu"]["selected"] && (mstart_list.back() - mstart_list[0]) < ref_is["rlu"]["selected"]) {
			int64_t local_support_sr = cluster[i].size();
			if (cl2) {
				boundary1.push_back(start_list[0]);
				boundary1.push_back(start_list.back());
				boundary2.push_back(mstart_list[0]);
				boundary2.push_back(mstart_list.back());
				support_sr.push_back(local_support_sr);
				selected_i = i;
			}
		} else {
			bpread[k].clear();
			continue;
		}
		if (0 != k) {
			break;
		}
		++k;
	}
}
void SplitReadSVCaller::sr_cluster_1_no_cl(vector<int64_t>& src_return, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& discord_sr) {
	auto& ref_is = options.is;
	vector<string> temp_cols;
	map<int64_t, vector<BamTools::BamAlignment>> cluster;
	int64_t count = 0;
	const char* delim_underscore = "_";
	for (auto& al : discord_sr) {
		string& readname = al.Name;
		int64_t seqid_num = al.RefID;
		int64_t start = al.Position;
		int64_t mseqid_num = al.MateRefID;
		int64_t mstart = al.MatePosition;
		if (seqid_num != mseqid_num) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(readname, delim_underscore, temp_cols);
		int64_t readlength = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
		int64_t topush = -1;
		for (int64_t i = count - 5; i <= count; ++i) {
			bool terminate_SR1 = false;
			if (cluster.end() != cluster.find(i) && !cluster[i].empty()) {
				auto& m_al = cluster[i].back();
				string& readnamet = m_al.Name;
				int64_t startt = m_al.Position;
				int64_t mstartt = m_al.MatePosition;
				castle::StringUtils::c_string_multi_split(readnamet, delim_underscore, temp_cols);
				int64_t readlengtht = boost::lexical_cast<int64_t>(temp_cols[temp_cols.size() - 1]);
				if (abs(abs(start - startt + mstart - mstartt) - abs(readlength - readlengtht)) <= 2 && abs(startt - start) < ref_is["rlu"]["selected"] && abs(mstartt - mstart) < ref_is["rlu"]["selected"]) {
					topush = i;
					terminate_SR1 = true;
					break;
				}
			}
			if (terminate_SR1) {
				break;
			}
		}
		if (-1 == topush) {
			++count;
			cluster[count].push_back(al);
		} else {
			cluster[topush].push_back(al);
		}
	}

	int64_t_value_desc_sortedset clustersize;
	for (int64_t i = 1; i <= count; ++i) {
		Int64Pair a_pair(i, cluster[i].size());
		clustersize.insert(a_pair);
	}

	int64_t k = 0;
	//cout << "sr_cluster_1_no_cl: " << count << "\n";
//	for (auto& cluster_entry : clustersize) {
//		int64_t i = cluster_entry.first;
//		//cout << "sr_cluster_1_no_cl: " << cluster_entry.first << ": " << cluster_entry.second << "\n";
//		for (auto& m_alt : cluster[i]) {
//			string& readname = m_alt.Name;
//			int64_t start = m_alt.Position;
//			int64_t mstart = m_alt.MatePosition;
//			//cout << "sr_cluster_1_no_cl: " << readname << ": " << start << "/" << mstart << "\n";
//		}
//	}

//	for (auto& cluster_entry : boost::adaptors::reverse(clustersize)) {
	for (auto& cluster_entry : clustersize) {
		vector<int64_t> start_list;
		vector<int64_t> mstart_list;
		vector<int64_t> isize_list;
		int64_t i = cluster_entry.first;
		for (auto& m_alt : cluster[i]) {
			string& readname = m_alt.Name;
			int64_t start = m_alt.Position;
			int64_t mstart = m_alt.MatePosition;
			start_list.push_back(start);
			mstart_list.push_back(mstart);
			bpread[k].push_back(readname);
		}
		sort(start_list.begin(), start_list.end());
		sort(mstart_list.begin(), mstart_list.end());
		//cout << "sr_cluster_1_no_cl: " << start_list[0] << "/" << start_list.back() << "/" << mstart_list[0] << "/" << mstart_list.back() << "\n";
		if ((start_list.back() - start_list[0]) < ref_is["rlu"]["selected"] && (mstart_list.back() - mstart_list[0]) < ref_is["rlu"]["selected"]) {
			int64_t local_support_sr = cluster[i].size();
			src_return.push_back(start_list[0]);
			src_return.push_back(start_list.back());
			src_return.push_back(mstart_list[0]);
			src_return.push_back(mstart_list.back());
			src_return.push_back(local_support_sr);
		} else {
			bpread[k].clear();
			continue;
		}
		if (0 != k) {
			break;
		}
		++k;
	}
}
} /* namespace meerkat */
