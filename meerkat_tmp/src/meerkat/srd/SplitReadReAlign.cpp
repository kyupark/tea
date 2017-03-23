/*
 * SplitReadReAlign.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#include "SplitReadReAlign.hpp"

namespace meerkat {

SplitReadReAlign::SplitReadReAlign() : n_blocks(1){
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

SplitReadReAlign::~SplitReadReAlign() {
}

void SplitReadReAlign::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
	n_cores = options.n_cores;
}


void SplitReadReAlign::align_split_reads() {
	map<string, map<int8_t, int32_t>> bp_weight;
	map<string, map<string, string>> intra_regions;
	map<string, map<string, map<string, string>>> inter_regions;
	boost::unordered_map<string, map<int8_t, string>> alg_mis;
	set<string> cluster_exist;
	collect_intra_regions(bp_weight, intra_regions, cluster_exist);
	collect_inter_regions(bp_weight, inter_regions, cluster_exist);
//	write_intra_regions(bp_weight, intra_regions, cluster_exist);
//	write_inter_regions(bp_weight, inter_regions, cluster_exist);
//	prepare_split_alignment();
//	prepare_split_alignment_single();
//	collect_misaligned_reads(alg_mis, bp_weight);
	collect_misaligned_reads_serial(alg_mis, bp_weight);
//	adjust_misaligned_reads(alg_mis);
	adjust_misaligned_reads_serial(alg_mis);
	create_misalignment_adjusted_bam();
}

void SplitReadReAlign::align_split_reads_alt() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.align_split_reads_alt");
	checker.start();

//	auto cut_sr = options.cut_sr;
//	const char* delim_tab = "\t";
//	auto& ref_is = options.is;
//	const char* delim_slash = "/";
//	const char* delim_colon = ":";
//	vector<string> data;
//	vector<string> cl;
////	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
//	string line;
//	string mpintra_outfile = options.prefix + ".mp.intra.out";
//	string mpinter_outfile = options.prefix + ".mp.inter.out";
//	string fafile = options.prefix + ".bp.fasta";
//	string bpinfofile = options.prefix + ".bp.info";
//	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//	string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
//	string sr_sai1 = options.prefix + ".sr.1.fq.gz.bak.sai";
//	string sr_sai2 = options.prefix + ".sr.2.fq.gz.bak.sai";
//	string sr_rawsam = options.prefix + ".sr.raw.sam";
//	string sr_sam = options.prefix + ".sr.sam";
//	string sr_bam = options.prefix + ".sr.bam";
//	string sr_sort = options.prefix + ".sr.sorted";
//	string sr_sortbam = options.prefix + ".sr.sorted.bam";
//	if(!options.working_dir.empty()) {
//		mpintra_outfile = options.working_prefix + ".mp.intra.out";
//		mpinter_outfile = options.working_prefix + ".mp.inter.out";
//		fafile = options.working_prefix + ".bp.fasta";
//		bpinfofile = options.working_prefix + ".bp.info";
//		fafile = options.working_prefix + ".bp.fasta";
//		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
//		sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
//		sr_rawsam = options.working_prefix + ".sr.raw.sam";
//		sr_sam = options.working_prefix + ".sr.sam";
//		sr_bam = options.working_prefix + ".sr.bam";
//		sr_sort = options.working_prefix + ".sr.sorted";
//		sr_sortbam = options.working_prefix + ".sr.sorted.bam";
//	}
////	auto sr_insertsize  = 2 * cut_sr + 100;
//
//	string n;
//	for (int64_t i = 0 ; i < 100 ; ++i) {
//		n += "N";
//	}

//	my ( $newline, $i, %bp_weight, %regions, %cluster_exist );

//#	%bp_weight: where is the real break point more likely to be located
//#	$bp_weight{regionname}[0]: left break point
//#	$bp_weight{regionname}[1]: right break point
//#	%regions, a list of existing regions
//#	$regions{$cluster_id}: region name, chr_p_p_chr_p_p_0/1, 0 same orientation, 1 opposite orientation
//#	%cluster_exist: the cluster with key id exists
//#	$cluster_exist{$cluster_id}: 1 exist, 0 not
//	vector<string> positiona;
//	vector<string> positionb;
	boost::unordered_map<string, map<int8_t, int32_t>> bp_weight;
	boost::unordered_map<string, map<string, string>> intra_regions;
	boost::unordered_map<string, map<string, map<string, string>>> inter_regions;
	boost::unordered_set<string> cluster_exist;
	collect_intra_regions_serial(bp_weight, intra_regions, cluster_exist);
	collect_inter_regions_serial(bp_weight, inter_regions, cluster_exist);

//	vector<int64_t> temp1;
//	vector<int64_t> temp2;
//	vector<string> data1;
//	vector<string> data2;

//	# merge overlapping regions
	write_intra_regions(bp_weight, intra_regions, cluster_exist);
	write_inter_regions(bp_weight, inter_regions, cluster_exist);

//	prepare_split_alignment_single();

//	#	adjust pair end mapping for mis-aligned reads
//	#	%alg_mis: mis-aligned reads
//	#	$alg_mis{readname}[0]: read id mis-aligned
//	#	$alg_mis{readname}[1]: read to be placed
	boost::unordered_map<string, map<int8_t, string>> alg_mis;
	collect_misaligned_reads_serial_alt(alg_mis, bp_weight);
	adjust_misaligned_reads_serial_alt(alg_mis);

////	vector<string> data;
//	vector<string> first_cols;
//	vector<string> xa_z_cols;
//	vector<string> alt_map_cols;
//	vector<string> p;
//
//	vector<string> match_cols;
////	const char* delim_tab = "\t";
//	const char* delim_underscore = "_";
//	const char* delim_semiconlon = ";";
//	const char* delim_comma = ",";
//	const string white_spaces = " \t\r\n";
//	string the_xa_z_pattern = "XA:Z:";
//	const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();

	create_misalignment_adjusted_bam();
//	system "$samtools_command faidx $fafile";
//	system "$samtools_command view -bt $faifile $sr_sam -o $sr_bam";
//	system "$samtools_command sort -@ $threads_bwa $sr_bam $sr_sort";
//	system "$samtools_command index $sr_sortbam";
	cout << checker;
}

void SplitReadReAlign::align_split_reads_par() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.align_split_reads_par");
	checker.start();
	boost::unordered_map<string, map<int8_t, int32_t>> bp_weight;
	boost::unordered_map<string, map<string, string>> intra_regions;
	boost::unordered_map<string, map<string, map<string, string>>> inter_regions;
	boost::unordered_set<string> cluster_exist;
	collect_intra_regions_serial(bp_weight, intra_regions, cluster_exist);
	collect_inter_regions_serial(bp_weight, inter_regions, cluster_exist);

//	# merge overlapping regions
	write_intra_regions(bp_weight, intra_regions, cluster_exist);
	write_inter_regions(bp_weight, inter_regions, cluster_exist);

	prepare_split_alignment();
	remove_and_merge_temporary_alignment_files();
//	#	adjust pair end mapping for mis-aligned reads
//	#	%alg_mis: mis-aligned reads
//	#	$alg_mis{readname}[0]: read id mis-aligned
//	#	$alg_mis{readname}[1]: read to be placed
	boost::unordered_map<string, map<int8_t, string>> alg_mis;
//	collect_misaligned_reads_serial_alt(alg_mis, bp_weight);
	collect_misaligned_reads_par(alg_mis, bp_weight);
//	adjust_misaligned_reads_serial_alt(alg_mis);
	adjust_misaligned_reads_par(alg_mis);

	create_misalignment_adjusted_bam();

	cout << checker;
}

void SplitReadReAlign::collect_intra_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions, boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.collect_intra_regions_serial");
	checker.start();
	auto& ref_is = options.is;
	auto cut_sr = options.cut_sr;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	const char* delim_colon = ":";
	vector<string> data;
	vector<string> cl;
//	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
	string line;
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}

	vector<string> positiona;
	vector<string> positionb;
	ifstream MPINTRAALG(mpintra_outfile, ios::binary);
	while (getline(MPINTRAALG, line, '\n')) {
		if(line.empty()) {
			continue;
		}
//		cout << line << "\n";
		castle::StringUtils::tokenize(line, delim_tab, data);
		if(data.empty()) {
			continue;
		}
//		#print "$data[1]\t$data[7]\n";
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		for(auto& a_cl : cl) {
			if(string::npos == a_cl.find("_0")) {
				a_cl = a_cl + "_0";
			}
		}
//		if(string::npos != data[1].find("73269")) {
//			cout << "[SplitReadReAlign.collect_intra_regions_serial] " << line << "\n";
//		}

		if ("del" == data[0]) {
			castle::StringUtils::c_string_multi_split(data[7], delim_colon, positiona);
			int64_t position1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t position2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t position3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t position4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			position1 = position1 - delta;
			position2 = position2 + delta;
			position3 = position3 - delta;
			position4 = position4 + delta;
			string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" + boost::lexical_cast<string>(position2) + "__" + data[3]
								+ "__" + boost::lexical_cast<string>(position3) + "__" + boost::lexical_cast<string>(position4) + "__0";
			bp_weight[name][0] = 1;
			bp_weight[name][1] = 4;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if (string::npos != data[0].find("inss")) {
			castle::StringUtils::c_string_multi_split(data[12], delim_colon, positiona);
			int64_t positiona1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t positiona2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t positiona3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t positiona4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			positiona1 = positiona1 - delta;
			positiona2 = positiona2 + delta;
			positiona3 = positiona3 - delta;
			positiona4 = positiona4 + delta;

			string namea = data[3] + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + data[3]
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;

			castle::StringUtils::c_string_multi_split(data[13], delim_colon, positionb);
			int64_t positionb1 = boost::lexical_cast<int64_t>(positionb[0]);
			int64_t positionb2 = boost::lexical_cast<int64_t>(positionb[1]);
			int64_t positionb3 = boost::lexical_cast<int64_t>(positionb[2]);
			int64_t positionb4 = boost::lexical_cast<int64_t>(positionb[3]);

			positionb1 = positionb1 - delta;
			positionb2 = positionb2 + delta;
			positionb3 = positionb3 - delta;
			positionb4 = positionb4 + delta;

			string nameb = data[3] + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + data[3]
					+ "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__0";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
				intra_regions[data[3]][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if (string::npos != data[0].find("inso")) {
			castle::StringUtils::c_string_multi_split(data[12], delim_colon, positiona);
			int64_t positiona1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t positiona2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t positiona3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t positiona4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			positiona1 = positiona1 - delta;
			positiona2 = positiona2 + delta;
			positiona3 = positiona3 - delta;
			positiona4 = positiona4 + delta;
			string namea = data[3] + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + data[3]
								+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";

			castle::StringUtils::c_string_multi_split(data[13], delim_colon, positionb);
			int64_t positionb1 = boost::lexical_cast<int64_t>(positionb[0]);
			int64_t positionb2 = boost::lexical_cast<int64_t>(positionb[1]);
			int64_t positionb3 = boost::lexical_cast<int64_t>(positionb[2]);
			int64_t positionb4 = boost::lexical_cast<int64_t>(positionb[3]);

			positionb1 = positionb1 - delta;
			positionb2 = positionb2 + delta;
			positionb3 = positionb3 - delta;
			positionb4 = positionb4 + delta;

			string nameb = data[3] + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + data[3]
								+ "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__1";
			if (string::npos != data[0].find("insod")) {
				bp_weight[namea][0] = 1;
				bp_weight[namea][1] = 3;
				bp_weight[nameb][0] = 2;
				bp_weight[nameb][1] = 4;
			} else {
				bp_weight[namea][0] = 2;
				bp_weight[namea][1] = 4;
				bp_weight[nameb][0] = 1;
				bp_weight[nameb][1] = 3;
			}
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
				intra_regions[data[3]][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if ("invers" == data[0] || "del_invers" == data[0]) {
			castle::StringUtils::c_string_multi_split(data[10], delim_colon, positiona);
			int64_t positiona1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t positiona2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t positiona3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t positiona4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			positiona1 = positiona1 - delta;
			positiona2 = positiona2 + delta;
			positiona3 = positiona3 - delta;
			positiona4 = positiona4 + delta;
			string namea = data[3] + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + data[3]
								+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";

			castle::StringUtils::c_string_multi_split(data[11], delim_colon, positionb);
			int64_t positionb1 = boost::lexical_cast<int64_t>(positionb[0]);
			int64_t positionb2 = boost::lexical_cast<int64_t>(positionb[1]);
			int64_t positionb3 = boost::lexical_cast<int64_t>(positionb[2]);
			int64_t positionb4 = boost::lexical_cast<int64_t>(positionb[3]);

			positionb1 = positionb1 - delta;
			positionb2 = positionb2 + delta;
			positionb3 = positionb3 - delta;
			positionb4 = positionb4 + delta;
			string nameb = data[3] + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + //
					data[3] + "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__1";
			if (!(covered(positiona1, positiona2, positionb1, positionb2) && covered(positiona3, positiona4, positionb3, positionb4))) {
				bp_weight[namea][0] = 1;
				bp_weight[namea][1] = 3;
				bp_weight[nameb][0] = 2;
				bp_weight[nameb][1] = 4;
			}
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
				intra_regions[data[3]][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if ("tandem_dup" == data[0]) {
			castle::StringUtils::c_string_multi_split(data[7], delim_colon, positiona);
			int64_t position1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t position2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t position3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t position4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			position1 = position1 - delta;
			position2 = position2 + delta;
			position3 = position3 - delta;
			position4 = position4 + delta;
			string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" + boost::lexical_cast<string>(position2) + "__" + data[3]
					+ "__" + boost::lexical_cast<string>(position3) + "__" + boost::lexical_cast<string>(position4) + "__0";
			bp_weight[name][0] = 2;
			bp_weight[name][1] = 3;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if ("invers_f" == data[0]) {
			castle::StringUtils::c_string_multi_split(data[7], delim_colon, positiona);
			int64_t position1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t position2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t position3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t position4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			position1 = position1 - delta;
			position2 = position2 + delta;
			position3 = position3 - delta;
			position4 = position4 + delta;
			string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" + boost::lexical_cast<string>(position2) + "__" + data[3]
					+ "__" + boost::lexical_cast<string>(position3) + "__" + boost::lexical_cast<string>(position4) + "__1";
			bp_weight[name][0] = 1;
			bp_weight[name][1] = 3;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if ("invers_r" == data[0]) {
			castle::StringUtils::c_string_multi_split(data[7], delim_colon, positiona);
			int64_t position1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t position2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t position3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t position4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			position1 = position1 - delta;
			position2 = position2 + delta;
			position3 = position3 - delta;
			position4 = position4 + delta;
			string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" + boost::lexical_cast<string>(position2) + "__" + data[3]
					+ "__" + boost::lexical_cast<string>(position3) + "__" + boost::lexical_cast<string>(position4) + "__1";
			bp_weight[name][0] = 2;
			bp_weight[name][1] = 4;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				intra_regions[data[3]][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		}
	}
	cout << checker;

}
void SplitReadReAlign::collect_intra_regions(map<string, map<int8_t, int32_t>>& bp_weight, map<string, map<string, string>>& regions,
		set<string>& cluster_exist) {
	cout << "[SplitReadReAlign.collect_intra_regions] start\n";
	auto& ref_is = options.is;
	vector<string> data;
	vector<string> cl;
	vector<string> positions;
	vector<string> positions_a;
	vector<string> positions_b;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	const char* delim_colon = ":";
	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
	string line;
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
	}
	ifstream MPINTRAALG(mpintra_outfile, ios::binary);
	while (getline(MPINTRAALG, line, '\n')) {
		if(line.empty()) {
			continue;
		}
		castle::StringUtils::tokenize(line, delim_tab, data);
		if(data.empty()) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		for(auto& a_cl : cl) {
			if(string::npos == a_cl.find("_0")) {
				a_cl += "_0";
			}
		}
//		const bool debug = data[1] == "14904_0" || data[1] == "4446_0/11983_0";
		const bool debug = false;
		if (debug) {
			cout << "[SplitReadReAlign.collect_intra_regions] line: " << line << "\n";
		}
//		const bool debug = false;
		int64_t positiona1 = 0;
		int64_t positiona2 = 0;
		int64_t positiona3 = 0;
		int64_t positiona4 = 0;

		int64_t positionb1 = 0;
		int64_t positionb2 = 0;
		int64_t positionb3 = 0;
		int64_t positionb4 = 0;
		string& type = data[0];
		string& ref_id = data[3];
		if("tandem_dup" == type || "invers_f" == type || "invers_r" == type || "transl_inter" == type || "del" == type) {
			castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_colon, positions_a);
			positiona1 = boost::lexical_cast<int64_t>(positions_a[0]);
			positiona2 = boost::lexical_cast<int64_t>(positions_a[1]);
			positiona3 = boost::lexical_cast<int64_t>(positions_a[2]);
			positiona4 = boost::lexical_cast<int64_t>(positions_a[3]);
//		} else if(string::npos != type.find("inssu") || string::npos != type.find("inssd") || string::npos != type.find("insod") || string::npos != type.find("insou") ||
//				string::npos != type.find("inso") || string::npos != type.find("inss") || string::npos != type.find("invers")) {
		} else {
			castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions_a);
			castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_colon, positions_b);
			positiona1 = boost::lexical_cast<int64_t>(positions_a[0]);
			positiona2 = boost::lexical_cast<int64_t>(positions_a[1]);
			positiona3 = boost::lexical_cast<int64_t>(positions_a[2]);
			positiona4 = boost::lexical_cast<int64_t>(positions_a[3]);

			positionb1 = boost::lexical_cast<int64_t>(positions_b[0]);
			positionb2 = boost::lexical_cast<int64_t>(positions_b[1]);
			positionb3 = boost::lexical_cast<int64_t>(positions_b[2]);
			positionb4 = boost::lexical_cast<int64_t>(positions_b[3]);
		}
		positiona1 -= delta;
		positiona2 += delta;
		positiona3 -= delta;
		positiona4 += delta;

		positionb1 -= delta;
		positionb2 += delta;
		positionb3 -= delta;
		positionb4 += delta;

		positiona1 = max(static_cast<int64_t>(0), positiona1);
		positiona2 = max(static_cast<int64_t>(0), positiona2);
		positiona3 = max(static_cast<int64_t>(0), positiona3);
		positiona4 = max(static_cast<int64_t>(0), positiona4);

		positionb1 = max(static_cast<int64_t>(0), positionb1);
		positionb2 = max(static_cast<int64_t>(0), positionb2);
		positionb3 = max(static_cast<int64_t>(0), positionb3);
		positionb4 = max(static_cast<int64_t>(0), positionb4);


		if ("del" == data[0]) {
			string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[name][0] = 1;
			bp_weight[name][1] = 4;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if (string::npos != data[0].find("inss")) {
			string namea = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;

			string nameb = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__0";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;

//			string namec = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
//					+ "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
//
//			bp_weight[namec][0] = 1;
//			bp_weight[namec][1] = 4;

			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_intra_regions] cl[0]: " << cl[0] << ", namea: " << namea << "\n";
				}
				regions[ref_id][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_intra_regions] cl[1]: " << cl[1] << ", nameb: " << nameb << "\n";
				}
				regions[ref_id][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}

		} else if (string::npos != data[0].find("inso")) {
			string namea = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";

			string nameb = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__1";

			if (string::npos != data[0].find("insod")) {
				bp_weight[namea][0] = 1;
				bp_weight[namea][1] = 3;
				bp_weight[nameb][0] = 2;
				bp_weight[nameb][1] = 4;
			} else {
				bp_weight[namea][0] = 2;
				bp_weight[namea][1] = 4;
				bp_weight[nameb][0] = 1;
				bp_weight[nameb][1] = 3;
			}
//			const bool debug = "del_insou" == data[0] && "514922_0/41657_0" == data[1];
//			if(debug) {
//				cout << "insou: here-0\n";
//			}
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
//				if(debug) {
//					cout << "insou: here-1: " << namea << "\n";
//				}
				regions[ref_id][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
//				if(debug) {
//					cout << "insou: here-2: " << nameb << "\n";
//				}
				regions[ref_id][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if ("invers" == data[0] || "del_invers" == data[0]) {
			string namea = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";
			string nameb = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + //
					ref_id + "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__1";
			if (!(covered(positiona1, positiona2, positionb1, positionb2) && covered(positiona3, positiona4, positionb3, positionb4))) {
				bp_weight[namea][0] = 1;
				bp_weight[namea][1] = 3;
				bp_weight[nameb][0] = 2;
				bp_weight[nameb][1] = 4;
			}
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if (cluster_exist.end() == cluster_exist.find(cl[1])) {
				regions[ref_id][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if ("tandem_dup" == data[0]) {
			string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[name][0] = 2;
			bp_weight[name][1] = 3;
			if(debug) {
				cout << "[SplitReadReAlign.collect_intra_regions] cl[0]: " << cl[0] << ", name: " << name << "\n";
			}
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if ("invers_f" == data[0]) {
			string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";
			bp_weight[name][0] = 1;
			bp_weight[name][1] = 3;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		} else if ("invers_r" == data[0]) {
			string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + ref_id
					+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__1";
			bp_weight[name][0] = 2;
			bp_weight[name][1] = 4;
			if (cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][cl[0]] = name;
				cluster_exist.insert(cl[0]);
			}
		}
	}

}

void SplitReadReAlign::collect_inter_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& inter_regions, boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
		checker.setTarget("SplitReadReAlign.collect_inter_regions_serial");
		checker.start();
		auto& ref_is = options.is;
		auto cut_sr = options.cut_sr;
		const char* delim_tab = "\t";
		const char* delim_slash = "/";
		const char* delim_colon = ":";
		vector<string> data;
		vector<string> cl;
	//	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
		string line;
//		string mpintra_outfile = options.prefix + ".mp.intra.out";
		string mpinter_outfile = options.prefix + ".mp.inter.out";
//		string fafile = options.prefix + ".bp.fasta";
//		string bpinfofile = options.prefix + ".bp.info";
//		string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//		string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
//		string sr_sai1 = options.prefix + ".sr.1.fq.gz.bak.sai";
//		string sr_sai2 = options.prefix + ".sr.2.fq.gz.bak.sai";
//		string sr_rawsam = options.prefix + ".sr.raw.sam";
//		string sr_sam = options.prefix + ".sr.sam";
//		string sr_bam = options.prefix + ".sr.bam";
//		string sr_sort = options.prefix + ".sr.sorted";
//		string sr_sortbam = options.prefix + ".sr.sorted.bam";
		if(!options.working_dir.empty()) {
//			mpintra_outfile = options.working_prefix + ".mp.intra.out";
			mpinter_outfile = options.working_prefix + ".mp.inter.out";
//			fafile = options.working_prefix + ".bp.fasta";
//			bpinfofile = options.working_prefix + ".bp.info";
//			fafile = options.working_prefix + ".bp.fasta";
//			sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
//			sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
//			sr_rawsam = options.working_prefix + ".sr.raw.sam";
//			sr_sam = options.working_prefix + ".sr.sam";
//			sr_bam = options.working_prefix + ".sr.bam";
//			sr_sort = options.working_prefix + ".sr.sorted";
//			sr_sortbam = options.working_prefix + ".sr.sorted.bam";
		}

	//	string n;
	//	for (int64_t i = 0 ; i < 100 ; ++i) {
	//		n += "N";
	//	}
		vector<string> positiona;
		vector<string> positionb;
	ifstream MPINTERALG(mpinter_outfile, ios::binary);
	while (getline(MPINTERALG, line, '\n')) {
		if(line.empty()) {
			continue;
		}
		castle::StringUtils::tokenize(line, delim_tab, data);
		if(data.empty()) {
			continue;
		}
//		const bool debug = ("inso" == data[0] && data[1] == "16624_0/16417_0");
//		if(string::npos != data[1].find("73269")) {
//			cout << "[SplitReadReAlign.collect_inter_regions_serial] " << line << "\n";
//		}
		const bool debug = false;
		if(debug) {
			cout << line << "\n";
		}

		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		for(auto& a_cl : cl) {
			if(string::npos == a_cl.find("_0")) {
				a_cl = a_cl + "_0";
			}
		}
		if (string::npos != data[0].find("inss")) {
			castle::StringUtils::c_string_multi_split(data[11], delim_colon, positiona);
			int64_t positiona1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t positiona2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t positiona3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t positiona4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			positiona1 = positiona1 - delta;
			positiona2 = positiona2 + delta;
			positiona3 = positiona3 - delta;
			positiona4 = positiona4 + delta;
			string namea = data[3] + "__" + boost::lexical_cast<string>(positiona1) + "__" + boost::lexical_cast<string>(positiona2) + "__" + data[7]
								+ "__" + boost::lexical_cast<string>(positiona3) + "__" + boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;

			castle::StringUtils::c_string_multi_split(data[12], delim_colon, positionb);
			int64_t positionb1 = boost::lexical_cast<int64_t>(positionb[0]);
			int64_t positionb2 = boost::lexical_cast<int64_t>(positionb[1]);
			int64_t positionb3 = boost::lexical_cast<int64_t>(positionb[2]);
			int64_t positionb4 = boost::lexical_cast<int64_t>(positionb[3]);

			positionb1 = positionb1 - delta;
			positionb2 = positionb2 + delta;
			positionb3 = positionb3 - delta;
			positionb4 = positionb4 + delta;
			string nameb = data[3] + "__" + boost::lexical_cast<string>(positionb1) + "__" + boost::lexical_cast<string>(positionb2) + "__" + //
					data[7] + "__" + boost::lexical_cast<string>(positionb3) + "__" + boost::lexical_cast<string>(positionb4) + "__0";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;
			if(cluster_exist.end() == cluster_exist.find(cl[0])) {
				inter_regions[data[3]][data[7]][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if(cluster_exist.end() == cluster_exist.find(cl[1])) {
				inter_regions[data[3]][data[7]][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if (string::npos != data[0].find("inso")) {
			castle::StringUtils::c_string_multi_split(data[11], delim_colon, positiona);
			int64_t positiona1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t positiona2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t positiona3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t positiona4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			positiona1 = positiona1 - delta;
			positiona2 = positiona2 + delta;
			positiona3 = positiona3 - delta;
			positiona4 = positiona4 + delta;
			string namea = data[3] + "__" + boost::lexical_cast<string>(positiona1) + "__"
					+ boost::lexical_cast<string>(positiona2) + "__"
					+ data[7] + "__" + boost::lexical_cast<string>(positiona3) + "__"
					+ boost::lexical_cast<string>(positiona4) + "__1";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;

			castle::StringUtils::c_string_multi_split(data[12], delim_colon, positionb);
			int64_t positionb1 = boost::lexical_cast<int64_t>(positionb[0]);
			int64_t positionb2 = boost::lexical_cast<int64_t>(positionb[1]);
			int64_t positionb3 = boost::lexical_cast<int64_t>(positionb[2]);
			int64_t positionb4 = boost::lexical_cast<int64_t>(positionb[3]);

			positionb1 = positionb1 - delta;
			positionb2 = positionb2 + delta;
			positionb3 = positionb3 - delta;
			positionb4 = positionb4 + delta;

			string nameb = data[3] + "__" + boost::lexical_cast<string>(positionb1) + "__" +
			boost::lexical_cast<string>(positionb2) + "__" + data[7] + "__" +
			boost::lexical_cast<string>(positionb3) + "__" +
			boost::lexical_cast<string>(positionb4) + "__1";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;
			if(cluster_exist.end() == cluster_exist.find(cl[0])) {
				inter_regions[data[3]][data[7]][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if(cluster_exist.end() == cluster_exist.find(cl[1])) {
				inter_regions[data[3]][data[7]][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if("transl_inter" == data[0]) {
			string& the_strand = data[5];
			string& the_mate_strand = data[8];
			castle::StringUtils::c_string_multi_split(data[9], delim_colon, positiona);
			int64_t position1 = boost::lexical_cast<int64_t>(positiona[0]);
			int64_t position2 = boost::lexical_cast<int64_t>(positiona[1]);
			int64_t position3 = boost::lexical_cast<int64_t>(positiona[2]);
			int64_t position4 = boost::lexical_cast<int64_t>(positiona[3]);
			int64_t delta = ref_is["rlu"]["selected"] - cut_sr + 10;
			position1 = position1 - delta;
			position2 = position2 + delta;
			position3 = position3 - delta;
			position4 = position4 + delta;
			if("1" == the_strand && "-1" == the_mate_strand) {
				string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" +
				boost::lexical_cast<string>(position2) + "__" + data[6] + "__" +
				boost::lexical_cast<string>(position3) + "__" +
				boost::lexical_cast<string>(position4) + "__0";
				bp_weight[name][0] = 1;
				bp_weight[name][1] = 4;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					inter_regions[data[3]][data[6]][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("-1" == the_strand && "1" == the_mate_strand) {
				string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" +
				boost::lexical_cast<string>(position2) + "__" + data[6] + "__" +
				boost::lexical_cast<string>(position3) + "__" +
				boost::lexical_cast<string>(position4) + "__0";
				bp_weight[name][0] = 2;
				bp_weight[name][1] = 3;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					inter_regions[data[3]][data[6]][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("1" == the_strand && "1" == the_mate_strand) {
				string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" +
				boost::lexical_cast<string>(position2) + "__" + data[6] + "__" +
				boost::lexical_cast<string>(position3) + "__" +
				boost::lexical_cast<string>(position4) + "__1";
				bp_weight[name][0] = 1;
				bp_weight[name][1] = 3;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					inter_regions[data[3]][data[6]][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("-1" == the_strand && "-1" == the_mate_strand) {
				string name = data[3] + "__" + boost::lexical_cast<string>(position1) + "__" +
				boost::lexical_cast<string>(position2) + "__" + data[6] + "__" +
				boost::lexical_cast<string>(position3) + "__" +
				boost::lexical_cast<string>(position4) + "__1";
				bp_weight[name][0] = 2;
				bp_weight[name][1] = 4;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					inter_regions[data[3]][data[6]][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			}
		}
	}
	cout << checker;
}

void SplitReadReAlign::collect_inter_regions(map<string, map<int8_t, int32_t>>& bp_weight, map<string, map<string, map<string, string>>>& regions, set<string>& cluster_exist) {
	cout << "[SplitReadReAlign.collect_inter_regions] start\n";
	auto& ref_is = options.is;
	vector<string> data;
	vector<string> cl;
	vector<string> positions;
	vector<string> positions_a;
	vector<string> positions_b;
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	const char* delim_colon = ":";
	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
	string line;
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	if(!options.working_dir.empty()) {
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
	}
	ifstream MPINTERALG(mpinter_outfile, ios::binary);
	while (getline(MPINTERALG, line, '\n')) {
		if(line.empty()) {
			continue;
		}
		castle::StringUtils::tokenize(line, delim_tab, data);
		if(data.empty()) {
			continue;
		}
//		const bool debug = ("inso" == data[0] && data[1] == "16624_0/16417_0");
		const bool debug = false;
		if(debug) {
			cout << line << "\n";
		}
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl);
		for(auto& a_cl : cl) {
			if(string::npos == a_cl.find("_0")) {
				a_cl += "_0";
			}
		}
		int64_t positiona1 = 0;
		int64_t positiona2 = 0;
		int64_t positiona3 = 0;
		int64_t positiona4 = 0;

		int64_t positionb1 = 0;
		int64_t positionb2 = 0;
		int64_t positionb3 = 0;
		int64_t positionb4 = 0;
		string& type = data[0];
		string& ref_id = data[3];
		string mate_ref_id = data[7];
		if("tandem_dup" == type || "invers_f" == type || "invers_r" == type || "transl_inter" == type || "del" == type) {
			castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_colon, positions_a);
			positiona1 = boost::lexical_cast<int64_t>(positions_a[0]);
			positiona2 = boost::lexical_cast<int64_t>(positions_a[1]);
			positiona3 = boost::lexical_cast<int64_t>(positions_a[2]);
			positiona4 = boost::lexical_cast<int64_t>(positions_a[3]);
//		} else if(string::npos != type.find("inssu") || string::npos != type.find("inssd") || string::npos != type.find("insod") || string::npos != type.find("insou") ||
//				string::npos != type.find("inso") || string::npos != type.find("inss") || string::npos != type.find("invers")) {
		} else {
			castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions_a);
			castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_colon, positions_b);
			positiona1 = boost::lexical_cast<int64_t>(positions_a[0]);
			positiona2 = boost::lexical_cast<int64_t>(positions_a[1]);
			positiona3 = boost::lexical_cast<int64_t>(positions_a[2]);
			positiona4 = boost::lexical_cast<int64_t>(positions_a[3]);

			positionb1 = boost::lexical_cast<int64_t>(positions_b[0]);
			positionb2 = boost::lexical_cast<int64_t>(positions_b[1]);
			positionb3 = boost::lexical_cast<int64_t>(positions_b[2]);
			positionb4 = boost::lexical_cast<int64_t>(positions_b[3]);
		}

		if("transl_inter" == type) {
			mate_ref_id = data[6];
		}

//		castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions_a);
//		castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_colon, positions_b);
//
//		int64_t positiona1 = boost::lexical_cast<int64_t>(positions_a[0]);
//		int64_t positiona2 = boost::lexical_cast<int64_t>(positions_a[1]);
//		int64_t positiona3 = boost::lexical_cast<int64_t>(positions_a[2]);
//		int64_t positiona4 = boost::lexical_cast<int64_t>(positions_a[3]);
		positiona1 -= delta;
		positiona2 += delta;
		positiona3 -= delta;
		positiona4 += delta;

//		int64_t positionb1 = boost::lexical_cast<int64_t>(positions_b[0]);
//		int64_t positionb2 = boost::lexical_cast<int64_t>(positions_b[1]);
//		int64_t positionb3 = boost::lexical_cast<int64_t>(positions_b[2]);
//		int64_t positionb4 = boost::lexical_cast<int64_t>(positions_b[3]);
		positionb1 -= delta;
		positionb2 += delta;
		positionb3 -= delta;
		positionb4 += delta;
		positiona1 = max(static_cast<int64_t>(0), positiona1);
		positiona2 = max(static_cast<int64_t>(0), positiona2);
		positiona3 = max(static_cast<int64_t>(0), positiona3);
		positiona4 = max(static_cast<int64_t>(0), positiona4);

		positionb1 = max(static_cast<int64_t>(0), positionb1);
		positionb2 = max(static_cast<int64_t>(0), positionb2);
		positionb3 = max(static_cast<int64_t>(0), positionb3);
		positionb4 = max(static_cast<int64_t>(0), positionb4);

		if (string::npos != data[0].find("inss")) {
			string namea = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" +
			boost::lexical_cast<string>(positiona2) + "__" +
			mate_ref_id + "__" +
			boost::lexical_cast<string>(positiona3) + "__" +
			boost::lexical_cast<string>(positiona4) + "__0";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;

			string nameb = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" +
			boost::lexical_cast<string>(positionb2) + "__" + mate_ref_id + "__" +
			boost::lexical_cast<string>(positionb3) + "__" +
			boost::lexical_cast<string>(positionb4) + "__0";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;
			if(cluster_exist.end() == cluster_exist.find(cl[0])) {
				regions[ref_id][mate_ref_id][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if(cluster_exist.end() == cluster_exist.find(cl[1])) {
				regions[ref_id][mate_ref_id][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			}
		} else if (string::npos != data[0].find("inso")) {
			string namea = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__"
					+ boost::lexical_cast<string>(positiona2) + "__"
					+ mate_ref_id + "__" + boost::lexical_cast<string>(positiona3) + "__"
					+ boost::lexical_cast<string>(positiona4) + "__1";
			bp_weight[namea][0] = 2;
			bp_weight[namea][1] = 3;
			string nameb = ref_id + "__" + boost::lexical_cast<string>(positionb1) + "__" +
			boost::lexical_cast<string>(positionb2) + "__" + mate_ref_id + "__" +
			boost::lexical_cast<string>(positionb3) + "__" +
			boost::lexical_cast<string>(positionb4) + "__1";
			bp_weight[nameb][0] = 1;
			bp_weight[nameb][1] = 4;
			if(debug) {
				cout << "inso here-0\n";
			}
			if(cluster_exist.end() == cluster_exist.find(cl[0])) {
				if(debug) {
					cout << "inso here-1\n";
				}
				regions[ref_id][mate_ref_id][cl[0]] = namea;
				cluster_exist.insert(cl[0]);
			}
			if(cluster_exist.end() == cluster_exist.find(cl[1])) {
				if(debug) {
					cout << "inso here-2\n";
				}
				regions[ref_id][mate_ref_id][cl[1]] = nameb;
				cluster_exist.insert(cl[1]);
			} else {
				if(debug) {
					cout << "inso here-3: << " << cl[1] << ": " << regions[ref_id][mate_ref_id][cl[1]] << "\n";
				}
			}
			if(debug) {
				cout << "inso here-4: << " << cl[0] << "/" << cl[1] << "\n";
			}
		} else if("transl_inter" == data[0]) {
			string& the_strand = data[5];
			string& the_mate_strand = data[8];
			if("1" == the_strand && "-1" == the_mate_strand) {
//				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions);
//castle::StringUtils::c_string_multi_split(data[9], delim_pos, positions);
//				int64_t position1 = boost::lexical_cast<int64_t>(positions[0]);
//				int64_t position2 = boost::lexical_cast<int64_t>(positions[1]);
//				int64_t position3 = boost::lexical_cast<int64_t>(positions[2]);
//				int64_t position4 = boost::lexical_cast<int64_t>(positions[3]);
//				position1 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position2 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position3 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position4 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);

				string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" +
				boost::lexical_cast<string>(positiona2) + "__" + mate_ref_id + "__" +
				boost::lexical_cast<string>(positiona3) + "__" +
				boost::lexical_cast<string>(positiona4) + "__0";
				bp_weight[name][0] = 1;
				bp_weight[name][1] = 4;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					regions[ref_id][mate_ref_id][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("-1" == the_strand && "1" == the_mate_strand) {
//				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions);
//				int64_t position1 = boost::lexical_cast<int64_t>(positions[0]);
//				int64_t position2 = boost::lexical_cast<int64_t>(positions[1]);
//				int64_t position3 = boost::lexical_cast<int64_t>(positions[2]);
//				int64_t position4 = boost::lexical_cast<int64_t>(positions[3]);
//				position1 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position2 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position3 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position4 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
				string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" +
				boost::lexical_cast<string>(positiona2) + "__" + mate_ref_id + "__" +
				boost::lexical_cast<string>(positiona3) + "__" +
				boost::lexical_cast<string>(positiona4) + "__0";
				bp_weight[name][0] = 2;
				bp_weight[name][1] = 3;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					regions[ref_id][mate_ref_id][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("1" == the_strand && "1" == the_mate_strand) {
//castle::StringUtils::c_string_multi_split(data[9], delim_pos, positions);
//				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions);
//				int64_t position1 = boost::lexical_cast<int64_t>(positions[0]);
//				int64_t position2 = boost::lexical_cast<int64_t>(positions[1]);
//				int64_t position3 = boost::lexical_cast<int64_t>(positions[2]);
//				int64_t position4 = boost::lexical_cast<int64_t>(positions[3]);
//				position1 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position2 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position3 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position4 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
				string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" +
				boost::lexical_cast<string>(positiona2) + "__" + mate_ref_id + "__" +
				boost::lexical_cast<string>(positiona3) + "__" +
				boost::lexical_cast<string>(positiona4) + "__1";
				bp_weight[name][0] = 1;
				bp_weight[name][1] = 3;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					regions[ref_id][mate_ref_id][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			} else if("-1" == the_strand && "-1" == the_mate_strand) {
//castle::StringUtils::c_string_multi_split(data[9], delim_pos, positions);
//				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_colon, positions);
//				int64_t position1 = boost::lexical_cast<int64_t>(positions[0]);
//				int64_t position2 = boost::lexical_cast<int64_t>(positions[1]);
//				int64_t position3 = boost::lexical_cast<int64_t>(positions[2]);
//				int64_t position4 = boost::lexical_cast<int64_t>(positions[3]);
//				position1 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position2 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position3 -= (ref_is["rlu"]["selected"] - options.cut_sr + 10);
//				position4 += (ref_is["rlu"]["selected"] - options.cut_sr + 10);
				string name = ref_id + "__" + boost::lexical_cast<string>(positiona1) + "__" +
				boost::lexical_cast<string>(positiona2) + "__" + mate_ref_id + "__" +
				boost::lexical_cast<string>(positiona3) + "__" +
				boost::lexical_cast<string>(positiona4) + "__1";
				bp_weight[name][0] = 2;
				bp_weight[name][1] = 4;
				if(cluster_exist.end() == cluster_exist.find(cl[0])) {
					regions[ref_id][mate_ref_id][cl[0]] = name;
					cluster_exist.insert(cl[0]);
				}
			}
		}
	}
}

void SplitReadReAlign::write_intra_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions, boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.write_intra_regions_serial");
	checker.start();
//	auto cut_sr = options.cut_sr;
	//	auto& ref_is = options.is;
	//	const char* delim_slash = "/";
	//	const char* delim_colon = ":";
		vector<string> data;
		vector<string> cl;
	//	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
		string line;
		string mpintra_outfile = options.prefix + ".mp.intra.out";
		string mpinter_outfile = options.prefix + ".mp.inter.out";
		string fafile = options.prefix + ".bp.fasta";
		string bpinfofile = options.prefix + ".bp.info";
		string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
		string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
		string sr_sai1 = options.prefix + ".sr.1.fq.gz.bak.sai";
		string sr_sai2 = options.prefix + ".sr.2.fq.gz.bak.sai";
//		string sr_rawsam = options.prefix + ".sr.raw.sam";
		string sr_sam = options.prefix + ".sr.sam";
		string sr_bam = options.prefix + ".sr.bam";
		string sr_sort = options.prefix + ".sr.sorted";
		string sr_sortbam = options.prefix + ".sr.sorted.bam";
		if(!options.working_dir.empty()) {
			mpintra_outfile = options.working_prefix + ".mp.intra.out";
			mpinter_outfile = options.working_prefix + ".mp.inter.out";
			fafile = options.working_prefix + ".bp.fasta";
			bpinfofile = options.working_prefix + ".bp.info";
			fafile = options.working_prefix + ".bp.fasta";
			sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
			sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
//			sr_rawsam = options.working_prefix + ".sr.raw.sam";
			sr_sam = options.working_prefix + ".sr.sam";
			sr_bam = options.working_prefix + ".sr.bam";
			sr_sort = options.working_prefix + ".sr.sorted";
			sr_sortbam = options.working_prefix + ".sr.sorted.bam";
		}

		string n;
		for (int64_t i = 0 ; i < 100 ; ++i) {
			n += "N";
		}

	vector<string> positiona;
	vector<string> positionb;
	const char* delim_double_underscore = "__";
	vector<int64_t> temp1;
	vector<int64_t> temp2;
	vector<int64_t> temp_merge;
	vector<string> data1;
	vector<string> data2;

	IndexedFasta fai(options.reference_path.c_str());
	ofstream FA(fafile, ios::binary);
	ofstream BPDT(bpinfofile, ios::binary);
	for(auto& chr_entry : intra_regions) {
		auto& chr = chr_entry.first;
		bool merge = true;
		while(merge) {
			merge = false;
			for(auto& key_entry: intra_regions[chr]) {
				auto& key = key_entry.first;
				castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
				if (data.size() < 5) {
					continue;
				}
				temp1.clear();
				temp1.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp1.push_back(boost::lexical_cast<int64_t>(data[2]));
				temp2.clear();
				temp2.push_back(boost::lexical_cast<int64_t>(data[4]));
				temp2.push_back(boost::lexical_cast<int64_t>(data[5]));
				if (covered(temp1[0], temp1[1], temp2[0], temp2[1])) {
					merge = true;
					temp_merge.clear();
					temp_merge.push_back(temp1[0]);
					temp_merge.push_back(temp1[2]);
					temp_merge.push_back(temp1[3]);
					temp_merge.push_back(temp1[4]);
					sort(temp_merge.begin(), temp_merge.end());
					intra_regions[chr][key] = data[0] + "__" + boost::lexical_cast<string>(temp_merge[0]) + "__" + boost::lexical_cast<string>(temp_merge.back());
//					#print STDERR "@data\n$regions{a}{$chr}{$key}\n";
				}
			}
//			bool next_ma = false;
			boost::unordered_set<string> removed_keys;
			for(auto& key1_entry : intra_regions[chr]) {
				 auto& key1 = key1_entry.first;
				 if(removed_keys.end() != removed_keys.find(key1)) {
					 continue;
				 }
				 castle::StringUtils::tokenize(key1_entry.second, delim_double_underscore, data1);
				 for(auto& key2_entry: intra_regions[chr]) {
					 auto& key2 = key2_entry.first;
					 if(removed_keys.end() != removed_keys.find(key2)) {
						 continue;
					 }
					 if(key1 == key2) {
						 continue;
					 }
					 castle::StringUtils::tokenize(key2_entry.second, delim_double_underscore, data2);
					 if(data1.size() > 4 && data2.size() > 4) {
						 temp1.clear();
						 temp2.clear();

						 temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
						 temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
						 temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
						 temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));

						 temp2.push_back(boost::lexical_cast<int64_t>(data1[4]));
						 temp2.push_back(boost::lexical_cast<int64_t>(data1[5]));
						 temp2.push_back(boost::lexical_cast<int64_t>(data2[4]));
						 temp2.push_back(boost::lexical_cast<int64_t>(data2[5]));

							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])
								&& covered(temp2[0], temp2[1], temp2[2], temp2[3])
							  ){
								merge = true;
								sort(temp1.begin(), temp1.end());
								sort(temp2.begin(), temp2.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;

								if (data1[6] == data2[6]) {
									intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
											+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" + boost::lexical_cast<string>(temp2[0])
											+ "__" + boost::lexical_cast<string>(temp2.back()) + "__" + data1[6];
								} else {
									intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
											+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" + boost::lexical_cast<string>(temp2[0])
											+ "__" + boost::lexical_cast<string>(temp2.back()) + "__10";
								}
//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next_ma = true;
								break;
//								next MA;
							}
						} else if (data1.size() > 4 && data2.size() < 4 && data2.size() > 2) {
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;
								intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
								+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" + data1[4]
								+ "__" + data1[5] + "__" + data1[6];
//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next_ma = true;
								break;
//								next MA;
							}
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[4]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[5]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;
								intra_regions[chr][newkey] = data1[0] + "__" + data1[1] + "__" + data1[2] + "__" + data1[3] + "__"
										+ boost::lexical_cast<string>(temp1[0]) + "__" + boost::lexical_cast<string>(temp1.back()) + "__" + data1[6];
//								delete $regions{a}{$chr}{$key1};
//								delete $regions{a}{$chr}{$key2};

//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next_ma = true;
								break;
//								next MA;
							}
						} else if (data2.size() > 4 && data1.size() < 4 && data1.size() > 2) {
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;
								intra_regions[chr][newkey] = data2[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
								+ boost::lexical_cast<string>(temp1.back()) + "__" + data2[3] + "__" + data2[4] + "__" + data2[5] + "__"
								+ data2[6];
//								delete $regions{a}{$chr}{$key1};
//								delete $regions{a}{$chr}{$key2};

//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next_ma = true;
								break;
//								next MA;
							}
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[4]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[5]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;
								intra_regions[chr][newkey] = data2[0] + "__" + data2[1] + "__" + data2[2] + "__" + data2[3] + "__"
										+ boost::lexical_cast<string>(temp1[0]) + "__" + boost::lexical_cast<string>(temp1.back()) + "__" + data2[6];
//								delete $regions{a}{$chr}{$key1};
//								delete $regions{a}{$chr}{$key2};

//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next_ma = true;
								break;
//								next MA;
							}
						} else if(3 == data1.size() && 3 == data2.size()) {
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
							if (covered(temp1[0], temp1[1], temp1[2], temp1[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								removed_keys.insert(key1);
								removed_keys.insert(key2);
								string newkey = key1 + "/" + key2;
								intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__" +  boost::lexical_cast<string>(temp1.back());
//								delete $regions{a}{$chr}{$key1};
//								delete $regions{a}{$chr}{$key2};
//								my $newkey = $key1 . '/' . $key2;
//								$regions{a}{$chr}{$newkey} =
//								  $data1[0] . '__' . $temp[0] . '__' . $temp[-1];
//								next_ma = true;
								break;
//					   #print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
//								next MA;
							}
						}
				 }
			}
		}
		for (auto key_entry : intra_regions[chr]) {
			auto& key = key_entry.first;
			castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
			temp1.clear();
			if(data.size() > 2) {
				temp1.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp1.push_back(boost::lexical_cast<int64_t>(data[2]));
			}
			if(data.size() > 5) {
				temp1.push_back(boost::lexical_cast<int64_t>(data[4]));
				temp1.push_back(boost::lexical_cast<int64_t>(data[5]));
				if(temp1[0] < 0 || temp1[1] < 0 || temp1[2] < 0 || temp1[3] < 0) {
					continue;
				}
			}
			string seq = fai.fetch_1_system(data[0], temp1[0], temp1[1]);
			if (data.size() > 3) {
				seq += n;
				if ("1" == data[6]) {
					seq += fai.fetch_1_system(data[3], temp1[3], temp1[2]);
				} else {
					seq += fai.fetch_1_system(data[3], temp1[2], temp1[3]);
				}
			}
			BPDT << (boost::format("%s\t%s\n") % key % key_entry.second).str();
			if (!seq.empty()) {
				FA << (boost::format(">%s\n%s\n") % key_entry.second % seq).str();
			}
		}
	}
	cout << checker;
}
void SplitReadReAlign::write_intra_regions(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions,
		boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.write_intra_regions");
	checker.start();
	string fafile = options.prefix + ".bp.fasta";
	string bpinfofile = options.prefix + ".bp.info";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
		bpinfofile = options.working_prefix + ".bp.info";
	}
//string faifile = options.prefix + ".bp.fasta.fai";
//int64_t sr_insertsize = 2 * options.cut_sr + 100;
//auto& ref_is = options.is;
	const char* delim_double_underscore = "__";
//# merge overlapping regions
	vector<string> data;
	vector<string> data1;
	vector<string> data2;
	vector<int64_t> temp;
	vector<int64_t> temp1;
	vector<int64_t> temp2;
	bool merge = true;
	map<string, IntervalEntryVector> local_intervals;
//map<string, IntervalEntryVector> local_intervals_second;
	ofstream FA(fafile, ios::binary);
	ofstream BPDT(bpinfofile, ios::binary);
	set<string> removed_entries;
	set<string> merged_no_matched_entries;
	set<string> merged_entries;
	string n(100, 'N');
	cout << (boost::format("[SplitReadReAlign.write_intra_regions] Ref: %s\n") % options.reference_path).str();
	IndexedFasta fai(options.reference_path.c_str());
	for (auto& chr_entry : intra_regions) {
		auto chr = chr_entry.first;
//cout << "create key chr: " << chr << "\n";
		for (auto& key_entry : chr_entry.second) {
			auto& key = key_entry.first;
			if (removed_entries.end() != removed_entries.find(key) || merged_no_matched_entries.end() != merged_no_matched_entries.find(key)) {
				continue;
			}
			castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
			if (data.size() < 5) {
				continue;
			}
//			cout << key << ":" << key_entry.second << "\n";
//			cout << data[1] << "/" << data[2] << "/" << data[4] << "/" << data[5] << "\n";

			temp1.clear();
			temp1.push_back(boost::lexical_cast<int64_t>(data[1]));
			temp1.push_back(boost::lexical_cast<int64_t>(data[2]));
			temp2.clear();
			temp2.push_back(boost::lexical_cast<int64_t>(data[4]));
			temp2.push_back(boost::lexical_cast<int64_t>(data[5]));
			sort(temp1.begin(), temp1.end());
			sort(temp2.begin(), temp2.end());
			if (covered(temp1[0], temp1[1], temp2[0], temp2[1])) {
				temp.clear();
				temp.push_back(temp1[0]);
				temp.push_back(temp1.back());
				temp.push_back(temp2[0]);
				temp.push_back(temp2.back());
				sort(temp.begin(), temp.end());
				intra_regions[chr][key] = data[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back());
				merged_entries.insert(key);
				IntervalEntryType additional_entry;
				additional_entry.p_id = key;
				IntervalEntry entry(temp[0], temp.back(), additional_entry);
				local_intervals[chr].push_back(entry);
				continue;
			}
			IntervalEntryType additional_entry;
			additional_entry.p_id = key;

			IntervalEntry entry(temp1[0], temp1.back(), additional_entry);
			IntervalEntry entry_second(temp2[0], temp2.back(), additional_entry);
			local_intervals[chr].push_back(entry);
			local_intervals[chr].push_back(entry_second);
//local_intervals_second[chr].push_back(entry_second);
//cout << "create key: " << key << "\n";
//cout << "create key entry: " << key_entry.second << "\n";
		}
	}
	cout << "[SplitReadReAlign.write_intra_regions] start merging\n";
	IntervalEntryVector results;
	results.reserve(10000);

	for (auto& chr_entry : intra_regions) {
		auto chr = chr_entry.first;
//cout << "Chr: " << chr << "(" << chr_entry.second.size() << ")"<< "\n";
//removed_entries.clear();
		merge = true;
		while (merge) {
			merge = false;
			for (auto& key_entry : intra_regions[chr]) {
				auto& key = key_entry.first;

				castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
				if (data.size() < 5) {
					continue;
				}
//				if(key == "20171_0/6276_0") {
//					cout << "prepare key_value: " << key_entry.second << "\n";
//				}
				temp.clear();
				temp.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp.push_back(boost::lexical_cast<int64_t>(data[2]));
				temp.push_back(boost::lexical_cast<int64_t>(data[4]));
				temp.push_back(boost::lexical_cast<int64_t>(data[5]));
				if (covered(temp[0], temp[1], temp[2], temp[3])) {
					merge = true;
					sort(temp.begin(), temp.end());
					intra_regions[chr][key] = data[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back());
					IntervalEntryType additional_entry;
					additional_entry.p_id = key;
					IntervalEntry entry(temp[0], temp.back(), additional_entry);
					local_intervals[chr].push_back(entry);
				}
			}
//cout << "Chr " << chr << " is initially merged\n";
			IntervalClusterTree local_interval_tree(local_intervals[chr]);
// MA
			for (auto& key1_entry : intra_regions[chr]) {
				auto& key1 = key1_entry.first;
				if (removed_entries.end() != removed_entries.find(key1) || merged_no_matched_entries.end() != merged_no_matched_entries.find(key1)) {
					continue;
				}
//const bool debug = true;
				const bool debug = false;
//const bool debug = string::npos != key1.find("20171") || string::npos != key1.find("6276");

if(debug) {
cout << "key1: " << key1 << ", key1 value: " << key1_entry.second << "\n";
}

				castle::StringUtils::tokenize(key1_entry.second, delim_double_underscore, data1);
				temp1.clear();
				results.clear();
//if(data1.size() > 6) {
//temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
//temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
//temp2.clear();
//temp2.push_back(boost::lexical_cast<int64_t>(data1[4]));
//temp2.push_back(boost::lexical_cast<int64_t>(data1[5]));
//sort(temp1.begin(), temp1.end());
//sort(temp2.begin(), temp2.end());
//local_interval_tree.find_overlap(temp1[0], temp1.back(), results);
//local_interval_tree.find_overlap(temp2[0], temp2.back(), results);
//} else {
				temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
				temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
				sort(temp1.begin(), temp1.end());
				local_interval_tree.find_overlap(temp1[0], temp1.back(), results);
//}

//if(debug) {
//cout << "key2 counts: " << results.size() << "\n";
//}
//cout << key1 << ": " << results.size() << "\n";
//for (auto key2_entry : intra_regions[chr]) {
				bool local_merged = false;
				for (auto& result_entry : results) {
					string key2 = result_entry.value.p_id;
//auto& key2 = key2_entry.first;

					if (key1 == key2) {
						continue;
					}
					if (merged_entries.end() != merged_entries.find(key2)) {
						if (string::npos != key1.find(key2) || string::npos != key2.find(key1)) {
							continue;
						}
					} else {
						if (removed_entries.end() != removed_entries.find(key2)) {
							continue;
						}
					}
					if(merged_no_matched_entries.end() != merged_no_matched_entries.find(key2)) {
						continue;
					}

//					if(string::npos != key1.find("/") && string::npos != key2.find("/")) {
//						continue;
//					}
					auto& key2_value = intra_regions[chr][key2];

					castle::StringUtils::tokenize(key2_value, delim_double_underscore, data2);
					if(3 == data2.size()) {
						int64_t first_value = boost::lexical_cast<int64_t>(data2[1]);
						int64_t second_value = boost::lexical_cast<int64_t>(data2[2]);
						if(-first_value == second_value) {
							removed_entries.insert(key2);
							continue;
						}
					}
					if(debug) {
					cout << "key2: " << key2 << ", key2 value: " << key2_value << "\n";
					}
					if (data1.size() > 6 && data2.size() > 6) {
						if(debug) {
cout << "case 1\n";
cout << data1.size() << "/" << data2.size() << "\n";
cout << "temp1: " << data1[1] << "/" << data1[2] << "/" << data2[1] << "/" << data2[2] << "\n";
cout << "temp2: " << data1[4] << "/" << data1[5] << "/" << data2[4] << "/" << data2[5] << "\n";

					}
						temp1.clear();
						temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));
						temp2.clear();
//cout << data1[4] << "/" << data1[5] << "/" << data2[4] << "/" << data2[5] << "\n";
						temp2.push_back(boost::lexical_cast<int64_t>(data1[4]));
						temp2.push_back(boost::lexical_cast<int64_t>(data1[5]));
						temp2.push_back(boost::lexical_cast<int64_t>(data2[4]));
						temp2.push_back(boost::lexical_cast<int64_t>(data2[5]));
						if (covered(temp1[0], temp1[1], temp1[2], temp1[3]) && covered(temp2[0], temp2[1], temp2[2], temp2[3])) {
							local_merged = true;
							merge = true;
							sort(temp1.begin(), temp1.end());
							sort(temp2.begin(), temp2.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;

							if (data1[6] == data2[6]) {
								if(debug) {
									cout << "case 1-a: " << newkey << "\n";
								}
								intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
										+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" + boost::lexical_cast<string>(temp2[0])
										+ "__" + boost::lexical_cast<string>(temp2.back()) + "__" + data1[6];
								IntervalEntryType additional_entry;
								additional_entry.p_id = newkey;
								IntervalEntry entry(temp1[0], temp1.back(), additional_entry);
								IntervalEntry entry_second(temp2[0], temp2.back(), additional_entry);
								local_intervals[chr].push_back(entry);
								local_intervals[chr].push_back(entry_second);
							} else {
								if(debug) {
									cout << "case 1-b: " << newkey << "\n";
								}
								intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
										+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" + boost::lexical_cast<string>(temp2[0])
										+ "__" + boost::lexical_cast<string>(temp2.back()) + "__10";
								IntervalEntryType additional_entry;
								additional_entry.p_id = newkey;
								IntervalEntry entry(temp1[0], temp1.back(), additional_entry);
								IntervalEntry entry_second(temp2[0], temp2.back(), additional_entry);
								local_intervals[chr].push_back(entry);
								local_intervals[chr].push_back(entry_second);
							}
							break;
						}
//						if(string::npos != key1.find("/") && !local_merged) {
//							if(debug) {
//								cout << "case 1 merged no matched: " << key1 << "\n";
//							}
//							merged_no_matched_entries.insert(key1);
//						}
					} else if (data1.size() > 6 && data2.size() > 1 && data2.size() < 4) {
						if(debug) {
cout << "case 2\n";
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__"
									+ boost::lexical_cast<string>(temp.back()) + "__" + data1[3] + "__" + data1[4] + "__" + data1[5] + "__"
									+ data1[6];
							if(debug) {
								cout << "case 2-a: " << newkey << "\n";
							}
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[4]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[5]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							if(debug) {
								cout << "case 2-b: " << newkey << "\n";
							}
							intra_regions[chr][newkey] = data1[0] + "__" + data1[1] + "__" + data1[2] + "__" + data1[3] + "__"
									+ boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back()) + "__" + data1[6];
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
//						if(string::npos != key1.find("/") && !local_merged) {
//							if(debug) {
//								cout << "case 2 merged no matched: " << key1 << "\n";
//							}
//							merged_no_matched_entries.insert(key1);
//						}
					} else if (data1.size() > 6 && data2.size() > 1 && data2.size() < 4) {
						if(debug) {
cout << "case 3\n";
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[3], temp[4])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							if(debug) {
								cout << "case 3-a: " << newkey << "\n";
							}
							intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__"
									+ boost::lexical_cast<string>(temp.back()) + "__" + data1[3] + "__" + data1[4] + "__" + data1[5] + "__"
									+ data1[6];
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[4]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[5]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							if(debug) {
								cout << "case 3-b: " << newkey << "\n";
							}
							intra_regions[chr][newkey] = data1[0] + "__" + data1[1] + "__" + data1[2] + "__" + data1[3] + "__"
									+ boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back()) + "__" + data1[6];
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
//						if(string::npos != key1.find("/") && !local_merged) {
//							if(debug) {
//								cout << "case 3 merged no matched: " << key1 << "\n";
//							}
//							merged_no_matched_entries.insert(key1);
//						}
					} else if (data1.size() < 4 && data1.size() > 1 && data2.size() > 6) {
						if(debug) {
cout << "case 4-1\n";
						}
						temp.clear();
//cout << data1.size() << "/" << data2.size() << "\n";
//cout << data1[0] << "/" << data1[1] << "/" << data1[2] << "/" << data2.size() << "\n";
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;

							intra_regions[chr][newkey] = data2[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__"
									+ boost::lexical_cast<string>(temp.back()) + "__" + data2[3] + "__" + data2[4] + "__" + data2[5] + "__"
									+ data2[6];
							if(debug) {
								cout << "case 4-a: " << newkey << " " << intra_regions[chr][newkey] << "\n";
							}
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[4]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[5]));

						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							if(debug) {
								cout << "case 4-b: " << newkey << "\n";
							}
							intra_regions[chr][newkey] = data2[0] + "__" + data2[1] + "__" + data2[2] + "__" + data2[3] + "__"
									+ boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back()) + "__" + data2[6];
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
//						if(string::npos != key1.find("/") && !local_merged) {
//							if(debug) {
//								cout << "case 4 merged no matched: " << key1 << "\n";
//							}
//							merged_no_matched_entries.insert(key1);
//						}
					} else if (data1.size() > 1 && data1.size() < 4 && data2.size() > 1 && data2.size() < 4) {
						if(debug) {
cout << "case 5\n";
						}
						temp.clear();
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
						if (covered(temp[0], temp[1], temp[2], temp[3])) {
							merge = true;
							local_merged = true;
							sort(temp.begin(), temp.end());
							removed_entries.insert(key1);
							removed_entries.insert(key2);
							string newkey = key1 + "/" + key2;
							if(debug) {
								cout << "case 5-a: " << newkey << "\n";
							}
							intra_regions[chr][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp[0]) + "__"
									+ boost::lexical_cast<string>(temp.back());
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp[0], temp.back(), additional_entry);
							local_intervals[chr].push_back(entry);
							break;
						}
//						if(string::npos != key1.find("/") && !local_merged) {
//							if(debug) {
//								cout << "case 5 merged no matched: " << key1 << "\n";
//							}
//							merged_no_matched_entries.insert(key1);
//						}
//temp.clear();
//temp.push_back(boost::lexical_cast<int64_t>(data1[4]));
//temp.push_back(boost::lexical_cast<int64_t>(data1[5]));
//temp.push_back(boost::lexical_cast<int64_t>(data2[1]));
//temp.push_back(boost::lexical_cast<int64_t>(data2[2]));
//if (covered(temp[0], temp[1], temp[2], temp[3])) {
//merge = true;
//sort(temp.begin(), temp.end());
//removed_entries.insert(key1);
//removed_entries.insert(key2);
//string newkey = key1 + "/" + key2;
//intra_regions[chr][newkey] = data2[0] + "__" + data2[1] + "__" + data2[2] + "__" + data2[3] + "__" + boost::lexical_cast<string>(temp[0]) + "__" + boost::lexical_cast<string>(temp.back()) + "__" + data2[6];
//next_ma = true;
//break;
//}
					}
				}
				if(string::npos != key1.find("/") && !local_merged) {
					if(debug) {
						cout << "merged no matched: " << key1 << " : " << intra_regions[chr][key1] << "\n";
					}
					merged_no_matched_entries.insert(key1);
				}
			}
		}
		for (auto& key_entry : intra_regions[chr]) {
			auto& key = key_entry.first;
//			const bool debug = string::npos != key.find("6276_0");
			const bool debug = false;
			if(debug) {
				cout << "last key entry (org): " << key << " " << key_entry.second << "\n";
			}
			if (removed_entries.end() != removed_entries.find(key)) {
				continue;
			}
			if(debug) {
				cout << "last key entry: " << key << " " << key_entry.second << "\n";
			}
			castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
			temp.clear();
			if (data.size() < 6) {
				temp.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp.push_back(boost::lexical_cast<int64_t>(data[2]));
				if (temp[0] < 0 || temp[1] < 0) {
					continue;
				}
				string seq = fai.fetch_1_system(data[0], temp[0], temp[1]);
				BPDT << (boost::format("%s\t%s\n") % key % key_entry.second).str();
				if (!seq.empty()) {
					FA << (boost::format(">%s\n%s\n") % key_entry.second % seq).str();
				}
				continue;
			}

			temp.push_back(boost::lexical_cast<int64_t>(data[1]));
			temp.push_back(boost::lexical_cast<int64_t>(data[2]));
			temp.push_back(boost::lexical_cast<int64_t>(data[4]));
			temp.push_back(boost::lexical_cast<int64_t>(data[5]));
//			if (temp[0] < 0 || temp[1] < 0 || temp[2] < 0 || temp[3] < 0) {
//				continue;
//			}
			temp[0] = max(static_cast<int64_t>(0), temp[0]);
			temp[1] = max(static_cast<int64_t>(0), temp[1]);
			temp[2] = max(static_cast<int64_t>(0), temp[2]);
			temp[3] = max(static_cast<int64_t>(0), temp[3]);
// ref_id, start => end
			string seq = fai.fetch_1_system(data[0], temp[0], temp[1]);
			if (data.size() > 2) {
				seq += n;
				bool is_reversed = ("1" == data[data.size() - 1]);
				if (temp[3] < temp[2]) {
					if(is_reversed) {
						seq += fai.fetch_1_system(data[3], temp[2], temp[3]);
					} else {
						seq += fai.fetch_1_system(data[3], temp[3], temp[2]);
					}
				} else {
					if(is_reversed) {
						seq += fai.fetch_1_system(data[3], temp[3], temp[2]);
					} else {
						seq += fai.fetch_1_system(data[3], temp[2], temp[3]);
					}
				}
			}
			BPDT << (boost::format("%s\t%s\n") % key % key_entry.second).str();
			if (!seq.empty()) {
				FA << (boost::format(">%s\n%s\n") % key_entry.second % seq).str();
			}
		}
	}
	cout << checker;
}
void SplitReadReAlign::write_inter_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& inter_regions, boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.write_inter_regions_serial");
	checker.start();
	string fafile = options.prefix + ".bp.fasta";
	string bpinfofile = options.prefix + ".bp.info";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
		bpinfofile = options.working_prefix + ".bp.info";
	}
//string faifile = options.prefix + ".bp.fasta.fai";
//int64_t sr_insertsize = 2 * options.cut_sr + 100;
//auto& ref_is = options.is;
	const char* delim_double_underscore = "__";
//# merge overlapping regions
	vector<string> data;
	vector<string> data1;
	vector<string> data2;
	vector<int64_t> temp;
	vector<int64_t> temp1;
	vector<int64_t> temp2;
	map<string, IntervalEntryVector> local_intervals;
//map<string, IntervalEntryVector> local_intervals_second;
	ofstream FA(fafile, ios::binary | ios::app);
	ofstream BPDT(bpinfofile, ios::binary | ios::app);
	set<string> removed_entries;
	set<string> merged_no_matched_entries;
	set<string> merged_entries;
	string n(100, 'N');
	cout << (boost::format("[SplitReadReAlign.write_intra_regions] Ref: %s\n") % options.reference_path).str();
	IndexedFasta fai(options.reference_path.c_str());
	for(auto& chr1_entry: inter_regions) {
			auto& chr1 = chr1_entry.first;
			for(auto& chr2_entry : inter_regions[chr1]) {
				auto& chr2 = chr2_entry.first;
				bool merge = true;
				while (merge) {
					merge = false;
					boost::unordered_set<string> removed_keys;
	//			  MI:
					for (auto& key1_entry : inter_regions[chr1][chr2]) {
						auto& key1 = key1_entry.first;
						if(removed_keys.end() != removed_keys.find(key1)) {
							 continue;
						 }
						castle::StringUtils::tokenize(key1_entry.second, delim_double_underscore, data1);
						for(auto key2_entry : inter_regions[chr1][chr2]) {
							auto& key2 = key2_entry.first;
							if(removed_keys.end() != removed_keys.find(key2)) {
								 continue;
							 }
							if(key1 == key2) {
								continue;
							}
							castle::StringUtils::tokenize(key2_entry.second, delim_double_underscore, data2);
							temp1.clear();
							temp1.push_back(boost::lexical_cast<int64_t>(data1[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data1[2]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));

							temp2.clear();
							temp2.push_back(boost::lexical_cast<int64_t>(data1[4]));
							temp2.push_back(boost::lexical_cast<int64_t>(data1[5]));
							temp2.push_back(boost::lexical_cast<int64_t>(data2[4]));
							temp2.push_back(boost::lexical_cast<int64_t>(data2[5]));

							if (covered(temp1[0], temp1[1], temp1[2], temp1[3]) && covered(temp2[0], temp2[1], temp2[2], temp2[3])) {
								merge = true;
								sort(temp1.begin(), temp1.end());
								sort(temp2.begin(), temp2.end());
								string newkey = key1 + "/" + key2;
								removed_keys.insert(key1);
								removed_keys.insert(key2);
	//							delete $regions{e}{$chr1}{$chr2}{$key1};
	//							delete $regions{e}{$chr1}{$chr2}{$key2};
	//							my $newkey = $key1 . '/' . $key2;

								if (data1[6] == data2[6]) {
									inter_regions[chr1][chr2][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
									+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" +
									boost::lexical_cast<string>(temp2[0]) + "__" + boost::lexical_cast<string>(temp2.back()) + "__"
									+ data1[6];
								} else {
									inter_regions[chr1][chr2][newkey] = data1[0] + "__"
									+ boost::lexical_cast<string>(temp1[0]) + "__"
									+ boost::lexical_cast<string>(temp1.back()) + "__"
									+ data1[3] + "__"
									+ boost::lexical_cast<string>(temp2[0]) + "__"
									+ boost::lexical_cast<string>(temp2.back()) + "__10";
								}

	//		   #print STDERR "@data1\n@data2\n$regions{e}{$chr1}{$chr2}{$newkey}\n";
	//							next MI;
								break;
							}
						}
					}
				}
				for (auto key_entry : inter_regions[chr1][chr2]) {
					auto& key = key_entry.first;
					castle::StringUtils::tokenize(key_entry.second, delim_double_underscore, data);
					temp1.clear();
					if(data.size() > 2) {
						temp1.push_back(boost::lexical_cast<int64_t>(data[1]));
						temp1.push_back(boost::lexical_cast<int64_t>(data[2]));
					}
					if(data.size() > 5) {
						temp1.push_back(boost::lexical_cast<int64_t>(data[4]));
						temp1.push_back(boost::lexical_cast<int64_t>(data[5]));
						if(temp1[0] < 0 || temp1[1] < 0 || temp1[2] < 0 || temp1[3] < 0) {
							continue;
						}
					}
					string seq = fai.fetch_1_system(data[0], temp1[0], temp1[1]);
					if (data.size() > 3) {
						seq += n;
						if ("1" == data[6]) {
							seq += fai.fetch_1_system(data[3], temp1[3], temp1[2]);
						} else {
							seq += fai.fetch_1_system(data[3], temp1[2], temp1[3]);
						}
					}
					BPDT << (boost::format("%s\t%s\n") % key % key_entry.second).str();
					if (!seq.empty()) {
						FA << (boost::format(">%s\n%s\n") % key_entry.second % seq).str();
					}
				}
			}
		}
	cout << checker;
}
void SplitReadReAlign::write_inter_regions(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& inter_regions, boost::unordered_set<string>& cluster_exist) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.write_inter_regions");
	checker.start();
	string fafile = options.prefix + ".bp.fasta";
	string bpinfofile = options.prefix + ".bp.info";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
		bpinfofile = options.working_prefix + ".bp.info";
	}
	const char* delim_double_underscore = "__";
	vector<string> data;
	vector<string> data1;
	vector<string> data2;
	vector<int64_t> temp;
	vector<int64_t> temp1;
	vector<int64_t> temp2;
	bool merge = true;
// chr1, chr2, Interval
	map<string, map<string, IntervalEntryVector>> local_intervals;
	for(auto& chr1_entry : inter_regions) {
		auto& chr1 = chr1_entry.first;
		for(auto& chr2_entry : chr1_entry.second) {
			auto& chr2 = chr2_entry.first;
			for(auto& cl_entry: chr2_entry.second) {
				castle::StringUtils::tokenize(cl_entry.second, delim_double_underscore, data);
					try {
				temp1.clear();
				temp1.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp1.push_back(boost::lexical_cast<int64_t>(data[2]));

				temp2.clear();
				temp2.push_back(boost::lexical_cast<int64_t>(data[4]));
				temp2.push_back(boost::lexical_cast<int64_t>(data[5]));
				} catch(exception& ex) {
					cout << "write_inter_regions: here-0: " << cl_entry.first << "/" << cl_entry.second << "\n";
					continue;
				}
				IntervalEntryType additional_entry;
				additional_entry.p_id = cl_entry.first;

				IntervalEntry entry(temp1[0], temp1.back(), additional_entry);
				IntervalEntry entry_second(temp2[0], temp2.back(), additional_entry);
				local_intervals[chr1][chr2].push_back(entry);
				local_intervals[chr1][chr2].push_back(entry_second);
			}
		}
	}
	IntervalEntryVector results;
	results.reserve(10000);
	ofstream FA(fafile, ios::binary | ios::app);
	ofstream BPDT(bpinfofile, ios::binary | ios::app);
	set<string> removed_entries;
	set<string> merged_keys;
	string n(100, 'N');
	IndexedFasta fai(options.reference_path.c_str());

	for(auto& chr1_entry : inter_regions) {
		auto& chr1 = chr1_entry.first;

		for(auto& chr2_entry : chr1_entry.second) {
			auto& chr2 = chr2_entry.first;
			IntervalClusterTree local_interval_tree(local_intervals[chr1][chr2]);
			merge = true;
			while (merge) {
				merge = false;
//  MI:
				for (auto& key1_entry : inter_regions[chr1][chr2]) {
					auto& key1 = key1_entry.first;
					string combined_key_1 = chr1 + "_" + chr2 + "_" + key1;
					if(removed_entries.end() != removed_entries.find(combined_key_1)) {
						continue;
					}

					castle::StringUtils::tokenize(key1_entry.second, delim_double_underscore, data1);
					results.clear();
					temp.clear();
					try {
						temp.push_back(boost::lexical_cast<int64_t>(data1[1]));
						temp.push_back(boost::lexical_cast<int64_t>(data1[2]));
					} catch (exception& ex) {
						cout << "write_inter_regions: here-1: " << combined_key_1 << ", " << key1_entry.second << "\n";
						continue;
					}
					sort(temp.begin(), temp.end());
					local_interval_tree.find_overlap(temp[0], temp[1], results);
//cout << results.size() << "/" << inter_regions[chr1][chr2].size() << "\n";
//for(auto& key2_entry : inter_regions[chr1][chr2])  {
					for(auto& overlap_entry : results) {
						auto& key2 = overlap_entry.value.p_id;
						auto& value = inter_regions[chr1][chr2][key2];
						if(key1 == key2) {
							continue;
						}
						string combined_key_2 = chr1 + "_" + chr2 + "_" + key2;
						if(removed_entries.end() != removed_entries.find(combined_key_2)) {
							continue;
						}

						castle::StringUtils::tokenize(value, delim_double_underscore, data2);
						temp1.clear();
						temp2.clear();
						temp1.push_back(temp[0]);
						temp1.push_back(temp[1]);
						try {
							temp1.push_back(boost::lexical_cast<int64_t>(data2[1]));
							temp1.push_back(boost::lexical_cast<int64_t>(data2[2]));

							temp2.push_back(boost::lexical_cast<int64_t>(data1[4]));
							temp2.push_back(boost::lexical_cast<int64_t>(data1[5]));
							temp2.push_back(boost::lexical_cast<int64_t>(data2[4]));
							temp2.push_back(boost::lexical_cast<int64_t>(data2[5]));
						} catch(exception& ex) {
							cout << "write_inter_regions: here-2: " << combined_key_2 << ", " << value << "\n";
							continue;
						}

						if (covered(temp1[0], temp1[1], temp1[2], temp1[3]) &&
						covered(temp2[0], temp2[1], temp2[2], temp2[3])) {
							merge = true;
							sort(temp1.begin(), temp1.end());
							sort(temp2.begin(), temp2.end());
							removed_entries.insert(combined_key_1);
							removed_entries.insert(combined_key_2);
							merged_keys.insert(key1);
							merged_keys.insert(key2);

							string newkey = key1 + "/" + key2;
							IntervalEntryType additional_entry;
							additional_entry.p_id = newkey;
							IntervalEntry entry(temp1[0], temp1.back(), additional_entry);
							IntervalEntry entry_second(temp2[0], temp2.back(), additional_entry);
							local_intervals[chr1][chr2].push_back(entry);
							local_intervals[chr1][chr2].push_back(entry_second);

							try {
							if (data1[6] == data2[6]) {
								inter_regions[chr1][chr2][newkey] = data1[0] + "__" + boost::lexical_cast<string>(temp1[0]) + "__"
								+ boost::lexical_cast<string>(temp1.back()) + "__" + data1[3] + "__" +
								boost::lexical_cast<string>(temp2[0]) + "__" + boost::lexical_cast<string>(temp2.back()) + "__"
								+ data1[6];
							} else {
								inter_regions[chr1][chr2][newkey] = data1[0] + "__"
								+ boost::lexical_cast<string>(temp1[0]) + "__"
								+ boost::lexical_cast<string>(temp1.back()) + "__"
								+ data1[3] + "__"
								+ boost::lexical_cast<string>(temp2[0]) + "__"
								+ boost::lexical_cast<string>(temp2.back()) + "__10";
							}
							} catch(exception& ex) {
								cout << "write_inter_regions: here-3: " << combined_key_2 << ", " << key1_entry.second << "\n";
							}
							break;
						}
					}
				}
			}
//cout << "Chr " << chr1 << " and Chr " << chr2 << " is merged\n";
			for (auto& cl_entry : inter_regions[chr1][chr2]) {
				auto& key = cl_entry.first;
				if(merged_keys.end() != merged_keys.find(key)) {
					continue;
				}
				auto& value = cl_entry.second;
//for(auto& value: cl_entry.second) {
				castle::StringUtils::tokenize(value, delim_double_underscore, data);
				temp.clear();
				if(data.size() < 6) {
					temp.push_back(boost::lexical_cast<int64_t>(data[1]));
					temp.push_back(boost::lexical_cast<int64_t>(data[2]));
					if (temp[0] < 0 || temp[1] < 0) {
						continue;
					}
					string seq = fai.fetch_1_system(data[0], temp[0], temp[1]);
					BPDT << (boost::format("%s\t%s\n") % key % value).str();
					if(!seq.empty()) {
						FA << (boost::format(">%s\n%s\n") % value % seq).str();
					}
					continue;
				}

				try {
				temp.push_back(boost::lexical_cast<int64_t>(data[1]));
				temp.push_back(boost::lexical_cast<int64_t>(data[2]));
				temp.push_back(boost::lexical_cast<int64_t>(data[4]));
				temp.push_back(boost::lexical_cast<int64_t>(data[5]));
				} catch(exception& ex) {
					cout << "write_inter_regions: here-4: " << key << ", " << value << "\n";
					continue;
				}
				if (temp[0] < 0 || temp[1] < 0 || temp[2] < 0 || temp[3] < 0 ) {
					continue;
				}
// ref_id, start => end
				string seq = fai.fetch_1_system(data[0], temp[0], temp[1]);
				if (data.size() > 2) {
					seq += n;
					bool is_reversed = ("1" == data[data.size() - 1]);
					if (temp[3] < temp[2]) {
						if(is_reversed) {
							seq += fai.fetch_1_system(data[3], temp[2], temp[3]);
						} else {
							seq += fai.fetch_1_system(data[3], temp[3], temp[2]);
						}
					} else {
						if(is_reversed) {
							seq += fai.fetch_1_system(data[3], temp[3], temp[2]);
						} else {
							seq += fai.fetch_1_system(data[3], temp[2], temp[3]);
						}
					}
//					if (temp[3] < temp[2]) {
//						seq += fai.fetch_1_system(data[3], temp[3], temp[2]);
//					} else {
//						seq += fai.fetch_1_system(data[3], temp[2], temp[3]);
//					}
				}
				BPDT << (boost::format("%s\t%s\n") % key % value).str();
				if(!seq.empty()) {
					FA << (boost::format(">%s\n%s\n") % value % seq).str();
				}
//}
			}
		}
	}
	cout << checker;
}

void SplitReadReAlign::prepare_split_alignment() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.prepare_split_alignment");
	checker.start();
	string fafile = options.prefix + ".bp.fasta";
	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
	string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
	string sr_rawsam = options.prefix + ".sr.raw.sam";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
		sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
	}

	string bwa_index_cmd = (boost::format("bwa index %s") % fafile).str();
	system(bwa_index_cmd.c_str());

	vector<function<void()> > tasks;

	int64_t sr_insertsize = 2 * options.cut_sr + 100;
	BWACaller bc;
	bc.set_n_cores(n_cores);
	n_blocks = bc.split_FASTQ_alt(sr_fq1, sr_fq2, 131072);
//	n_blocks = bc.split_FASTQ_alt(sr_fq1, sr_fq2);
	string aln_param = (boost::format("-l %d -k 1") % options.cut_sr).str();
	string sampe_param = (boost::format("-a %d -P -N 20") % sr_insertsize).str();

	//			string bwa_aln_cmd_1 = (boost::format("bwa aln %s -l %d -k 1 -t 1 %s > %s 2>%s") % fafile % options.cut_sr % sr_fq_block_1
	//					% sr_sai1 % bwa_err_file).str();
	//			string bwa_sampe_cmd = (boost::format("bwa sampe -a %d -P -N 20 %s %s %s %s %s >%s 2>>%s") % sr_insertsize % fafile % sr_sai1 % sr_sai2 % sr_fq_block_1 % sr_fq_block_2
	//					% sr_rawsam % bwa_err_file).str();

	bc.collect_align_tasks_alt(tasks, sr_rawsam, fafile, aln_param, sampe_param, sr_fq1, sr_fq2, n_blocks);
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id] {
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			int64_t sr_insertsize = 2 * options.cut_sr + 100;
//			string sr_fq_block_1 = options.prefix + ".sr.1.fq." + str_block_id;
//			string sr_fq_block_2 = options.prefix + ".sr.2.fq." + str_block_id;
//			string bwa_err_file = options.prefix + ".bwa.err." + str_block_id;
//			string sr_sai1 = options.prefix + ".sr.1.sai." + str_block_id;
//			string sr_sai2 = options.prefix + ".sr.2.sai." + str_block_id;
//			string sr_rawsam = options.prefix + ".sr.raw.sam." + str_block_id;
//
//			if(!options.working_dir.empty()) {
//				sr_fq_block_1 = options.working_prefix + ".sr.1.fq." + str_block_id;
//				sr_fq_block_2 = options.working_prefix + ".sr.2.fq." + str_block_id;
//				bwa_err_file = options.working_prefix + ".bwa.err." + str_block_id;
//				sr_sai1 = options.working_prefix + ".sr.1.sai." + str_block_id;
//				sr_sai2 = options.working_prefix + ".sr.2.sai." + str_block_id;
//				sr_rawsam = options.working_prefix + ".sr.raw.sam." + str_block_id;
//			}
//
//			string bwa_aln_cmd_1 = (boost::format("bwa aln %s -l %d -k 1 -t 1 %s > %s 2>%s") % fafile % options.cut_sr % sr_fq_block_1
//					% sr_sai1 % bwa_err_file).str();
//			string bwa_aln_cmd_2 = (boost::format("bwa aln %s -l %d -k 1 -t 1 %s > %s 2>>%s") % fafile % options.cut_sr % sr_fq_block_2
//					% sr_sai2 % bwa_err_file).str();
//			string bwa_sampe_cmd = (boost::format("bwa sampe -a %d -P -N 20 %s %s %s %s %s >%s 2>>%s") % sr_insertsize % fafile % sr_sai1 % sr_sai2 % sr_fq_block_1 % sr_fq_block_2
//					% sr_rawsam % bwa_err_file).str();
//			system(bwa_aln_cmd_1.c_str());
//			system(bwa_aln_cmd_2.c_str());
//			system(bwa_sampe_cmd.c_str());
//		});
//	}

//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id] {
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			int64_t sr_insertsize = 2 * options.cut_sr + 100;
//			string sr_fq_block_1 = options.prefix + ".sr.1.fq." + str_block_id;
//			string sr_fq_block_2 = options.prefix + ".sr.2.fq." + str_block_id;
//			string bwa_err_file = options.prefix + ".bwa.err." + str_block_id;
//			string sr_sai1 = options.prefix + ".sr.1.sai." + str_block_id;
//			string sr_sai2 = options.prefix + ".sr.2.sai." + str_block_id;
//			string sr_rawsam = options.prefix + ".sr.raw.sam." + str_block_id;
//
//			if(!options.working_dir.empty()) {
//				sr_fq_block_1 = options.working_prefix + ".sr.1.fq." + str_block_id;
//				sr_fq_block_2 = options.working_prefix + ".sr.2.fq." + str_block_id;
//				bwa_err_file = options.working_prefix + ".bwa.err." + str_block_id;
//				sr_sai1 = options.working_prefix + ".sr.1.sai." + str_block_id;
//				sr_sai2 = options.working_prefix + ".sr.2.sai." + str_block_id;
//				sr_rawsam = options.working_prefix + ".sr.raw.sam." + str_block_id;
//			}
//
//			string bwa_aln_cmd_1 = (boost::format("bwa aln %s -l %d -k 1 -t 1 %s > %s 2>%s") % fafile % options.cut_sr % sr_fq_block_1
//					% sr_sai1 % bwa_err_file).str();
//			string bwa_aln_cmd_2 = (boost::format("bwa aln %s -l %d -k 1 -t 1 %s > %s 2>>%s") % fafile % options.cut_sr % sr_fq_block_2
//					% sr_sai2 % bwa_err_file).str();
//			string bwa_sampe_cmd = (boost::format("bwa sampe -a %d -P -N 20 %s %s %s %s %s >%s 2>>%s") % sr_insertsize % fafile % sr_sai1 % sr_sai2 % sr_fq_block_1 % sr_fq_block_2
//					% sr_rawsam % bwa_err_file).str();
//			system(bwa_aln_cmd_1.c_str());
//			system(bwa_aln_cmd_2.c_str());
//			system(bwa_sampe_cmd.c_str());
//		});
//	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void SplitReadReAlign::prepare_split_alignment_single() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.prepare_split_alignment_serial");
	checker.start();
	string fafile = options.prefix + ".bp.fasta";
	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
	string sr_sai1 = options.prefix + ".sr.1.fq.gz.bak.sai";
	string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
	string sr_sai2 = options.prefix + ".sr.2.fq.gz.bak.sai";
	string sr_rawsam = options.prefix + ".sr.raw.sam";
	string bwa_err_file = options.prefix + ".sr.bwa.err";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
		sr_sai1 = options.working_prefix + ".sr.1.fq.gz.bak.sai";
		sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
		sr_sai2 = options.working_prefix + ".sr.2.fq.gz.bak.sai";
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
		bwa_err_file = options.working_prefix + ".sr.bwa.err";
	}

	int64_t sr_insertsize = 2 * options.cut_sr + 100;

	string bwa_index_cmd = (boost::format("bwa index %s") % fafile).str();
	system(bwa_index_cmd.c_str());

	string bwa_aln_cmd_1 = (boost::format("bwa aln %s -l %d -k 1 -t %d %s > %s 2>%s") % fafile % options.cut_sr % n_cores % sr_fq1
		% sr_sai1 % bwa_err_file).str();
	string bwa_aln_cmd_2 = (boost::format("bwa aln %s -l %d -k 1 -t %d %s > %s 2>>%s") % fafile % options.cut_sr % n_cores % sr_fq2
		% sr_sai2 % bwa_err_file).str();

	string bwa_sampe_cmd = (boost::format("bwa sampe -a %d -P -N 20 %s %s %s %s %s >%s 2>>%s") % sr_insertsize % fafile % sr_sai1 % sr_sai2 % sr_fq1 % sr_fq2
			% sr_rawsam % bwa_err_file).str();
	system(bwa_aln_cmd_1.c_str());
	system(bwa_aln_cmd_2.c_str());
	system(bwa_sampe_cmd.c_str());

	cout << checker;
}

void SplitReadReAlign::remove_and_merge_temporary_alignment_files() {
	vector<function<void()> > tasks;
	vector<string> removal_file_names;
	vector<string> merge_sam_file_names;
	vector<string> merge_bwa_err_file_names;
//	n_blocks = 187;
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string local_sr_fq_1 = options.prefix + ".sr.1.fq.gz.bak." + str_block_id;
		string local_sr_fq_2 = options.prefix + ".sr.2.fq.gz.bak." + str_block_id;
		string local_sr_raw_sai_1 = options.prefix + ".sr.1.fq.gz.bak.sai." + str_block_id;
		string local_sr_raw_sai_2 = options.prefix + ".sr.2.fq.gz.bak.sai." + str_block_id;
		string local_sr_raw_sam = options.prefix + ".sr.raw.sam." + str_block_id;
		string local_sr_raw_tmp_sam = options.prefix + ".sr.raw.tmp.sam." + str_block_id;
//		string local_sr_sam = options.prefix + ".sr.sam." + str_block_id;
		string local_bwa_err_1 = options.prefix + ".sr.1.fq.gz.bak.bwa.err." + str_block_id;
		string local_bwa_err_2 = options.prefix + ".sr.2.fq.gz.bak.bwa.err." + str_block_id;
		if(!options.working_dir.empty()) {
			local_sr_fq_1 = options.working_prefix + ".sr.1.fq.gz.bak." + str_block_id;
			local_sr_fq_2 = options.working_prefix + ".sr.2.fq.gz.bak." + str_block_id;
			local_sr_raw_sai_1 = options.working_prefix + ".sr.1.fq.gz.bak.sai." + str_block_id;
			local_sr_raw_sai_2 = options.working_prefix + ".sr.2.fq.gz.bak.sai." + str_block_id;
			local_sr_raw_sam = options.working_prefix + ".sr.raw.sam." + str_block_id;
			local_sr_raw_tmp_sam = options.working_prefix + ".sr.raw.tmp.sam." + str_block_id;
//			local_sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
			local_bwa_err_1 = options.working_prefix + ".sr.1.fq.gz.bak.bwa.err." + str_block_id;
			local_bwa_err_2 = options.working_prefix + ".sr.2.fq.gz.bak.bwa.err." + str_block_id;
		}

		merge_sam_file_names.push_back(local_sr_raw_tmp_sam);
		merge_bwa_err_file_names.push_back(local_bwa_err_1);
		merge_bwa_err_file_names.push_back(local_bwa_err_2);

		removal_file_names.push_back(local_sr_fq_1);
		removal_file_names.push_back(local_sr_fq_2);
		removal_file_names.push_back(local_sr_raw_sai_1);
		removal_file_names.push_back(local_sr_raw_sai_2);
		removal_file_names.push_back(local_sr_raw_tmp_sam);
		removal_file_names.push_back(local_sr_raw_sam);
	}

//	string sr_bam = options.prefix + ".sr.bam";
//	if(!options.working_dir.empty()) {
//		sr_bam = options.working_prefix + ".sr.bam";
//	}
//	string sambamba_cmd = (boost::format("sambamba merge -l 1 -t %d %s %s") % n_cores % sr_bam
//			% (castle::StringUtils::join(merge_bam_file_names, " "))).str();
//	system(sambamba_cmd.c_str());
	// remove SAM headers
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id]{
			string str_block_id = boost::lexical_cast<string>(block_id);
			string local_sr_raw_sam = options.prefix + ".sr.raw.sam." + str_block_id;
			string local_sr_out_sam = options.prefix + ".sr.raw.tmp.sam." + str_block_id;
			if(!options.working_dir.empty()) {
				local_sr_raw_sam = options.working_prefix + ".sr.raw.sam." + str_block_id;
				local_sr_out_sam = options.working_prefix + ".sr.raw.tmp.sam." + str_block_id;
			}
			{
				ifstream in_MDSAM(local_sr_raw_sam, ios::binary);
				ofstream out_MDSAM(local_sr_out_sam, ios::binary);
				string line;
				while (getline(in_MDSAM, line, '\n')) {
					if ('@' == line[0]) {
						if(0 != block_id) {
							continue;
						}
					}
					out_MDSAM << line << "\n";
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string out_sr_raw_sam = options.prefix + ".sr.raw.sam";
	string out_bwa_err_name = options.prefix + ".bwa.err";
	if(!options.working_dir.empty()) {
		out_sr_raw_sam = options.working_prefix + ".sr.raw.sam";
		out_bwa_err_name = options.working_prefix + ".bwa.err";
	}
	castle::IOUtils::plain_file_merge(out_sr_raw_sam, merge_sam_file_names, n_cores, true);
	castle::IOUtils::plain_file_merge(out_bwa_err_name, merge_bwa_err_file_names, n_cores, true);

	castle::IOUtils::remove_files(removal_file_names, n_cores);
}

void SplitReadReAlign::collect_misaligned_reads(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const map<string, map<int8_t, int32_t>>& bp_weight) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.collect_misaligned_reads");
	checker.start();
	vector<function<void()> > tasks;
//	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//	if(!options.working_dir.empty()) {
//		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
//	}
//#%alg_mis: mis-aligned reads
//#$alg_mis{readname}[0]: read id mis-aligned
//#$alg_mis{readname}[1]: read to be placed
//	n_blocks = 187;

	vector<map<string, map<int8_t, string>>> alg_mis_lists(n_blocks);
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
//cout << str_block_id << "\n";
				string sr_rawsam = options.prefix + ".sr.raw.sam." + str_block_id;
				if(!options.working_dir.empty()) {
					sr_rawsam = options.working_prefix + ".sr.raw.sam." + str_block_id;
				}
				auto& local_alg_mis = alg_mis_lists[block_id];
				vector<string> data;
				vector<string> first_cols;
				vector<string> xa_z_cols;
				vector<string> alt_map_cols;
				vector<string> p;
				vector<string> match_cols;
				const char* delim_tab = "\t";
				const char* delim_underscore = "_";
				const char* delim_semiconlon = ";";
				const char* delim_comma = ",";
				const string white_spaces = " \t\r\n";
				string the_xa_z_pattern = "XA:Z:";
				const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();
				string line;
				ifstream RAWSAM(sr_rawsam, ios::binary);

				while (getline(RAWSAM, line, '\n')) {
					if('@' == line[0]) {
						continue;
					}
					castle::StringUtils::c_string_multi_split(line, delim_tab, data);
					castle::StringUtils::c_string_multi_split(data[0], delim_underscore, first_cols);
//					const bool debug = "ST-E00104:502:HFJN5CCXX:3:2202:14286:70645_85" == data[0];
					const bool debug = false;
					if(debug) {
						cout << line << "\n";
					}
					string readlength_str = first_cols[first_cols.size() - 1];
					int64_t readlength = boost::lexical_cast<int64_t>(readlength_str);
//cout << line << "\n";
					auto the_xa_z_pos = line.find(the_xa_z_pattern);

					if ("=" != data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads] here-0\n";
						}
						string sr_alt = line.substr(the_xa_z_pos + the_xa_z_pattern_size);

						castle::StringUtils::c_string_multi_split(sr_alt, delim_semiconlon, xa_z_cols);
						string select_chr;
						int64_t select_p = -1;
						uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
						for(string& an_alt : xa_z_cols) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads] here-1: " << an_alt << "/" << data[6] << "\n";
							}
							auto the_mate_ref_name_pos = an_alt.find(data[6]);
							if (string::npos != the_mate_ref_name_pos) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads] here-2\n";
								}
								auto& current_match = an_alt;
								if(local_alg_mis[data[0]].end() == local_alg_mis[data[0]].find(1)) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-3\n";
									}
									local_alg_mis[data[0]][1] = current_match;
									//    is first mate?
									if (sam_flag & 0x40) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-4\n";
										}
										local_alg_mis[data[0]][0] = "1";
									}
									//# second in mate
									else if (sam_flag & 0x80) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-5\n";
										}
										local_alg_mis[data[0]][0] = "2";
									}
									castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);

									select_chr = alt_map_cols[0];
									select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-6: " << select_chr << "/" << select_p << "\n";
									}
									continue;
								}

								//TODO: check here
//								castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);
//								select_chr = alt_map_cols[0];
//								select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
								//TODO: check here
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads] here-7: " << select_chr << "/" << select_p << "\n";
								}
								castle::StringUtils::c_string_multi_split(select_chr, delim_underscore, p);
								vector<int64_t> temp_p;
								temp_p.push_back(boost::lexical_cast<int64_t>(p[1]));
								temp_p.push_back(boost::lexical_cast<int64_t>(p[2]));
								int64_t left_size = temp_p[1] - temp_p[0] + 1;
								if (p.size() > 5) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-6\n";
									}
									p[3] = p[4];
									p[4] = p[5];
									temp_p.push_back(boost::lexical_cast<int64_t>(p[3]));
									temp_p.push_back(boost::lexical_cast<int64_t>(p[4]));
									temp_p[3] = temp_p[3] - temp_p[2] + left_size + 100;
									temp_p[2] = left_size + 100;
								}

								temp_p[1] -= temp_p[0];
								temp_p[0] = 1;
								castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
								string current_p_str = match_cols[1].substr(1);
								int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
								int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);

								if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-7\n";
									}
									// first in mate
									if (sam_flag & 0x40) {
										local_alg_mis[data[0]][0] = "1";
									}
									// second in mate
									else if (sam_flag & 0x80) {
										local_alg_mis[data[0]][0] = "2";
									}
									local_alg_mis[data[0]][1] = current_match;
									select_p = current_p;
									break;
								}
								auto bp_weight_itr = bp_weight.find(data[6]);
								if(bp_weight.end() != bp_weight_itr) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-8\n";
									}
									auto bp_weight_first_itr = bp_weight_itr->second.find(0);
									if(bp_weight_itr->second.end() != bp_weight_first_itr) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-9\n";
										}
										int32_t bp_weight_index_1 = bp_weight_first_itr->second - 1;
										int64_t bp_weight_value_1 = temp_p[bp_weight_index_1];
										if (select_p < left_size && current_p < left_size) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads] here-10\n";
											}
											if (abs(current_p - bp_weight_value_1) < abs(select_p - bp_weight_value_1)) {
												if(debug) {
													cout << "[SplitReadReAlign.collect_misaligned_reads] here-11\n";
												}
												// first in mate
												if (sam_flag & 0x40) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-12\n";
													}
													local_alg_mis[data[0]][0] = "1";
												}
												// second in mate
												else if (sam_flag & 0x80) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-13\n";
													}
													local_alg_mis[data[0]][0] = "2";
												}
												local_alg_mis[data[0]][1] = current_match;
												select_p = current_p;
											}
										}
									}
									auto bp_weight_second_itr = bp_weight_itr->second.find(1);
									if(bp_weight_itr->second.end() != bp_weight_second_itr) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-14\n";
										}
										int32_t bp_weight_index_2 = bp_weight_second_itr->second - 1;
										int64_t bp_weight_value_2 = temp_p[bp_weight_index_2];
										if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads] here-15\n";
											}
											if (abs(current_p - bp_weight_value_2) < abs(select_p - bp_weight_value_2)) {
												if(debug) {
													cout << "[SplitReadReAlign.collect_misaligned_reads] here-16\n";
												}
												// first in mate
												if (sam_flag & 0x40) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-17\n";
													}
													local_alg_mis[data[0]][0] = "1";
												} // second in mate
												else if (sam_flag & 0x80) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-18\n";
													}
													local_alg_mis[data[0]][0] = "2";
												}
												local_alg_mis[data[0]][1] = current_match;
												select_p = current_p;
											}
										}
									}
								}
							}
						}
					} else if ("=" == data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
						string sr_alt_str = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads] here-19: " << sr_alt_str << "\n";
						}
						castle::StringUtils::c_string_multi_split(sr_alt_str, delim_semiconlon, xa_z_cols);
						int64_t select_p = boost::lexical_cast<int64_t>(data[3]);
						uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
						for(string& an_alt : xa_z_cols) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads] here-20: " << an_alt << "\n";
							}
							auto the_ref_name_pos = an_alt.find(data[2]);
							if (string::npos != the_ref_name_pos) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads] here-21\n";
								}
								auto current_match = an_alt;
								castle::StringUtils::c_string_multi_split(data[2], delim_underscore, p);
								vector<int64_t> temp_p;
								temp_p.push_back(boost::lexical_cast<int64_t>(p[1]));
								temp_p.push_back(boost::lexical_cast<int64_t>(p[2]));
								int64_t left_size = temp_p[1] - temp_p[0] + 1;
								if (p.size() > 5) {
									p[3] = p[4];
									p[4] = p[5];
									temp_p.push_back(boost::lexical_cast<int64_t>(p[3]));
									temp_p.push_back(boost::lexical_cast<int64_t>(p[4]));
									temp_p[3] = temp_p[3] - temp_p[2] + left_size + 100;
									temp_p[2] = left_size + 100;
								}
								temp_p[1] -= temp_p[0];
								temp_p[0] = 1;
								if (p.size() > 5) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-22-a: temp_p: " << temp_p[0] << "/" << temp_p[1] << "/" << temp_p[2] << "/" << temp_p[3] << "\n";
									}
								}

								castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
								string current_p_str = match_cols[1].substr(1);
								int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
								int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);

								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads] here-22-b: " <<
											select_p << "<-?" << current_p << "/" << mate_pos << "\n";
								}

								if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-23\n";
									}
									// first in mate
									if (sam_flag & 0x40) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-24\n";
										}
										local_alg_mis[data[0]][0] = "1";
									} // second in mate
									else if (sam_flag & 0x80) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-25\n";
										}
										local_alg_mis[data[0]][0] = "2";
									}
									local_alg_mis[data[0]][1] = current_match;
									select_p = current_p;
									break;
								}
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads] here-26-a: " << data[2] << "\n";
								}
								auto bp_weight_itr = bp_weight.find(data[2]);
								if(bp_weight.end() != bp_weight_itr) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads] here-26-b\n";
									}
									auto bp_weight_first_itr = bp_weight_itr->second.find(0);
									if(bp_weight_itr->second.end() != bp_weight_first_itr) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-27\n";
										}
										int32_t bp_weight_index_1 = bp_weight_first_itr->second - 1;
										int64_t bp_weight_value_1 = temp_p[bp_weight_index_1];
										if (select_p < left_size && current_p < left_size) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads] here-28\n";
											}
											if(abs(current_p - bp_weight_value_1) < abs(select_p - bp_weight_value_1)) {
												// # first in mate
												if(sam_flag & 0x40) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-30\n";
													}
													local_alg_mis[data[0]][0] = "1";
												}
												// # second in mate
												else if (sam_flag & 0x80) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-31\n";
													}
													local_alg_mis[data[0]][0] = "2";
												}
												local_alg_mis[data[0]][1] = current_match;
												select_p = current_p;
												if(debug) {
													cout << "[SplitReadReAlign.collect_misaligned_reads] here-29 selected: " << select_p << ", " << current_match << "\n";
												}
											}
										}
									}
									auto bp_weight_second_itr = bp_weight_itr->second.find(1);
									if(bp_weight_itr->second.end() != bp_weight_second_itr) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads] here-32: " << select_p << "<-?" << current_p << "/" << left_size << "\n";
										}
										int32_t bp_weight_index_2 = bp_weight_second_itr->second - 1;
										int64_t bp_weight_value_2 = temp_p[bp_weight_index_2];
										if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads] here-33\n";
											}
											if (abs(current_p - bp_weight_value_2) < abs(select_p - bp_weight_value_2)) {
												if(debug) {
													cout << "[SplitReadReAlign.collect_misaligned_reads] here-34\n";
												}
												// # first in mate
												if (sam_flag & 0x40) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-35\n";
													}
													local_alg_mis[data[0]][0] = "1";
												} // second in mate
												else if (sam_flag & 0x80) {
													if(debug) {
														cout << "[SplitReadReAlign.collect_misaligned_reads] here-36\n";
													}
													local_alg_mis[data[0]][0] = "2";
												}
												local_alg_mis[data[0]][1] = current_match;
												select_p = current_p;
											}
										}
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
		auto& local_alg_mis = alg_mis_lists[block_id];
		for (auto& first_entry : local_alg_mis) {
			auto& first_key = first_entry.first;
			for (auto& second_entry : first_entry.second) {
				alg_mis[first_key][second_entry.first] = second_entry.second;
			}
		}
	}
	cout << checker;
}

void SplitReadReAlign::collect_misaligned_reads_serial(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const map<string, map<int8_t, int32_t>>& bp_weight) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.collect_misaligned_reads_serial");
	checker.start();
//	vector<function<void()> > tasks;
//	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//	if(!options.working_dir.empty()) {
//		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
//	}
//#%alg_mis: mis-aligned reads
//#$alg_mis{readname}[0]: read id mis-aligned
//#$alg_mis{readname}[1]: read to be placed
//	n_blocks = 187;

//	vector<map<string, map<int8_t, string>>> alg_mis_lists(n_blocks);
//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id] {
//			string str_block_id = boost::lexical_cast<string>(block_id);
//cout << str_block_id << "\n";
	int64_t n_processed = 0;
	int64_t n_processed_after_debug = 0;
		string sr_rawsam = options.prefix + ".sr.raw.sam";
		if(!options.working_dir.empty()) {
			sr_rawsam = options.working_prefix + ".sr.raw.sam";
		}
		auto& local_alg_mis = alg_mis;
		vector<string> data;
		vector<string> first_cols;
		vector<string> xa_z_cols;
		vector<string> alt_map_cols;
		vector<string> p;
		vector<string> match_cols;
		const char* delim_tab = "\t";
		const char* delim_underscore = "_";
		const char* delim_semiconlon = ";";
		const char* delim_comma = ",";
		const string white_spaces = " \t\r\n";
		string the_xa_z_pattern = "XA:Z:";
		const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();
		string line;
		ifstream RAWSAM(sr_rawsam, ios::binary);

		while (getline(RAWSAM, line, '\n')) {
			++n_processed;
			if(0 == (n_processed & 1048575)) {
				cout << (boost::format("[SplitReadReAlign.collect_misaligned_reads_serial] processed: %d\n") % n_processed).str();
			}
			if('@' == line[0]) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			castle::StringUtils::c_string_multi_split(data[0], delim_underscore, first_cols);
			const bool debug = "ST-E00104:502:HFJN5CCXX:4:1206:6421:54453_151" == data[0];
//			const bool debug = false;
			if(debug) {
				cout << line << "\n";
				n_processed_after_debug = 1;
			}
			if(n_processed_after_debug > 0) {
				++n_processed_after_debug;
				if(n_processed_after_debug > 10) {
					break;
				}
			}
			string readlength_str = first_cols[first_cols.size() - 1];
			int64_t readlength = boost::lexical_cast<int64_t>(readlength_str);
//cout << line << "\n";
			auto the_xa_z_pos = line.find(the_xa_z_pattern);

			if ("=" != data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-0\n";
				}
				string sr_alt = line.substr(the_xa_z_pos + the_xa_z_pattern_size);

				castle::StringUtils::c_string_multi_split(sr_alt, delim_semiconlon, xa_z_cols);
				string select_chr;
				int64_t select_p = -1;
				uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
				for(string& an_alt : xa_z_cols) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-1: " << an_alt << "/" << data[6] << "\n";
					}
					auto the_mate_ref_name_pos = an_alt.find(data[6]);
					if (string::npos != the_mate_ref_name_pos) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-2\n";
						}
						auto& current_match = an_alt;
						if(local_alg_mis[data[0]].end() == local_alg_mis[data[0]].find(1)) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-3\n";
							}
							local_alg_mis[data[0]][1] = current_match;
							//    is first mate?
							if (sam_flag & 0x40) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-4\n";
								}
								local_alg_mis[data[0]][0] = "1";
							}
							//# second in mate
							else if (sam_flag & 0x80) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-5\n";
								}
								local_alg_mis[data[0]][0] = "2";
							}
							castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);

							select_chr = alt_map_cols[0];
							select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-6: " << select_chr << "/" << select_p << "\n";
							}
							continue;
						}

						castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);
						select_chr = alt_map_cols[0];
						select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-7: " << select_chr << "/" << select_p << "\n";
						}
						castle::StringUtils::c_string_multi_split(select_chr, delim_underscore, p);
						vector<int64_t> temp_p;
						temp_p.push_back(boost::lexical_cast<int64_t>(p[1]));
						temp_p.push_back(boost::lexical_cast<int64_t>(p[2]));
						int64_t left_size = temp_p[1] - temp_p[0] + 1;
						if (p.size() > 5) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-6\n";
							}
							p[3] = p[4];
							p[4] = p[5];
							temp_p.push_back(boost::lexical_cast<int64_t>(p[3]));
							temp_p.push_back(boost::lexical_cast<int64_t>(p[4]));
							temp_p[3] = temp_p[3] - temp_p[2] + left_size + 100;
							temp_p[2] = left_size + 100;
						}

						temp_p[1] -= temp_p[0];
						temp_p[0] = 1;
						castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
						string current_p_str = match_cols[1].substr(1);
						int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
						int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);

						if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-7\n";
							}
							// first in mate
							if (sam_flag & 0x40) {
								local_alg_mis[data[0]][0] = "1";
							}
							// second in mate
							else if (sam_flag & 0x80) {
								local_alg_mis[data[0]][0] = "2";
							}
							local_alg_mis[data[0]][1] = current_match;
							select_p = current_p;
							break;
						}
						auto bp_weight_itr = bp_weight.find(data[6]);
						if(bp_weight.end() != bp_weight_itr) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-8\n";
							}
							auto bp_weight_first_itr = bp_weight_itr->second.find(0);
							if(bp_weight_itr->second.end() != bp_weight_first_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-9\n";
								}
								int32_t bp_weight_index_1 = bp_weight_first_itr->second - 1;
								int64_t bp_weight_value_1 = temp_p[bp_weight_index_1];
								if (select_p < left_size && current_p < left_size) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-10\n";
									}
									if (abs(current_p - bp_weight_value_1) < abs(select_p - bp_weight_value_1)) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-11\n";
										}
										// first in mate
										if (sam_flag & 0x40) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-12\n";
											}
											local_alg_mis[data[0]][0] = "1";
										}
										// second in mate
										else if (sam_flag & 0x80) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-13\n";
											}
											local_alg_mis[data[0]][0] = "2";
										}
										local_alg_mis[data[0]][1] = current_match;
										select_p = current_p;
									}
								}
							}
							auto bp_weight_second_itr = bp_weight_itr->second.find(1);
							if(bp_weight_itr->second.end() != bp_weight_second_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-14\n";
								}
								int32_t bp_weight_index_2 = bp_weight_second_itr->second - 1;
								int64_t bp_weight_value_2 = temp_p[bp_weight_index_2];
								if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-15\n";
									}
									if (abs(current_p - bp_weight_value_2) < abs(select_p - bp_weight_value_2)) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-16\n";
										}
										// first in mate
										if (sam_flag & 0x40) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-17\n";
											}
											local_alg_mis[data[0]][0] = "1";
										} // second in mate
										else if (sam_flag & 0x80) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-18\n";
											}
											local_alg_mis[data[0]][0] = "2";
										}
										local_alg_mis[data[0]][1] = current_match;
										select_p = current_p;
									}
								}
							}
						}
					}
				}
			} else if ("=" == data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
				string sr_alt_str = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-19: " << sr_alt_str << "\n";
				}
				castle::StringUtils::c_string_multi_split(sr_alt_str, delim_semiconlon, xa_z_cols);
				int64_t select_p = boost::lexical_cast<int64_t>(data[3]);
				uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
				for(string& an_alt : xa_z_cols) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-20: " << an_alt << "\n";
					}
					auto the_ref_name_pos = an_alt.find(data[2]);
					if (string::npos != the_ref_name_pos) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-21\n";
						}
						auto current_match = an_alt;
						castle::StringUtils::c_string_multi_split(data[2], delim_underscore, p);
						vector<int64_t> temp_p;
						temp_p.push_back(boost::lexical_cast<int64_t>(p[1]));
						temp_p.push_back(boost::lexical_cast<int64_t>(p[2]));
						int64_t left_size = temp_p[1] - temp_p[0] + 1;
						if (p.size() > 5) {
							p[3] = p[4];
							p[4] = p[5];
							temp_p.push_back(boost::lexical_cast<int64_t>(p[3]));
							temp_p.push_back(boost::lexical_cast<int64_t>(p[4]));
							temp_p[3] = temp_p[3] - temp_p[2] + left_size + 100;
							temp_p[2] = left_size + 100;
						}
						temp_p[1] -= temp_p[0];
						temp_p[0] = 1;
						if (p.size() > 5) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-22-a: temp_p: " << temp_p[0] << "/" << temp_p[1] << "/" << temp_p[2] << "/" << temp_p[3] << "\n";
							}
						}

						castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
						string current_p_str = match_cols[1].substr(1);
						int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
						int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);

						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-22-b: " <<
									select_p << "<-?" << current_p << "/" << mate_pos << "\n";
						}

						if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-23\n";
							}
							// first in mate
							if (sam_flag & 0x40) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-24\n";
								}
								local_alg_mis[data[0]][0] = "1";
							} // second in mate
							else if (sam_flag & 0x80) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-25\n";
								}
								local_alg_mis[data[0]][0] = "2";
							}
							local_alg_mis[data[0]][1] = current_match;
							select_p = current_p;
							break;
						}
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-26-a: " << data[2] << "\n";
						}
						auto bp_weight_itr = bp_weight.find(data[2]);
						if(bp_weight.end() != bp_weight_itr) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-26-b\n";
							}
							auto bp_weight_first_itr = bp_weight_itr->second.find(0);
							if(bp_weight_itr->second.end() != bp_weight_first_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-27\n";
								}
								int32_t bp_weight_index_1 = bp_weight_first_itr->second - 1;
								int64_t bp_weight_value_1 = temp_p[bp_weight_index_1];
								if (select_p < left_size && current_p < left_size) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-28\n";
									}
									if(abs(current_p - bp_weight_value_1) < abs(select_p - bp_weight_value_1)) {
										// # first in mate
										if(sam_flag & 0x40) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-30\n";
											}
											local_alg_mis[data[0]][0] = "1";
										}
										// # second in mate
										else if (sam_flag & 0x80) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-31\n";
											}
											local_alg_mis[data[0]][0] = "2";
										}
										local_alg_mis[data[0]][1] = current_match;
										select_p = current_p;
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-29 selected: " << select_p << ", " << current_match << "\n";
										}
									}
								}
							}
							auto bp_weight_second_itr = bp_weight_itr->second.find(1);
							if(bp_weight_itr->second.end() != bp_weight_second_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-32: " << select_p << "<-?" << current_p << "/" << left_size << "\n";
								}
								int32_t bp_weight_index_2 = bp_weight_second_itr->second - 1;
								int64_t bp_weight_value_2 = temp_p[bp_weight_index_2];
								if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-33\n";
									}
									if (abs(current_p - bp_weight_value_2) < abs(select_p - bp_weight_value_2)) {
										if(debug) {
											cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-34\n";
										}
										// # first in mate
										if (sam_flag & 0x40) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-35\n";
											}
											local_alg_mis[data[0]][0] = "1";
										} // second in mate
										else if (sam_flag & 0x80) {
											if(debug) {
												cout << "[SplitReadReAlign.collect_misaligned_reads_serial] here-36\n";
											}
											local_alg_mis[data[0]][0] = "2";
										}
										local_alg_mis[data[0]][1] = current_match;
										select_p = current_p;
									}
								}
							}
						}
					}
				}
			}

		}
//			});
//	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		auto& local_alg_mis = alg_mis_lists[block_id];
//		for (auto& first_entry : local_alg_mis) {
//			auto& first_key = first_entry.first;
//			for (auto& second_entry : first_entry.second) {
//				alg_mis[first_key][second_entry.first] = second_entry.second;
//			}
//		}
//	}
	cout << checker;
}

void SplitReadReAlign::collect_misaligned_reads_serial_alt(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.collect_misaligned_reads_serial_alt");
	checker.start();

	int64_t n_processed = 0;
	int64_t n_processed_after_debug = 0;
	string sr_rawsam = options.prefix + ".sr.raw.sam";
	if(!options.working_dir.empty()) {
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
	}
//	auto& local_alg_mis = alg_mis;
	vector<string> data;
	vector<string> first_cols;
	vector<string> xa_z_cols;
	vector<string> alt_map_cols;
	vector<string> p;
	vector<string> match_cols;
	const char* delim_tab = "\t";
	const char* delim_underscore = "_";
	const char* delim_semiconlon = ";";
	const char* delim_comma = ",";
	const string white_spaces = " \t\r\n";
	string the_xa_z_pattern = "XA:Z:";
	const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();
	string line;

	ifstream RAWSAM(sr_rawsam, ios::binary);
	while (getline(RAWSAM, line, '\n')) {
		++n_processed;
		if(0 == (n_processed & 1048575)) {
			cout << (boost::format("[SplitReadReAlign.collect_misaligned_reads_serial_alt] processed: %d\n") % n_processed).str();
		}
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);

//		const bool debug = "ST-E00104:502:HFJN5CCXX:4:1206:6421:54453_151" == data[0];
		const bool debug = false;
		if(debug) {
			cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] " << line << "\n";
			n_processed_after_debug = 1;
		}
		if(n_processed_after_debug > 0) {
			++n_processed_after_debug;
			if(n_processed_after_debug > 10) {
				break;
			}
		}
//		cout << line << "\n";

		castle::StringUtils::c_string_multi_split(data[0], delim_underscore, first_cols);
		string readlength_str = first_cols[first_cols.size() - 1];
		int64_t readlength = boost::lexical_cast<int64_t>(readlength_str);
		auto the_xa_z_pos = line.find(the_xa_z_pattern);
		if ("=" != data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
			if(debug) {
				cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-0\n";
			}
			string sr_alt = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
			castle::StringUtils::c_string_multi_split(sr_alt, delim_semiconlon, xa_z_cols);
			string select_chr;
			int64_t select_p = -1;
			uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
			for(string& an_alt : xa_z_cols) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-1: " << an_alt << "\n";
				}
//					#print "$newline\n";
				auto the_mate_ref_name_pos = an_alt.find(data[6]);
				if(string::npos == the_mate_ref_name_pos) {
					continue;
				}
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-2\n";
				}
				auto& current_match = an_alt;
				auto alt_map_itr = alg_mis.find(data[0]);
				if(alg_mis.end() == alt_map_itr) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-3\n";
					}
					alg_mis[data[0]][1] = current_match;
					//    is first mate?
					if (sam_flag & 0x40) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-4\n";
						}
						alg_mis[data[0]][0] = "1";
					}
					//# second in mate
					if (sam_flag & 0x80) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-5\n";
						}
						alg_mis[data[0]][0] = "2";
					}
					castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);

					select_chr = alt_map_cols[0];
					select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
					continue;
				}

				castle::StringUtils::c_string_multi_split(select_chr, delim_underscore, p);
				vector<int64_t> value_p;
				// dummy value
				value_p.push_back(0);
				if (p.size() > 5) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-6\n";
					}
					p[3] = p[4];
					p[4] = p[5];
				}
				if(p.size() > 4) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-7\n";
					}
					value_p.push_back(boost::lexical_cast<int64_t>(p[1]));
					value_p.push_back(boost::lexical_cast<int64_t>(p[2]));
					value_p.push_back(boost::lexical_cast<int64_t>(p[3]));
					value_p.push_back(boost::lexical_cast<int64_t>(p[4]));
				}
				int64_t left_size = value_p[2] - value_p[1] + 1;
				if(value_p.size() > 4) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-8\n";
					}
					value_p[4] = value_p[4] - value_p[3] + left_size + 100;
				}
				if(value_p.size() > 3) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-9\n";
					}
					value_p[3] = left_size + 100;
				}
				value_p[2] = value_p[2] - value_p[1];
				value_p[1] = 1;

				castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
				string current_p_str = match_cols[1].substr(1);
				int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
				int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);
				if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-10\n";
					}
					// first in mate
					if (sam_flag & 0x40) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-11\n";
						}
						alg_mis[data[0]][0] = "1";
					}
					// second in mate
					else if (sam_flag & 0x80) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-12\n";
						}
						alg_mis[data[0]][0] = "2";
					}
					alg_mis[data[0]][1] = current_match;
					select_p = current_p;
					break;
				}
				if (select_p < left_size && current_p < left_size) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-13\n";
					}
					int64_t cur_p_idx = 0;
					auto bp_itr = bp_weight.find(data[6]);
					if(bp_weight.end() != bp_itr) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-14\n";
						}
						auto weight_itr = bp_itr->second.find(0);
						if(bp_itr->second.end() != weight_itr) {
							cur_p_idx = weight_itr->second;
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-15: cur_p_idx " << cur_p_idx << "\n";
							}
						}
					}
					if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-16: sel_p: " << select_p << ", cur_p: " << current_p << ", p[cur_p]: " << value_p[cur_p_idx] << "\n";
						}
						// first in mate
						if (sam_flag & 0x40) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-17\n";
							}
							alg_mis[data[0]][0] = "1";
						}
						// second in mate
						if (sam_flag & 0x80) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-18\n";
							}
							alg_mis[data[0]][0] = "2";
						}
						alg_mis[data[0]][1] = current_match;
						select_p = current_p;
					}
				}
				if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-19\n";
					}
					int64_t cur_p_idx = 0;
					auto bp_itr = bp_weight.find(data[6]);
					if(bp_weight.end() != bp_itr) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-20\n";
						}
						auto weight_itr = bp_itr->second.find(1);
						if(bp_itr->second.end() != weight_itr) {
							cur_p_idx = weight_itr->second;
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-21: " << cur_p_idx << "\n";
							}
						}
					}
					if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-22\n";
						}
						// first in mate
						if (sam_flag & 0x40) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-23\n";
							}
							alg_mis[data[0]][0] = "1";
						}
						// second in mate
						if (sam_flag & 0x80) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-24\n";
							}
							alg_mis[data[0]][0] = "2";
						}
						alg_mis[data[0]][1] = current_match;
						select_p = current_p;
					}
				}
			}
//	#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
	} else if ("=" == data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
		if(debug) {
			cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-25\n";
		}
		string sr_alt_str = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
		castle::StringUtils::c_string_multi_split(sr_alt_str, delim_semiconlon, xa_z_cols);
		int64_t select_p = boost::lexical_cast<int64_t>(data[3]);
		uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
		for(string& an_alt : xa_z_cols) {
			if(debug) {
				cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-26: " << an_alt << "\n";
			}
			auto the_ref_name_pos = an_alt.find(data[2]);
			if (string::npos == the_ref_name_pos) {
				continue;
			}
			if(debug) {
				cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-27\n";
			}
//					#print "$newline\n";
			auto current_match = an_alt;
			castle::StringUtils::c_string_multi_split(data[2], delim_underscore, p);
			vector<int64_t> value_p;
			// dummy value
			value_p.push_back(0);
			if (p.size() > 5) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-28\n";
				}
				p[3] = p[4];
				p[4] = p[5];
			}
			if(p.size() > 4) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-29\n";
				}
				value_p.push_back(boost::lexical_cast<int64_t>(p[1]));
				value_p.push_back(boost::lexical_cast<int64_t>(p[2]));
				value_p.push_back(boost::lexical_cast<int64_t>(p[3]));
				value_p.push_back(boost::lexical_cast<int64_t>(p[4]));
			}
//					if (p.size() > 5) {
//						p[3] = p[4];
//						p[4] = p[5];
//					}
			int64_t left_size = value_p[2] - value_p[1] + 1;
			if(value_p.size() > 4) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-30\n";
				}
				value_p[4] = value_p[4] - value_p[3] + left_size + 100;
			}
			if(value_p.size() > 3) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-31\n";
				}
				value_p[3] = left_size + 100;
			}
			value_p[2] = value_p[2] - value_p[1];
			value_p[1] = 1;
//					int64_t left_size = p[2] - p[1] + 1;
//					p[4] = p[4] - p[3] + left_size + 100;
//					p[3] = left_size + 100;
//					p[2] = p[2] - p[1];
//					p[1] = 1;
			castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
			string current_p_str = match_cols[1].substr(1);
			int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
			int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);
			if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-32\n";
				}
				// first in mate
				if (sam_flag & 0x40) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-33\n";
					}
					alg_mis[data[0]][0] = "1";
				} // second in mate
				if (sam_flag & 0x80) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-34\n";
					}
					alg_mis[data[0]][0] = "2";
				}
				alg_mis[data[0]][1] = current_match;
				select_p = current_p;
				break;
				}
				if (select_p < left_size && current_p < left_size) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-35\n";
					}
					int64_t cur_p_idx = 0;
					auto bp_itr = bp_weight.find(data[2]);
					if(bp_weight.end() != bp_itr) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-36\n";
						}
						auto weight_itr = bp_itr->second.find(0);
						if(bp_itr->second.end() != weight_itr) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-37\n";
							}
							cur_p_idx = weight_itr->second;
						}
					}
					if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx]) ) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-38\n";
						}
						// # first in mate
						if(sam_flag & 0x40) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-39\n";
							}
							alg_mis[data[0]][0] = "1";
						}
						// # second in mate
						if (sam_flag & 0x80) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-40\n";
							}
							alg_mis[data[0]][0] = "2";
						}
						alg_mis[data[0]][1] = current_match;
						select_p = current_p;
					}
				}
				if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-41\n";
					}
					int64_t cur_p_idx = 0;
					auto bp_itr = bp_weight.find(data[2]);
					if(bp_weight.end() != bp_itr) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-42\n";
						}
						auto weight_itr = bp_itr->second.find(1);
						if(bp_itr->second.end() != weight_itr) {
							cur_p_idx = weight_itr->second;
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-43: " << cur_p_idx << "\n";
							}
						}
					}
					if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-44\n";
						}
//#print "$data[2]\t$data[0]\t$select_p\t$current_p\t@p\t$p[$bp_weight{$data[2]}[1]]\n";
						// # first in mate
						if (sam_flag & 0x40) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-45\n";
							}
							alg_mis[data[0]][0] = "1";
						} // second in mate
						if (sam_flag & 0x80) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_serial_alt] here-46\n";
							}
							alg_mis[data[0]][0] = "2";
						}
						alg_mis[data[0]][1] = current_match;
						select_p = current_p;
					}
				}
			}
//	#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
		}
	}
	cout << (boost::format("[SplitReadReAlign.collect_misaligned_reads_serial_alt] # entries: %d\n") % alg_mis.size()).str();
	cout << checker;
}

void SplitReadReAlign::collect_misaligned_reads_par(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.collect_misaligned_reads_par");
	checker.start();


	string sr_rawsam = options.prefix + ".sr.raw.sam";
	if(!options.working_dir.empty()) {
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
	}

	int64_t max_index = castle::IOUtils::get_file_size(sr_rawsam);
	int64_t chunk_size = max_index / (double) n_cores;
	int64_t local_n_blocks = max_index / (double) chunk_size;

	vector<uint64_t> skip_points;

	vector<boost::unordered_map<string, map<int8_t, string>>> alg_mis_lists(local_n_blocks);

//	castle::IOUtils::find_skip_points(skip_points, sr_rawsam, chunk_size, max_index, local_n_blocks, n_cores);
	{
		skip_points.resize(local_n_blocks + 1);

		vector<function<void()> > tasks;
		for (uint32_t id = 0; id < local_n_blocks - 1; ++id) {
			tasks.push_back([&, id, local_n_blocks, chunk_size, max_index] {
				ifstream input_stream(sr_rawsam, ios::binary);
				int64_t cur_pos = castle::IOUtils::get_next_end_pos(id, local_n_blocks - 1, chunk_size, max_index, 0);
				input_stream.seekg(cur_pos, ios::beg);
				string line;
				vector<string> data;
				const char * delim_tab = "\t";
				getline(input_stream, line, '\n');
				cur_pos += line.size() + 1;
				string prev_id;
				string cur_id;
				while(getline(input_stream, line, '\n')) {
					castle::StringUtils::c_string_multi_split(line, delim_tab, data);
					cur_id = data[0];
					if(!prev_id.empty() && cur_id != prev_id) {
						break;
					}
					cur_pos += line.size() + 1;
					prev_id = cur_id;
				}
				skip_points[id + 1] = cur_pos;
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		skip_points[0] = 0;
		skip_points[local_n_blocks] = max_index;
	}
	vector<function<void()> > tasks;

	for (uint32_t block_id = 0; block_id < local_n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
//			int64_t n_processed = 0;
//			int64_t n_processed_after_debug = 0;
			int64_t cur_pos = skip_points[block_id];
			int64_t next_pos = skip_points[block_id + 1];
			auto& local_alg_mis = alg_mis_lists[block_id];

			vector<string> data;
			vector<string> first_cols;
			vector<string> xa_z_cols;
			vector<string> alt_map_cols;
			vector<string> p;
			vector<string> match_cols;
			const char* delim_tab = "\t";
			const char* delim_underscore = "_";
			const char* delim_semiconlon = ";";
			const char* delim_comma = ",";
			const string white_spaces = " \t\r\n";
			string the_xa_z_pattern = "XA:Z:";
			const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();
			string line;
			ifstream in(sr_rawsam, ios::binary);
			if(0 != block_id) {
				in.seekg(cur_pos, ios::beg);
			}
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;

				if('@' == line[0]) {
					continue;
				}
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);

		//		const bool debug = "ST-E00104:502:HFJN5CCXX:4:1206:6421:54453_151" == data[0];
				const bool debug = false;
				if(debug) {
					cout << "[SplitReadReAlign.collect_misaligned_reads_par] " << line << "\n";
//					n_processed_after_debug = 1;
				}
//				if(n_processed_after_debug > 0) {
//					++n_processed_after_debug;
//					if(n_processed_after_debug > 10) {
//						break;
//					}
//				}
		//		cout << line << "\n";

				castle::StringUtils::c_string_multi_split(data[0], delim_underscore, first_cols);
				string readlength_str = first_cols[first_cols.size() - 1];
				int64_t readlength = boost::lexical_cast<int64_t>(readlength_str);
				auto the_xa_z_pos = line.find(the_xa_z_pattern);
				if ("=" != data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-0\n";
					}
					string sr_alt = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
					castle::StringUtils::c_string_multi_split(sr_alt, delim_semiconlon, xa_z_cols);
					string select_chr;
					int64_t select_p = -1;
					uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
					for(string& an_alt : xa_z_cols) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-1: " << an_alt << "\n";
						}
		//					#print "$newline\n";
						auto the_mate_ref_name_pos = an_alt.find(data[6]);
						if (string::npos == the_mate_ref_name_pos) {
							continue;
						}
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-2\n";
						}
						auto& current_match = an_alt;
						if(local_alg_mis.end() == local_alg_mis.find(data[0])) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-3\n";
							}
							local_alg_mis[data[0]][1] = current_match;
							//    is first mate?
							if (sam_flag & 0x40) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-4\n";
								}
								local_alg_mis[data[0]][0] = "1";
							}
							//# second in mate
							if (sam_flag & 0x80) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-5\n";
								}
								local_alg_mis[data[0]][0] = "2";
							}
							castle::StringUtils::c_string_multi_split(an_alt, delim_comma, alt_map_cols);

							select_chr = alt_map_cols[0];
							select_p = boost::lexical_cast<int64_t>(alt_map_cols[1].substr(1));
							continue;
						}
						castle::StringUtils::c_string_multi_split(select_chr, delim_underscore, p);
						vector<int64_t> value_p;
						// dummy value
						value_p.push_back(0);
						if (p.size() > 5) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-6\n";
							}
							p[3] = p[4];
							p[4] = p[5];
						}
						if(p.size() > 4) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-7\n";
							}
							value_p.push_back(boost::lexical_cast<int64_t>(p[1]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[2]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[3]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[4]));
						}
						int64_t left_size = value_p[2] - value_p[1] + 1;
						if(value_p.size() > 4) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-8\n";
							}
							value_p[4] = value_p[4] - value_p[3] + left_size + 100;
						}
						if(value_p.size() > 3) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-9\n";
							}
							value_p[3] = left_size + 100;
						}
						value_p[2] = value_p[2] - value_p[1];
						value_p[1] = 1;

						castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
						string current_p_str = match_cols[1].substr(1);
						int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
						int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);
						if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-10\n";
							}
							// first in mate
							if (sam_flag & 0x40) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-11\n";
								}
								local_alg_mis[data[0]][0] = "1";
							}
							// second in mate
							else if (sam_flag & 0x80) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-12\n";
								}
								local_alg_mis[data[0]][0] = "2";
							}
							local_alg_mis[data[0]][1] = current_match;
							select_p = current_p;
							break;
						}
						if (select_p < left_size && current_p < left_size) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-13\n";
							}
							int64_t cur_p_idx = 0;
							auto bp_itr = bp_weight.find(data[6]);
							if(bp_weight.end() != bp_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-14\n";
								}
								auto weight_itr = bp_itr->second.find(0);
								if(bp_itr->second.end() != weight_itr) {
									cur_p_idx = weight_itr->second;
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-15: cur_p_idx " << cur_p_idx << "\n";
									}
								}
							}
							if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-16: sel_p: " << select_p << ", cur_p: " << current_p << ", p[cur_p]: " << value_p[cur_p_idx] << "\n";
								}
								// first in mate
								if (sam_flag & 0x40) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-17\n";
									}
									local_alg_mis[data[0]][0] = "1";
								}
								// second in mate
								if (sam_flag & 0x80) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-18\n";
									}
									local_alg_mis[data[0]][0] = "2";
								}
								local_alg_mis[data[0]][1] = current_match;
								select_p = current_p;
							}
						}
						if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-19\n";
							}
							int64_t cur_p_idx = 0;
							auto bp_itr = bp_weight.find(data[6]);
							if(bp_weight.end() != bp_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-20\n";
								}
								auto weight_itr = bp_itr->second.find(1);
								if(bp_itr->second.end() != weight_itr) {
									cur_p_idx = weight_itr->second;
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-21: " << cur_p_idx << "\n";
									}
								}
							}
							if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-22\n";
								}
								// first in mate
								if (sam_flag & 0x40) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-23\n";
									}
									local_alg_mis[data[0]][0] = "1";
								}
								// second in mate
								if (sam_flag & 0x80) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-24\n";
									}
									local_alg_mis[data[0]][0] = "2";
								}
								local_alg_mis[data[0]][1] = current_match;
								select_p = current_p;
							}
						}
					}
				//	#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
				} else if ("=" == data[6] && data.size() > 11 && string::npos != the_xa_z_pos) {
					if(debug) {
						cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-25\n";
					}
					string sr_alt_str = line.substr(the_xa_z_pos + the_xa_z_pattern_size);
					castle::StringUtils::c_string_multi_split(sr_alt_str, delim_semiconlon, xa_z_cols);
					int64_t select_p = boost::lexical_cast<int64_t>(data[3]);
					uint64_t sam_flag = boost::lexical_cast<uint64_t>(data[1]);
					for(string& an_alt : xa_z_cols) {
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-26: " << an_alt << "\n";
						}
						auto the_ref_name_pos = an_alt.find(data[2]);
						if (string::npos == the_ref_name_pos) {
							continue;
						}
						if(debug) {
							cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-27\n";
						}
		//					#print "$newline\n";
						auto current_match = an_alt;
						castle::StringUtils::c_string_multi_split(data[2], delim_underscore, p);
						vector<int64_t> value_p;
						// dummy value
						value_p.push_back(0);
						if (p.size() > 5) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-28\n";
							}
							p[3] = p[4];
							p[4] = p[5];
						}
						if(p.size() > 4) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-29\n";
							}
							value_p.push_back(boost::lexical_cast<int64_t>(p[1]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[2]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[3]));
							value_p.push_back(boost::lexical_cast<int64_t>(p[4]));
						}
		//					if (p.size() > 5) {
		//						p[3] = p[4];
		//						p[4] = p[5];
		//					}
						int64_t left_size = value_p[2] - value_p[1] + 1;
						if(value_p.size() > 4) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-30\n";
							}
							value_p[4] = value_p[4] - value_p[3] + left_size + 100;
						}
						if(value_p.size() > 3) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-31\n";
							}
							value_p[3] = left_size + 100;
						}
						value_p[2] = value_p[2] - value_p[1];
						value_p[1] = 1;
		//					int64_t left_size = p[2] - p[1] + 1;
		//					p[4] = p[4] - p[3] + left_size + 100;
		//					p[3] = left_size + 100;
		//					p[2] = p[2] - p[1];
		//					p[1] = 1;
						castle::StringUtils::c_string_multi_split(current_match, delim_comma, match_cols);
						string current_p_str = match_cols[1].substr(1);
						int64_t current_p = boost::lexical_cast<int64_t>(current_p_str);
						int64_t mate_pos = boost::lexical_cast<int64_t>(data[7]);
						if (abs(abs(current_p - mate_pos) - (readlength - options.cut_sr)) < 2) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-32\n";
							}
							// first in mate
							if (sam_flag & 0x40) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-33\n";
								}
								local_alg_mis[data[0]][0] = "1";
							} // second in mate
							if (sam_flag & 0x80) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-34\n";
								}
								local_alg_mis[data[0]][0] = "2";
							}
							local_alg_mis[data[0]][1] = current_match;
							select_p = current_p;
							break;
						}
						if (select_p < left_size && current_p < left_size) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-35\n";
							}
							int64_t cur_p_idx = 0;
							auto bp_itr = bp_weight.find(data[2]);
							if(bp_weight.end() != bp_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-36\n";
								}
								auto weight_itr = bp_itr->second.find(0);
								if(bp_itr->second.end() != weight_itr) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-37\n";
									}
									cur_p_idx = weight_itr->second;
								}
							}
							if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx]) ) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-38\n";
								}
								// # first in mate
								if(sam_flag & 0x40) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-39\n";
									}
									local_alg_mis[data[0]][0] = "1";
								}
								// # second in mate
								if (sam_flag & 0x80) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-40\n";
									}
									local_alg_mis[data[0]][0] = "2";
								}
								local_alg_mis[data[0]][1] = current_match;
								select_p = current_p;
							}
						}
						if (select_p > (left_size + 100) && current_p > (left_size + 100)) {
							if(debug) {
								cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-41\n";
							}
							int64_t cur_p_idx = 0;
							auto bp_itr = bp_weight.find(data[2]);
							if(bp_weight.end() != bp_itr) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-42\n";
								}
								auto weight_itr = bp_itr->second.find(1);
								if(bp_itr->second.end() != weight_itr) {
									cur_p_idx = weight_itr->second;
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-43: " << cur_p_idx << "\n";
									}
								}
							}
							if (abs(current_p - value_p[cur_p_idx]) < abs(select_p - value_p[cur_p_idx])) {
								if(debug) {
									cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-44\n";
								}
	//#print "$data[2]\t$data[0]\t$select_p\t$current_p\t@p\t$p[$bp_weight{$data[2]}[1]]\n";
								// # first in mate
								if (sam_flag & 0x40) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-45\n";
									}
									local_alg_mis[data[0]][0] = "1";
								} // second in mate
								if (sam_flag & 0x80) {
									if(debug) {
										cout << "[SplitReadReAlign.collect_misaligned_reads_par] here-46\n";
									}
									local_alg_mis[data[0]][0] = "2";
								}
								local_alg_mis[data[0]][1] = current_match;
								select_p = current_p;
							}
						}
					}
		//	#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	//vector<boost::unordered_map<string, map<int8_t, string>>> alg_mis_lists(local_n_blocks);
	// should overwrite previous values
	for (uint32_t block_id = 0; block_id < local_n_blocks; ++block_id) {
		auto& local_alg_mis = alg_mis_lists[block_id];
		for(auto& read_name_entry: local_alg_mis) {
			auto& a_read_name = read_name_entry.first;
			if(alg_mis.end() != alg_mis.find(a_read_name)) {
				cout << "exists: " << a_read_name << "\n";
			}
			for(auto& match_entry : read_name_entry.second) {
				int64_t the_value_id = match_entry.first;
				// the_value_id[0]: pair number
				// the_value_id[1]: alternative match
				auto& the_value = match_entry.second;
				alg_mis[a_read_name][the_value_id] = the_value;
			}
		}
	}

	cout << (boost::format("[SplitReadReAlign.collect_misaligned_reads_par] # entries: %d\n") % alg_mis.size()).str();
	cout << checker;
}

void SplitReadReAlign::adjust_misaligned_reads(boost::unordered_map<string, map<int8_t, string>>& alg_mis) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.adjust_misaligned_reads");
	checker.start();
	vector<function<void()> > tasks;
//	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//	int64_t file_size = castle::IOUtils::get_file_size(sr_fq1);
//	int64_t n_blocks = file_size / (double) BLOCK_SIZE;
	tasks.push_back([&]{
		string fafile = options.prefix + ".bp.fasta";
		if(!options.working_dir.empty()) {
			fafile = options.working_prefix + ".bp.fasta";
		}
		string samtools_faidx_cmd = (boost::format("samtools faidx %s") % fafile).str();
		system(samtools_faidx_cmd.c_str());
	});
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			vector<string> data;
			vector<string> data1;
			vector<string> alt_map_cols;
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			string line;

			string sr_rawsam = options.prefix + ".sr.raw.sam." + str_block_id;
			string sr_sam = options.prefix + ".sr.sam." + str_block_id;
			if(!options.working_dir.empty()) {
				sr_rawsam = options.working_prefix + ".sr.raw.sam." + str_block_id;
				sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
			}
			ifstream RAWSAM(sr_rawsam, ios::binary);
			ofstream MDSAM(sr_sam, ios::binary);
			while (getline(RAWSAM, line, '\n')) {
				if('@' == line[0]) {
					if(0 == block_id) {
						MDSAM << line << "\n";
					}
					continue;
				}


				castle::StringUtils::tokenize(line, delim_tab, data);
//				const bool debug = "ST-E00104:502:HFJN5CCXX:3:2202:14286:70645_85" == data[0];
				const bool debug = false;
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads] " << line << "\n";
				}
				int64_t sam_flag = boost::lexical_cast<int64_t>(data[1]);
				vector<int64_t> temp_values;

				if (alg_mis.end() != alg_mis.find(data[0])) {
					if(((sam_flag & 0x40) && "1" == alg_mis[data[0]][0]) || ((sam_flag & 0x80) && "2" == alg_mis[data[0]][0])) {
						if ((sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("+")) {
							sam_flag -= 16;
						} else if (!(sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("-")) {
							sam_flag += 16;
						}
						if(debug) {
							cout << "[SplitReadReAlign.adjust_misaligned_reads] data[0]: " << data[0] << "\n";
						}
						const auto an_alt_map = alg_mis[data[0]];
						const auto an_alt = an_alt_map.find(1);
						if (an_alt_map.end() != an_alt) {
							data[6] = "=";
							castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
							data[2] = alt_map_cols[0];
							data[3] = alt_map_cols[1].substr(1);
							temp_values.push_back(boost::lexical_cast<int64_t>(data[3]));
							temp_values.push_back(boost::lexical_cast<int64_t>(data[7]));
							data[8] = boost::lexical_cast<string>(temp_values[1] - temp_values[0]);
							if(debug) {
								cout << "[SplitReadReAlign.adjust_misaligned_reads] temp_values[0]: " << temp_values[0] << ", temp_values[1]: " << temp_values[1] << "\n";
							}
						}
						if(!(sam_flag & 0x2)) {
							sam_flag += 2;
						}
					}
//# modify mate
				else {
//#print "$data[0]\t2\t$data[1]\t$alg_mis{$data[0]}[0]\n";
					if ((sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("+")) {
						sam_flag -= 32;
					} else if (!(sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("-")) {
						sam_flag += 32;
					}
					const auto& an_alt_map = alg_mis[data[0]];
					const auto& an_alt = an_alt_map.find(1);
					if (an_alt_map.end() != an_alt) {
						data[6] = "=";
						castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
						data[7] = alt_map_cols[1].substr(1);
						temp_values.push_back(boost::lexical_cast<int64_t>(data[3]));
						temp_values.push_back(boost::lexical_cast<int64_t>(data[7]));
						data[8] = boost::lexical_cast<string>(temp_values[1] - temp_values[0]);
					}
					if(sam_flag & 0x2) {
						sam_flag -= 2;
					}
				}
				data[1] = boost::lexical_cast<string>(sam_flag);
				string the_result = castle::StringUtils::join(data, "\t");
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads] the result: " << the_result << "\n";
				}
				MDSAM << the_result << "\n";
			} else {
				MDSAM << line << "\n";
			}
		}
	});
	}

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	vector<string> removal_file_names;
	vector<string> merge_sam_file_names;
	vector<string> merge_bwa_err_file_names;
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string local_sr_fq_1 = options.prefix + ".sr.1.fq.gz.bak." + str_block_id;
		string local_sr_fq_2 = options.prefix + ".sr.2.fq.gz.bak." + str_block_id;
		string local_sr_raw_sai_1 = options.prefix + ".sr.1.fq.gz.bak.sai." + str_block_id;
		string local_sr_raw_sai_2 = options.prefix + ".sr.2.fq.gz.bak.sai." + str_block_id;
		string local_sr_raw_sam = options.prefix + ".sr.raw.sam." + str_block_id;
		string local_sr_sam = options.prefix + ".sr.sam." + str_block_id;
		string local_bwa_err_1 = options.prefix + ".sr.1.fq.gz.bak.bwa.err." + str_block_id;
		string local_bwa_err_2 = options.prefix + ".sr.2.fq.gz.bak.bwa.err." + str_block_id;
		if(!options.working_dir.empty()) {
			local_sr_fq_1 = options.working_prefix + ".sr.1.fq.gz.bak." + str_block_id;
			local_sr_fq_2 = options.working_prefix + ".sr.2.fq.gz.bak." + str_block_id;
			local_sr_raw_sai_1 = options.working_prefix + ".sr.1.fq.gz.bak.sai." + str_block_id;
			local_sr_raw_sai_2 = options.working_prefix + ".sr.2.fq.gz.bak.sai." + str_block_id;
			local_sr_raw_sam = options.working_prefix + ".sr.raw.sam." + str_block_id;
			local_sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
			local_bwa_err_1 = options.working_prefix + ".sr.1.fq.gz.bak.bwa.err." + str_block_id;
			local_bwa_err_2 = options.working_prefix + ".sr.2.fq.gz.bak.bwa.err." + str_block_id;
		}

		merge_sam_file_names.push_back(local_sr_sam);
		merge_bwa_err_file_names.push_back(local_bwa_err_1);
		merge_bwa_err_file_names.push_back(local_bwa_err_2);

		removal_file_names.push_back(local_sr_fq_1);
		removal_file_names.push_back(local_sr_fq_2);
		removal_file_names.push_back(local_sr_raw_sai_1);
		removal_file_names.push_back(local_sr_raw_sai_2);
//		removal_file_names.push_back(local_sr_raw_sam);
//		removal_file_names.push_back(local_sr_sam);
	}
//	vector<string> merge_bam_file_names;
//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		string str_block_id = boost::lexical_cast<string>(block_id);
//		string local_sr_bam = options.prefix + ".sr.bam." + str_block_id;
//		string local_sr_sorted_bam = options.prefix + ".sr.sorted.bam." + str_block_id;
//		if(!options.working_dir.empty()) {
//			local_sr_bam = options.working_prefix + ".sr.bam." + str_block_id;
//			local_sr_sorted_bam = options.working_prefix + ".sr.sorted.bam." + str_block_id;
//		}
//		removal_file_names.push_back(local_sr_bam);
//		removal_file_names.push_back(local_sr_sorted_bam);
//		merge_bam_file_names.push_back(local_sr_sorted_bam);
//	}
	// convert SAM to BAM
//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id]{
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			string local_sr_sam = options.prefix + ".sr.sam." + str_block_id;
//			string local_sr_bam = options.prefix + ".sr.bam." + str_block_id;
//			string local_sr_bam_sorted = options.prefix + ".sr.sorted.bam." + str_block_id;
//			string faifile = options.prefix + ".bp.fasta.fai";
//			if(!options.working_dir.empty()) {
//				local_sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
//				local_sr_bam = options.working_prefix + ".sr.bam." + str_block_id;
//				local_sr_bam_sorted = options.working_prefix + ".sr.sorted.bam." + str_block_id;
//				faifile = options.working_prefix + ".bp.fasta.fai";
//			}
//			string samtools_sam_to_bam_cmd = (boost::format("samtools view -bt %s %s -o %s") % faifile % local_sr_sam % local_sr_bam).str();
//			system(samtools_sam_to_bam_cmd.c_str());
//			string sambamba_sort_cmd = (boost::format("sambamba sort -o %s -l 1 -t 1 %s") % local_sr_bam_sorted % local_sr_bam).str();
//			system(sambamba_sort_cmd.c_str());
//		});
//		}
//
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	string sr_bam = options.prefix + ".sr.bam";
	if(!options.working_dir.empty()) {
		sr_bam = options.working_prefix + ".sr.bam";
	}
//	string sambamba_cmd = (boost::format("sambamba merge -l 1 -t %d %s %s") % n_cores % sr_bam
//			% (castle::StringUtils::join(merge_bam_file_names, " "))).str();
//	system(sambamba_cmd.c_str());
	// remove SAM headers
//	for (int64_t block_id = 1; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id]{
//			string str_block_id = boost::lexical_cast<string>(block_id);
//			string local_sr_sam = options.prefix + ".sr.sam." + str_block_id;
//			string local_sr_out_sam = options.prefix + ".sr.sam.tmp." + str_block_id;
//			if(!options.working_dir.empty()) {
//				local_sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
//				local_sr_out_sam = options.working_prefix + ".sr.sam.tmp." + str_block_id;
//			}
//			{
//				ifstream in_MDSAM(local_sr_sam, ios::binary);
//				ofstream out_MDSAM(local_sr_out_sam, ios::binary);
//				string line;
//				while (getline(in_MDSAM, line, '\n')) {
//					if ('@' == line[0]) {
//						continue;
//					}
//					out_MDSAM << line << "\n";
//				}
//			}
//			boost::filesystem::copy_file(local_sr_out_sam, local_sr_sam, boost::filesystem::copy_option::overwrite_if_exists);
//			boost::filesystem::remove(local_sr_out_sam);
//		});
//	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string out_sr_sam = options.prefix + ".sr.sam";
	string out_bwa_err_name = options.prefix + ".bwa.err";
	if(!options.working_dir.empty()) {
		out_sr_sam = options.working_prefix + ".sr.sam";
		out_bwa_err_name = options.working_prefix + ".bwa.err";
	}
	castle::IOUtils::plain_file_merge(out_sr_sam, merge_sam_file_names, n_cores, false);
	castle::IOUtils::plain_file_merge(out_bwa_err_name, merge_bwa_err_file_names, n_cores, false);

	castle::IOUtils::remove_files(removal_file_names, n_cores);

	cout << checker;
}

void SplitReadReAlign::adjust_misaligned_reads_serial(boost::unordered_map<string, map<int8_t, string>>& alg_mis) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.adjust_misaligned_reads_serial");
	checker.start();
//	vector<function<void()> > tasks;
//	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
//	int64_t file_size = castle::IOUtils::get_file_size(sr_fq1);
//	int64_t n_blocks = file_size / (double) BLOCK_SIZE;
	string fafile = options.prefix + ".bp.fasta";
	if(!options.working_dir.empty()) {
		fafile = options.working_prefix + ".bp.fasta";
	}
	string samtools_faidx_cmd = (boost::format("samtools faidx %s") % fafile).str();
	system(samtools_faidx_cmd.c_str());

//	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
//		tasks.push_back([&, block_id] {
//			string str_block_id = boost::lexical_cast<string>(block_id);
			vector<string> data;
			vector<string> data1;
			vector<string> alt_map_cols;
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			string line;

			string sr_rawsam = options.prefix + ".sr.raw.sam";
			string sr_sam = options.prefix + ".sr.sam";
			if(!options.working_dir.empty()) {
				sr_rawsam = options.working_prefix + ".sr.raw.sam";
				sr_sam = options.working_prefix + ".sr.sam";
			}
			int64_t n_processed = 0;
			int64_t n_processed_after_debug = 0;
			ifstream RAWSAM(sr_rawsam, ios::binary);
			ofstream MDSAM(sr_sam, ios::binary);
			while (getline(RAWSAM, line, '\n')) {
				++n_processed;
				if(0 == (n_processed & 1048575)) {
					cout << (boost::format("[SplitReadReAlign.adjust_misaligned_reads_serial] processed: %d\n") % n_processed).str();
				}
				if('@' == line[0]) {
					MDSAM << line << "\n";
					continue;
				}

				castle::StringUtils::tokenize(line, delim_tab, data);
//				const bool debug = "ST-E00104:502:HFJN5CCXX:6:2103:7395:60993_2_151" == data[0];
				const bool debug = false;
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads_serial] " << line << "\n";
					n_processed_after_debug = 1;
				}
				if(n_processed_after_debug > 0) {
					++n_processed_after_debug;
					if(n_processed_after_debug > 10) {
						break;
					}
				}
				int64_t sam_flag = boost::lexical_cast<int64_t>(data[1]);
				vector<int64_t> temp_values;

				if (alg_mis.end() != alg_mis.find(data[0])) {
					if(((sam_flag & 0x40) && "1" == alg_mis[data[0]][0]) || ((sam_flag & 0x80) && "2" == alg_mis[data[0]][0])) {
						if ((sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("+")) {
							sam_flag -= 16;
						} else if (!(sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("-")) {
							sam_flag += 16;
						}
						if(debug) {
							cout << "[SplitReadReAlign.adjust_misaligned_reads_serial] data[0]: " << data[0] << "\n";
						}
						const auto an_alt_map = alg_mis[data[0]];
						const auto an_alt = an_alt_map.find(1);
						if (an_alt_map.end() != an_alt) {
							data[6] = "=";
							castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
							data[2] = alt_map_cols[0];
							data[3] = alt_map_cols[1].substr(1);
							temp_values.push_back(boost::lexical_cast<int64_t>(data[3]));
							temp_values.push_back(boost::lexical_cast<int64_t>(data[7]));
							data[8] = boost::lexical_cast<string>(temp_values[1] - temp_values[0]);
							if(debug) {
								cout << "[SplitReadReAlign.adjust_misaligned_reads_serial] temp_values[0]: " << temp_values[0] << ", temp_values[1]: " << temp_values[1] << "\n";
							}
						}
						if(!(sam_flag & 0x2)) {
							sam_flag += 2;
						}
					}
//# modify mate
				else {
//#print "$data[0]\t2\t$data[1]\t$alg_mis{$data[0]}[0]\n";
					if ((sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("+")) {
						sam_flag -= 32;
					} else if (!(sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("-")) {
						sam_flag += 32;
					}
					const auto& an_alt_map = alg_mis[data[0]];
					const auto& an_alt = an_alt_map.find(1);
					if (an_alt_map.end() != an_alt) {
						data[6] = "=";
						castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
						data[7] = alt_map_cols[1].substr(1);
						temp_values.push_back(boost::lexical_cast<int64_t>(data[3]));
						temp_values.push_back(boost::lexical_cast<int64_t>(data[7]));
						data[8] = boost::lexical_cast<string>(temp_values[1] - temp_values[0]);
					}
					if(sam_flag & 0x2) {
						sam_flag -= 2;
					}
				}
				data[1] = boost::lexical_cast<string>(sam_flag);
				string the_result = castle::StringUtils::join(data, "\t");
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads_serial] the result: " << the_result << "\n";
				}
				MDSAM << the_result << "\n";
			} else {
				MDSAM << line << "\n";
			}
		}
//	});
//	}

	string sr_bam = options.prefix + ".sr.bam";
	if(!options.working_dir.empty()) {
		sr_bam = options.working_prefix + ".sr.bam";
	}
//	string sambamba_cmd = (boost::format("sambamba merge -l 1 -t %d %s %s") % n_cores % sr_bam
//			% (castle::StringUtils::join(merge_bam_file_names, " "))).str();
//	system(sambamba_cmd.c_str());
//	string out_sr_sam = options.prefix + ".sr.sam";
//	string out_bwa_err_name = options.prefix + ".bwa.err";
//	if(!options.working_dir.empty()) {
//		out_sr_sam = options.working_prefix + ".sr.sam";
//		out_bwa_err_name = options.working_prefix + ".bwa.err";
//	}
//	castle::IOUtils::plain_file_merge(out_sr_sam, merge_sam_file_names, n_cores, false);
//	castle::IOUtils::plain_file_merge(out_bwa_err_name, merge_bwa_err_file_names, n_cores, false);
//
//	castle::IOUtils::remove_files(removal_file_names, n_cores);

	cout << checker;

}

void SplitReadReAlign::adjust_misaligned_reads_serial_alt(boost::unordered_map<string, map<int8_t, string>>& alg_mis) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.adjust_misaligned_reads_serial_alt");
	checker.start();
//	auto cut_sr = options.cut_sr;
	const char* delim_tab = "\t";
//	auto& ref_is = options.is;
//	const char* delim_slash = "/";
//	const char* delim_colon = ":";
	vector<string> data;
	vector<string> cl;
//	const int64_t delta = ref_is["rlu"]["selected"] - options.cut_sr + 10;
	string line;
	string mpintra_outfile = options.prefix + ".mp.intra.out";
	string mpinter_outfile = options.prefix + ".mp.inter.out";
	string fafile = options.prefix + ".bp.fasta";
	string bpinfofile = options.prefix + ".bp.info";
	string sr_fq1 = options.prefix + ".sr.1.fq.gz.bak";
	string sr_fq2 = options.prefix + ".sr.2.fq.gz.bak";
	string sr_sai1 = options.prefix + ".sr.1.fq.gz.bak.sai";
	string sr_sai2 = options.prefix + ".sr.2.fq.gz.bak.sai";
	string sr_rawsam = options.prefix + ".sr.raw.sam";
	string sr_sam = options.prefix + ".sr.sam";
	string sr_bam = options.prefix + ".sr.bam";
	string sr_sort = options.prefix + ".sr.sorted";
	string sr_sortbam = options.prefix + ".sr.sorted.bam";
	if(!options.working_dir.empty()) {
		mpintra_outfile = options.working_prefix + ".mp.intra.out";
		mpinter_outfile = options.working_prefix + ".mp.inter.out";
		fafile = options.working_prefix + ".bp.fasta";
		bpinfofile = options.working_prefix + ".bp.info";
		fafile = options.working_prefix + ".bp.fasta";
		sr_fq1 = options.working_prefix + ".sr.1.fq.gz.bak";
		sr_fq2 = options.working_prefix + ".sr.2.fq.gz.bak";
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
		sr_sam = options.working_prefix + ".sr.sam";
		sr_bam = options.working_prefix + ".sr.bam";
		sr_sort = options.working_prefix + ".sr.sorted";
		sr_sortbam = options.working_prefix + ".sr.sorted.bam";
	}
//	auto sr_insertsize  = 2 * cut_sr + 100;

//	string n;
//	for (int64_t i = 0 ; i < 100 ; ++i) {
//		n += "N";
//	}

	vector<string> first_cols;
	vector<string> xa_z_cols;
	vector<string> alt_map_cols;
	vector<string> p;

//	vector<string> match_cols;
//	const char* delim_tab = "\t";
//	const char* delim_underscore = "_";
//	const char* delim_semiconlon = ";";
	const char* delim_comma = ",";
//	const string white_spaces = " \t\r\n";
//	string the_xa_z_pattern = "XA:Z:";
//	const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();

//	my ( $newline, $i, %bp_weight, %regions, %cluster_exist );

//#	%bp_weight: where is the real break point more likely to be located
//#	$bp_weight{regionname}[0]: left break point
//#	$bp_weight{regionname}[1]: right break point
//#	%regions, a list of existing regions
//#	$regions{$cluster_id}: region name, chr_p_p_chr_p_p_0/1, 0 same orientation, 1 opposite orientation
//#	%cluster_exist: the cluster with key id exists
//#	$cluster_exist{$cluster_id}: 1 exist, 0 not
	vector<string> positiona;
	vector<string> positionb;
	int64_t n_processed = 0;
	int64_t n_processed_after_debug = 0;

	ifstream RAWSAM(sr_rawsam, ios::binary);
	ofstream MDSAM(sr_sam, ios::binary);
	while (getline(RAWSAM, line, '\n')) {
		++n_processed;
		if(0 == (n_processed & 1048575)) {
			cout << (boost::format("[SplitReadReAlign.adjust_misaligned_reads_serial_alt] processed: %d\n") % n_processed).str();
		}
		if('@' == line[0]) {
			MDSAM << line << "\n";
			continue;
		}
		castle::StringUtils::tokenize(line, delim_tab, data);

//		const bool debug = "ST-E00104:502:HFJN5CCXX:4:1206:6421:54453_151" == data[0];
		const bool debug = false;
		if(debug) {
			cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] " << line << "\n";
			n_processed_after_debug = 1;
		}
		if(n_processed_after_debug > 0) {
			++n_processed_after_debug;
			if(n_processed_after_debug > 10) {
				break;
			}
		}

//		vector<string> data1 = data;

		if (alg_mis.end() != alg_mis.find(data[0])) {
			if(debug) {
				cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-0\n";
			}
			int64_t sam_flag = boost::lexical_cast<int64_t>(data[1]);
//			#print "$newline\n";
//			# modify read
			if(((sam_flag & 0x40) && "1" == alg_mis[data[0]][0]) || ((sam_flag & 0x80) && "2" == alg_mis[data[0]][0])) {
//#print "$data[0]\t1\t$data[1]\t$alg_mis{$data[0]}[0]\n";
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-1\n";
				}
				if ((sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("+")) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-2\n";
					}
					sam_flag -= 16;
				} else if (!(sam_flag & 0x10) && string::npos != alg_mis[data[0]][1].find("-")) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-3\n";
					}
					sam_flag += 16;
				}
				const auto an_alt_map = alg_mis[data[0]];
				const auto an_alt = an_alt_map.find(1);
//				if ( $alg_mis{ $data[0] }[1] =~ /(\S{0,100}?),(\S{0,100}?),(\S{0,10}?),(\S{0,10})/ )
				if (an_alt_map.end() != an_alt) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-4: " << an_alt->second << "\n";
					}
//#print "$1\t$2\t$3\t$4\n";
					data[6] = "=";
					castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
					data[2] = alt_map_cols[0];
					data[3] = alt_map_cols[1].substr(1);
					data[8] = boost::lexical_cast<string>(boost::lexical_cast<int64_t>(data[7]) - boost::lexical_cast<int64_t>(data[3]));
				}
				if(!(sam_flag & 0x2)) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-5\n";
					}
					sam_flag += 2;
				}
//			   #print "$data1[0]\t$data1[1]\t$data1[6]\t$data1[7]\t$data1[8]\n";
			}
//			# modify mate
			else {
				if(debug) {
					cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-6\n";
				}
//				#print "$data[0]\t2\t$data[1]\t$alg_mis{$data[0]}[0]\n";
				if ((sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("+")) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-7\n";
					}
					sam_flag -= 32;
				} else if (!(sam_flag & 0x20) && string::npos != alg_mis[data[0]][1].find("-")) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-8\n";
					}
					sam_flag += 32;
				}
				const auto& an_alt_map = alg_mis[data[0]];
				const auto& an_alt = an_alt_map.find(1);
				if (an_alt_map.end() != an_alt) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-9: " << an_alt->second << "\n";
					}
					data[6] = "=";
					castle::StringUtils::c_string_multi_split(an_alt->second, delim_comma, alt_map_cols);
					data[7] = alt_map_cols[1].substr(1);
					data[8] = boost::lexical_cast<string>(boost::lexical_cast<int64_t>(data[7]) - boost::lexical_cast<int64_t>(data[3]));
				}
				if(sam_flag & 0x2) {
					if(debug) {
						cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-10\n";
					}
					sam_flag -= 2;
				}
//				if ( $alg_mis{ $data[0] }[1] =~
//					/(\S{0,100}?),(\S{0,100}?),(\S{0,10}?),(\S{0,10})/ )
//				{
//
//					#print "$data[0]\t$1\t$2\t$3\t$4\n";
//					$data1[6] = '=';
//					$data1[7] = substr( $2, 1 );
//					$data1[8] = $data1[7] - $data1[3];
//				}
//				$data1[1] -= 2 if ( $data[1] & 0x2 );

//			   #print "$data1[0]\t$data1[1]\t$data1[6]\t$data1[7]\t$data1[8]\n";
			}
			data[1] = boost::lexical_cast<string>(sam_flag);
			string the_result = castle::StringUtils::join(data, "\t");
			if(debug) {
				cout << "[SplitReadReAlign.adjust_misaligned_reads_serial_alt] here-11: " << the_result << "\n";
			}
			MDSAM << the_result << "\n";
		} else {
			MDSAM << line << "\n";
		}
	}
	cout << checker;
}

void SplitReadReAlign::adjust_misaligned_reads_par(const boost::unordered_map<string, map<int8_t, string>>& alg_mis) {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.adjust_misaligned_reads_par");
	checker.start();

	string sr_rawsam = options.prefix + ".sr.raw.sam";
	if(!options.working_dir.empty()) {
		sr_rawsam = options.working_prefix + ".sr.raw.sam";
	}

	int64_t max_index = castle::IOUtils::get_file_size(sr_rawsam);
	int64_t chunk_size = max_index / (double) n_cores;
	int64_t local_n_blocks = max_index / (double) chunk_size;

	vector<uint64_t> skip_points;

	vector<string> out_filenames(local_n_blocks);
//	castle::IOUtils::find_skip_points(skip_points, sr_rawsam, chunk_size, max_index, local_n_blocks, n_cores);
	{
		skip_points.resize(local_n_blocks + 1);

		vector<function<void()> > tasks;
		for (uint32_t id = 0; id < local_n_blocks - 1; ++id) {
			tasks.push_back([&, id, local_n_blocks, chunk_size, max_index] {
				ifstream input_stream(sr_rawsam, ios::binary);
				int64_t cur_pos = castle::IOUtils::get_next_end_pos(id, local_n_blocks - 1, chunk_size, max_index, 0);
				input_stream.seekg(cur_pos, ios::beg);
				string line;
				vector<string> data;
				const char * delim_tab = "\t";
				getline(input_stream, line, '\n');
				cur_pos += line.size() + 1;
				string prev_id;
				string cur_id;
				while(getline(input_stream, line, '\n')) {
					castle::StringUtils::c_string_multi_split(line, delim_tab, data);
					cur_id = data[0];
					if(!prev_id.empty() && cur_id != prev_id) {
						break;
					}
					cur_pos += line.size() + 1;
					prev_id = cur_id;
				}
				skip_points[id + 1] = cur_pos;
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		skip_points[0] = 0;
		skip_points[local_n_blocks] = max_index;
	}
	vector<function<void()> > tasks;

	for (uint32_t block_id = 0; block_id < local_n_blocks; ++block_id) {
		tasks.push_back([&, block_id] {
//			int64_t n_processed = 0;
//			int64_t n_processed_after_debug = 0;
			int64_t cur_pos = skip_points[block_id];
			int64_t next_pos = skip_points[block_id + 1];
//			auto& local_alg_mis = alg_mis_lists[block_id];

			vector<string> data;
			vector<string> first_cols;
			vector<string> xa_z_cols;
			vector<string> alt_map_cols;
			vector<string> p;
			vector<string> match_cols;
			const char* delim_tab = "\t";
//			const char* delim_underscore = "_";
//			const char* delim_semiconlon = ";";
			const char* delim_comma = ",";
			const string white_spaces = " \t\r\n";
			string the_xa_z_pattern = "XA:Z:";
//			const int64_t the_xa_z_pattern_size = the_xa_z_pattern.size();
			string line;
			string str_block_id = boost::lexical_cast<string>(block_id);
			string sr_sam = options.prefix + ".sr.sam." + str_block_id;
			if(!options.working_dir.empty()) {
				sr_sam = options.working_prefix + ".sr.sam." + str_block_id;
			}
			out_filenames[block_id] = sr_sam;

			ifstream RAWSAM(sr_rawsam, ios::binary);
			ofstream MDSAM(sr_sam, ios::binary);
			if(0 != block_id) {
				RAWSAM.seekg(cur_pos, ios::beg);
			}
			while(getline(RAWSAM, line, '\n')) {
				cur_pos += line.size() + 1;
				if('@' == line[0]) {
					if(0 == block_id) {
						MDSAM << line << "\n";
					}
					continue;
				}
				castle::StringUtils::tokenize(line, delim_tab, data);
				auto first_entry_itr = alg_mis.find(data[0]);
				if (alg_mis.end() != first_entry_itr) {
					int64_t sam_flag = boost::lexical_cast<int64_t>(data[1]);
		//			#print "$newline\n";
		//			# modify read
					const auto an_alt_map = first_entry_itr->second;
					auto pair_number_itr = an_alt_map.find(0);
					auto an_alt_itr = an_alt_map.find(1);
					if(an_alt_map.end() != pair_number_itr && an_alt_map.end() != an_alt_itr) {
						if(((sam_flag & 0x40) && "1" == pair_number_itr->second) || ((sam_flag & 0x80) && "2" == pair_number_itr->second)) {
			//#print "$data[0]\t1\t$data[1]\t$alg_mis{$data[0]}[0]\n";
							if ((sam_flag & 0x10) && string::npos != an_alt_itr->second.find("+")) {
								sam_flag -= 16;
							} else if (!(sam_flag & 0x10) && string::npos != an_alt_itr->second.find("-")) {
								sam_flag += 16;
							}
//							const auto an_alt_map = alg_mis[data[0]];
//							const auto an_alt = an_alt_map.find(1);
			//				if ( $alg_mis{ $data[0] }[1] =~ /(\S{0,100}?),(\S{0,100}?),(\S{0,10}?),(\S{0,10})/ )
							if (an_alt_map.end() != an_alt_itr) {
			//#print "$1\t$2\t$3\t$4\n";
								data[6] = "=";
								castle::StringUtils::c_string_multi_split(an_alt_itr->second, delim_comma, alt_map_cols);
								data[2] = alt_map_cols[0];
								data[3] = alt_map_cols[1].substr(1);
								data[8] = boost::lexical_cast<string>(boost::lexical_cast<int64_t>(data[7]) - boost::lexical_cast<int64_t>(data[3]));
							}
							if(!(sam_flag & 0x2)) {
								sam_flag += 2;
							}
			//			   #print "$data1[0]\t$data1[1]\t$data1[6]\t$data1[7]\t$data1[8]\n";
						}
			//			# modify mate
						else {
			//				#print "$data[0]\t2\t$data[1]\t$alg_mis{$data[0]}[0]\n";
							if ((sam_flag & 0x20) && string::npos != an_alt_itr->second.find("+")) {
								sam_flag -= 32;
							} else if (!(sam_flag & 0x20) && string::npos != an_alt_itr->second.find("-")) {
								sam_flag += 32;
							}
//							const auto& an_alt_map = alg_mis[data[0]];
//							const auto& an_alt = an_alt_map.find(1);
							if (an_alt_map.end() != an_alt_itr) {
								data[6] = "=";
								castle::StringUtils::c_string_multi_split(an_alt_itr->second, delim_comma, alt_map_cols);
								data[7] = alt_map_cols[1].substr(1);
								data[8] = boost::lexical_cast<string>(boost::lexical_cast<int64_t>(data[7]) - boost::lexical_cast<int64_t>(data[3]));
							}
							if(sam_flag & 0x2) {
								sam_flag -= 2;
							}
						}
					}
					data[1] = boost::lexical_cast<string>(sam_flag);
					string the_result = castle::StringUtils::join(data, "\t");
					MDSAM << the_result << "\n";
				} else {
					MDSAM << line << "\n";
				}
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string sr_sam = options.prefix + ".sr.sam";
	if(!options.working_dir.empty()) {
		sr_sam = options.working_prefix + ".sr.sam";
	}

	castle::IOUtils::plain_file_merge(sr_sam, out_filenames, n_cores, true);
	cout << checker;
}

void SplitReadReAlign::create_misalignment_adjusted_bam() {
	castle::TimeChecker checker;
	checker.setTarget("SplitReadReAlign.create_misalignment_adjusted_bam");
	checker.start();
//	string fafile = options.prefix + ".bp.fasta";
//	string faifile = options.prefix + ".bp.fasta.fai";
	string sr_sam = options.prefix + ".sr.sam";
	string sr_bam = options.prefix + ".sr.bam";
//	string sr_sort = options.prefix + ".sr.sorted";
	string sr_sortbam = options.prefix + ".sr.sorted.bam";
	if(!options.working_dir.empty()) {
//		fafile = options.working_prefix + ".bp.fasta";
//		faifile = options.working_prefix + ".bp.fasta.fai";
		sr_sam = options.working_prefix + ".sr.sam";
		sr_bam = options.working_prefix + ".sr.bam";
//		sr_sort = options.working_prefix + ".sr.sorted";
		sr_sortbam = options.working_prefix + ".sr.sorted.bam";
	}
//	string samtools_faidx_cmd = (boost::format("samtools faidx %s") % fafile).str();
//	string samtools_sam_to_bam_cmd = (boost::format("samtools view -bt %s -o %s %s") % faifile % sr_bam % sr_sam).str();
//	string sam_to_bam_cmd = (boost::format("sambamba view -S -f bam -l 1 -t %d -o %s %s") % n_cores % sr_bam % sr_sam).str();
	string sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % sr_bam % sr_sam).str();
	system(sam_to_bam_cmd.c_str());

//	string samtools_sort_cmd = (boost::format("samtools sort -@ %d %s %s") % n_cores % sr_bam % sr_sort).str();
	string samtools_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % sr_sortbam % sr_bam).str();
//	system("sambamba sort --tmpdir=$is_actual_dir -N -l 1 -t $threads_bwa -o $disc_sort $disc_bam");
//	string samtools_index_cmd = (boost::format("samtools index %s") % sr_sortbam).str();

//	system(samtools_faidx_cmd.c_str());
	system(samtools_sort_cmd.c_str());
//	system(samtools_index_cmd.c_str());
	cout << checker;
}


} /* namespace meerkat */
