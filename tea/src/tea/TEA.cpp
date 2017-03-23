/*
 * TEA.cpp
 *
 *  Created on: Aug 17, 2016
 *      Author: el174
 */

#include "TEA.hpp"

namespace tea {

TEA::TEA() :
		n_cores(1), qenc(0), min_qual('#') {
}

TEA::~TEA() {
}

void TEA::set_option_parser(const TEAOptionParser& an_option) {
	options = an_option;
	n_cores = options.n_cores;
}
void TEA::preprocess() {
	format_isize();
	int32_t start_step = options.sub_name_map[options.start_step];
	if("all" == options.sub_module) {
		switch (start_step) {
			case 0:
				create_disc_FASTQs();
			case 1:
				generate_ra_bams();
			case 2:
				generate_ram_files();
			case 3:
				generate_cbam_files();
			case 4:
				run_rid();
			default:
				break;
		}
	} else {
		int32_t selected_step = options.sub_name_map[options.sub_module];
		switch (selected_step) {
			case 0:
				create_disc_FASTQs();
				break;
			case 1:
				generate_ra_bams();
				break;
			case 2:
				generate_ram_files();
				break;
			case 3:
				generate_cbam_files();
				break;
			case 4:
				run_rid();
				break;
			default:
				break;
		}
	}
}

void TEA::preprocess_v() {
	format_isize();
	string cl_sorted_bam_path = options.prefix + ".cl.sorted.bam";
	string cl_sorted_bai_path = options.prefix + ".cl.sorted.bam.bai";
	string cl_sorted_bni_path = options.prefix + ".cl.sorted.bam.bni";
	if (!options.working_dir.empty()) {
		cl_sorted_bam_path = options.working_prefix + ".cl.sorted.bam";
		cl_sorted_bai_path = options.working_prefix + ".cl.sorted.bam.bai";
		cl_sorted_bni_path = options.working_prefix + ".cl.sorted.bam.bni";
	}
	generate_um_bams(cl_sorted_bam_path, cl_sorted_bai_path, cl_sorted_bni_path);
	create_um_FASTQs();
	generate_va_bams();
	generate_cbam_files();
	generate_vam_files();
}

void TEA::preprocess_u() {
	options.rasym = "um";
	format_isize();
	string cl_sorted_bam_path = options.prefix + ".cl.sorted.bam";
	string cl_sorted_bai_path = options.prefix + ".cl.sorted.bam.bai";
	string cl_sorted_bni_path = options.prefix + ".cl.sorted.bam.bni";
	string input_bam_path = options.prefix + ".bam";
	string input_bai_path = options.prefix + ".bam.bai";
	string input_bni_path = options.prefix + ".bam.bni";

	if (!options.working_dir.empty()) {
		cl_sorted_bam_path = options.working_prefix + ".cl.sorted.bam";
		cl_sorted_bai_path = options.working_prefix + ".cl.sorted.bam.bai";
		cl_sorted_bni_path = options.working_prefix + ".cl.sorted.bam.bni";
		input_bni_path = options.working_prefix + ".bam.bni";
	}
	if(boost::filesystem::exists(cl_sorted_bam_path)) {
		generate_um_bams(cl_sorted_bam_path, cl_sorted_bai_path, cl_sorted_bni_path);
	} else if(boost::filesystem::exists(input_bam_path)) {
		generate_um_bams(input_bam_path, input_bai_path, input_bni_path);
	} else {
		cout << (boost::format("[TEA.preprocess_u] Neither %s nor %s exists\n") % cl_sorted_bam_path % input_bam_path).str();
		exit(1);
	}
	generate_cbam_files();
}

//# rid can detect both endogeneous and exogeneous insertions
//# rasym needs to be set up to represent different sequence libraries
//# (example)
//# ra: TE sequence libraries
//# va: virus sequence libraries
//# um: no sequence library was used but identifying the clusters of unmapped reads
void TEA::run_rid() {
	string all_assembly_files = options.prefix + ".asssemblies.target";
	if (!options.working_dir.empty()) {
		all_assembly_files = options.working_prefix + ".asssemblies.target";
	}
	if(!options.is_force && boost::filesystem::exists(all_assembly_files)) {
		return;
	}
	set<string> ref_support;
	ref_support.insert("hg18");
	ref_support.insert("hg19");
	ref_support.insert("ponAbe2");
	ref_support.insert("panTro3");
	ref_support.insert("rheMac2");
	string ref = options.ref;
	if (ref_support.end() == ref_support.find(ref)) {
		cout << (boost::format("[TEA.run_rid] The reference %s is not supported\n") % ref).str();
		exit(0);
	}

	boost::unordered_map<string, pair<string, string>> rannot;
	load_repeat_annotation(rannot);

//	boost::unordered_map<string, pair<string, string>> vannot;
//	if("va" == options.rasym) {
//		load_virus_annotation()
//	}

	set<string> chrl;
	boost::unordered_map<string, vector<pair<int64_t, int64_t>>> gap_annot;
	boost::unordered_map<string, RefRepeatIntervalVector> ril_annot_alt;
	boost::unordered_map<string, GeneIntervalVector> gene_annot;

	bool out_chrl = true;
	bool out_gap = false;

	load_ref_annotation(chrl, rannot, ril_annot_alt, gap_annot, gene_annot, out_chrl, out_gap);

//	for(auto& an_itr : ril_annot_alt) {
//		cout << (boost::format("[TEA.run_rid] entry: %s\n") % an_itr.first).str();
//	}

	string cbam_file;
	if (!options.no_clipped && options.is_sampe) {
		cbam_file = options.prefix + ".softclips.consd.bam";
	}

	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string contig_dir = options.prefix + "/assembly_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		if (!options.no_clipped && options.is_sampe) {
			cbam_file = options.working_prefix + ".softclips.consd.bam";
		}
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		contig_dir = options.working_prefix + "/assembly_" + options.rasym + "m";
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	if (!options.no_clipped && options.is_sampe) {
		if (!boost::filesystem::exists(cbam_file)) {
			cout << (boost::format("[TEA.run_rid] there is no clipped bam file: %s\n") % cbam_file).str();
			exit(1);
		}
	}
	if (!boost::filesystem::exists(cl_dir)) {
		boost::filesystem::create_directories(cl_dir);
	}
	if (!boost::filesystem::exists(contig_dir)) {
		boost::filesystem::create_directories(contig_dir);
	}

	boost::unordered_map<string, int32_t> rl;
	load_read_length(rl);
	boost::unordered_map<string, boost::unordered_map<string, double>> is;
	load_insert_size(is, rl);

	auto& the_rep_is = is["all"];
	double rmasker_filter_margin = 500;

	cout << (boost::format("[TEA.run_rid] fragment: %d (mu: %d, sd: %d), intra.gap: %d, inter.gap: %d, ins.margin: %d, rmasker.filter.margin: %d\n") % the_rep_is["fr"] % the_rep_is["mu"] % the_rep_is["sd"] % the_rep_is["intra_gap"] % the_rep_is["inter_gap"] % the_rep_is["ins_margin"]
					% rmasker_filter_margin).str();

	const bool rm_dup = false;

	boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>> ram;
	load_ram(ram, rannot, rm_dup);
	for (auto& c_entry : ram) {
		auto& c = c_entry.first;
		// create positive and negative ram storage for the parallelization
		if (ram[c][1].size() > 0) {
			;
		}
		if (ram[c][-1].size() > 0) {
			;
		}
	}
	string the_first_stat_file = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		the_first_stat_file = options.working_prefix + ".firststat";
	}
	if (!boost::filesystem::exists(the_first_stat_file)) {
		cout << "[TEA.run_rid] ERROR:: could not open " << the_first_stat_file << "\n";
		exit(0);
	}

	{
		string line;
		ifstream in(the_first_stat_file);
		for (int64_t line_id = 0; line_id < 5; ++line_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int32_t>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}

//	string c = "22";
	vector<string> chr_names;
	for(auto& c_entry : ram) {
		auto c = c_entry.first;
		chr_names.push_back(c);
	}
	sort(chr_names.begin(), chr_names.end(), [&](const string& lhs, const string& rhs)->bool{
		string local_lhs = lhs;
		string local_rhs = rhs;
		boost::replace_all(local_lhs, "chr", "");
		boost::replace_all(local_rhs, "chr", "");
		int64_t lhs_v = numeric_limits<int64_t>::max();
		int64_t rhs_v = lhs_v;
		try {
			lhs_v = boost::lexical_cast<int64_t>(local_lhs);
		} catch(exception& ex) {}
		try {
			rhs_v = boost::lexical_cast<int64_t>(local_rhs);
		} catch(exception& ex) {}
		if(lhs_v < rhs_v) {
			return true;
		} else if(lhs_v > rhs_v) {
			return false;
		}
		return local_lhs < local_rhs;
	});
	vector<function<void()> > tasks;
	bool headless = false;
	for (auto c: chr_names) {
//		auto c = c_entry.first;
		tasks.push_back([&, c, headless] {
			RAMIntervalVector p_cl;
			RAMIntervalVector n_cl;
			if (ram[c][1].size() > 0) {
				get_cluster_alt(c, p_cl, ram[c][1], rannot, 1, is["all"]["intra_gap"]);
			}
			if (ram[c][-1].size() > 0) {
				get_cluster_alt(c, n_cl, ram[c][-1], rannot, -1, is["all"]["intra_gap"]);
			}
			map<int64_t, int64_t> pm_cl;
			vector<int64_t> unpaired_pidx;
			vector<int64_t> unpaired_nidx;

			if (!p_cl.empty() && !n_cl.empty()) {
				bool stringent_pair = ("um" == options.rasym);
				for (uint64_t r_id = 0; r_id < p_cl.size(); ++r_id) {
					p_cl[r_id].value.global_cluster_id = r_id;
				}
				for (uint64_t r_id = 0; r_id < n_cl.size(); ++r_id) {
					n_cl[r_id].value.global_cluster_id = r_id;
				}
				pair_cluster_alt(pm_cl, p_cl, n_cl, is["all"]["inter_gap"], rl["all"], stringent_pair);
			}
			boost::unordered_set<int64_t> positive_paired;
			boost::unordered_set<int64_t> negative_paired;

			for (auto an_entry : pm_cl) {
				auto& positive_entry = p_cl[an_entry.first];
				auto& negative_entry = n_cl[an_entry.second];
				positive_paired.insert(positive_entry.value.global_cluster_id);
				negative_paired.insert(negative_entry.value.global_cluster_id);
			}
			boost::unordered_set<int64_t> positive_only;
			boost::unordered_set<int64_t> negative_only;

			for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
				if (positive_paired.end() != positive_paired.find(r_id)) {
					continue;
				}
				positive_only.insert(r_id);
			}
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
				if (negative_paired.end() != negative_paired.find(r_id)) {
					continue;
				}
				negative_only.insert(r_id);
			}
			output_raw_file(c, cl_prefix, p_cl, n_cl, pm_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], headless);
			int64_t gene_margin = 2;
			count_clipped(ril_annot_alt, gene_annot, c, cl_prefix, contig_dir, pm_cl, p_cl, n_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], rmasker_filter_margin, gene_margin, headless);
		});
		headless = true;
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	output_mate_fa(ram);

	vector<string> cluster_raw_files;
	vector<string> clipped_files;
	vector<string> germline_files;
	vector<string> cluster_files;
	vector<string> germline_contig_files;

	vector<string> removal_files;
	for (auto c : chr_names) {
//		auto c = c_entry.first;
		string tmp_chr_name(c);
		if(string::npos == tmp_chr_name.find("chr")) {
			tmp_chr_name = "chr" + c;
		}
		string p_mate_rname_file = cl_prefix + "." + tmp_chr_name + ".p.mate.rname";
		string n_mate_rname_file = cl_prefix + "." + tmp_chr_name + ".n.mate.rname";
		string p_clipped_fname_file = cl_prefix + "." + tmp_chr_name + ".p.clipped.fname";
		string n_clipped_fname_file = cl_prefix + "." + tmp_chr_name + ".n.clipped.fname";
		removal_files.push_back(p_mate_rname_file);
		removal_files.push_back(n_mate_rname_file);
		removal_files.push_back(p_clipped_fname_file);
		removal_files.push_back(n_clipped_fname_file);

		string cluster_raw_file = cl_prefix + "." + tmp_chr_name + ".cluster.raw";
		string clipped_file = cl_prefix + "." + tmp_chr_name + ".clipped";
		string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";
		string cluster_file = cl_prefix + "." + tmp_chr_name + ".cluster";
		string germline_contig_file = cl_prefix + "." + tmp_chr_name + ".germline.contig";
		cluster_raw_files.push_back(cluster_raw_file);
		clipped_files.push_back(clipped_file);
		germline_files.push_back(germline_file);
		cluster_files.push_back(cluster_file);
		germline_contig_files.push_back(germline_contig_file);
	}
	string out_cluster_raw_file = cl_prefix + ".cluster.raw";
	string out_clipped_file = cl_prefix + ".clipped";
	string out_germline_file = cl_prefix + ".germline";
	string out_cluster_file = cl_prefix + ".cluster";
	string out_germline_contig_file = cl_prefix + ".germline.contig";

	castle::IOUtils::plain_file_merge(out_cluster_raw_file, cluster_raw_files, n_cores, true);
	castle::IOUtils::plain_file_merge(out_clipped_file, clipped_files, n_cores, true);
	castle::IOUtils::plain_file_merge(out_germline_file, germline_files, n_cores, true);
	castle::IOUtils::plain_file_merge(out_cluster_file, cluster_files, n_cores, true);
	castle::IOUtils::plain_file_merge(out_germline_contig_file, germline_contig_files, n_cores, true);
	castle::IOUtils::remove_files(removal_files, n_cores);
}

void TEA::run_vid() {
	set<string> ref_support;
	ref_support.insert("hg18");
	ref_support.insert("hg19");
	ref_support.insert("ponAbe2");
	ref_support.insert("panTro3");
	ref_support.insert("rheMac2");
	string ref = options.ref;
	if (ref_support.end() == ref_support.find(ref)) {
		cout << (boost::format("The reference %s is not supported\n") % ref).str();
		exit(0);
	}

	boost::unordered_map<string, pair<string, string>> rannot;
	load_repeat_annotation(rannot);

	map<int64_t, string> vannot;
	if("va" == options.rasym) {
		load_virus_annotation(vannot);
	}

	set<string> chrl;
	boost::unordered_map<string, vector<pair<int64_t, int64_t>>> gap_annot;
	boost::unordered_map<string, RefRepeatIntervalVector> ril_annot_alt;
	boost::unordered_map<string, GeneIntervalVector> gene_annot;

	bool out_chrl = true;
	bool out_gap = false;

	load_ref_annotation(chrl, rannot, ril_annot_alt, gap_annot, gene_annot, out_chrl, out_gap);

	string cbam_file;
	if (!options.no_clipped && options.is_sampe) {
		cbam_file = options.prefix + ".softclips.consd.bam";
	}

	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string contig_dir = options.prefix + "/assembly_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		if (!options.no_clipped && options.is_sampe) {
			cbam_file = options.working_prefix + ".softclips.consd.bam";
		}
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		contig_dir = options.prefix + "/assembly_" + options.rasym + "m";
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	if (!options.no_clipped && options.is_sampe) {
		if (!boost::filesystem::exists(cbam_file)) {
			cout << (boost::format("there is no clipped bam file: %s\n") % cbam_file).str();
			exit(1);
		}
	}
	if (!boost::filesystem::exists(cl_dir)) {
		boost::filesystem::create_directories(cl_dir);
	}
	if (!boost::filesystem::exists(contig_dir)) {
		boost::filesystem::create_directories(contig_dir);
	}

	boost::unordered_map<string, int32_t> rl;
	load_read_length(rl);
	boost::unordered_map<string, boost::unordered_map<string, double>> is;
	load_insert_size(is, rl);

	auto& the_rep_is = is["all"];
	double rmasker_filter_margin = 500;

	cout << (boost::format("[TEA.run_vid] fragment: %d (mu: %d, sd: %d), intra.gap: %d, inter.gap: %d, ins.margin: %d, rmasker.filter.margin: %d\n") % the_rep_is["fr"] % the_rep_is["mu"] % the_rep_is["sd"] % the_rep_is["intra_gap"] % the_rep_is["inter_gap"] % the_rep_is["ins_margin"] %
			rmasker_filter_margin).str();

	const bool rm_dup = false;

	boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>> ram;
	load_ram(ram, rannot, rm_dup);
	for (auto& c_entry : ram) {
		auto& c = c_entry.first;
		// create positive and negative ram storage for the parallelization
		if (ram[c][1].size() > 0) {
			;
		}
		if (ram[c][-1].size() > 0) {
			;
		}
	}
	string the_first_stat_file = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		the_first_stat_file = options.working_prefix + ".firststat";
	}
	if (!boost::filesystem::exists(the_first_stat_file)) {
		cout << "ERROR:: could not open " << the_first_stat_file << "\n";
		exit(0);
	}

	{
		string line;
		ifstream in(the_first_stat_file);
		for (int64_t line_id = 0; line_id < 5; ++line_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int32_t>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}

//	string c = "22";
	vector<function<void()> > tasks;
	bool headless = false;
	for (auto& c_entry : ram) {
		auto c = c_entry.first;
		tasks.push_back([&, c, headless] {
			RAMIntervalVector p_cl;
			RAMIntervalVector n_cl;
			if (ram[c][1].size() > 0) {
				get_cluster_alt(c, p_cl, ram[c][1], rannot, 1, is["all"]["intra_gap"]);
			}
			if (ram[c][-1].size() > 0) {
				get_cluster_alt(c, n_cl, ram[c][-1], rannot, -1, is["all"]["intra_gap"]);
			}
			map<int64_t, int64_t> pm_cl;
			vector<int64_t> unpaired_pidx;
			vector<int64_t> unpaired_nidx;

			if (!p_cl.empty() && !n_cl.empty()) {
				bool stringent_pair = ("um" == options.rasym);
				for (uint64_t r_id = 0; r_id < p_cl.size(); ++r_id) {
					p_cl[r_id].value.global_cluster_id = r_id;
				}
				for (uint64_t r_id = 0; r_id < n_cl.size(); ++r_id) {
					n_cl[r_id].value.global_cluster_id = r_id;
				}
				pair_cluster_alt(pm_cl, p_cl, n_cl, is["all"]["inter_gap"], rl["all"], stringent_pair);
			}
			boost::unordered_set<int64_t> positive_paired;
			boost::unordered_set<int64_t> negative_paired;

			for (auto an_entry : pm_cl) {
				auto& positive_entry = p_cl[an_entry.first];
				auto& negative_entry = n_cl[an_entry.second];
				positive_paired.insert(positive_entry.value.global_cluster_id);
				negative_paired.insert(negative_entry.value.global_cluster_id);
			}
			boost::unordered_set<int64_t> positive_only;
			boost::unordered_set<int64_t> negative_only;
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
				if (positive_paired.end() != positive_paired.find(r_id)) {
					continue;
				}
				positive_only.insert(r_id);
			}
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
				if (negative_paired.end() != negative_paired.find(r_id)) {
					continue;
				}
				negative_only.insert(r_id);
			}
			output_raw_file(c, cl_prefix, p_cl, n_cl, pm_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], headless);
			int64_t gene_margin = 2;
			count_clipped_v(ril_annot_alt, vannot, c, cl_prefix, contig_dir, pm_cl, p_cl, n_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], rmasker_filter_margin, gene_margin);
		});
		headless = true;
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	output_mate_fa_v(ram);
}

void TEA::run_uid() {
	if("all" != options.sub_module) {
		return;
	}
	string all_assembly_files = options.prefix + ".asssemblies.target";
	if (!options.working_dir.empty()) {
		all_assembly_files = options.working_prefix + ".asssemblies.target";
	}
	set<string> ref_support;
	ref_support.insert("hg18");
	ref_support.insert("hg19");
	ref_support.insert("ponAbe2");
	ref_support.insert("panTro3");
	ref_support.insert("rheMac2");
	string ref = options.ref;
	if (ref_support.end() == ref_support.find(ref)) {
		cout << (boost::format("The reference %s is not supported\n") % ref).str();
		exit(0);
	}

	boost::unordered_map<string, pair<string, string>> rannot;
	load_repeat_annotation(rannot);

//	boost::unordered_map<string, pair<string, string>> vannot;
//	if("va" == options.rasym) {
//		load_virus_annotation()
//	}

	set<string> chrl;
	boost::unordered_map<string, vector<pair<int64_t, int64_t>>> gap_annot;
	boost::unordered_map<string, RefRepeatIntervalVector> ril_annot_alt;
	boost::unordered_map<string, GeneIntervalVector> gene_annot;

	bool out_chrl = true;
	bool out_gap = false;

	load_ref_annotation(chrl, rannot, ril_annot_alt, gap_annot, gene_annot, out_chrl, out_gap);

//	for(auto& an_itr : ril_annot_alt) {
//		cout << (boost::format("[TEA.run_rid] entry: %s\n") % an_itr.first).str();
//	}

	string cbam_file;
	if (!options.no_clipped && options.is_sampe) {
		cbam_file = options.prefix + ".softclips.consd.bam";
	}

	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string contig_dir = options.prefix + "/assembly_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		if (!options.no_clipped && options.is_sampe) {
			cbam_file = options.working_prefix + ".softclips.consd.bam";
		}
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		contig_dir = options.prefix + "/assembly_" + options.rasym + "m";
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	if (!options.no_clipped && options.is_sampe) {
		if (!boost::filesystem::exists(cbam_file)) {
			cout << (boost::format("there is no clipped bam file: %s\n") % cbam_file).str();
			exit(1);
		}
	}
	if (!boost::filesystem::exists(cl_dir)) {
		boost::filesystem::create_directories(cl_dir);
	}
	if (!boost::filesystem::exists(contig_dir)) {
		boost::filesystem::create_directories(contig_dir);
	}

	boost::unordered_map<string, int32_t> rl;
	load_read_length(rl);
	boost::unordered_map<string, boost::unordered_map<string, double>> is;
	load_insert_size(is, rl);

	auto& the_rep_is = is["all"];
	double rmasker_filter_margin = 500;

	cout << (boost::format("[TEA.run_uid] fragment: %d (mu: %d, sd: %d), intra.gap: %d, inter.gap: %d, ins.margin: %d, rmasker.filter.margin: %d\n") % the_rep_is["fr"] % the_rep_is["mu"] % the_rep_is["sd"] % the_rep_is["intra_gap"] % the_rep_is["inter_gap"] % the_rep_is["ins_margin"]
					% rmasker_filter_margin).str();

	const bool rm_dup = false;

	boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>> ram;
	load_ram(ram, rannot, rm_dup);
	for (auto& c_entry : ram) {
		auto& c = c_entry.first;
		// create positive and negative ram storage for the parallelization
		if (ram[c][1].size() > 0) {
			;
		}
		if (ram[c][-1].size() > 0) {
			;
		}
	}
	string the_first_stat_file = options.prefix + ".firststat";
	if (!options.working_dir.empty()) {
		the_first_stat_file = options.working_prefix + ".firststat";
	}
	if (!boost::filesystem::exists(the_first_stat_file)) {
		cout << "ERROR:: could not open " << the_first_stat_file << "\n";
		exit(0);
	}

	{
		string line;
		ifstream in(the_first_stat_file);
		for (int64_t line_id = 0; line_id < 5; ++line_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int32_t>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}

//	string c = "22";
	vector<function<void()> > tasks;
	bool headless = false;
	for (auto& c_entry : ram) {
		auto c = c_entry.first;
		tasks.push_back([&, c, headless] {
			RAMIntervalVector p_cl;
			RAMIntervalVector n_cl;
			if (ram[c][1].size() > 0) {
				get_cluster_alt(c, p_cl, ram[c][1], rannot, 1, is["all"]["intra_gap"]);
			}
			if (ram[c][-1].size() > 0) {
				get_cluster_alt(c, n_cl, ram[c][-1], rannot, -1, is["all"]["intra_gap"]);
			}
			map<int64_t, int64_t> pm_cl;
			vector<int64_t> unpaired_pidx;
			vector<int64_t> unpaired_nidx;

			if (!p_cl.empty() && !n_cl.empty()) {
				bool stringent_pair = ("um" == options.rasym);
				for (uint64_t r_id = 0; r_id < p_cl.size(); ++r_id) {
					p_cl[r_id].value.global_cluster_id = r_id;
				}
				for (uint64_t r_id = 0; r_id < n_cl.size(); ++r_id) {
					n_cl[r_id].value.global_cluster_id = r_id;
				}
				pair_cluster_alt(pm_cl, p_cl, n_cl, is["all"]["inter_gap"], rl["all"], stringent_pair);
			}
			boost::unordered_set<int64_t> positive_paired;
			boost::unordered_set<int64_t> negative_paired;

			for (auto an_entry : pm_cl) {
				auto& positive_entry = p_cl[an_entry.first];
				auto& negative_entry = n_cl[an_entry.second];
				positive_paired.insert(positive_entry.value.global_cluster_id);
				negative_paired.insert(negative_entry.value.global_cluster_id);
			}
			boost::unordered_set<int64_t> positive_only;
			boost::unordered_set<int64_t> negative_only;

			for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
				if (positive_paired.end() != positive_paired.find(r_id)) {
					continue;
				}
				positive_only.insert(r_id);
			}
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
				if (negative_paired.end() != negative_paired.find(r_id)) {
					continue;
				}
				negative_only.insert(r_id);
			}
//			cout << (boost::format("[TEA.run_uid] chr: %s, paired: %s, pos_only: %s, neg_only: %s\n") % c % pm_cl.size() % positive_only.size() % negative_only.size()).str();
			output_raw_file(c, cl_prefix, p_cl, n_cl, pm_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], headless);
			int64_t gene_margin = 2;
			count_clipped(ril_annot_alt, gene_annot, c, cl_prefix, contig_dir, pm_cl, p_cl, n_cl, positive_only, negative_only, rl["all"], the_rep_is["fr"], rmasker_filter_margin, gene_margin, headless);
		});
		headless = true;
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	output_mate_fa(ram);
}

void TEA::run_transduction() {
	options.rasym = "um";
	set<string> chrs;
	load_chr(chrs);
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string transduction_dir = options.prefix + "/transduction_" + options.rasym + "m";
	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		transduction_dir = options.working_prefix + "/transduction_" + options.rasym + "m";
	}
	if(!boost::filesystem::exists(transduction_dir)) {
		boost::filesystem::create_directories(transduction_dir);
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	string transduction_prefix = transduction_dir + "/" + naive_prefix;

	cout << (boost::format("[TEA.run_transduction] copy %sm files to transduction folder and create the discord.bed file\n") % options.rasym).str();
	vector<function<void()> > tasks;

	string in_discord = options.prefix + ".discord";
	string out_discord = transduction_prefix + ".discord.bed";
	string in_cluster = options.prefix + ".clusters";
	if (!options.working_dir.empty()) {
		in_discord = options.working_prefix + ".discord";
		in_cluster = options.working_prefix + ".clusters";
	}
	create_discord_bed(out_discord, in_discord);
	vector<boost::unordered_set<string>> o_lists(chrs.size());
	vector<boost::unordered_map<string, int32_t>> cluster_entry_lists(chrs.size());
	vector<boost::unordered_map<string, int32_t>> dup_cnts_lists(chrs.size());
	boost::unordered_map<string, int64_t> chr_ids;
	int64_t chr_id = 0;
	for(auto& chr: chrs) {
		chr_ids[chr] = chr_id;
		++chr_id;
	}

	const bool debug = false;
	for(auto& chr: chrs) {
		tasks.push_back([&, debug]{
			string in_germline = cl_prefix + "." + chr + ".germline";
			string in_clusterraw = cl_prefix + "." + chr + ".cluster.raw";

			string out_getrmline = transduction_prefix + "." + chr + ".germline";
			string out_clusterraw = transduction_prefix + "." + chr + ".cluster.raw";
			boost::filesystem::copy_file(in_germline, out_getrmline, boost::filesystem::copy_option::overwrite_if_exists);
			boost::filesystem::copy_file(in_clusterraw, out_clusterraw, boost::filesystem::copy_option::overwrite_if_exists);

			string out_getrmline_transduction = transduction_prefix + "." + chr + ".germline.transduction";
			if(debug) {
				cout << (boost::format("[TEA.run_transduction] create germline transduction: %s\n") % chr).str();
			}
			create_germline_transduction(out_getrmline_transduction, out_getrmline);
			string out_intersect_bed = transduction_prefix + "." + chr + ".intersected.bed";
			string bedtools_cmd = (boost::format("bedtools pairtobed -type xor -S -a %s -b %s > %s") % out_discord % out_getrmline_transduction % out_intersect_bed).str();
			if(debug) {
				cout << bedtools_cmd << "\n";
			}
			system(bedtools_cmd.c_str());

			boost::unordered_set<string> ram_id;
			auto& cluster_entries = cluster_entry_lists[chr_ids[chr]];
			if(debug) {
				cout << (boost::format("[TEA.run_transduction] read_ram_ids: %s\n") % chr).str();
			}
			read_ram_ids(ram_id, out_intersect_bed, cluster_entries);

			auto& o = o_lists[chr_ids[chr]];

			{
				string in_clusterraw = transduction_prefix + "." + chr + ".cluster.raw";
				string out_umm_tmp1 = transduction_prefix + "." + chr + ".umm.tmp.1";
				if(debug) {
					cout << (boost::format("[TEA.run_transduction] create .umm.tmp (%s) from .cluster.raw (%s) %s\n") % out_umm_tmp1 % in_clusterraw % chr).str();
				}
				create_umm_tmp_from_cluster_raw(out_umm_tmp1, in_clusterraw, o, ram_id);
			}
			if(debug) {
				cout << (boost::format("[TEA.run_transduction] count dup counts (%s)\n") % chr).str();
			}
			// count seq_ids {seq_id, counts}
			auto& local_dup_cnt = dup_cnts_lists[chr_ids[chr]];
			{
				string line;
				vector<string> fso;
				const char* delim_tab = "\t";
				// skips the first description line
				string out_umm_tmp1 = transduction_prefix + "." + chr + ".umm.tmp.1";

				ifstream in(out_umm_tmp1, ios::binary);
				while(getline(in, line, '\n')) {
					castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
				    ++local_dup_cnt[fso[0]];
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	if(debug) {
		cout << "[TEA.run_transduction] read ram read ids\n";
	}
	boost::unordered_set<string> ram_read_ids;
	read_ram_read_ids(ram_read_ids);

	boost::unordered_set<string> o;
	for(auto& a_set : o_lists) {
		o.insert(a_set.begin(), a_set.end());
	}

	boost::unordered_map<string, int32_t> cluster_entries;
	for(auto& a_map : cluster_entry_lists) {
		for(auto& the_itr : a_map) {
			cluster_entries[the_itr.first] = the_itr.second;
		}
	}

	boost::unordered_map<string, int32_t> dup_cnt;
	for(auto& dup_cnt_map : dup_cnts_lists) {
		for(auto& dup_entry : dup_cnt_map) {
			auto& seq_id = dup_entry.first;
			auto& dup_n = dup_entry.second;
			dup_cnt[seq_id] += dup_n;
		}
	}

	string out_umm_tmp2 = transduction_prefix + ".umm.tmp.2";
	cout << (boost::format("[TEA.run_transduction] create %s from .clusters %s\n") % out_umm_tmp2 % in_cluster).str();
	create_umm_tmp_from_cluster(out_umm_tmp2, in_cluster, o, cluster_entries, ram_read_ids);
	{
		string line;
		vector<string> fso;
		const char* delim_tab = "\t";
		// skips the first description line
		ifstream in(out_umm_tmp2, ios::binary);
		while(getline(in, line, '\n')) {
			castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
			++dup_cnt[fso[0]];
		}
	}
	// write non-duplicated entries to umm
	for(auto& chr: chrs) {
		tasks.push_back([&, debug]{
			string in_umm_tmp = transduction_prefix + "." + chr + ".umm.tmp.1";
			string out_umm = transduction_prefix + "." + chr + ".umm.1";
			if(debug) {
				cout << (boost::format("[TEA.run_transduction] write non dup umm(%s) from umm_tmp(%s)\n") % out_umm % in_umm_tmp).str();
			}
			write_non_dup_umm(out_umm, in_umm_tmp, dup_cnt);
		});
	}
	tasks.push_back([&]{
		string out_umm2 = transduction_prefix + ".umm.2";
		if(debug) {
			cout << (boost::format("[TEA.run_transduction] write non dup umm(%s) from umm_tmp(%s)\n") % out_umm2 % out_umm_tmp2).str();
		}
		write_non_dup_umm(out_umm2, out_umm_tmp2, dup_cnt);
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void TEA::create_germline_transduction(const string& out_path, const string& in_path) {
	string line;
	boost::unordered_set<string> o;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream in(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	// skips the first description line
	getline(in, line, '\n');
	while(getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		for(auto& a_col : fso) {
			castle::StringUtils::trim(a_col);
			boost::replace_all(a_col, "\r", "");
			boost::replace_all(a_col, " ", "");
		}
		if (string::npos != fso[1].find("chr")) {
			boost::replace_all(fso[1], "chr", "");
		}
		string name2 = fso[2] + '@' + fso[3];
		string name = fso[6] + '@' + fso[7];
		if(o.end() == o.find(name)) {
			if(string::npos != fso[1].find("GL")) {
				continue;
			}
			o.insert(name);
			if("0" != fso[22]) {
				int64_t p_ram_end = boost::lexical_cast<int64_t>(fso[22]);
				int64_t s = p_ram_end - 1000;
				int64_t e = p_ram_end + 1000;
				string strand = "+";
				string a_line = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\n") % fso[1] % s % e % name % name2 % strand).str();
				g << a_line;
			}
			if("0" != fso[23]) {
				int64_t n_ram_start = boost::lexical_cast<int64_t>(fso[23]);
				int64_t s = n_ram_start - 1000;
				int64_t e = n_ram_start + 1000;
				string strand = "-";
				string a_line = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\n") % fso[1] % s % e % name % name2 % strand).str();
				g << a_line;
			}
		}
	}
}

void TEA::create_discord_bed(const string& out_path, const string& in_path) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.create_discord_bed");
	checker.start();
	vector<function<void()> > tasks;
	int64_t max_pos = castle::IOUtils::get_file_size(in_path);
	int64_t chunk_size = max_pos / (double)n_cores;
	int64_t n_blocks = max_pos / (double)chunk_size;
	vector<uint64_t> skip_points;
	castle::IOUtils::find_skip_points(skip_points, in_path, chunk_size, max_pos, n_blocks, n_cores);
	vector<string> out_files(n_blocks);

	for(int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, block_id]{
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = skip_points[block_id];
			int64_t end_pos = skip_points[block_id + 1];

			vector<string> fso;
			const char* delim_tab = "\t";

			string line;
			string out_file_name = out_path + "." + str_block_id;
			ifstream in(in_path, ios::binary);
			in.seekg(cur_pos, ios::beg);
			out_files[block_id] = out_file_name;
			ofstream out(out_file_name, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				if(cur_pos > end_pos) {
					break;
				}
				castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
				if((string::npos != fso[3].find("GL")) || (string::npos != fso[6].find("GL"))) {
					continue;
				}
				int64_t cluster_cnt = boost::lexical_cast<int64_t>(fso[2]);
				if (cluster_cnt < 3) {
					continue;
				}
				for(auto& a_col : fso) {
					castle::StringUtils::trim(a_col);
					boost::replace_all(a_col, "\r", "");
					boost::replace_all(a_col, " ", "");
				}

				int64_t first_pos = boost::lexical_cast<int64_t>(fso[4]);
				int64_t start = first_pos - 1000;
				int64_t end = first_pos + 1000;

				int64_t second_pos = boost::lexical_cast<int64_t>(fso[7]);
				int64_t start2 = second_pos - 1000;
				int64_t end2 = second_pos + 1000;

				auto& first_strand = fso[5];
				if ("1" == first_strand) {
					first_strand = "+";
				} else {
					first_strand = "-";
				}

				auto& second_strand = fso[8];
				if ("1" == second_strand) {
					second_strand = "+";
				} else {
					second_strand = "-";
				}

				auto& first_chr = fso[3];
				auto& second_chr = fso[6];
				if (first_chr == second_chr) {
					int64_t diff = abs(second_pos - first_pos);
					if (diff > 0 && diff <= 10000) {
						continue;
					}
				}

				auto& cluster_p_id = fso[0];
				auto& cluster_s_id = fso[1];
				string a_line = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % first_chr % start % end % second_chr % start2 % end2 % cluster_p_id % cluster_s_id % first_strand % second_strand).str();
				out << a_line;
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	castle::IOUtils::plain_file_merge(out_path, out_files, n_cores, true);
	cout << checker;
}

void TEA::read_ram_read_ids(boost::unordered_set<string>& ram_read_ids) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.read_ram_read_ids");
	checker.start();
	string in_ram = options.prefix + ".ram";
	if (!options.working_dir.empty()) {
		in_ram = options.working_prefix + ".ram";
	}

	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	ifstream in(in_ram, ios::binary);

	while(getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		ram_read_ids.insert(fso[0]);
	}
	cout << checker;
}
void TEA::read_ram_ids(boost::unordered_set<string>& ram_ids, const string& in_path, boost::unordered_map<string, int32_t>& cluster_entries) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream in(in_path);

	while (getline(in, line, '\n')){
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& ram_id = fso[14];
		auto& ram_strand = fso[15];
		string a_new_ram_id = ram_id + "@" + ram_strand;
	    ram_ids.insert(a_new_ram_id);
	    int64_t on = 0;
	    auto& chr1 = fso[0];
	    auto& chr2 = fso[3];
	    auto& center_chr = fso[10];
	    auto& strand1 = fso[8];
		auto& center_strand = fso[15];
		auto& strand2 = fso[9];
	    if (chr1 == chr2) {
	    	int64_t start1 = boost::lexical_cast<int64_t>(fso[1]);
	    	int64_t center = boost::lexical_cast<int64_t>(fso[11]);
	    	int64_t start2 = boost::lexical_cast<int64_t>(fso[4]);

	        if (abs(start1 - center) < abs(start2 - center)) {
	            if (center_strand != strand1) {
	                on = 1;
	            }
	        } else {
	        	if (center_strand != strand2) {
	                on = 2;
	        	}
	        }
	    } else {
	        if (center_chr == chr1 && center_strand != strand1) {
	            on = 1;
	        }
	        if (center_chr == chr2 && center_strand != strand2) {
	            on = 2;
	        }
	    }
	    if(0 == on) {
	    	continue;
	    }
		auto& cluster_pid = fso[6];
		auto& cluster_sid = fso[7];
		string cluster_id = cluster_pid + "@" + cluster_sid;
		cluster_entries[cluster_id] = on;
	}
}

void TEA::create_umm_tmp_from_cluster_raw(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_set<string>& ram_id) {
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	string line;
	vector<string> fso;
	vector<string> read_names;
	vector<string> positions;

	const char* delim_tab = "\t";
	const char* delim_comma = ",";
	// skips the first description line
	getline(f, line, '\n');
	while(getline(f, line, '\n')) {
		castle::StringUtils::tokenize(line, delim_tab, fso);
		string temp = fso[1] + "@" + fso[2];
		if (string::npos != fso[0].find("chr")) {
			boost::replace_all(fso[0], "chr", "");
		}

		int64_t n_ram = boost::lexical_cast<int64_t>(fso[7]);

		if (n_ram < 3) {
			continue;
		}
		string temp_id_plus = temp + "@+";
		string temp_id_minus = temp + "@-";

		if(ram_id.end() != ram_id.find(temp_id_plus)) {
			auto& chromosome = fso[0];
			auto& read_names_plus = fso[20];
			castle::StringUtils::c_string_multi_split(read_names_plus, delim_comma, read_names);
			auto& positions_plus = fso[14];
			castle::StringUtils::c_string_multi_split(positions_plus, delim_comma, positions);
			for(uint64_t e_id = 0; e_id < read_names.size(); ++e_id) {
				string rename = read_names[e_id] + "@" + chromosome + "@" + boost::replace_all_copy(positions[e_id], "-", "");
				if(o.end() == o.find(rename)) {
					o.insert(rename);
					string a_line = (boost::format("%s\t%s\t%s\tx\n") % read_names[e_id] % chromosome % positions[e_id]).str();
					g << a_line;
				}
			}
		} else if(ram_id.end() != ram_id.find(temp_id_minus)) {
			auto& chromosome = fso[0];
			auto& read_names_minus = fso[21];
			castle::StringUtils::c_string_multi_split(read_names_minus, delim_comma, read_names);
			auto& positions_minus = fso[15];
			castle::StringUtils::c_string_multi_split(positions_minus, delim_comma, positions);
			for(uint64_t e_id = 0; e_id < read_names.size(); ++e_id) {
				string rename = read_names[e_id] + "@" + chromosome + "@" + boost::replace_all_copy(positions[e_id], "-", "");
				if(o.end() == o.find(rename)) {
					o.insert(rename);
					string a_line = (boost::format("%s\t%s\t%s\tx\n") % read_names[e_id] % chromosome % positions[e_id]).str();
					g << a_line;
				}
			}
		}
	}
}

void TEA::create_umm_tmp_from_cluster(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_map<string, int32_t> cluster_entries, const boost::unordered_set<string>& ram_read_ids) {
	if(castle::IOUtils::get_file_size(out_path) > 0) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.create_umm_tmp_from_cluster");
	checker.start();
	string line;
	vector<string> fso;
	vector<string> read_names;
	vector<string> positions;

	const char* delim_tab = "\t";
	ifstream f(in_path, ios::binary);
	// skips the first description line
	getline(f, line, '\n');
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& seq_id = fso[5];
		if(ram_read_ids.end() == ram_read_ids.find(seq_id)) {
			continue;
		}
		auto& cluster_pid = fso[0];
		auto& cluster_sid = fso[1];
		string cluster_id = cluster_pid + "@" + cluster_sid;
		auto cluster_itr = cluster_entries.find(cluster_id);
		if(cluster_entries.end() == cluster_itr) {
			continue;
		}
		int32_t cluster_type = cluster_itr->second;
		auto& chromosome1 = fso[7];
		auto& strand1 = fso[8];
		auto& pos1 = fso[9];
		int64_t rlen1 = boost::lexical_cast<int64_t>(fso[10]);
		auto& chromosome2 = fso[11];
		auto& strand2 = fso[12];
		auto& pos2 = fso[13];
		int64_t rlen2 = boost::lexical_cast<int64_t>(fso[14]);

		if (1 == cluster_type) {
			if (rlen1 >= rlen2) {
				boost::replace_all(strand1, "1", "");
				string rename = seq_id + "@" + chromosome1 + "@" + pos1;
				if (o.end() == o.find(rename)) {
					string a_line = (boost::format("%s\t%s\t%s%s\tx\n") % seq_id % chromosome1 % strand1 % pos1).str();
					g << a_line;
					o.insert(rename);
				}
			}
		} else if(2 == cluster_type) {
			if (rlen2 >= rlen1) {
				boost::replace_all(strand2, "1", "");
				string rename = seq_id + "@" + chromosome2 + "@" + pos2;
				if (o.end() == o.find(rename)) {
					string a_line = (boost::format("%s\t%s\t%s%s\tx\n") % seq_id % chromosome2 % strand2 % pos2).str();
					g << a_line;
					o.insert(rename);
				}
			}
		}
	}
	cout << checker;
}

void TEA::write_non_dup_umm(const string& out_path, const string& in_path, const boost::unordered_map<string, int32_t>& dup_cnt_map) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& seq_id = fso[0];
		auto dup_entry_itr = dup_cnt_map.find(seq_id);
	    if (dup_cnt_map.end() == dup_entry_itr) {
	    	continue;
	    }
	    auto dup_n = dup_entry_itr->second;
	    if(1 == dup_n) {
	        g << line << "\n";
	    }
	}
}

void TEA::run_transduction_contig() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.run_transduction_contig");
	checker.start();
	set<string> chrs;
	load_chr(chrs);
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string transduction_dir = options.prefix + "/transduction_" + options.rasym + "m";
	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		transduction_dir = options.working_prefix + "/transduction_" + options.rasym + "m";
	}
	if(!boost::filesystem::exists(transduction_dir)) {
		boost::filesystem::create_directories(transduction_dir);
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	string transduction_prefix = transduction_dir + "/" + naive_prefix;

	vector<function<void()> > tasks;
	boost::unordered_map<string, int64_t> chr_ids;

	int64_t chr_id = 0;
	for(auto& chr: chrs) {
		chr_ids[chr] = chr_id;
		++chr_id;
	}
	// create two ram files
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			string in_germline_contig = cl_prefix + "." + chr + ".germline.contig";
			string out_getrmline_contig_two_ram = transduction_prefix + "." + chr + ".germline.contig.two_ram";
			create_contig_two_ram(out_getrmline_contig_two_ram, in_germline_contig);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	// create insertion sequences
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			string in_two_ram = transduction_prefix + "." + chr + ".germline.contig.two_ram";
			string out_two_ram_fa = transduction_prefix + "." + chr + ".germline.contig.two_ram.fa";
			string out_two_only_ram_fa = transduction_prefix + "." + chr + ".germline.contig.two_only_ram.fa";
			create_fa_from_germline_contig(out_two_ram_fa, out_two_only_ram_fa, in_two_ram);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	cout << "[TEA.run_transduction_contig] align sequences\n";
//	boost::mutex print_mutex;
	// align insertion sequences
	for(auto& chr: chrs) {
		// two_ram against repeat ref
		tasks.push_back([&]{
			string in_two_ram_fa = transduction_prefix + "." + chr + ".germline.contig.two_ram.fa";
			string out_repeat_aln_sam = transduction_prefix + "." + chr + ".germline.contig.repeat.aln.sam";
//			if(castle::IOUtils::get_file_size(out_repeat_aln_sam) > 0) {
//				return;
//			}
			string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.repeat_reference % in_two_ram_fa % out_repeat_aln_sam % options.repeat_reference % in_two_ram_fa).str();
//			{
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << aln_cmd << "\n";
//			}
			system(aln_cmd.c_str());
		});
		// two_only_ram against repeat ref
		tasks.push_back([&]{
			string in_two_only_ram_fa = transduction_prefix + "." + chr + ".germline.contig.two_only_ram.fa";
			string out_repeat_aln_ram_sam = transduction_prefix + "." + chr + ".germline.contig.repeat.aln.ram.sam";
//			if(castle::IOUtils::get_file_size(out_repeat_aln_ram_sam) > 0) {
//				return;
//			}
			string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.repeat_reference % in_two_only_ram_fa % out_repeat_aln_ram_sam % options.repeat_reference % in_two_only_ram_fa).str();
//			{
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << aln_cmd << "\n";
//			}
			system(aln_cmd.c_str());
		});
		// two_ram against ref
		tasks.push_back([&]{
			string in_two_ram_fa = transduction_prefix + "." + chr + ".germline.contig.two_ram.fa";
			string out_ref_aln_sam = transduction_prefix + "." + chr + ".germline.contig.ref.aln.sam";
//			if(castle::IOUtils::get_file_size(out_ref_aln_sam) > 0) {
//				return;
//			}
			string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.human_reference % in_two_ram_fa % out_ref_aln_sam % options.human_reference % in_two_ram_fa).str();
//			{
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << aln_cmd << "\n";
//			}
			system(aln_cmd.c_str());
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	for(auto& chr: chrs) {
		tasks.push_back([&]{
			boost::unordered_map<int64_t, string> two_ram_map;
			string in_two_ram = transduction_prefix + "." + chr + ".germline.contig.two_ram";
			collect_two_ram_map(two_ram_map, in_two_ram);
			string in_ref_aln_sam = transduction_prefix + "." + chr + ".germline.contig.ref.aln.sam";
			boost::unordered_map<string, string> ref_aligned_map;
			set<string> ref_selected_seq_id;
			collect_two_ram_seq_id_ref_set(ref_aligned_map, ref_selected_seq_id, in_ref_aln_sam, two_ram_map);
			string in_repeat_aln_sam = transduction_prefix + "." + chr + ".germline.contig.repeat.aln.sam";
			set<string> repeat_selected_seq_id;
			boost::unordered_set<string> repeat_two_ram_id;
			collect_aln_sam_repeat(repeat_selected_seq_id, repeat_two_ram_id, in_repeat_aln_sam);
			string in_repeat_aln_ram_sam = transduction_prefix + "." + chr + ".germline.contig.repeat.aln.ram.sam";
			boost::unordered_map<string, string> repeat_ram_aligned_map;
			collect_aln_ram_sam_repeat(repeat_ram_aligned_map, in_repeat_aln_ram_sam);

			set<string> candidate;
			set_difference(ref_selected_seq_id.begin(), ref_selected_seq_id.end(), repeat_selected_seq_id.begin(), repeat_selected_seq_id.end(),
				inserter(candidate, candidate.begin()));
			set<string> gold;
			{
				string line;
				vector<string> temp;
				const char* delim_underscore = "_";

				for(auto& gene : candidate) {
					castle::StringUtils::c_string_multi_split(gene, delim_underscore, temp);
					auto& two_ram_line_id = temp[0];
				    if (repeat_two_ram_id.end() == repeat_two_ram_id.find(two_ram_line_id)) {
				        gold.insert(temp[0]);
				    }
				}
			}
//			const bool debug = "chr12" == chr;
//			if(debug) {
//				cout << (boost::format("[TEA.run_transduction_contig] Chr: %s, # ram_map: %s\n") % chr % repeat_ram_aligned_map.size());
//				for(auto ram_map_itr : repeat_ram_aligned_map) {
//					cout << (boost::format("[TEA.run_transduction_contig] ram_map: %s, %s\n") % ram_map_itr.first % ram_map_itr.second);
//				}
////				for(auto ref_id : ref_selected_seq_id) {
////					cout << (boost::format("[TEA.run_transduction_contig] Ref: %s\n") % ref_id);
////				}
////				for(auto repeat_id : repeat_selected_seq_id) {
////					cout << (boost::format("[TEA.run_transduction_contig] Repeat: %s\n") % repeat_id);
////				}
//				for(auto gold_id : gold) {
//					cout << (boost::format("[TEA.run_transduction_contig] Gold: %s\n") % gold_id);
//				}
//			}
//			cout << (boost::format("[TEA.run_transduction_contig] # Gold: %s\n") % gold.size()).str();
			string out_transduction = transduction_prefix + "." + chr + ".germline.contig.transduction";
			create_germline_transduction(out_transduction, gold, repeat_ram_aligned_map, ref_aligned_map, in_two_ram);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void TEA::create_contig_two_ram(const string& out_path, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// germline file
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& pram = fso[12];
		auto& nram = fso[13];
	    if ("0" != pram && "0" != nram) {
	    	auto& pclipped = fso[35];
	    	auto& nclipped = fso[36];
	    	auto& prammate = fso[37];
	    	auto& nrammate = fso[38];
	    	castle::StringUtils::trim(pclipped);
	    	castle::StringUtils::trim(nclipped);
	    	castle::StringUtils::trim(prammate);
	    	castle::StringUtils::trim(nrammate);
	        if (pclipped.size() > 0 && nclipped.size() > 0 && prammate.size() > 0 && nrammate.size() > 0) {
	            g << line << "\n";
	        }
	    }
	}
}

void TEA::create_fa_from_germline_contig(const string& out_path_1, const string& out_path_2, const string& in_path) {
//	string poly_a_str(5, 'A');
//	string poly_t_str(5, 'T');

//	const bool debug = "/dev/shm/pfg050_cancer/cluster_ram/pfg050_cancer.chr1.germline.contig.two_ram";
	const bool debug = false;

	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	int32_t min_polyAT = 5;
	int64_t line_id = 1;
	// two_ram file
	ifstream f(in_path, ios::binary);
	ofstream g(out_path_1, ios::binary);
	ofstream g1(out_path_2, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);

		auto& polyA_tag = fso[33];
		auto& polyT_tag = fso[34];
		auto& pclipped = fso[35];
		auto& nclipped = fso[36];
		auto& prammate = fso[37];
		auto& nrammate = fso[38];

		if(debug) {
			cout << line << "\n";
			cout << polyA_tag << "/" << polyT_tag << "/" << pclipped << "\n";
		}
		// polyA
		if("polyA" == polyA_tag && "-" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped.substr(0, n_pos) << "\n";
				}
			} else {
				if(pclipped.size() >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped << "\n";
				}
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate << "\n";
				}
			}
			if (nclipped.size() >= 15) {
				g << ">" << line_id << "_nclipped\n";
				g << nclipped << "\n";
			}
			if (prammate.size() >= 15) {
				g << ">" << line_id << "_pram\n";
				g << prammate << "\n";
			}
		}
		// polyT
		if("-" == polyA_tag && "polyT" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped.substr(n_pos) << "\n";
				}
			} else {
				if(nclipped.size() >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped << "\n";
				}
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate << "\n";
				}
			}
			if (pclipped.size() >= 15) {
				g << ">" << line_id << "_pclipped\n";
				g << pclipped << "\n";
			}
			if (nrammate.size() >= 15) {
				g << ">" << line_id << "_nram\n";
				g << nrammate << "\n";
			}
		}

		// output nram to g1
		{
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g1 << ">" << line_id << "_nram\n";
					g1 << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g1 << ">" << line_id << "_nram\n";
					g1 << nrammate << "\n";
				}
			}
		}
		// output pram to g1
		{
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_len;
				if(seq_size >= 15) {
					g1 << ">" << line_id << "_pram\n";
					g1 << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g1 << ">" << line_id << "_pram\n";
					g1 << prammate << "\n";
				}
			}
		}
		++line_id;
	}
}

int64_t TEA::rfind_the_last_of(int64_t& n_len, const string& str, const char chr) {
	n_len = 0;
	int64_t the_pos = str.size();
	--the_pos;
	for(; the_pos >= 0; --the_pos, ++n_len) {
		if(chr != str[the_pos]) {
			return the_pos;
		}
	}
	return the_pos;
}

int64_t TEA::find_the_last_of(int64_t& n_len, const string& str, const char chr) {
	n_len = 0;
	int64_t the_pos = 0;
	int64_t max_pos = str.size();
	for(; the_pos < max_pos; ++the_pos, ++n_len) {
		if(chr != str[the_pos]) {
			return the_pos;
		}
	}
	return -1;
}

void TEA::collect_two_ram_map(boost::unordered_map<int64_t, string>& two_ram_map, const string& in_path) {
	int64_t line_id = 1;
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& chr_name = fso[1];
		boost::replace_all(chr_name, "chr", "");
		auto& start_pos = fso[2];
		auto& end_pos = fso[3];
		string tmp = chr_name + "@" + start_pos + "@" + end_pos;
		two_ram_map[line_id] = tmp;
		++line_id;
	}
}

void TEA::collect_two_ram_seq_id_ref_set(boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map) {
	string line;
	vector<string> fso;
	vector<string> cols;
	vector<string> out;

	const char* delim_tab = "\t";
	const char* delim_underscore = "_";
	const char* delim_ampersand = "&";

	// read from .germline.contig.ref.aln.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
	    if ("*" != fso[2]) {
	        if ("0" == fso[1]) {
	            fso[1] = "+";
	        } else {
	            fso[1] = "-";
	        }
	        auto& seq_id = fso[0];
	        auto& chr_name = fso[2];
	        auto& start_pos = fso[3];
	        auto& cigar_str = fso[5];
	        auto& cigar_flag = fso[1];
	        string entry_value = chr_name + "," + start_pos + "," + cigar_str + "," + cigar_flag;
	        aligned_map[seq_id] = entry_value;
	    }
	    int64_t mapq = boost::lexical_cast<int64_t>(fso[4]);
	    if (mapq >= 25) {
	    	castle::StringUtils::c_string_multi_split(fso[0], delim_underscore, cols);
	    	int64_t cur_line_id = boost::lexical_cast<int64_t>(cols[0]);
	    	auto two_ram_itr = two_ram_map.find(cur_line_id);
	    	if(two_ram_map.end() == two_ram_itr) {
	    		continue;
	    	}
	    	auto& two_ram_entry_value = two_ram_itr->second;
	    	castle::StringUtils::c_string_multi_split(two_ram_entry_value, delim_ampersand, out);

	        bool on = false;
	        // is two_ram_id equal?
	        if (fso[2] == out[0]) {
	        	int64_t sam_start_pos = boost::lexical_cast<int64_t>(fso[3]);
	        	int64_t two_ram_start_pos = boost::lexical_cast<int64_t>(out[1]);
	        	int64_t two_ram_end_pos = boost::lexical_cast<int64_t>(out[2]);
	        	if ((sam_start_pos >= (two_ram_start_pos - 5000)) && sam_start_pos <= (two_ram_end_pos + 5000)) {
	                on = true;
	            }
	        }

	        if (!on) {
	        	two_ram_seq_id_set.insert(fso[0]);
	        }
	    }
	}
}

void TEA::collect_aln_sam_repeat(set<string>& repeat_selected_seq_id, boost::unordered_set<string>& repeat_two_ram_id, const string& in_path) {
	string line;
	vector<string> fso;
	vector<string> temp;
	const char* delim_tab = "\t";
	const char* delim_underscore = "_";

	// read from .germline.contig.repeat.aln.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		int64_t mapq = boost::lexical_cast<int64_t>(fso[4]);
		auto& ref_name = fso[2];
		if(37 == mapq && "RICKSHA" != ref_name) {
			auto& query_seq_id = fso[0];
			repeat_selected_seq_id.insert(query_seq_id);
			castle::StringUtils::c_string_multi_split(query_seq_id, delim_underscore, temp);
			repeat_two_ram_id.insert(temp[0]);
		}
	}
}

void TEA::collect_aln_ram_sam_repeat(boost::unordered_map<string, string>& repeat_ram_aligned_map, const string& in_path) {
	string line;
	vector<string> fso;
	vector<string> temp;
	const char* delim_tab = "\t";

	// read from .germline.contig.repeat.aln.ram.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);

		auto& ref_name = fso[2];

		if("*" == ref_name){
			continue;
		}
		auto& query_seq_id = fso[0];
		repeat_ram_aligned_map[query_seq_id] = ref_name;
	}
}

void TEA::create_germline_transduction(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& repeat_ram_aligned_map, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .germline.contig.two_ram
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.transduction
	ofstream g(out_path, ios::binary);
	int64_t i = 1;

//	const bool is_chr = string::npos != in_path.find("chr12");
	while(getline(f, line, '\n')) {
		string str_i = boost::lexical_cast<string>(i);
		if(gold.end() == gold.find(str_i)) {
			++i;
			continue;
		}
//		const bool debug = is_chr && i == 100;
//		if(debug) {
//			cout << "[TEA.create_germline_transduction] " << line << "\n";
//		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		int64_t on = 0;
		auto& polyA_tag = fso[33];
		auto& polyT_tag = fso[34];
		set<string> rf;

		// polyA
		if("polyA" == polyA_tag && "-" == polyT_tag) {
			string the_pram_id = str_i + "_pram";
			string the_nclipped_id = str_i + "_nclipped";
			auto pram_itr = repeat_ram_aligned_map.find(the_pram_id);
			if (repeat_ram_aligned_map.end() != pram_itr) {
				rf.insert(pram_itr->second);
			}
			auto nclipped_itr = repeat_ram_aligned_map.find(the_nclipped_id);
			if (repeat_ram_aligned_map.end() != nclipped_itr) {
				rf.insert(nclipped_itr->second);
			}
			if(0 != rf.size()) {
				on = 1;
			}
		}
		// polyT
		if("-" == polyA_tag && "polyT" == polyT_tag) {
			string the_nram_id = str_i + "_nram";
			string the_pclipped_id = str_i + "_pclipped";
			auto nram_itr = repeat_ram_aligned_map.find(the_nram_id);
			if (repeat_ram_aligned_map.end() != nram_itr) {
				rf.insert(nram_itr->second);
			}
			auto pclipped_itr = repeat_ram_aligned_map.find(the_pclipped_id);
			if (repeat_ram_aligned_map.end() != pclipped_itr) {
				rf.insert(pclipped_itr->second);
			}
			if(0 != rf.size()) {
				on = 1;
			}
		}

        if (1 == on) {
            g << fso[0];
            uint64_t tt = 1;
            while (tt <= 7) {
                g << "\t" << fso[tt];
                ++tt;
            }
            g << "\t";
            for (auto& re : rf) {
                g << re <<  ",";
            }
            g << "\t";
            for (auto& re : rf) {
				g << re <<  ",";
			}
            tt = 10;
            while (tt < fso.size()) {
                g << "\t" << fso[tt];
                ++tt;
            }

            string pc = str_i + "_pclipped";
            auto pclipped_itr = ref_aligned_map.find(pc);
            if (repeat_ram_aligned_map.end() == pclipped_itr) {
            	g << "\tna";
            } else {
            	g << "\t" << pclipped_itr->second;
            }

            string nc = str_i + "_nclipped";
            auto nclipped_itr = ref_aligned_map.find(nc);
            if (repeat_ram_aligned_map.end() == nclipped_itr) {
            	g << "\tna";
            } else {
            	g << "\t" << nclipped_itr->second;
            }

            string pram = str_i + "_pram";
            auto pram_itr = ref_aligned_map.find(pram);
			if (repeat_ram_aligned_map.end() == pram_itr) {
				g << "\tna";
			} else {
				g << "\t" << pram_itr->second;
			}

            string nram = str_i + "_nram";
            auto nram_itr = ref_aligned_map.find(nram);
			if (repeat_ram_aligned_map.end() == nram_itr) {
				g << "\tna";
			} else {
				g << "\t" << nram_itr->second;
			}
			g << "\n";
        }
        ++i;
	}
}

void TEA::run_orphan() {
	options.rasym = "um";
	set<string> chrs;
	load_chr(chrs);
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string orphan_dir = options.prefix + "/orphan_" + options.rasym + "m";
	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		orphan_dir = options.working_prefix + "/orphan_" + options.rasym + "m";
	}
	if(!boost::filesystem::exists(orphan_dir)) {
		boost::filesystem::create_directories(orphan_dir);
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	string orphan_prefix = orphan_dir + "/" + naive_prefix;

	cout << (boost::format("[TEA.run_transduction] copy %sm files to orphan folder and create the discord.bed file\n") % options.rasym).str();
	vector<function<void()> > tasks;

	string in_discord = options.prefix + ".discord";
	string out_discord = orphan_prefix + ".discord.bed";
	string in_cluster = options.prefix + ".clusters";
	if (!options.working_dir.empty()) {
		in_discord = options.working_prefix + ".discord";
		in_cluster = options.working_prefix + ".clusters";
	}
	create_discord_bed(out_discord, in_discord);

	string out_discord_intersected = orphan_prefix + ".discord.bed.intersected";
	string bedtools_cmd = (boost::format("bedtools pairtopair -is -rdn -a %s -b %s > %s") % out_discord % out_discord % out_discord_intersected).str();
	cout << (boost::format("[TEA.run_transduction] %s\n") % bedtools_cmd).str();
	system(bedtools_cmd.c_str());

	set<string> oo;
	set<string> gene;

	cout << "[TEA.run_transduction] collect gene sets\n";

	collect_gene_sets(oo, gene, out_discord_intersected);

	cout << "[TEA.run_transduction] collect rescue sets\n";
	set<string> rescue;
	collect_rescue_sets(rescue, gene, out_discord_intersected);

	set<string> ooo;
	set_difference(oo.begin(), oo.end(), gene.begin(), gene.end(), inserter(ooo, ooo.begin()));
	ooo.insert(rescue.begin(), rescue.end());

	cout << "[TEA.run_transduction] create bed intersected filtered\n";
	string out_discord_bed_intersected_filtered = orphan_prefix + ".discord.bed.intersected.filtered";
	create_bed_filtered_intersected(out_discord_bed_intersected_filtered, ooo, out_discord_intersected);

	cout << "[TEA.run_transduction] collect cluster id pos map\n";
	boost::unordered_map<string, string> cluster_id_pos_map;
	collect_cluster_id_pos_map(cluster_id_pos_map, in_cluster);

	cout << "[TEA.run_transduction] create bed intersected filtered insertion\n";
	string out_bed_intersected_filtered_insertion = orphan_prefix + ".discord.bed.intersected.filtered.insertion";
	create_intersected_filtered_insertion(out_bed_intersected_filtered_insertion, cluster_id_pos_map, out_discord_bed_intersected_filtered);

	boost::unordered_map<string, string> cluster_id_pair_map;
	collect_cluster_id_pair_map(cluster_id_pair_map, out_bed_intersected_filtered_insertion);

	string out_umm = orphan_prefix + ".umm";
	create_orphan_umm_from_cluster(out_umm, cluster_id_pair_map, in_cluster);
}

void TEA::collect_gene_sets(set<string>& oo, set<string>& gene, const string& in_path) {
	set<string> o;
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		vector<int64_t> a(2);
		a[0] = boost::lexical_cast<int64_t>(fso[6]);
	    a[1] = boost::lexical_cast<int64_t>(fso[16]);
	    sort(a.begin(), a.end());
	    string temp = (boost::format("%s@%s") % a[0] % a[1]).str();
	    if(o.end() != o.find(temp)) {
	    	continue;
	    }
		o.insert(temp);

		if (oo.end() != oo.find(fso[6])) {
			gene.insert(fso[6]);
		}
		if (oo.end() != oo.find(fso[16])) {
			gene.insert(fso[16]);
		}
		oo.insert(fso[6]);
		oo.insert(fso[16]);
	}
}

void TEA::collect_rescue_sets(set<string>& rescue, set<string>& gene, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	for(auto& val : gene) {
		ifstream f(in_path, ios::binary);
		set<string> direc;
		uint64_t c = 0;
		while(getline(f, line, '\n')) {
			castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
			if (fso[6] == val) {
				string a_val = (boost::format("%s@%s@%s@%s") % fso[8] % fso[9] % fso[18] % fso[19]).str();
				direc.insert(a_val);
				++c;
			}
			if (fso[16] == val) {
				string a_val = (boost::format("%s@%s@%s@%s") % fso[18] % fso[19] % fso[8] % fso[9]).str();
				direc.insert(a_val);
				++c;
			}
		}
		if (direc.size() == (c/2)) {
			rescue.insert(val);
		}
	}
}

void TEA::create_bed_filtered_intersected(const string& out_path, const set<string>& ooo, const string& in_path) {
	set<string> o;
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		vector<int64_t> a(2);
		a[0] = boost::lexical_cast<int64_t>(fso[6]);
		a[1] = boost::lexical_cast<int64_t>(fso[16]);
		sort(a.begin(), a.end());
		string temp = (boost::format("%s@%s") % a[0] % a[1]).str();
		if(o.end() != o.find(temp)) {
			continue;
		}
		o.insert(temp);

		if (ooo.end() != ooo.find(fso[6]) && ooo.end() != ooo.find(fso[16])) {
			g << line << "\n";
		}
	}
}

void TEA::collect_cluster_id_pos_map(boost::unordered_map<string, string>& cluster_id_pos_map, const string& in_path) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.collect_cluster_id_pos_map");
	checker.start();

	vector<function<void()> > tasks;

	const int64_t BLOCK_SIZE = 4 * 1024 * 1024;
	const int64_t ref_file_size = castle::IOUtils::get_file_size(in_path);
	int64_t n_blocks = (ref_file_size / (double) BLOCK_SIZE) + 1;
	cout << (boost::format("[TEA.collect_cluster_id_pos_map] Ref. File size: %d\n") % ref_file_size).str();
	cout << (boost::format("[TEA.collect_cluster_id_pos_map] # blocks: %d\n") % n_blocks).str();
	vector<int64_t> block_boundary;
	block_boundary.resize(n_blocks);
	block_boundary[0] = 0;
	block_boundary[n_blocks - 1] = numeric_limits<int64_t>::max();
	for (int64_t block_id = 1; block_id < n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id, BLOCK_SIZE] {
			int64_t cur_boundary_pos = block_id * BLOCK_SIZE;
			string line;
			ifstream in(in_path, ios::binary);
			in.seekg(cur_boundary_pos, ios::beg);
			getline(in, line, '\n');
			cur_boundary_pos += line.size() + 1;
			string previous_pid;
			while(getline(in, line, '\n')) {
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
	vector<boost::unordered_map<string, string>> cluster_id_pos_map_list(n_blocks - 1);

	for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_boundary_pos = block_boundary[block_id];
			int64_t cur_pos = cur_boundary_pos;
			int64_t next_boundary_pos = block_boundary[block_id + 1];
			string line;
			vector<string> fso;
			vector<int64_t> a;
			vector<int64_t> b;
			string prev_id;
			const char* delim_tab = "\t";
			auto& local_cluster_id_pos_map = cluster_id_pos_map_list[block_id];
			ifstream in(in_path, ios::binary);
			in.seekg(cur_boundary_pos, ios::beg);
			while (getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
				int64_t f_pos1 = boost::lexical_cast<int64_t>(fso[9]);
				int64_t f_pos2 = boost::lexical_cast<int64_t>(fso[13]);
				a.push_back(f_pos1);
				b.push_back(f_pos2);

				if(prev_id.empty()) {
					prev_id == fso[0];
				}
				if(prev_id != fso[0]) {
					if(a.size() > 0 && b.size() > 0) {
						auto a_biggest = max_element(a.begin(), a.end());
						auto b_biggest = max_element(b.begin(), b.end());
						local_cluster_id_pos_map[(boost::format("%s@%s") % fso[0] % fso[1]).str()] = (boost::format("%s@%s") % *a_biggest % *b_biggest).str();
					}
					a.clear();
					b.clear();
					prev_id = fso[0];
				}
				if(cur_pos >= next_boundary_pos) {
					break;
				}
			}
			if(a.size() > 0 && b.size() > 0) {
				auto a_biggest = max_element(a.begin(), a.end());
				auto b_biggest = max_element(b.begin(), b.end());
				local_cluster_id_pos_map[(boost::format("%s@%s") % fso[0] % fso[1]).str()] = (boost::format("%s@%s") % *a_biggest % *b_biggest).str();
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
		auto& local_cluster_id_pos_map = cluster_id_pos_map_list[block_id];
		for(auto& an_entry : local_cluster_id_pos_map) {
			cluster_id_pos_map[an_entry.first] = an_entry.second;
		}
	}
	cout << checker;
}

void TEA::create_intersected_filtered_insertion(const string& out_path, const boost::unordered_map<string, string> cluster_id_pos_map, const string& in_path) {
	string line;
	vector<string> fso;
	vector<string> temp1;
	vector<string> temp2;
	const char* delim_tab = "\t";
	const char* delim_ampersand = "@";
	// from .discord.bed.intersected.filtered
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		string cl_id1 = (boost::format("%s@%s") % fso[6] % fso[7]).str();
		string cl_id2 = (boost::format("%s@%s") % fso[16] % fso[17]).str();
		auto itr1 = cluster_id_pos_map.find(cl_id1);
		auto itr2 = cluster_id_pos_map.find(cl_id2);
		if(itr1 == cluster_id_pos_map.end() || itr2 == cluster_id_pos_map.end()) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(itr1->second, delim_ampersand, temp1);
		castle::StringUtils::c_string_multi_split(itr2->second, delim_ampersand, temp2);

		int64_t pos1 = boost::lexical_cast<int64_t>(temp1[0]);
		int64_t pos2 = boost::lexical_cast<int64_t>(temp2[0]);
		int64_t pos3 = boost::lexical_cast<int64_t>(temp1[1]);
		int64_t pos4 = boost::lexical_cast<int64_t>(temp2[1]);
		string fir;
		if (pos1 < pos2) {
			fir = fso[8] + fso[18];
		} else {
			fir = fso[18] + fso[8];
		}
		string sec;

		if (pos3 < pos4) {
			sec = fso[9] + fso[19];
		} else {
			sec = fso[19] + fso[9];
		}
		if ("+-" == fir && "-+" == sec) {
			g << line << "\t1\tno\n";
		}
		if ("-+" == fir && "+-" == sec) {
			g << line << "\t2\tno\n";
		}

		//inversion
		if ("+-" == fir && ("--" == sec || "++" == sec)) {
			g << line << "\t1\tinv\n";
		}
		if (("--" == fir || "++" == fir) && "+-" == sec) {
			g << line << "\t2\tinv\n";
		}
	}
}

void TEA::collect_cluster_id_pair_map(boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// from .discord.bed.intersected.filtered.insertion
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		cluster_id_pair_map[fso[6]] = fso[20];
		cluster_id_pair_map[fso[16]] = fso[20];
	}
}

void TEA::create_orphan_umm_from_cluster(const string& out_path, const boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.create_orphan_umm_from_cluster");
	checker.start();

	vector<function<void()> > tasks;

	const int64_t BLOCK_SIZE = 4 * 1024 * 1024;
	const int64_t ref_file_size = castle::IOUtils::get_file_size(in_path);
	int64_t n_blocks = (ref_file_size / (double) BLOCK_SIZE) + 1;
	vector<int64_t> block_boundary;
	block_boundary.resize(n_blocks);
	block_boundary[0] = 0;
	block_boundary[n_blocks - 1] = numeric_limits<int64_t>::max();
	for (int64_t block_id = 1; block_id < n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id, BLOCK_SIZE] {
			int64_t cur_boundary_pos = block_id * BLOCK_SIZE;
			string line;
			ifstream in(in_path, ios::binary);
			in.seekg(cur_boundary_pos, ios::beg);
			getline(in, line, '\n');
			cur_boundary_pos += line.size() + 1;
			string previous_pid;
			while(getline(in, line, '\n')) {
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
	vector<boost::unordered_map<string, string>> cluster_id_pos_map_list(n_blocks - 1);

	vector<string> out_files(n_blocks - 1);
	for (int64_t block_id = 0; block_id < n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);

			int64_t cur_boundary_pos = block_boundary[block_id];
			int64_t cur_pos = cur_boundary_pos;
			int64_t next_boundary_pos = block_boundary[block_id + 1];
			string line;
			vector<string> fso;
			vector<int64_t> a;
			vector<int64_t> b;
			string prev_id;
			const char* delim_tab = "\t";
			boost::unordered_set<string> o;

			string out_filename = out_path + "." + str_block_id;
			out_files[block_id] = out_filename;

			ifstream in(in_path, ios::binary);
			in.seekg(cur_boundary_pos, ios::beg);
			ofstream g(out_filename, ios::binary);
			while (getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
				auto& p_id = fso[0];
				auto id_pair_itr = cluster_id_pair_map.find(p_id);
				if (cluster_id_pair_map.end() == id_pair_itr) {
					if(cur_pos >= next_boundary_pos) {
						break;
					}
					continue;
				}

				auto& strand1 = fso[8];
				auto& strand2 = fso[12];
				boost::replace_all(strand1, "1", "");
				boost::replace_all(strand2, "1", "");
				auto& seq_id = fso[5];
				int64_t rlen1 = boost::lexical_cast<int64_t>(fso[10]);
				int64_t rlen2 = boost::lexical_cast<int64_t>(fso[14]);

				if(o.end() == o.find(seq_id)) {
					if(rlen1 >= rlen2 && "1" == id_pair_itr->second) {
						string a_line = (boost::format("%s\t%s\t%s%s\tx\n") % seq_id % fso[7] % fso[8] % fso[9]).str();
						g << a_line;
						o.insert(fso[5]);
					}
					if(rlen1 <= rlen2 && "2" == id_pair_itr->second) {
						string a_line = (boost::format("%s\t%s\t%s%s\tx\n") % seq_id % fso[11] % fso[12] % fso[13]).str();
						g << a_line;
						o.insert(fso[5]);
					}
				}
				if(cur_pos >= next_boundary_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge(out_path, out_files, n_cores, true);
	cout << checker;
}

void TEA::run_orphan_contig() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.run_orphan_contig");
	checker.start();
	set<string> chrs;
	load_chr(chrs);
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string orphan_dir = options.prefix + "/orphan_" + options.rasym + "m";
	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		orphan_dir = options.working_prefix + "/orphan_" + options.rasym + "m";
	}
	if(!boost::filesystem::exists(orphan_dir)) {
		boost::filesystem::create_directories(orphan_dir);
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	string orphan_prefix = orphan_dir + "/" + naive_prefix;

	vector<function<void()> > tasks;
	boost::unordered_map<string, int64_t> chr_ids;

	int64_t chr_id = 0;
	for(auto& chr: chrs) {
		chr_ids[chr] = chr_id;
		++chr_id;
	}
	// create two ram files
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			string in_germline_contig = cl_prefix + "." + chr + ".germline.contig";
			string out_getrmline_contig_two_ram = orphan_prefix + "." + chr + ".germline.contig.two_ram";
			create_contig_two_ram(out_getrmline_contig_two_ram, in_germline_contig);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	// create insertion sequences
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			string in_two_ram = orphan_prefix + "." + chr + ".germline.contig.two_ram";
			string out_two_ram_fa = orphan_prefix + "." + chr + ".germline.contig.two_ram.fa";
			create_orphan_fa_from_germline_contig(out_two_ram_fa, in_two_ram);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	cout << "[TEA.run_transduction_contig] align sequences\n";
	// align insertion sequences
	for(auto& chr: chrs) {
		// two_ram against repeat ref
		tasks.push_back([&]{
			string in_two_ram_fa = orphan_prefix + "." + chr + ".germline.contig.two_ram.fa";
			string out_repeat_aln_sam = orphan_prefix + "." + chr + ".germline.contig.repeat.aln.sam";
//			if(castle::IOUtils::get_file_size(out_repeat_aln_sam) > 0) {
//				return;
//			}
			string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.repeat_reference % in_two_ram_fa % out_repeat_aln_sam % options.repeat_reference % in_two_ram_fa).str();
//			{
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << aln_cmd << "\n";
//			}
			system(aln_cmd.c_str());
		});
		// two_ram against ref
		tasks.push_back([&]{
			string in_two_ram_fa = orphan_prefix + "." + chr + ".germline.contig.two_ram.fa";
			string out_ref_aln_sam = orphan_prefix + "." + chr + ".germline.contig.ref.aln.sam";
//			if(castle::IOUtils::get_file_size(out_ref_aln_sam) > 0) {
//				return;
//			}
			string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.human_reference % in_two_ram_fa % out_ref_aln_sam % options.human_reference % in_two_ram_fa).str();
//			{
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << aln_cmd << "\n";
//			}
			system(aln_cmd.c_str());
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			boost::unordered_map<int64_t, string> two_ram_map;
			string in_two_ram = orphan_prefix + "." + chr + ".germline.contig.two_ram";
			collect_two_ram_map(two_ram_map, in_two_ram);
			string in_ref_aln_sam = orphan_prefix + "." + chr + ".germline.contig.ref.aln.sam";
			boost::unordered_map<string, string> ref_aligned_map;
			set<string> ref_selected_seq_id;
			collect_two_ram_seq_id_ref_set(ref_aligned_map, ref_selected_seq_id, in_ref_aln_sam, two_ram_map);
			string in_repeat_aln_sam = orphan_prefix + "." + chr + ".germline.contig.repeat.aln.sam";
			set<string> repeat_selected_seq_id;
			boost::unordered_set<string> repeat_two_ram_id;
			collect_aln_sam_repeat(repeat_selected_seq_id, repeat_two_ram_id, in_repeat_aln_sam);

			set<string> candidate;
			set_difference(ref_selected_seq_id.begin(), ref_selected_seq_id.end(), repeat_selected_seq_id.begin(), repeat_selected_seq_id.end(),
				inserter(candidate, candidate.begin()));
			set<string> gold;
			{
				string line;
				vector<string> temp;
				const char* delim_underscore = "_";

				for(auto& gene : candidate) {
					castle::StringUtils::c_string_multi_split(gene, delim_underscore, temp);
					auto& two_ram_line_id = temp[0];
					if (repeat_two_ram_id.end() == repeat_two_ram_id.find(two_ram_line_id)) {
						gold.insert(temp[0]);
					}
				}
			}
			string out_transduction = orphan_prefix + "." + chr + ".germline.contig.orphan";
			create_germline_orphan(out_transduction, gold, ref_aligned_map, in_two_ram);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}
void TEA::create_orphan_fa_from_germline_contig(const string& out_path, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	int32_t min_polyAT = 5;
	int64_t line_id = 1;
	// two_ram file
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);

		auto& polyA_tag = fso[33];
		auto& polyT_tag = fso[34];
		auto& pclipped = fso[35];
		auto& nclipped = fso[36];
		auto& prammate = fso[37];
		auto& nrammate = fso[38];

		// polyA
		if("polyA" == polyA_tag && "-" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped.substr(0, n_pos) << "\n";
				}
			} else {
				if(pclipped.size() >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped << "\n";
				}
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate << "\n";
				}
			}
		}
		// polyT
		if("-" == polyA_tag && "polyT" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped.substr(n_pos) << "\n";
				}
			} else {
				if(nclipped.size() >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped << "\n";
				}
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate << "\n";
				}
			}
		}
		// Neither poly
		if("-" == polyA_tag && "-" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped.substr(0, n_pos) << "\n";
				}
			} else {
				if(pclipped.size() >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped << "\n";
				}
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate << "\n";
				}
			}

			n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped.substr(n_pos) << "\n";
				}
			} else {
				if(nclipped.size() >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped << "\n";
				}
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate << "\n";
				}
			}
		}
		++line_id;
	}
}


void TEA::create_germline_orphan(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .germline.contig.two_ram
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.transduction
	ofstream g(out_path, ios::binary);
	int64_t i = 1;

//	const bool is_chr = string::npos != in_path.find("chr12");
	while(getline(f, line, '\n')) {
		string str_i = boost::lexical_cast<string>(i);
		if(gold.end() == gold.find(str_i)) {
			++i;
			continue;
		}
//		const bool debug = is_chr && i == 100;
//		if(debug) {
//			cout << "[TEA.create_germline_transduction] " << line << "\n";
//		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		if(fso.size() == 40) {
			g << line;
		}
		if(fso.size() == 39) {
			g << line << "\tNA";
		}
		string pc = str_i + "_pclipped";
		auto pclipped_itr = ref_aligned_map.find(pc);
		if (ref_aligned_map.end() == pclipped_itr) {
			g << "\tna";
		} else {
			g << "\t" << pclipped_itr->second;
		}

		string nc = str_i + "_nclipped";
		auto nclipped_itr = ref_aligned_map.find(nc);
		if (ref_aligned_map.end() == nclipped_itr) {
			g << "\tna";
		} else {
			g << "\t" << nclipped_itr->second;
		}

		string pram = str_i + "_pram";
		auto pram_itr = ref_aligned_map.find(pram);
		if (ref_aligned_map.end() == pram_itr) {
			g << "\tna";
		} else {
			g << "\t" << pram_itr->second;
		}

		string nram = str_i + "_nram";
		auto nram_itr = ref_aligned_map.find(nram);
		if (ref_aligned_map.end() == nram_itr) {
			g << "\tna";
		} else {
			g << "\t" << nram_itr->second;
		}
		g << "\n";
        ++i;
	}
}

void TEA::clean() {
	if(!options.is_cleaning) {
		return;
	}
	string prefix = options.prefix;
	if(!options.working_dir.empty()) {
		prefix = options.working_prefix;
	}
	string naive_prefix = options.naive_prefix;
	boost::filesystem::path root_path(prefix);
	root_path = boost::filesystem::absolute(root_path.parent_path());
	string including_regex = (boost::format("%s.*") % naive_prefix).str();

	cout << (boost::format("[TEA.clean] pattern: %s\n") % including_regex).str();
	cout << (boost::format("[TEA.clean] path: %s\n") % root_path.string()).str();

	vector<string> removal_paths;
//	castle::IOUtils::selective_search(root_path, including_regex, removal_paths);
//
//	vector<string> excluding_set;
//	excluding_set.push_back(".firststat");
//	excluding_set.push_back(".bam");
//	excluding_set.push_back(".bam.bai");
//	excluding_set.push_back(".isize");
//	excluding_set.push_back(".rl");
//	excluding_set.push_back(".ram");
//	excluding_set.push_back(".ram.bam");
//	excluding_set.push_back(".ram.bam.bai");
//	excluding_set.push_back(".softclips.consd.bam.bai");
//	excluding_set.push_back(".softclips.consd.bam");
//	excluding_set.push_back(".softclips.consd.cpos");
//	excluding_set.push_back(".pre.log");
//	excluding_set.push_back(".dre.log");
//	excluding_set.push_back(".pdf");
//	excluding_set.push_back(".cl.sorted.bam");
//	excluding_set.push_back(".cl.sorted.bam.bai");
//	excluding_set.push_back(".disc.num.bam");
//	excluding_set.push_back(".disc.num.bam.bai");
//	excluding_set.push_back(".cl.sorted.disc.bam");
//	excluding_set.push_back(".cl.sorted.disc.bam.bai");
//	excluding_set.push_back(".mapped_um.bam");
//	excluding_set.push_back(".mapped_um.bam.bai");
//	excluding_set.push_back("*.germline");
//	excluding_set.push_back("*.contig");
//	excluding_set.push_back("*.cluster");
//	excluding_set.push_back(".clusters");
//	excluding_set.push_back(".discord");
//	excluding_set.push_back("*.cluster.raw");
//	excluding_set.push_back("*.clipped");
//
//	removal_paths.erase(remove_if(removal_paths.begin(), removal_paths.end(),
//	        [&](const string& o) {
//		bool found_match = false;
//		for(auto a_pattern : excluding_set) {
//			if('*' == a_pattern[0]) {
//				a_pattern = a_pattern.substr(1);
//				if(boost::ends_with(o, a_pattern)) {
//					found_match = true;
//					break;
//				}
//			} else {
//				a_pattern = naive_prefix + a_pattern;
//				if(boost::ends_with(o, a_pattern)) {
//					found_match = true;
//					break;
//				}
//			}
//
//		}
//		return found_match;
//	}), removal_paths.end());

	string contig_dir = prefix + "/assembly_" + options.rasym + "m";
	string original_bam = options.prefix + ".bam";

	if(boost::filesystem::exists(contig_dir)) {
		removal_paths.push_back(contig_dir);
	}

	vector<string> read_groups;
	map<string, int64_t> read_groups_reverse_index;
	map<string, int64_t> ref_reverse_index;
	{
		BamTools::BamReader local_reader;
		string an_index_path;
		get_bai_index_path(original_bam, an_index_path);
		if (!local_reader.Open(original_bam, an_index_path)) {
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
		string cl_fq1 = prefix + "/none_1.fq";
		if (boost::filesystem::exists(cl_fq1)) {
			string last_key = "none";
			read_groups_reverse_index[last_key] = read_groups.size();
			read_groups.push_back(last_key);
		}
		local_reader.Close();
	}

	for (auto& rgname : read_groups) {
		string cl_fq1 = prefix + "/" + rgname + "_1.fq";
		if(boost::filesystem::exists(cl_fq1)) {
			removal_paths.push_back(cl_fq1);
		}
		string cl_fq2 = prefix + "/" + rgname + "_2.fq";
		if(boost::filesystem::exists(cl_fq2)) {
			removal_paths.push_back(cl_fq2);
		}

		for (int64_t block_id = 0;; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string cl_fq1 = prefix + "/" + rgname + "_1." + str_block_id + ".fq";
			string cl_fq2 = prefix + "/" + rgname + "_2." + str_block_id + ".fq";
			string cl_sai1 = prefix + "/" + rgname + "_1.fq.sai." + str_block_id;
			string cl_sai2 = prefix + "/" + rgname + "_2.fq.sai." + str_block_id;

			string cl_tmp_sam = prefix + "/" + rgname + ".tmp.sam." + str_block_id;
			string cl_sam = prefix + "/" + rgname + ".sam." + str_block_id;
			string bwa_err = prefix + "/" + rgname + ".fq.bwa.err." + str_block_id;
			string cl_fq1_bwa_err = prefix + "/" + rgname + "_1.fq.bwa.err." + str_block_id;
			string cl_fq2_bwa_err = prefix + "/" + rgname + "_2.fq.bwa.err." + str_block_id;

			bool has_file = false;
			if(boost::filesystem::exists(cl_fq1)) {
				has_file = true;
				removal_paths.push_back(cl_fq1);
			}
			if(boost::filesystem::exists(cl_fq2)) {
				has_file = true;
				removal_paths.push_back(cl_fq2);
			}
			if(boost::filesystem::exists(cl_sai1)) {
				has_file = true;
				removal_paths.push_back(cl_sai1);
			}
			if(boost::filesystem::exists(cl_sai2)) {
				has_file = true;
				removal_paths.push_back(cl_sai2);
			}
			if(boost::filesystem::exists(cl_tmp_sam)) {
				has_file = true;
				removal_paths.push_back(cl_tmp_sam);
			}
			if(boost::filesystem::exists(bwa_err)) {
				has_file = true;
				removal_paths.push_back(bwa_err);
			}
			if(boost::filesystem::exists(cl_fq1_bwa_err)) {
				has_file = true;
				removal_paths.push_back(cl_fq1_bwa_err);
			}
			if(boost::filesystem::exists(cl_fq2_bwa_err)) {
				has_file = true;
				removal_paths.push_back(cl_fq2_bwa_err);
			}
			if(boost::filesystem::exists(cl_sam)) {
				has_file = true;
				removal_paths.push_back(cl_sam);
			}
			if(!has_file) {
				break;
			}
		}
	}
	string tmp_dir = prefix + "/tmp";
	if(boost::filesystem::exists(tmp_dir)) {
		removal_paths.push_back(tmp_dir);
	}
	string cl_sorted = prefix + ".cl.sorted";
	if(boost::filesystem::exists(cl_sorted)) {
		removal_paths.push_back(cl_sorted);
	}

	string tmp_name = prefix + ".assemblies.log";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".assemblies.target";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".blacklist.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".blacklist.gz.bak";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.fq.gz.bwa.err";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.fq.gz.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.ra.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.ra.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}

	tmp_name = prefix + ".cl.disc_2.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.fq.gz.bwa.err";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.fq.gz.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.ra.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.ra.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}

	tmp_name = prefix + ".cl.disc.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}

	tmp_name = prefix + ".cl.dup.bam.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.dup.bam.tmp.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.disc.num.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.disc.num.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.disc.num.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.disc.sorted.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.dup.num.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.dup.num.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.sorted.dup.sorted.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.fq.gz.bwa.err";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.fq.gz.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.ra.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.ra.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.ra.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_1.ra.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.ra.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_1.ra.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.fq.gz.bwa.err";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.fq.gz.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.ra.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.ra.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.ra.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc_2.ra.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.ra.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc_2.ra.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}

	tmp_name = prefix + ".disc.bam.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.bam.tmp.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc.bam.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".cl.disc.bam.tmp.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.num.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.num.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.num.bam.bni";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".disc.sorted.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".dup.bam.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".dup.bam.tmp.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".dup.num.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".dup.num.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".dup.sorted.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".insert.r";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_sc.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_sc.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_sc.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_sc.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_um.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_um.bam.bai";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_um.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".mapped_um.tmp.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".merged.ram.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".ram.raw.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".unmapped.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".unmapped.fq.gz.bak";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".unmapped.rdist";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".softclips.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".softclips.fq.gz.bak";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".sr.1.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".sr.1.fq.gz.bak";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".sr.2.fq.gz";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".sr.2.fq.gz.bak";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".softclips.rdist";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".softclips.consd.raw.sam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	tmp_name = prefix + ".softclips.consd.raw.bam";
	if(boost::filesystem::exists(tmp_name)) {
		removal_paths.push_back(tmp_name);
	}
	cout << (boost::format("[TEA.clean] %s entries will be removed\n") % removal_paths.size()).str();
//	for(auto& a_path: removal_paths) {
//		cout << a_path << "\n";
//	}
	castle::IOUtils::remove_files(removal_paths, n_cores);
}

void TEA::post_process() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.post_process");
	checker.start();
	set<string> chrs;
	load_chr(chrs);
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string transduction_dir = options.prefix + "/transduction_" + options.rasym + "m";
	string orphan_dir = options.prefix + "/orphan_" + options.rasym + "m";
	string germline_tmp_dir = options.prefix + "/germline_tmp_" + options.rasym + "m";
	string germline_result_dir = options.prefix + "/germline_result_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		transduction_dir = options.working_prefix + "/transduction_" + options.rasym + "m";
		orphan_dir = options.working_prefix + "/orphan_" + options.rasym + "m";
		germline_tmp_dir = options.working_prefix + "/germline_tmp_" + options.rasym + "m";
		germline_result_dir = options.working_prefix + "/germline_result_" + options.rasym + "m";
	}

	if(!boost::filesystem::exists(transduction_dir)) {
		cout << "[TEA.post_process] Your transduction directory is missing\n";
		return;
	}
	if(!boost::filesystem::exists(orphan_dir)) {
		cout << "[TEA.post_process] Your orphan directory is missing\n";
		return;
	}
	if(!boost::filesystem::exists(germline_tmp_dir)) {
		boost::filesystem::create_directories(germline_tmp_dir);
	}
	if(!boost::filesystem::exists(germline_result_dir)) {
		boost::filesystem::create_directories(germline_result_dir);
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;
	string orphan_prefix = orphan_dir + "/" + naive_prefix;
	string transduction_prefix = transduction_dir + "/" + naive_prefix;
	string germline_tmp_prefix = germline_tmp_dir + "/" + naive_prefix;
	string germline_result_prefix = germline_result_dir + "/" + naive_prefix;

	vector<function<void()> > tasks;
	boost::unordered_map<string, int64_t> chr_ids;

	int64_t chr_id = 0;
	for(auto& chr: chrs) {
		chr_ids[chr] = chr_id;
		++chr_id;
	}
	// create orphan list
	for(auto& chr: chrs) {
		tasks.push_back([&]{
			string in_contig_orphan = orphan_prefix + "." + chr + ".germline.contig.orphan";
			string out_contig_orphan_list = germline_result_prefix + "." + chr + ".germline.contig.orphan.list";
			set<string> orphan;
			// create orphan.list
			create_orphan_list(out_contig_orphan_list, orphan, in_contig_orphan);
			// create transduction.filtered
			set<string> trans;
			string in_contig_transduction = transduction_prefix + "." + chr + ".germline.contig.transduction";
			string out_contig_transduction_filtered = germline_tmp_prefix + "." + chr + ".germline.contig.transduction.filtered";
			create_transduction_filtered(out_contig_transduction_filtered, trans, in_contig_transduction, orphan);
			set<string> allset;
			allset.insert(orphan.begin(), orphan.end());
			allset.insert(trans.begin(), trans.end());

			// create contig.filtered.fa
			string in_germline_contig = cl_prefix + "." + chr + ".germline.contig";
			string out_getrmline_contig_filtered_fa = germline_tmp_prefix + "." + chr + ".germline.contig.filtered.fa";
			create_contig_filtered_fa(out_getrmline_contig_filtered_fa, in_germline_contig);

			// align filtered fa against repeat sequences
			string out_repeat_aln_sam = germline_tmp_prefix + "." + chr + ".germline.contig.filtered.fa.aln.sam";
			if(castle::IOUtils::get_file_size(out_repeat_aln_sam) <= 0) {
				string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.repeat_reference % out_getrmline_contig_filtered_fa % out_repeat_aln_sam % options.repeat_reference % out_getrmline_contig_filtered_fa).str();
				system(aln_cmd.c_str());
			}
			// collect o (repeat selected seq id(e.g. 5_nram_ALR_ALPHA -> '5') from sam file
			set<string> o;
			collect_aln_repeat_selected_seq_id(o, out_repeat_aln_sam);
			// create germline contig tmp
			string out_germline_contig_tmp = germline_tmp_prefix + "." + chr + ".germline.contig.tmp";
			create_germline_contig_tmp(out_germline_contig_tmp, in_germline_contig, o);

			// create germline contig tmp tmp
			string out_germline_contig_tmp_tmp = germline_tmp_prefix + "." + chr + ".germline.contig.tmp.tmp";
			create_germline_contig_tmp_tmp(out_germline_contig_tmp_tmp, out_germline_contig_tmp, allset);

			// refine contig tmp.tmp
			string out_germline_contig_tmp_tmp_refined = germline_tmp_prefix + "." + chr + ".germline.contig.tmp.tmp.refined";
			refine_germline_contig_tmp_tmp(out_germline_contig_tmp_tmp_refined, out_germline_contig_tmp_tmp);

			string out_germlime_contig_tmp_tmp_filtered_fa = germline_tmp_prefix + "." + chr + ".germline.contig.tmp.tmp.filtered.fa";
			create_germline_contig_tmp_tmp_filtered_fa(out_germlime_contig_tmp_tmp_filtered_fa, out_germline_contig_tmp_tmp_refined);

			// align tmp.tmp.filtered fa against repeat sequence
			string out_repeat_tmp_tmp_repeat_aln_sam = germline_tmp_prefix + "." + chr + ".germline.contig.tmp.tmp.filtered.fa.repeat.aln.sam";
			if(castle::IOUtils::get_file_size(out_repeat_tmp_tmp_repeat_aln_sam) <= 0) {
				string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.repeat_reference % out_germlime_contig_tmp_tmp_filtered_fa
						% out_repeat_tmp_tmp_repeat_aln_sam % options.repeat_reference % out_germlime_contig_tmp_tmp_filtered_fa).str();
				system(aln_cmd.c_str());
			}
			// align tmp.tmp.filtered fa against ref sequence
			string out_repeat_tmp_tmp_ref_aln_sam = germline_tmp_prefix + "." + chr + ".germline.contig.tmp.tmp.filtered.fa.ref.aln.sam";
			if(castle::IOUtils::get_file_size(out_repeat_tmp_tmp_ref_aln_sam) <= 0) {
				string aln_cmd = (boost::format("bwa aln -t 1 %s %s | bwa samse -f %s %s - %s") % options.human_reference % out_germlime_contig_tmp_tmp_filtered_fa
						% out_repeat_tmp_tmp_ref_aln_sam % options.human_reference % out_germlime_contig_tmp_tmp_filtered_fa).str();
				system(aln_cmd.c_str());
			}
			// read from .tmp.tmp.refined to list aa and bb (refined_map)
			boost::unordered_map<int64_t, string> refined_map;
			collect_two_ram_map(refined_map, out_germline_contig_tmp_tmp_refined);

			// read from .ref.aln.sam to list a and b, and set o
			set<string> rname;
			boost::unordered_map<string, string> ref_aligned_map;
			// this is set o
			set<string> ref_selected_seq_id;
			collect_refined_aln_sam_ref(rname, ref_aligned_map, ref_selected_seq_id, out_repeat_tmp_tmp_ref_aln_sam, refined_map);

			set<string> rrname;
			set<string> oo;
			set<string> ooo;
			collect_refined_aln_sam_repeat(rrname, oo, ooo, out_repeat_tmp_tmp_repeat_aln_sam, ref_selected_seq_id);
			set<string> candidate;
			set_difference(oo.begin(), oo.end(), ooo.begin(), ooo.end(), inserter(candidate, candidate.begin()));
			set<string> gold;
			{
				string line;
				vector<string> temp;
				const char* delim_underscore = "_";

				for(auto& gene : candidate) {
					castle::StringUtils::c_string_multi_split(gene, delim_underscore, temp);
					auto& two_ram_line_id = temp[0];
					if (ooo.end() == ooo.find(two_ram_line_id)) {
						gold.insert(temp[0]);
					}
				}
			}
			set<string> overlap;
			set_intersection(rname.begin(), rname.end(), rrname.begin(), rrname.end(), inserter(overlap, overlap.begin()));
			// ref_aligned_map <= dt

//			cout << (boost::format("[TEA.post_process] %s, o: %s, oo: %s, ooo: %s, cand: %s, gold: %s\n") % chr % ref_selected_seq_id.size() % oo.size() % ooo.size() % candidate.size() % gold.size()).str();
			string out_short_transduction_list = germline_tmp_prefix + "." + chr + ".germline.contig.short.transduction.list";
			create_short_transduction_list(out_short_transduction_list, candidate, o, overlap, ref_aligned_map, out_germline_contig_tmp_tmp_refined);
			vector<string> transduction_list_files;
			transduction_list_files.push_back(out_contig_transduction_filtered);
			transduction_list_files.push_back(out_short_transduction_list);
			string out_transduction_list = germline_result_prefix + "." + chr + ".germline.contig.transduction.list";
			castle::IOUtils::plain_file_merge_serial(out_transduction_list, transduction_list_files, 1, false);
			// create post contig list
			o.clear();
			collect_transduction_set(o, out_short_transduction_list);
			string out_germline_contig_list = germline_result_prefix + "." + chr + ".germline.contig.list";
			create_post_contig_list(out_germline_contig_list, o, out_germline_contig_tmp_tmp_refined);
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void TEA::create_orphan_list(const string& out_path, set<string>& orphan, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .germline.contig.orphan
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.orphan.list
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		bool on = false;
		if (boost::starts_with(fso[33], "poly") || boost::starts_with(fso[34], "poly")) {
			uint64_t i = 40;
			while (i < fso.size()) {
				if ("na" != fso[i]) {
					on = true;
					break;
				}
				++i;
			}
		} else {
			if ("na" != fso[42] && "na" != fso[43]) {
				on = true;
			}
			if (("na" != fso[42] || "na" != fso[43]) && ("na" != fso[40] || "na" != fso[41])) {
				on = true;
			}
		}

		if (on) {
			string tmp = fso[1] + "@" + fso[6] + "@" + fso[7];
			orphan.insert(tmp);
			g << line << "\n";
		}
	}
}

void TEA::create_transduction_filtered(const string& out_path, set<string>& trans, const string& in_path, const set<string>& orphan) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .germline.contig.transduction
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.transduction.filtered
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		int32_t on = 0;
		uint64_t i = 40;
		while (i < fso.size()) {
			if ("na" != fso[i]) {
				++on;
			}
			++i;
		}

		if (4 != on) {
			string tmp = fso[1] + "@" + fso[6] + "@" + fso[7];
			if (boost::starts_with(fso[33], "poly") || boost::starts_with(fso[34], "poly")) {
				if (orphan.end() == orphan.find(tmp)) {
					trans.insert(tmp);
					g << line << "\n";
				}
			}
		}
	}
}

void TEA::create_contig_filtered_fa(const string& out_path, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	int32_t min_polyAT = 5;
	int64_t line_id = 1;
	// two_ram file
	ifstream f(in_path, ios::binary);
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);

//		auto& polyA_tag = fso[33];
//		auto& polyT_tag = fso[34];
//		auto& pclipped = fso[35];
//		auto& nclipped = fso[36];
		auto& prammate = fso[37];
		auto& nrammate = fso[38];

		// polyA
		if(nrammate.size() > 0) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_nram_" << fso[8] << "\n";
					g << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g << ">" << line_id << "_nram_" << fso[8] << "\n";
					g << nrammate << "\n";
				}
			}
		}
		// polyT
		if(prammate.size() > 0) {
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_pram_" << fso[8] << "\n";
					g << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g << ">" << line_id << "_pram_" << fso[8] << "\n";
					g << prammate << "\n";
				}
			}
		}
		++line_id;
	}
}

void TEA::collect_aln_repeat_selected_seq_id(set<string>& repeat_selected_seq_id, const string& in_path) {
	string line;
	vector<string> fso;
	vector<string> temp;
	const char* delim_tab = "\t";
	const char* delim_underscore = "_";

	// read from .germline.contig.repeat.aln.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& query_seq_id = fso[0];
		auto& repeat_seq_id = fso[2];
		boost::replace_all(repeat_seq_id, "/", "_");
		if(string::npos == query_seq_id.find(repeat_seq_id)) {
			continue;
		}

		castle::StringUtils::c_string_multi_split(query_seq_id, delim_underscore, temp);
		repeat_selected_seq_id.insert(temp[0]);
	}
}

void TEA::create_germline_contig_tmp(const string& out_path, const string& in_path, const set<string>& o) {
	string line;
	vector<string> fso;
	// read from .germline.contig
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.tmp
	ofstream g(out_path, ios::binary);
	int64_t i = 1;

//	const bool is_chr = string::npos != in_path.find("chr12");
	while(getline(f, line, '\n')) {
		string str_i = boost::lexical_cast<string>(i);
		if(o.end() == o.find(str_i)) {
			++i;
			continue;
		}
		g << line << "\n";
		++i;
	}
}

void TEA::create_germline_contig_tmp_tmp(const string& out_path, const string& in_path, const set<string>& allset) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .germline.contig.tmp
	ifstream f(in_path, ios::binary);
	// write to .germline.contig.tmp.tmp
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		string tmp = fso[1] + "@" + fso[6] + "@" + fso[7];
		auto& pram = fso[12];
		auto& nram = fso[13];
		if (allset.end() == allset.find(tmp)) {
			if("0" != pram && "0" != nram) {
				g << line << "\n";
			}
		}
	}
}

void TEA::refine_germline_contig_tmp_tmp(const string& out_path, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	int32_t min_polyAT = 5;
	// read from .tmp.tmp
	ifstream f(in_path, ios::binary);
	// read from .tmp.tmp.refined
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& polyA_tag = fso[33];
		auto& polyT_tag = fso[34];
		auto& pclipped = fso[35];
		auto& nclipped = fso[36];
		auto& prammate = fso[37];
		auto& nrammate = fso[38];

		// polyA
		if("polyA" == polyA_tag && "-" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(n_pos - 1 > 0 && 'A' == pclipped[n_pos - 1]) {
				pclipped[n_pos] = 'A';
			}
			n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				pclipped = pclipped.substr(0, n_pos);
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(n_pos - 1 > 0 && 'A' == nrammate[n_pos - 1]) {
				nrammate[n_pos] = 'A';
			}
			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				nrammate = nrammate.substr(0, n_pos);
			}
		}
		// polyT
		else if("-" == polyA_tag && "polyT" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(n_pos + 1 < static_cast<int64_t>(nclipped.size()) && 'T' == nclipped[n_pos + 1]) {
				nclipped[n_pos] = 'T';
			}
			n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 0) {
					nclipped = nclipped.substr(n_pos);
				}
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(n_pos + 1 < static_cast<int64_t>(prammate.size()) && 'T' == prammate[n_pos + 1]) {
				prammate[n_pos] = 'T';
			}
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 0) {
					prammate = prammate.substr(n_pos);
				}
			}
		} else {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(n_pos - 1 > 0 && 'A' == pclipped[n_pos - 1]) {
				pclipped[n_pos] = 'A';
			}
			n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				pclipped = pclipped.substr(0, n_pos);
				polyA_tag = "polyA";
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(n_pos - 1 > 0 && 'A' == nrammate[n_pos - 1]) {
				nrammate[n_pos] = 'A';
			}
			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				nrammate = nrammate.substr(0, n_pos);
			}
			n_len = 0;
			n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(n_pos + 1 < static_cast<int64_t>(nclipped.size()) && 'T' == nclipped[n_pos + 1]) {
				nclipped[n_pos] = 'T';
			}
			n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 0) {
					nclipped = nclipped.substr(n_pos);
				}
				polyT_tag = "polyT";
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(n_pos + 1 < static_cast<int64_t>(prammate.size()) && 'T' == prammate[n_pos + 1]) {
				prammate[n_pos] = 'T';
			}
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 0) {
					prammate = prammate.substr(n_pos);
				}
			}
		}

		g << fso[0];
		int64_t col_i = 1;
		while(col_i <= 32) {
			g << "\t" << fso[col_i];
			++col_i;
		}
		g << "\t" << polyA_tag << "\t" << polyT_tag << "\t" << pclipped << "\t" << nclipped << "\t" << prammate << "\t" << nrammate << "\t" << (fso.size() > 39 ? fso[39] : "") << "\t" << (fso.size() > 40 ? fso[40] : "") << "\n";
	}
}

void TEA::create_germline_contig_tmp_tmp_filtered_fa(const string& out_path, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";

	int32_t min_polyAT = 5;
	int64_t line_id = 1;
	// read from .tmp.tmp.refined
	ifstream f(in_path, ios::binary);
	// write to .tmp.tmp.filtered.fa
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::tokenize(line, delim_tab, fso);

		if(fso.size() <= 38) {
			cout << in_path << "\n";
			cout << line << "\n";
			exit(0);
		}
		auto& polyA_tag = fso[33];
		auto& polyT_tag = fso[34];
		auto& pclipped = fso[35];
		auto& nclipped = fso[36];
		auto& prammate = fso[37];
		auto& nrammate = fso[38];

		// polyA
		if("polyA" == polyA_tag && "-" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = rfind_the_last_of(n_len, pclipped, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped.substr(0, n_pos) << "\n";
				}
			} else {
				if(pclipped.size() >= 15) {
					g << ">" << line_id << "_pclipped\n";
					g << pclipped << "\n";
				}
			}

			n_pos = rfind_the_last_of(n_len, nrammate, 'A');
			if(-1 != n_pos && n_len >= min_polyAT) {
				if(n_pos >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate.substr(0, n_pos) << "\n";
				}
			} else {
				if(nrammate.size() >= 15) {
					g << ">" << line_id << "_nram\n";
					g << nrammate << "\n";
				}
			}
		}
		// polyT
		if("-" == polyA_tag && "polyT" == polyT_tag) {
			int64_t n_len = 0;
			int64_t n_pos = find_the_last_of(n_len, nclipped, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = nclipped.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped.substr(n_pos) << "\n";
				}
			} else {
				if(nclipped.size() >= 15) {
					g << ">" << line_id << "_nclipped\n";
					g << nclipped << "\n";
				}
			}

			n_pos = find_the_last_of(n_len, prammate, 'T');
			if(-1 != n_pos && n_len >= min_polyAT) {
				int64_t seq_size = prammate.size() - n_pos;
				if(seq_size >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate.substr(n_pos) << "\n";
				}
			} else {
				if(prammate.size() >= 15) {
					g << ">" << line_id << "_pram\n";
					g << prammate << "\n";
				}
			}
		}
		++line_id;
	}
}

void TEA::collect_refined_aln_sam_ref(set<string>& rname, boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map) {
	string line;
	vector<string> fso;
	vector<string> cols;
	vector<string> out;

	const char* delim_tab = "\t";
	const char* delim_underscore = "_";
	const char* delim_ampersand = "&";

	// read from .germline.contig.ref.aln.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		if ("*" != fso[2]) {
			if ("0" == fso[1]) {
				fso[1] = "+";
			} else {
				fso[1] = "-";
			}
			auto& seq_id = fso[0];
			auto& chr_name = fso[2];
			auto& start_pos = fso[3];
			auto& cigar_str = fso[5];
			auto& cigar_flag = fso[1];
			string entry_value = chr_name + "," + start_pos + "," + cigar_str + "," + cigar_flag;
			aligned_map[seq_id] = entry_value;
		}
		int64_t mapq = boost::lexical_cast<int64_t>(fso[4]);
		if (mapq >= 25) {
			rname.insert(fso[0]);
			castle::StringUtils::c_string_multi_split(fso[0], delim_underscore, cols);
			int64_t cur_line_id = boost::lexical_cast<int64_t>(cols[0]);
			auto two_ram_itr = two_ram_map.find(cur_line_id);
			if(two_ram_map.end() == two_ram_itr) {
				continue;
			}
			auto& two_ram_entry_value = two_ram_itr->second;
			castle::StringUtils::c_string_multi_split(two_ram_entry_value, delim_ampersand, out);

			bool on = false;
			// is two_ram_id equal?
			if (fso[2] == out[0]) {
				int64_t sam_start_pos = boost::lexical_cast<int64_t>(fso[3]);
				int64_t two_ram_start_pos = boost::lexical_cast<int64_t>(out[1]);
				int64_t two_ram_end_pos = boost::lexical_cast<int64_t>(out[2]);
				if ((sam_start_pos >= (two_ram_start_pos - 5000)) && sam_start_pos <= (two_ram_end_pos + 5000)) {
					on = true;
				}
			}

			if (!on) {
				two_ram_seq_id_set.insert(fso[0]);
			}
		}
	}
}
void TEA::collect_refined_aln_sam_repeat(set<string>& rrname, set<string>& oo, set<string>& ooo, const string& in_path, const set<string>& o) {
	string line;
	vector<string> fso;
	vector<string> temp;
	const char* delim_tab = "\t";
	const char* delim_underscore = "_";

	// read from .germline.contig.tmp.tmp.filtered.fa.repeat.aln.sam
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		if('@' == line[0]) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& query_seq_id = fso[0];
		int32_t sam_flag = boost::lexical_cast<int64_t>(fso[1]);
		auto& ref_name = fso[2];
		int64_t mapq = boost::lexical_cast<int64_t>(fso[4]);

		if(o.end() != o.find(query_seq_id)) {
			if (4 == sam_flag || (4 != sam_flag && (0 != mapq && 37 != mapq)) || ("RICKSHA" == ref_name)) {
				rrname.insert(query_seq_id);
				castle::StringUtils::c_string_multi_split(query_seq_id, delim_underscore, temp);
				oo.insert(temp[0]);
			}
		}
		if(37 == mapq && "RICKSHA" != ref_name) {
			if(string::npos != query_seq_id.find("clipped")) {
				castle::StringUtils::c_string_multi_split(query_seq_id, delim_underscore, temp);
				ooo.insert(temp[0]);
			}
		}
	}
}

void TEA::create_short_transduction_list(const string& out_path, const set<string>& candidate, const set<string>& o, const set<string>& overlap, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path) {
	string line;
	// read from .tmp.tmp.refined
	ifstream f(in_path, ios::binary);
	// write to .short.transduction.list
	ofstream g(out_path, ios::binary);
	int64_t i = 1;

//	const bool is_chr = string::npos != in_path.find("chr12");
	while(getline(f, line, '\n')) {
		string str_i = boost::lexical_cast<string>(i);
		if(candidate.end() == candidate.find(str_i)) {
			++i;
			continue;
		}
		g << line;
		// p clipped
		string pc = str_i + "_pclipped";
		if(o.end() != o.find(pc) && overlap.end() != overlap.find(pc)) {
			auto dt_val = ref_aligned_map.find(pc);
			if(dt_val != ref_aligned_map.end()) {
				g << "\t" << dt_val->second;
			} else {
				g << "\tna";
			}
		} else {
			g << "\tna";
		}

		//n clipped
		string nc = str_i + "_nclipped";
		if(o.end() != o.find(nc) && overlap.end() != overlap.find(nc)) {
			auto dt_val = ref_aligned_map.find(nc);
			if(dt_val != ref_aligned_map.end()) {
				g << "\t" << dt_val->second;
			} else {
				g << "\tna";
			}
		} else {
			g << "\tna";
		}
		// pram
		string pram = str_i + "_pram";
		if(o.end() != o.find(pram) && overlap.end() != overlap.find(pram)) {
			auto dt_val = ref_aligned_map.find(pram);
			if(dt_val != ref_aligned_map.end()) {
				g << "\t" << dt_val->second;
			} else {
				g << "\tna";
			}
		} else {
			g << "\tna";
		}
		// nram
		string nram = str_i + "_nram";
		if(o.end() != o.find(nram) && overlap.end() != overlap.find(nram)) {
			auto dt_val = ref_aligned_map.find(nram);
			if(dt_val != ref_aligned_map.end()) {
				g << "\t" << dt_val->second;
			} else {
				g << "\tna";
			}
		} else {
			g << "\tna";
		}
		g << "\n";
		++i;
	}
}

void TEA::collect_transduction_set(set<string>& o, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .short.transuction.list
	ifstream f(in_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, fso);
		auto& chr_name = fso[1];
//		boost::replace_all(chr_name, "chr", "");
		auto& pbp = fso[6];
		auto& nbp = fso[7];
		string tmp = chr_name + "@" + pbp + "@" + nbp;
		o.insert(tmp);
	}
}

void TEA::create_post_contig_list(const string& out_path, const set<string>& o, const string& in_path) {
	string line;
	vector<string> fso;
	const char* delim_tab = "\t";
	// read from .tmp.tmp.refined
	ifstream f(in_path, ios::binary);
	// write to .contig.list
	ofstream g(out_path, ios::binary);
	while(getline(f, line, '\n')) {
		castle::StringUtils::tokenize(line, delim_tab, fso);
		auto& chr_name = fso[1];
	//		boost::replace_all(chr_name, "chr", "");
		auto& pbp = fso[6];
		auto& nbp = fso[7];
		string tmp = chr_name + "@" + pbp + "@" + nbp;
		auto& pram = fso[12];
		auto& nram = fso[13];
		if(o.end() != o.find(tmp) || "0" == pram || "0" == nram) {
			continue;
		}
		g << line << "\n";
	}
}

void TEA::output_raw_file(const string& chr, const string& cl_prefix, const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl, const map<int64_t, int64_t>& pm_cl, const boost::unordered_set<int64_t>& positive_only, const boost::unordered_set<int64_t>& negative_only, const int64_t read_length,
		const int64_t fragment_size, const bool headless) {

	string header_cl_raw = "chr\ts\te\tsize\trep.repeat\tfamily\tclass\tram\tram1\tram2\ts1\te1\ts2\te2\tpos1\tpos2\trep.repeat1\trep.repeat2\trepeat.name1\trepeat.name2\t"
			"rname1\trname2\n";
	string tmp_chr_name(chr);
	if (string::npos == tmp_chr_name.find("chr")) {
		tmp_chr_name = "chr" + chr;
	}
	string cl_raw_file = cl_prefix + "." + tmp_chr_name + ".cluster.raw";
	ofstream out_cl_raw(cl_raw_file, ios::binary);

	if(!headless) {
		out_cl_raw << header_cl_raw;
	}
	{
		for (auto an_entry : pm_cl) {
			auto& positive_entry = p_cl[an_entry.first];
			auto& negative_entry = n_cl[an_entry.second];
			int64_t the_ram_boundary_start = positive_entry.start;
			int64_t the_ram_boundary_end = negative_entry.stop + read_length;
			int64_t ram_boundary_size = the_ram_boundary_end - the_ram_boundary_start + 1;

			string p_rep_repeat = castle::StringUtils::join(positive_entry.value.rep_repeat, ",");
			string n_rep_repeat = castle::StringUtils::join(negative_entry.value.rep_repeat, ",");

			set<string> rep_family;
			rep_family.insert(positive_entry.value.family.begin(), positive_entry.value.family.end());
			string repeat_family = castle::StringUtils::join(rep_family, ",");
			set<string> rep_class;
			rep_class.insert(positive_entry.value.repeat_class.begin(), positive_entry.value.repeat_class.end());
			string repeat_class = castle::StringUtils::join(rep_class, ",");

			int64_t ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();

			int64_t pram = positive_entry.value.pos.size();
			int64_t nram = negative_entry.value.pos.size();
			int64_t pram_start = positive_entry.value.pos[0];
			int64_t pram_end = positive_entry.value.pos.back();
			int64_t nram_start = negative_entry.value.pos[0];
			int64_t nram_end = negative_entry.value.pos.back();

			string pram_pos = castle::StringUtils::join(positive_entry.value.pos, ",");
			vector<int64_t> negative_pos(negative_entry.value.pos);
			for (auto& a_val : negative_pos) {
				a_val = -a_val;
			}
			string nram_pos = castle::StringUtils::join(negative_pos, ",");

			string p_repeat_name = castle::StringUtils::join(positive_entry.value.repeat_name, ",");
			string n_repeat_name = castle::StringUtils::join(negative_entry.value.repeat_name, ",");

			string p_read_name = castle::StringUtils::join(positive_entry.value.rname, ",");
			string n_read_name = castle::StringUtils::join(negative_entry.value.rname, ",");

			out_cl_raw << tmp_chr_name << "\t" << the_ram_boundary_start << "\t" << the_ram_boundary_end << "\t" << ram_boundary_size << "\t" << p_rep_repeat << "\t" << repeat_family << "\t" << repeat_class << "\t" << ram << "\t" << pram << "\t" << nram << "\t" << pram_start << "\t" << pram_end
					<< "\t" << nram_start << "\t" << nram_end << "\t" << pram_pos << "\t" << nram_pos << "\t" << p_rep_repeat << "\t" << n_rep_repeat << "\t" << p_repeat_name << "\t" << n_repeat_name << "\t" << p_read_name << "\t" << n_read_name << "\n";
		}
	}
	{
//		string cl_raw_file = cl_prefix + "." + tmp_chr_name + ".cluster.p.raw";
//		ofstream out_cl_raw(cl_raw_file, ios::binary);
//		out_cl_raw << header_cl_raw;
		for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
			if (positive_only.end() == positive_only.find(r_id)) {
				continue;
			}
			auto& positive_entry = p_cl[r_id];
			int64_t the_ram_boundary_start = positive_entry.start;
			int64_t the_ram_boundary_end = positive_entry.stop + read_length + fragment_size;
			int64_t ram_boundary_size = the_ram_boundary_end - the_ram_boundary_start + 1;

			string p_rep_repeat = castle::StringUtils::join(positive_entry.value.rep_repeat, ",");

			set<string> rep_family;
			rep_family.insert(positive_entry.value.family.begin(), positive_entry.value.family.end());
			string repeat_family = castle::StringUtils::join(rep_family, ",");
			set<string> rep_class;
			rep_class.insert(positive_entry.value.repeat_class.begin(), positive_entry.value.repeat_class.end());
			string repeat_class = castle::StringUtils::join(rep_class, ",");

			int64_t ram = positive_entry.value.pos.size();

			int64_t pram = positive_entry.value.pos.size();
			int64_t pram_start = positive_entry.value.pos[0];
			int64_t pram_end = positive_entry.value.pos.back();

			string pram_pos = castle::StringUtils::join(positive_entry.value.pos, ",");

			string p_repeat_name = castle::StringUtils::join(positive_entry.value.repeat_name, ",");

			string p_read_name = castle::StringUtils::join(positive_entry.value.rname, ",");

			out_cl_raw << tmp_chr_name << "\t" << the_ram_boundary_start << "\t" << the_ram_boundary_end << "\t" << ram_boundary_size << "\t" << p_rep_repeat << "\t" << repeat_family << "\t" << repeat_class << "\t" << ram << "\t" << pram << "\t0\t" << pram_start << "\t" << pram_end << "\t0\t0\t"
					<< pram_pos << "\t\t" << p_rep_repeat << "\t\t" << p_repeat_name << "\t\t" << p_read_name << "\t\n";
		}
	}
	{
//		string cl_raw_file = cl_prefix + "." + tmp_chr_name + ".cluster.n.raw";
//		ofstream out_cl_raw(cl_raw_file, ios::binary);
//		out_cl_raw << header_cl_raw;
		for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
			if (negative_only.end() == negative_only.find(r_id)) {
				continue;
			}
			auto& negative_entry = n_cl[r_id];
			int64_t the_ram_boundary_start = negative_entry.start - fragment_size;
			int64_t the_ram_boundary_end = negative_entry.stop + read_length;
			int64_t ram_boundary_size = the_ram_boundary_end - the_ram_boundary_start + 1;

			string n_rep_repeat = castle::StringUtils::join(negative_entry.value.rep_repeat, ",");

			set<string> rep_family;
			rep_family.insert(negative_entry.value.family.begin(), negative_entry.value.family.end());
			string repeat_family = castle::StringUtils::join(rep_family, ",");
			set<string> rep_class;
			rep_class.insert(negative_entry.value.repeat_class.begin(), negative_entry.value.repeat_class.end());
			string repeat_class = castle::StringUtils::join(rep_class, ",");

			int64_t ram = negative_entry.value.pos.size();

			int64_t nram = negative_entry.value.pos.size();
			int64_t nram_start = negative_entry.value.pos[0];
			int64_t nram_end = negative_entry.value.pos.back();

			vector<int64_t> negative_pos(negative_entry.value.pos);
			for (auto& a_val : negative_pos) {
				a_val = -a_val;
			}
			string nram_pos = castle::StringUtils::join(negative_pos, ",");

			string n_repeat_name = castle::StringUtils::join(negative_entry.value.repeat_name, ",");

			string n_read_name = castle::StringUtils::join(negative_entry.value.rname, ",");

			out_cl_raw << tmp_chr_name << "\t" << the_ram_boundary_start << "\t" << the_ram_boundary_end << "\t" << ram_boundary_size << "\t" << n_rep_repeat << "\t" << repeat_family << "\t" << repeat_class << "\t" << ram << "\t0\t" << nram << "\t0\t0\t" << nram_start << "\t" << nram_end << "\t\t"
					<< nram_pos << "\t\t" << n_rep_repeat << "\t\t" << n_repeat_name << "\t\t" << n_read_name << "\n";
		}
	}
}

void TEA::BAM_to_FASTQ_serial(const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.BAM_to_FASTQ_serial");
	checker.start();

	string a_path(input_BAM_name);
	string an_index_path(a_path);
	an_index_path += ".bai";
	vector<function<void()> > tasks;
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		return;
	}

	BamTools::BamAlignment local_alignment_entry;
	boost::unordered_map<string, pair<string, string>> unprocessed_pair;
	auto& local_pair = unprocessed_pair;
	string out_1_name = disc_1_FASTQ_name;
	string out_2_name = disc_2_FASTQ_name;
	ofstream out_disc_1(out_1_name, ios::binary);
	ofstream out_disc_2(out_2_name, ios::binary);
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//		const bool debug = "FCD1JLLACXX:7:2315:9464:54561#AAAAAGATmu1" == local_alignment_entry.Name;
		const bool debug = false;
		if (local_alignment_entry.IsReverseStrand()) {
			local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
		}
		string an_entry = (boost::format("@%s\n%s\n+\n%s\n") % local_alignment_entry.Name % local_alignment_entry.QueryBases % local_alignment_entry.Qualities).str();
		auto the_pair_itr = local_pair.find(local_alignment_entry.Name);
		if (local_pair.end() == the_pair_itr) {
			if (debug) {
				cout << "here-0\n";
			}
			if (local_alignment_entry.IsSecondMate()) {
				if (debug) {
					cout << "here-1\n";
				}
				local_pair[local_alignment_entry.Name].second = an_entry;
			} else {
				if (debug) {
					cout << "here-2\n";
				}
				local_pair[local_alignment_entry.Name].first = an_entry;
			}
		} else {
			if (debug) {
				cout << "here-3\n";
			}
			if (local_alignment_entry.IsSecondMate()) {
				if (debug) {
					cout << "here-4\n";
				}
				the_pair_itr->second.second = an_entry;
			} else {
				if (debug) {
					cout << "here-5\n";
				}
				the_pair_itr->second.first = an_entry;
			}
			if (!the_pair_itr->second.first.empty() && !the_pair_itr->second.second.empty()) {
				out_disc_1 << the_pair_itr->second.first;
				out_disc_2 << the_pair_itr->second.second;
				local_pair.erase(the_pair_itr);
			}
		}
	}

	local_reader.Close();

	cout << "[TEA.BAM_to_FASTQ_serial] gather scattered information\n";
	string the_orphan_file_name = orphan_FASTQ_name;
	string the_disc_1_file_name = disc_1_FASTQ_name;
	string the_disc_2_file_name = disc_2_FASTQ_name;

	tasks.push_back([&] {
		ofstream out_disc_1(the_disc_1_file_name, ios::binary | ios::app);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_disc_1 << the_pair_entry.second.first;
			}
		}
	});
	tasks.push_back([&] {
		ofstream out_disc_2(the_disc_2_file_name, ios::binary | ios::app);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_disc_2 << the_pair_entry.second.second;
			}
		}
	});

	tasks.push_back([&] {
		ofstream out_orphan(the_orphan_file_name, ios::binary);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				continue;
			}
			if(the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_orphan << the_pair_entry.second.second;
			} else if(the_pair_entry.second.second.empty() && !the_pair_entry.second.first.empty()) {
				out_orphan << the_pair_entry.second.first;
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}
void TEA::BAM_to_FASTQ(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name) {
	int64_t size_block = 8192000;
//	vector<meerkat::BlockBoundary> fixed_size_blocks;
//	collect_boundaries(fixed_size_blocks, input_BAM_name, size_block);
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, a_path, a_bai_path, a_bni_path, size_block);
	_BAM_to_FASTQ(local_unmapped_included_blocks, a_path, orphan_FASTQ_name, disc_1_FASTQ_name, disc_2_FASTQ_name);
}

void TEA::_BAM_to_FASTQ(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.BAM_to_FASTQ");
	checker.start();

	string a_path(input_BAM_name);
	string an_index_path(a_path);
	an_index_path += ".bai";

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

// only for debugging
	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	if (verbose) {
//		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
//			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
//		}
//	}
	// map<read name, pair<first read, second read>>
	vector<boost::unordered_map<string, pair<string, string>>> pair_lists(calculated_n_blocks - 1);
	vector<string> disc_1_filenames(calculated_n_blocks - 1);
	vector<string> disc_2_filenames(calculated_n_blocks - 1);
//	vector<string> read_groups;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;

			BamAlignment local_alignment_entry;
//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;

			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;

//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.BAM_to_FASTQ] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.BAM_to_FASTQ] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			auto& local_pair = pair_lists[block_id];
			string out_1_name = disc_1_FASTQ_name + "." + str_block_id;
			string out_2_name = disc_2_FASTQ_name + "." + str_block_id;
			disc_1_filenames[block_id] = out_1_name;
			disc_2_filenames[block_id] = out_2_name;
			ofstream out_disc_1(out_1_name, ios::binary);
			ofstream out_disc_2(out_2_name, ios::binary);
			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				if(verbose && 0 == num_total) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (first-wrong) %s\n")
								% block_id % a_block_boundary_str).str();
					} else {
						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (first) %s\n")
								% block_id % a_block_boundary_str).str();
					}
				}
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
				++num_total;
				if(local_alignment_entry.IsReverseStrand()) {
					local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
				}
				string an_entry = (boost::format("@%s\n%s\n+\n%s\n") % local_alignment_entry.Name %
						local_alignment_entry.QueryBases % local_alignment_entry.Qualities).str();
				auto the_pair_itr = local_pair.find(local_alignment_entry.Name);
				if(local_pair.end() == the_pair_itr) {
					if(local_alignment_entry.IsSecondMate()) {
						local_pair[local_alignment_entry.Name].second = an_entry;
					} else {
						local_pair[local_alignment_entry.Name].first = an_entry;
					}
				} else {
					if(local_alignment_entry.IsSecondMate()) {
						the_pair_itr->second.second = an_entry;
					} else {
						the_pair_itr->second.first = an_entry;
					}
					if(!the_pair_itr->second.first.empty() && !the_pair_itr->second.second.empty()) {
						out_disc_1 << the_pair_itr->second.first;
						out_disc_2 << the_pair_itr->second.second;
						local_pair.erase(the_pair_itr);
					}
				}
			}

			local_reader.Close();
//				done_vector[block_id] = 'D';

//				if(verbose) {
//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//							% local_alignment_entry.AlignmentFlag).str();
//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
//								% block_id % a_block_boundary_str).str();
//					} else {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//				} else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
			});
	}
// sometimes the unaligned reads are in the last portion of BAM file, hence
// changing the order.
	if (tasks.size() > 0) {
		swap(tasks[0], tasks.back());
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	cout << "[TEA.BAM_to_FASTQ] gather scattered information\n";
	boost::unordered_map<string, pair<string, string>> unprocessed_pair;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		auto& local_pair = pair_lists[block_id];
		for (auto& an_entry : local_pair) {
			auto& a_read_name = an_entry.first;
			auto& the_pair = an_entry.second;
			if (unprocessed_pair[a_read_name].first.empty() && !the_pair.first.empty()) {
				unprocessed_pair[a_read_name].first = the_pair.first;
			}
			if (unprocessed_pair[a_read_name].second.empty() && !the_pair.second.empty()) {
				unprocessed_pair[a_read_name].second = the_pair.second;
			}
		}
//		local_pair.clear();
	}
	string the_orphan_file_name = orphan_FASTQ_name;
	string the_disc_1_file_name = disc_1_FASTQ_name;
	string the_disc_2_file_name = disc_2_FASTQ_name;
	if (0 < actual_blocks.size()) {
		the_disc_1_file_name = disc_1_filenames[0];
		the_disc_2_file_name = disc_2_filenames[0];
	}

	tasks.push_back([&] {
		ofstream out_disc_1(the_disc_1_file_name, ios::binary | ios::app);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_disc_1 << the_pair_entry.second.first;
			}
		}
	});
	tasks.push_back([&] {
		ofstream out_disc_2(the_disc_2_file_name, ios::binary | ios::app);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_disc_2 << the_pair_entry.second.second;
			}
		}
	});

	tasks.push_back([&] {
		ofstream out_orphan(the_orphan_file_name, ios::binary);
		for(auto the_pair_entry: unprocessed_pair) {
			if(!the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				continue;
			}
			if(the_pair_entry.second.first.empty() && !the_pair_entry.second.second.empty()) {
				out_orphan << the_pair_entry.second.second;
			} else if(the_pair_entry.second.second.empty() && !the_pair_entry.second.first.empty()) {
				out_orphan << the_pair_entry.second.first;
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge_serial(disc_1_FASTQ_name, disc_1_filenames, n_cores, true);
	castle::IOUtils::plain_file_merge_serial(disc_2_FASTQ_name, disc_2_filenames, n_cores, true);
	cout << checker;
}

void TEA::MEMBAM_to_FASTQ(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.MEMBAM_to_FASTQ");
	checker.start();
	cout << (boost::format("[TEA.MEMBAM_to_FASTQ] input: %s\n") % input_BAM_name).str();
	cout << (boost::format("[TEA.MEMBAM_to_FASTQ] output-1: %s\n") % disc_1_FASTQ_name).str();
	cout << (boost::format("[TEA.MEMBAM_to_FASTQ] output-2: %s\n") % disc_2_FASTQ_name).str();
	string a_path_index(input_BAM_name + ".bfi");
	vector<int64_t> block_boundary;
	if (boost::filesystem::exists(a_path_index)) {
		castle::ParallelRunner::load_bfi_index(block_boundary, a_path_index);
	} else {
		const int64_t N_BLOCK_ENTRIES = 262144;
		castle::ParallelRunner::create_bfi_index(block_boundary, input_BAM_name, a_path_index, N_BLOCK_ENTRIES);
	}
	_MEMBAM_to_FASTQ(block_boundary, input_BAM_name, disc_1_FASTQ_name, disc_2_FASTQ_name);
	cout << checker;
}
void TEA::_MEMBAM_to_FASTQ(vector<int64_t>& block_boundary, const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name)  {
	string a_path(input_BAM_name);
	BamReader serial_reader;
	string a_path_index(a_path + ".bfi");
	if (!serial_reader.Open(a_path)) {
		return;
	}
	vector<function<void()> > tasks;
	int64_t n_block_boundaries = block_boundary.size();
	vector<string> output_files_1(n_block_boundaries - 1);
	vector<string> output_files_2(n_block_boundaries - 1);

	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path)) {
				return;
			}
			const RefVector& m_ref = local_reader.GetReferenceData();
			boost::unordered_map<string, int64_t> ref_reverse_index;
			for (uint64_t ref_id = 0; ref_id < m_ref.size(); ++ref_id) {
				auto& a_ref = m_ref[ref_id];
				ref_reverse_index[a_ref.RefName] = ref_id;
			}
			auto& local_data = local_reader.GetBGZF();
			if(0 != block_boundary[block_id]) {
				local_data.Seek(block_boundary[block_id]);
			}
			int64_t next_block_pos = block_boundary[block_id + 1];
			int64_t cur_pos = -1;
			string str_block_id = boost::lexical_cast<string>(block_id);

			vector<BamAlignment> pair1_alns;
			vector<BamAlignment> pair2_alns;
			vector<string> pair1_seqs;
			vector<string> pair2_seqs;
			vector<string> pair1_quals;
			vector<string> pair2_quals;

			string prev_read_name;

			BamAlignment al;

			string out_path_1(disc_1_FASTQ_name + "." + str_block_id);
			string out_path_2(disc_2_FASTQ_name + "." + str_block_id);

			ofstream out_disc_1(out_path_1, ios::binary);
			ofstream out_disc_2(out_path_2, ios::binary);
			output_files_1[block_id] = out_path_1;
			output_files_2[block_id] = out_path_2;
//			const bool debug = (0 == block_id);
			const bool debug = false;
			while (local_reader.LoadNextAlignmentCore(al)) {
				auto& cur_name = al.Name;

				if(prev_read_name != cur_name) {
//					const bool debug = string::npos != cur_name.find("HSQ700642:191:D18JJACXX:1:2312:4453:42941");
		//			for(auto& an_aln : pair1_alns) {
		//				bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.back().Type;
		//				if(an_aln.IsReverseStrand() && has_H) {
		////					debug = true;
		//					break;
		//				}
		//			}
		//			for(auto& an_aln : pair2_alns) {
		//				bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.back().Type;
		//				if(an_aln.IsReverseStrand() && has_H) {
		////					debug = true;
		//					break;
		//				}
		//			}
					if(debug) {
						for(auto& an_aln : pair1_alns) {
							cout << "[TEA.MEMBAM_to_FASTQ] 1-1 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
						}
						for(auto& an_aln : pair2_alns) {
							cout << "[TEA.MEMBAM_to_FASTQ] 1-2 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
						}
						cout << "\n";
					}
					if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
						auto original_aln1 = pair1_alns;
						auto original_aln2 = pair2_alns;
						try{
							remove_entry_enclosed_with_large_H(pair1_alns);
							remove_entry_enclosed_with_large_H(pair2_alns);
							convert_H_to_S(pair1_alns);
							convert_H_to_S(pair2_alns);
							if(debug) {
								for(auto& an_aln : pair1_alns) {
									cout << "[TEA.MEMBAM_to_FASTQ] 2-1 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
								}
								for(auto& an_aln : pair2_alns) {
									cout << "[TEA.MEMBAM_to_FASTQ] 2-2 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
								}
								cout << "\n";
							}
							split_query_to_segments(pair1_seqs, pair1_quals, pair1_alns);
							split_query_to_segments(pair2_seqs, pair2_quals, pair2_alns);
						} catch(exception& ex) {
							cout << ex.what() << "\n";
							cout << cur_name << "\n";
							for(auto& an_aln : original_aln1) {
								cout << "[TEA.MEMBAM_to_FASTQ] 3-1 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
							for(auto& an_aln : original_aln2) {
								cout << "[TEA.MEMBAM_to_FASTQ] 3-2 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
							for(auto& an_aln : pair1_alns) {
								cout << "[TEA.MEMBAM_to_FASTQ] 4-1 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
							for(auto& an_aln : pair2_alns) {
								cout << "[TEA.MEMBAM_to_FASTQ] 4-2 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
						}
						// break a query bases of pair2 alns. to each segment
						if(debug) {
							for(auto& an_aln : pair1_alns) {
								cout << "[TEA.MEMBAM_to_FASTQ] 5-1 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
							for(auto& an_aln : pair2_alns) {
								cout << "[TEA.MEMBAM_to_FASTQ] 5-2 " << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
							}
							cout << "\n";

							for(auto& a_seq : pair1_seqs) {
								cout << a_seq << ", ";
							}
							cout << "\n";
							for(auto& a_seq : pair2_seqs) {
								cout << a_seq << ",";
							}
							cout << "\n";
						}
						// write to files
						int64_t entry_id = 1;
						for(uint64_t s_id1 = 0; s_id1 < pair1_seqs.size(); ++s_id1) {
							auto& a_seq = pair1_seqs[s_id1];
							auto& a_qual = pair1_quals[s_id1];
							for(uint64_t s_id2 = 0; s_id2 < pair2_seqs.size(); ++s_id2) {
								string pair1_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % prev_read_name % entry_id
										% a_seq % a_qual).str();
								auto& a_seq2 = pair2_seqs[s_id2];
								auto& a_qual2 = pair2_quals[s_id2];
								string pair2_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % prev_read_name % entry_id
									% a_seq2 % a_qual2).str();
								out_disc_1 << pair1_entry;
								out_disc_2 << pair2_entry;
								++entry_id;
							}
						}
					}
					prev_read_name = cur_name;
					pair1_alns.clear();
					pair2_alns.clear();
					pair1_seqs.clear();
					pair2_seqs.clear();
					pair1_quals.clear();
					pair2_quals.clear();
				}
				if("(null)" != cur_name) {
					if(al.IsFirstMate()) {
						pair1_alns.push_back(al);
					} else if(al.IsSecondMate()) {
						pair2_alns.push_back(al);
					}
				}
				cur_pos = local_data.Tell();
				if(cur_pos >= next_block_pos) {
					break;
				}
			}
			if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
				auto original_aln1 = pair1_alns;
				auto original_aln2 = pair2_alns;
				try{
					remove_entry_enclosed_with_large_H(pair1_alns);
					remove_entry_enclosed_with_large_H(pair2_alns);
					convert_H_to_S(pair1_alns);
					convert_H_to_S(pair2_alns);
					split_query_to_segments(pair1_seqs, pair1_quals, pair1_alns);
					split_query_to_segments(pair2_seqs, pair2_quals, pair2_alns);
				} catch(exception& ex) {
					cout << ex.what() << "\n";
					for(auto& an_aln : original_aln1) {
						cout << "1(B)" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : original_aln2) {
						cout << "2(B)" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair1_alns) {
						cout << "1(A)" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair2_alns) {
						cout << "2(A)" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
				}
				// write to files
				int64_t entry_id = 1;
				for(uint64_t s_id1 = 0; s_id1 < pair1_seqs.size(); ++s_id1) {
					auto& a_seq = pair1_seqs[s_id1];
					auto& a_qual = pair1_quals[s_id1];
					for(uint64_t s_id2 = 0; s_id2 < pair2_seqs.size(); ++s_id2) {
						string pair1_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % prev_read_name % entry_id
								% a_seq % a_qual).str();
						auto& a_seq2 = pair2_seqs[s_id2];
						auto& a_qual2 = pair2_quals[s_id2];
						string pair2_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % prev_read_name % entry_id
							% a_seq2 % a_qual2).str();
						out_disc_1 << pair1_entry;
						out_disc_2 << pair2_entry;
						++entry_id;
					}
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge(disc_1_FASTQ_name, output_files_1, n_cores, true);
	castle::IOUtils::plain_file_merge(disc_2_FASTQ_name, output_files_2, n_cores, true);
}
void TEA::_MEMBAM_to_FASTQ_serial(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name) {
	BamTools::BamReader local_reader;
	if (!local_reader.Open(input_BAM_name)) {
		return;
	}
	auto& m_ref = local_reader.GetReferenceData();
	BamAlignment local_alignment_entry;
	ofstream out_disc_1(disc_1_FASTQ_name, ios::binary);
	ofstream out_disc_2(disc_2_FASTQ_name, ios::binary);
	vector<BamAlignment> pair1_alns;
	vector<BamAlignment> pair2_alns;
	vector<string> pair1_seqs;
	vector<string> pair2_seqs;
	vector<string> pair1_quals;
	vector<string> pair2_quals;

	string prev_read_name;
	int64_t n_counts = 0;
//	const bool debug = true;
	const bool debug = false;
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		if(debug) {
			cout << "[TEA.MEMBAM_to_FASTQ_serial] " << BamTools::BamWriter::GetSAMAlignment(local_alignment_entry, m_ref) << "\n";
		}
		auto& cur_name = local_alignment_entry.Name;
		++n_counts;
//		cout << n_counts << "\n";
		if (0 == (n_counts & 1048575)) {
			cout << (boost::format("[TEA.MEMBAM_to_FASTQ] %d processed\n") % n_counts).str();
		}
		if(prev_read_name != cur_name) {

//			for(auto& an_aln : pair1_alns) {
//				bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.back().Type;
//				if(an_aln.IsReverseStrand() && has_H) {
////					debug = true;
//					break;
//				}
//			}
//			for(auto& an_aln : pair2_alns) {
//				bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.back().Type;
//				if(an_aln.IsReverseStrand() && has_H) {
////					debug = true;
//					break;
//				}
//			}

			if(pair1_alns.size() > 0 && pair1_alns.size() > 0) {
				if(debug) {
					for(auto& an_aln : pair1_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 1-1" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair2_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 1-2" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					cout << "\n";
				}
				auto original_aln1 = pair1_alns;
				auto original_aln2 = pair2_alns;
				try{
					remove_entry_enclosed_with_large_H(pair1_alns);
					remove_entry_enclosed_with_large_H(pair2_alns);
					convert_H_to_S(pair1_alns);
					convert_H_to_S(pair2_alns);
					if(debug) {
						for(auto& an_aln : pair1_alns) {
							cout << "[TEA.MEMBAM_to_FASTQ_serial] 2-1" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
						}
						for(auto& an_aln : pair2_alns) {
							cout << "[TEA.MEMBAM_to_FASTQ_serial] 2-2" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
						}
						cout << "\n";
					}
					split_query_to_segments(pair1_seqs, pair1_quals, pair1_alns);
					split_query_to_segments(pair2_seqs, pair2_quals, pair2_alns);
				} catch(exception& ex) {
					cout << ex.what() << "\n";
					cout << cur_name << "\n";
					for(auto& an_aln : original_aln1) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 3-1" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : original_aln2) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 3-2" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair1_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 4-1" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair2_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 4-2" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
				}
				// break a query bases of pair2 alns. to each segment
				if(debug) {
					for(auto& an_aln : pair1_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 5-1" << (an_aln.IsReverseStrand() ? "R:" : "F:") << BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					for(auto& an_aln : pair2_alns) {
						cout << "[TEA.MEMBAM_to_FASTQ_serial] 5-2" << (an_aln.IsReverseStrand() ? "R:" : "F:") <<BamTools::BamWriter::GetSAMAlignment(an_aln, m_ref) << "\n";
					}
					cout << "\n";

					for(auto& a_seq : pair1_seqs) {
						cout << a_seq << ", ";
					}
					cout << "\n";
					for(auto& a_seq : pair2_seqs) {
						cout << a_seq << ",";
					}
					cout << "\n";
				}
				// write to files
				int64_t entry_id = 1;
				for(uint64_t s_id1 = 0; s_id1 < pair1_seqs.size(); ++s_id1) {
					auto& a_seq = pair1_seqs[s_id1];
					auto& a_qual = pair1_quals[s_id1];
					for(uint64_t s_id2 = 0; s_id2 < pair2_seqs.size(); ++s_id2) {
						string pair1_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % local_alignment_entry.Name % entry_id
								% a_seq % a_qual).str();
						auto& a_seq2 = pair2_seqs[s_id2];
						auto& a_qual2 = pair2_quals[s_id2];
						string pair2_entry = (boost::format("@%s:m%d\n%s\n+\n%s\n") % local_alignment_entry.Name % entry_id
							% a_seq2 % a_qual2).str();
						if(debug) {
							cout << "[TEA.MEMBAM_to_FASTQ_serial] 6-1: " << pair1_entry;
							cout << "[TEA.MEMBAM_to_FASTQ_serial] 6-2: " << pair2_entry;
						}
						out_disc_1 << pair1_entry;
						out_disc_2 << pair2_entry;
						++entry_id;
					}
				}
			}
			prev_read_name = cur_name;
			pair1_alns.clear();
			pair2_alns.clear();
			pair1_seqs.clear();
			pair2_seqs.clear();
		}
		if("(null)" != cur_name) {
			if(local_alignment_entry.IsFirstMate()) {
				pair1_alns.push_back(local_alignment_entry);
			} else if(local_alignment_entry.IsSecondMate()) {
				pair2_alns.push_back(local_alignment_entry);
			}
		}
	}
}

void TEA::remove_entry_enclosed_with_large_H(vector<BamAlignment>& alns) {
	uint64_t n_H1 = 0;
	uint64_t n_any_H = 0;
	// check if all the entries have 'H' CIGAR string
	for (auto& an_aln : alns) {
		bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.front().Type && 'H' == an_aln.CigarData.back().Type;
		bool has_any_H = an_aln.CigarData.size() > 0 && ('H' == an_aln.CigarData.front().Type || 'H' == an_aln.CigarData.back().Type);
		if(has_H) {
			++n_H1;
		}
		if(has_any_H) {
			++n_any_H;
		}
	}
	//	 all entries have 'H' CIGAR
	if(alns.size() == n_any_H) {
		int64_t max_aln_len = 0;
		for (uint64_t aln_id = 0; aln_id < alns.size(); ++aln_id) {
			auto& an_aln = alns[aln_id];
			int64_t n_matches = 0;
			for(auto& a_cigar : an_aln.CigarData) {
				if('M' == a_cigar.Type || 'I' == a_cigar.Type) {
					n_matches += a_cigar.Length;
				}
			}
			if(n_matches > max_aln_len) {
				max_aln_len = n_matches;
			}
		}
		for (uint64_t aln_id = 0; aln_id < alns.size();) {
			auto& an_aln = alns[aln_id];
			int64_t n_matches = 0;
			for(auto& a_cigar : an_aln.CigarData) {
				if('M' == a_cigar.Type || 'I' == a_cigar.Type) {
					n_matches += a_cigar.Length;
				}
			}
			if(n_matches != max_aln_len) {
				alns.erase(alns.begin() + aln_id);
				continue;
			}
			++aln_id;
		}
	} else if(n_H1 > 1) {
		int64_t max_aln_len = 0;
		for (uint64_t aln_id = 0; aln_id < alns.size(); ++aln_id) {
			auto& an_aln = alns[aln_id];
			bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.front().Type && 'H' == an_aln.CigarData.back().Type;
			if(!has_H) {
				continue;
			}
			int64_t n_matches = 0;
			for(auto& a_cigar : an_aln.CigarData) {
				if('M' == a_cigar.Type || 'I' == a_cigar.Type) {
					n_matches += a_cigar.Length;
				}
			}
			if(n_matches > max_aln_len) {
				max_aln_len = n_matches;
			}
		}
		for (uint64_t aln_id = 0; aln_id < alns.size();) {
			auto& an_aln = alns[aln_id];
			bool has_H = an_aln.CigarData.size() > 0 && 'H' == an_aln.CigarData.front().Type && 'H' == an_aln.CigarData.back().Type;
			if(!has_H) {
				++aln_id;
				continue;
			}
			int64_t n_matches = 0;
			for(auto& a_cigar : an_aln.CigarData) {
				if('M' == a_cigar.Type || 'I' == a_cigar.Type) {
					n_matches += a_cigar.Length;
				}
			}
			if(n_matches != max_aln_len) {
				alns.erase(alns.begin() + aln_id);
				continue;
			}
			++aln_id;
		}
	}
}

void TEA::convert_H_to_S(vector<BamAlignment>& alns) {
	if(1 == alns.size()) {
		return;
	}
	uint64_t max_len = 0;
	string max_seq;
	string max_qual;
	for (auto& an_aln : alns) {
		if(an_aln.QueryBases.size() > max_len) {
			max_len = an_aln.QueryBases.size();
			max_seq = an_aln.QueryBases;
			max_qual = an_aln.Qualities;
			if(an_aln.IsReverseStrand()) {
				max_seq = castle::StringUtils::get_reverse_complement(max_seq);
				max_qual = castle::StringUtils::get_reverse_string(max_qual);
			}
		}
	}
	for(auto& an_aln : alns) {
		bool has_H = false;
		for(auto& a_cigar : an_aln.CigarData) {
			if('H' == a_cigar.Type) {
				has_H = true;
				break;
			}
		}
		if(has_H) {
			if(an_aln.QueryBases.size() != max_len) {
				for(auto& a_cigar : an_aln.CigarData) {
					if('H' == a_cigar.Type) {
						a_cigar.Type = 'S';
					}
				}
				an_aln.QueryBases = max_seq;
				an_aln.Qualities = max_qual;
				if(an_aln.IsReverseStrand()) {
//					if(an_aln.CigarData.size() > 0 && 'S' == an_aln.CigarData.back().Type) {
//						an_aln.Position -= an_aln.CigarData.back().Length;
//					}
					an_aln.QueryBases = castle::StringUtils::get_reverse_complement(an_aln.QueryBases);
					an_aln.Qualities = castle::StringUtils::get_reverse_string(an_aln.Qualities);
				}
//				else {
//					if(an_aln.CigarData.size() > 0 && 'S' == an_aln.CigarData.front().Type) {
//						an_aln.Position -= an_aln.CigarData.front().Length;
//					}
//				}
			}
		}
	}
}

void TEA::fill_H_to_S(BamAlignment& aln, const AlnSeqQualEntry& aln_seq_entry) {
	int64_t expected_len = 0;
	for(auto& cigar : aln.CigarData) {
		switch (cigar.Type) {
			case 'H':
			case 'S':
			case 'I':
			case 'M':
				expected_len += cigar.Length;
				break;
			case 'D':
				expected_len -= cigar.Length;
				break;
			default:
				break;
		}
	}
	if(aln.IsFirstMate()) {
		if(expected_len == static_cast<int64_t>(aln_seq_entry.seq1.size())) {
			aln.QueryBases = aln_seq_entry.seq1;
			aln.Qualities = aln_seq_entry.qual1;
			for(auto& cigar : aln.CigarData) {
				switch (cigar.Type) {
					case 'H':
						cigar.Type = 'S';
						break;
					default:
						break;
				}
			}
		}
	} else if(aln.IsSecondMate()) {
		if(expected_len == static_cast<int64_t>(aln_seq_entry.seq2.size())) {
			aln.QueryBases = aln_seq_entry.seq2;
			aln.Qualities = aln_seq_entry.qual2;
			for(auto& cigar : aln.CigarData) {
				switch (cigar.Type) {
					case 'H':
						cigar.Type = 'S';
						break;
					default:
						break;
				}
			}
		}
	}
}

void TEA::split_query_to_segments(vector<string>& seqs, vector<string>& quals, vector<BamAlignment>& alns) {
	const bool debug = false;
	for(auto& an_aln : alns) {
		auto& a_seq = an_aln.QueryBases;
		auto& a_qual = an_aln.Qualities;
		if(!an_aln.IsMapped()) {
			seqs.push_back(a_seq);
			quals.push_back(a_qual);
			continue;
		}
		auto& cigars = an_aln.CigarData;
		if(debug) {
			for(auto& a_cigar : cigars) {
				cout << a_cigar.Length << a_cigar.Type;
			}
			cout << "\n";
		}

		for (int64_t c_id = 0; c_id < (static_cast<int64_t>(cigars.size()) - 1);) {
			auto& cur_cigar = cigars[c_id];
			auto& next_cigar = cigars[c_id + 1];
			if('S' == cur_cigar.Type) {
				++c_id;
				continue;
			}
			if ('D' == cur_cigar.Type || 'H' == cur_cigar.Type) {
				cigars.erase(cigars.begin() + c_id);
				continue;
			}
			if (('M' == cur_cigar.Type && 'M' == next_cigar.Type) || ('M' == cur_cigar.Type && 'I' == next_cigar.Type)) {
				cur_cigar.Length += next_cigar.Length;
				cigars.erase(cigars.begin() + (c_id + 1));
				continue;
			}
			++c_id;
		}
		if(debug) {
			for(auto& a_cigar : cigars) {
				cout << a_cigar.Length << a_cigar.Type;
			}
			cout << "\n";
		}
		if('H' == cigars.back().Type) {
			cigars.erase(cigars.begin() + (cigars.size() - 1));
		}
		uint64_t prev_cut_pos = 0;
		for(auto& a_cigar : an_aln.CigarData) {
			string a_segment_seq = a_seq.substr(prev_cut_pos, a_cigar.Length);
			string a_segment_qual = a_qual.substr(prev_cut_pos, a_cigar.Length);
			if(debug) {
				cout << a_cigar.Length << a_cigar.Type << "\n";
				cout << a_segment_seq << "\n";
			}
			if(static_cast<int64_t>(a_segment_seq.size()) >= options.min_clipped_len) {
				if(an_aln.IsReverseStrand()) {
					a_segment_seq = castle::StringUtils::get_reverse_complement(a_segment_seq);
					a_segment_qual = castle::StringUtils::get_reverse_string(a_segment_qual);
				}
				seqs.push_back(a_segment_seq);
				quals.push_back(a_segment_qual);
			}
			prev_cut_pos += a_cigar.Length;
		}
	}
}

void TEA::collect_boundaries_un(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block) {
//	string a_path(input_BAM_name);
//	string an_index_path;
//	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, a_bai_path)) {
		std::cout << "[TEA.collect_boundaries_un] ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.collect_boundaries_un");
	checker.start();
	BamTools::BamAlignment al;
	fixed_size_blocks.clear();
	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	cout << (boost::format("[TEA.collect_boundaries_un] # refs: %d\n") % n_refs).str();
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
		} else {
			if (a_ref_data.RefLength > size_block) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
				++boundary_id;
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			++boundary_id;
		}
	}

	int64_t calculated_n_blocks = boundary_positions.size();

	cout << "[TEA.collect_boundaries_un] total ref. size: " << size_total_ref << "\n";

	fixed_size_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	cout << "[TEA.collect_boundaries_un] collect boundary positions\n";
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
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
						cout << (boost::format("[TEA.collect_boundaries_un] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						cout << (boost::format("[TEA.collect_boundaries_un] unaligned: %s %d:%d\n")
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
		collect_boundaries_alt(unmapped_offsets, a_path, a_bai_path, a_bni_path);
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
	if (fixed_size_blocks.size() > 0 && 0 == fixed_size_blocks[0].offset) {
		fixed_size_blocks.erase(fixed_size_blocks.begin());
	}
	fixed_size_blocks.erase(unique(fixed_size_blocks.begin(), fixed_size_blocks.end()), fixed_size_blocks.end());
	unmapped_included_blocks = fixed_size_blocks;
	meerkat::BlockBoundary a_block_boundary;
	a_block_boundary.read_name = "last";
	a_block_boundary.ref_id = n_refs;
	a_block_boundary.offset = numeric_limits<int64_t>::max();
	a_block_boundary.pos = 0;
	a_block_boundary.jump_pos = -1;
	fixed_size_blocks.push_back(a_block_boundary);
	for (auto& an_offset : unmapped_offsets) {
		if (-1 == an_offset.ref_id) {
			meerkat::BlockBoundary a_boundary;
			a_boundary.offset = an_offset.offset;
			a_boundary.ref_id = an_offset.ref_id;
			a_boundary.pos = an_offset.position;
			unmapped_included_blocks.push_back(a_boundary);
		}
	}

	cout << (boost::format("[TEA.collect_boundaries_un] actual # blocks: %d\n") % fixed_size_blocks.size()).str();
	independent_blocks = fixed_size_blocks;
	cout << checker;
}

void TEA::collect_boundaries_pos(vector<meerkat::BlockBoundary>& fixed_size_blocks, vector<meerkat::BlockBoundary>& unmapped_included_blocks, vector<meerkat::BlockBoundary>& independent_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block) {
	BamTools::BamReader reader;
	if (!reader.Open(a_path, a_bai_path)) {
		std::cout << "[TEA.collect_boundaries_pos] ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.collect_boundaries_pos");
	checker.start();
	BamTools::BamAlignment al;
	fixed_size_blocks.clear();
	const BamTools::RefVector& a_ref_vector = reader.GetReferenceData();
	int64_t n_refs = a_ref_vector.size();
	cout << (boost::format("[TEA.collect_boundaries_pos] # refs: %d\n") % n_refs).str();
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
		} else {
			if (a_ref_data.RefLength > size_block) {
				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
				++boundary_id;
			}

			++ref_id;
			if (ref_id >= n_refs) {
				break;
			}
			boundary_positions.push_back(make_pair(ref_id, 0));
			++boundary_id;
		}
	}

	int64_t calculated_n_blocks = boundary_positions.size();

	cout << "[TEA.collect_boundaries_pos] total ref. size: " << size_total_ref << "\n";

	fixed_size_blocks.resize(calculated_n_blocks);
	vector<function<void()> > tasks;
	cout << "[TEA.collect_boundaries_pos] collect boundary positions\n";
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
				return;
			}
			auto& m_bgzf = local_reader.GetBGZF();
			BamTools::BamAlignment local_alignment_entry;
			bool success = local_reader.Jump(boundary_positions[block_id].first, boundary_positions[block_id].second);
			if(success) {
				string the_pos;
				string old_pos;
				do {
					int64_t cur_offset = m_bgzf.Tell();
					success = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
					if(success) {
						the_pos = (boost::format("%s:%s") % local_alignment_entry.RefID % local_alignment_entry.Position).str();
						fixed_size_blocks[block_id].read_name = local_alignment_entry.Name;
						fixed_size_blocks[block_id].ref_id = local_alignment_entry.RefID;
						fixed_size_blocks[block_id].offset = cur_offset;
						fixed_size_blocks[block_id].pos = local_alignment_entry.Position;
						fixed_size_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
						fixed_size_blocks[block_id].jump_pos = boundary_positions[block_id].second;
						if(!old_pos.empty() && old_pos != the_pos) {
							break;
						}
					}
					else {
						cout << (boost::format("[TEA.collect_boundaries_pos] Failed: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
						break;
					}
					old_pos = the_pos;
					if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
						cout << (boost::format("[TEA.collect_boundaries_pos] unaligned: %s %d:%d\n")
								% local_alignment_entry.Name
								% local_alignment_entry.RefID
								% local_alignment_entry.Position).str();
					}
				}while(true);
			}
			local_reader.Close();
		});
	}
	vector<BlockOffset> unmapped_offsets;
	tasks.push_back([&] {
		collect_boundaries_alt(unmapped_offsets, a_path, a_bai_path, a_bni_path);
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
	if (fixed_size_blocks.size() > 0 && 0 == fixed_size_blocks[0].offset) {
		fixed_size_blocks.erase(fixed_size_blocks.begin());
	}
	fixed_size_blocks.erase(unique(fixed_size_blocks.begin(), fixed_size_blocks.end()), fixed_size_blocks.end());
	unmapped_included_blocks = fixed_size_blocks;
	meerkat::BlockBoundary a_block_boundary;
	a_block_boundary.read_name = "last";
	a_block_boundary.ref_id = n_refs;
	a_block_boundary.offset = numeric_limits<int64_t>::max();
	a_block_boundary.pos = 0;
	a_block_boundary.jump_pos = -1;
	fixed_size_blocks.push_back(a_block_boundary);
	for (auto& an_offset : unmapped_offsets) {
		if (-1 == an_offset.ref_id) {
			meerkat::BlockBoundary a_boundary;
			a_boundary.offset = an_offset.offset;
			a_boundary.ref_id = an_offset.ref_id;
			a_boundary.pos = an_offset.position;
			unmapped_included_blocks.push_back(a_boundary);
		}
	}

	cout << (boost::format("[TEA.collect_boundaries_pos] actual # blocks: %d\n") % fixed_size_blocks.size()).str();
	independent_blocks = fixed_size_blocks;
	cout << checker;
}

void TEA::collect_boundaries_alt(vector<BlockOffset>& offset_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path) {
//	string a_path(input_BAM_name);
//	string bni_index_path;
//	get_bni_index_path(a_path, bni_index_path);
	if (!boost::filesystem::exists(a_bni_path)) {
		create_bni_even_index(a_path, a_bai_path, a_bni_path);
	}
	offset_blocks.clear();
	cout << (boost::format("[TEA.collect_boundaries_alt] bni: %s\n") % a_bni_path).str();
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in(a_bni_path, ios::binary);
	while (getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		BamTools::BlockOffset an_offset;
		an_offset.offset = boost::lexical_cast<int64_t>(data[0]);
		an_offset.ref_id = boost::lexical_cast<int32_t>(data[1]);
		an_offset.position = boost::lexical_cast<int64_t>(data[2]);
		offset_blocks.push_back(an_offset);
	}
	cout << (boost::format("[TEA.collect_boundaries_alt] # blocks: %d\n") % offset_blocks.size()).str();
}

void TEA::create_bni_even_index(const string& a_path, const string& a_bai_path, const string& a_bni_path) {
//	get_bni_index_path(a_path, bni_index_path);
	if (boost::filesystem::exists(a_bni_path)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.create_bni_even_index");
	checker.start();


	cout << (boost::format("[TEA.create_bni_even_index] path: %s\n") % a_path).str();
	cout << (boost::format("[TEA.create_bni_even_index] bai: %s\n") % a_bai_path).str();
	cout << (boost::format("[TEA.create_bni_even_index] bni: %s\n") % a_bni_path).str();
	BamTools::BamReader serial_reader;
	if (!serial_reader.Open(a_path, a_bai_path)) {
		cout << checker;
		return;
	}

	vector<BamTools::BlockOffset> block_boundary;
	serial_reader.GetBlockOffsets(block_boundary);
	ofstream out(a_bni_path, ios::binary);
	for (uint64_t block_id = 0; block_id < block_boundary.size(); ++block_id) {
		out << block_boundary[block_id].offset << "\t" << block_boundary[block_id].ref_id << "\t" << block_boundary[block_id].position << "\n";
	}
	cout << (boost::format("[TEA.create_bni_even_index] %d entries were indexed\n") % block_boundary.size()).str();
	cout << checker;
}

//void TEA::collect_boundaries(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& input_BAM_name, int64_t size_block, bool verbose) {
//	string a_path(input_BAM_name);
//	string an_index_path(a_path);
//	an_index_path += ".bai";
//	BamTools::BamReader reader;
//	if (!reader.Open(a_path, an_index_path)) {
//		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
//		exit(1);
//	}
//	castle::TimeChecker checker;
//	checker.setTarget("TEA.collect_boundaries");
//	if (verbose) {
//		checker.start();
//	} else {
//		checker.start_without_output();
//	}
//	BamAlignment al;
////	bool has_index = reader.HasIndex();
////	if (!has_index) {
////		if (verbose) {
////			cout << "[TEA.collect_boundaries] no index is found so it is generating.";
////		}
////		reader.CreateIndex();
////	}
//	fixed_size_blocks.clear();
//	const RefVector& a_ref_vector = reader.GetReferenceData();
//	int64_t n_refs = a_ref_vector.size();
//	if (verbose) {
//		cout << (boost::format("[TEA.collect_boundaries] # refs: %d\n") % n_refs).str();
//	}
//	int64_t size_total_ref = 0;
//	for (int32_t ref_id = 0; ref_id < n_refs; ++ref_id) {
//		auto& a_ref_data = a_ref_vector[ref_id];
//		//cout << a_ref_data.RefName << ": " << a_ref_data.RefLength << "\n";
//		size_total_ref += a_ref_data.RefLength;
//	}
//
//	int64_t estimated_n_blocks = size_total_ref / (double) size_block;
//	++estimated_n_blocks;
//
//	if (estimated_n_blocks < n_cores) {
//		size_block = (size_total_ref / (double) n_cores);
//		estimated_n_blocks = size_total_ref / (double) size_block;
//		++estimated_n_blocks;
//	}
//
//	// calculate the positions of boundary entries.
//	vector<pair<int32_t, int32_t>> boundary_positions;
//	boundary_positions.push_back(make_pair(0, 0));
//	int64_t boundary_id = 0;
//	int32_t ref_id = 0;
//	int64_t n_remaining_bases = a_ref_vector[ref_id].RefLength;
//	while (n_remaining_bases >= 0) {
//		auto& a_ref_data = a_ref_vector[ref_id];
//		int64_t last_base_pos = boundary_positions[boundary_id].second;
//		n_remaining_bases = a_ref_data.RefLength - last_base_pos;
//		if (n_remaining_bases >= size_block) {
//			boundary_positions.push_back(make_pair(ref_id, last_base_pos + size_block));
//			++boundary_id;
//		} else {
//			if (a_ref_data.RefLength > size_block) {
//				boundary_positions.push_back(make_pair(ref_id, min(static_cast<int32_t>(last_base_pos + size_block), a_ref_data.RefLength)));
//				++boundary_id;
//			}
//
//			++ref_id;
//			if (ref_id >= n_refs) {
//				break;
//			}
//			boundary_positions.push_back(make_pair(ref_id, 0));
//			++boundary_id;
//		}
//	}
//
//	//for(uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
//	//auto& a_ref_data = a_ref_vector[ref_id];
//	//boundary_positions.push_back(make_pair(ref_id, a_ref_data.RefLength));
//	//}
//
//	int64_t calculated_n_blocks = boundary_positions.size();
//
//	if (verbose) {
//		cout << "[TEA.collect_boundaries] estimated # blocks: " << estimated_n_blocks << "\n";
//		cout << "[TEA.collect_boundaries] calculated # blocks: " << calculated_n_blocks << "\n";
//		cout << "[TEA.collect_boundaries] total ref. size: " << size_total_ref << "\n";
//	}
//
//	fixed_size_blocks.resize(calculated_n_blocks);
//	vector<function<void()> > tasks;
//	if (verbose) {
//		cout << "[TEA.collect_boundaries] collect boundary positions\n";
//	}
//	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
//		tasks.push_back([&, block_id] {
//			BamTools::BamReader local_reader;
//			if (!local_reader.Open(a_path, an_index_path)) {
//				return;
//			}
//
//			BamAlignment local_alignment_entry;
//			auto& m_bgzf = local_reader.GetBGZF();
//			bool jump_success = local_reader.Jump(boundary_positions[block_id].first, boundary_positions[block_id].second);
//
//			if(jump_success) {
//				bool align_succes = false;
//				int64_t cur_offset = m_bgzf.Tell();
//				do {
//					cur_offset = m_bgzf.Tell();
//					align_succes = local_reader.LoadNextAlignmentWithName(local_alignment_entry);
//					if(align_succes) {
////						const bool debug = boundary_positions[block_id].first == 2 && "DHFC08P1:370:C1784ACXX:1:1101:6090:32147" == local_alignment_entry.Name;
//
//				fixed_size_blocks[block_id].read_name = local_alignment_entry.Name;
//				fixed_size_blocks[block_id].ref_id = local_alignment_entry.RefID;
//				fixed_size_blocks[block_id].offset = cur_offset;
//				fixed_size_blocks[block_id].pos = local_alignment_entry.Position;
//				fixed_size_blocks[block_id].aln_flag = local_alignment_entry.AlignmentFlag;
//				fixed_size_blocks[block_id].jump_pos = boundary_positions[block_id].second;
////						if(debug) {
////								cout << (boost::format("[TEA.collect_boundaries] block-%d: %s->%s %s %s %s\n") % block_id %
////										boundary_positions[block_id].first % local_alignment_entry.RefID % local_alignment_entry.Name
////										% local_alignment_entry.Position
////									% boundary_positions[block_id].second).str();
////						}
//				break;
//			}
//			else {
//				if(verbose) {
//					cout << (boost::format("[TEA.collect_boundaries] Failed: %s %d:%d\n")
//							% local_alignment_entry.Name
//							% local_alignment_entry.RefID
//							% local_alignment_entry.Position).str();
//				}
//				break;
//			}
//			if(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position) {
//				if(verbose) {
//					cout << (boost::format("[TEA.collect_boundaries] unaligned: %s %d:%d\n")
//							% local_alignment_entry.Name
//							% local_alignment_entry.RefID
//							% local_alignment_entry.Position).str();
//				}
//			}
//		}while(-1 == local_alignment_entry.RefID || -1 == local_alignment_entry.Position);
//	}
//	local_reader.Close();
//});
//	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
////	for(int64_t block_id = 0; block_id < static_cast<int64_t>(fixed_size_blocks.size()); ++block_id) {
////		auto& a_boundary = fixed_size_blocks[block_id];
////		if(2 != a_boundary.ref_id) {
////			continue;
////		}
////		cout << (boost::format("[TEA.collect_boundaries] block-%d(before): Name(%s), Ref(%d), Pos(%d), Jump(%d)\n") % block_id
////				% a_boundary.read_name % a_boundary.ref_id
////				% a_boundary.pos % a_boundary.jump_pos).str();
////	}
//
//	fixed_size_blocks.erase(unique(fixed_size_blocks.begin(), fixed_size_blocks.end()), fixed_size_blocks.end());
////	{
////		string removed(fixed_size_blocks.size(), 'F');
////		for (uint64_t block_id = 0; block_id < fixed_size_blocks.size(); ++block_id) {
////			tasks.push_back([&, block_id] {
////				BamTools::BamReader local_reader;
////				if (!local_reader.Open(a_path, an_index_path)) {
////					return;
////				}
////				auto& o = fixed_size_blocks[block_id];
////				if(!local_reader.Jump(o.ref_id, o.jump_pos)) {
////					removed[block_id] = 'T';
////				}
////				local_reader.Close();
////			});
////		}
////		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
////		vector<meerkat::BlockBoundary> tmp_fixed_block;
////		tmp_fixed_block.reserve(fixed_size_blocks.size());
////		for (uint64_t block_id = 0; block_id < fixed_size_blocks.size(); ++block_id) {
////			if ('F' == removed[block_id]) {
////				tmp_fixed_block.push_back(fixed_size_blocks[block_id]);
////			}
////		}
////		fixed_size_blocks.swap(tmp_fixed_block);
////	}
////	sort(fixed_size_blocks.begin(), fixed_size_blocks.end());
//	meerkat::BlockBoundary a_block_boundary;
//	a_block_boundary.read_name = "last";
//	a_block_boundary.ref_id = n_refs;
//	a_block_boundary.pos = 0;
//	a_block_boundary.jump_pos = -1;
//	fixed_size_blocks.push_back(a_block_boundary);
//
////	cout << (boost::format("[TEA.collect_boundaries] actual # blocks: %d\n") % fixed_size_blocks.size()).str();
////	for(int64_t block_id = 0; block_id < static_cast<int64_t>(fixed_size_blocks.size()); ++block_id) {
////		auto& a_boundary = fixed_size_blocks[block_id];
////		if(2 != a_boundary.ref_id) {
////			continue;
////		}
////		cout << (boost::format("[TEA.collect_boundaries] block-%d: Name(%s), Ref(%d), Pos(%d), Jump(%d)\n") % block_id
////				% a_boundary.read_name % a_boundary.ref_id
////				% a_boundary.pos % a_boundary.jump_pos).str();
////	}
//	if (verbose) {
//		cout << checker;
//	}
//}

void TEA::format_isize() {
	string isize_file = options.prefix + ".isize";
	string rl_file = options.prefix + ".rl";
	if (!options.working_dir.empty()) {
		isize_file = options.working_prefix + ".isize";
		rl_file = options.working_prefix + ".rl";
	}
	if (!options.is_force && boost::filesystem::exists(isize_file) && boost::filesystem::exists(rl_file)) {
		return;
	}
	vector<string> cols;
	const char* delim_white_space = ":\t";
	map<string, int32_t> rl;
	map<string, double> is_mean;
	map<string, double> is_median;
	map<string, double> is_sd;
	string line;
	ofstream O1(isize_file, ios::binary);
	ofstream O2(rl_file, ios::binary);
	ifstream in(options.isinfoname, ios::binary);
	while (getline(in, line, '\n')) {
		if (string::npos != line.find("Read length")) {
			castle::StringUtils::c_string_multi_split(line, delim_white_space, cols);
			string& rg = cols.back();
			getline(in, line, '\n');
			int32_t val = boost::lexical_cast<int32_t>(line);
			rl[rg] = val;
			//#print "wrote $val into rl hash with key $rg\n";
		} else if (string::npos != line.find("Mean")) {
			castle::StringUtils::c_string_multi_split(line, delim_white_space, cols);
			string& rg = cols.back();
			getline(in, line, '\n');
			double val = boost::lexical_cast<double>(line);
			is_mean[rg] = val;
//			#print "wrote $val into is_mean hash with key $rg\n";
		} else if (string::npos != line.find("Median")) {
			castle::StringUtils::c_string_multi_split(line, delim_white_space, cols);
			string& rg = cols.back();
			getline(in, line, '\n');
			double val = boost::lexical_cast<double>(line);
			is_median[rg] = val;
//			#print "wrote $val into is_median hash with key $rg\n";
		} else if (string::npos != line.find("Standard deviation")) {
			castle::StringUtils::c_string_multi_split(line, delim_white_space, cols);
			string& rg = cols.back();
			getline(in, line, '\n');
			double val = boost::lexical_cast<double>(line);
			is_sd[rg] = val;
//			#print "wrote $val into is_sd hash with key $rg\n";
		}
	}
//		# representative values (mean) for all read groups
	vector<double> is_mean_v;
	is_mean_v.reserve(is_mean.size());
	for (auto e : is_mean) {
		is_mean_v.push_back(e.second);
	}
	vector<double> is_median_v;
	is_median_v.reserve(is_median.size());
	for (auto e : is_median) {
		is_median_v.push_back(e.second);
	}

	vector<double> is_sd_v;
	is_sd_v.reserve(is_sd.size());
	for (auto e : is_sd) {
		is_sd_v.push_back(e.second);
	}

	vector<double> rl_v;
	rl_v.reserve(rl.size());
	for (auto e : rl) {
		rl_v.push_back(e.second);
	}

	int64_t all_is_mean = castle::MathUtils::mean(is_mean_v);
//		int64_t all_is_median = castle::MathUtils::mean(is_median_v);
	int64_t all_is_sd = castle::MathUtils::mean(is_sd_v);
	int64_t all_rl = castle::MathUtils::mean(rl_v);

	for (auto& an_entry : is_mean) {
		auto& rg = an_entry.first;
		O1 << (boost::format("rg%s\t%d\t%d\n") % rg % static_cast<int64_t>(is_mean[rg]) % static_cast<int64_t>(is_sd[rg])).str();
	}
	O1 << (boost::format("all\t%d\t%d\n") % all_is_mean % all_is_sd).str();

	for (auto& an_entry : rl) {
		O2 << (boost::format("rg%s\t%s\n") % an_entry.first % an_entry.second).str();
	}
	O2 << (boost::format("all\t%d\n") % all_rl).str();
}
void TEA::create_disc_FASTQs() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.create_disc_FASTQs");
	checker.start();
	string the_prefix = options.prefix;
	if (!options.working_dir.empty()) {
		the_prefix = options.working_prefix;
	}

	if(options.is_sampe) {
		string disc_file_name = the_prefix + ".disc.num.bam";
		string disc_bai_name = the_prefix + ".disc.num.bam.bai";
		string disc_bni_name = the_prefix + ".disc.num.bam.bni";
		string disc_out_file_name = the_prefix + ".disc.fq.gz";
		string disc_1_out_file_name = the_prefix + ".disc_1.fq.gz";
		string disc_2_out_file_name = the_prefix + ".disc_2.fq.gz";

		string cl_disc_file_name = the_prefix + ".cl.sorted.disc.num.bam";
		string cl_disc_bai_name = the_prefix + ".cl.sorted.disc.num.bam.bai";
		string cl_disc_bni_name = the_prefix + ".cl.sorted.disc.num.bam.bni";

		string cl_disc_out_file_name = the_prefix + ".cl.disc.fq.gz";
		string cl_disc_1_out_file_name = the_prefix + ".cl.disc_1.fq.gz";
		string cl_disc_2_out_file_name = the_prefix + ".cl.disc_2.fq.gz";
		if (options.is_force || (!boost::filesystem::exists(disc_1_out_file_name) && !boost::filesystem::exists(disc_2_out_file_name))) {
			BAM_to_FASTQ(disc_file_name, disc_bai_name, disc_bni_name, disc_out_file_name, disc_1_out_file_name, disc_2_out_file_name);
		}
		if (options.is_force || (!boost::filesystem::exists(cl_disc_1_out_file_name) && !boost::filesystem::exists(cl_disc_2_out_file_name))) {
			BAM_to_FASTQ(cl_disc_file_name, cl_disc_bai_name, cl_disc_bni_name, cl_disc_out_file_name, cl_disc_1_out_file_name, cl_disc_2_out_file_name);
		}
	} else {
		string disc_file_name = the_prefix + ".disc.sorted.bam";
		string disc_out_file_name = the_prefix + ".disc.fq.gz";
		string disc_1_out_file_name = the_prefix + ".disc_1.fq.gz";
		string disc_2_out_file_name = the_prefix + ".disc_2.fq.gz";

		string cl_disc_file_name = the_prefix + ".cl.sorted.disc.sorted.bam";
		string cl_disc_out_file_name = the_prefix + ".cl.disc.fq.gz";
		string cl_disc_1_out_file_name = the_prefix + ".cl.disc_1.fq.gz";
		string cl_disc_2_out_file_name = the_prefix + ".cl.disc_2.fq.gz";
		if (options.is_force || (!boost::filesystem::exists(disc_1_out_file_name) && !boost::filesystem::exists(disc_2_out_file_name))) {
			MEMBAM_to_FASTQ(disc_file_name, disc_1_out_file_name, disc_2_out_file_name);
		}

		if (options.is_force || (!boost::filesystem::exists(cl_disc_1_out_file_name) && !boost::filesystem::exists(cl_disc_2_out_file_name))) {
			MEMBAM_to_FASTQ(cl_disc_file_name, cl_disc_1_out_file_name, cl_disc_2_out_file_name);
		}
	}
	cout << checker;
}

void TEA::create_um_FASTQs() {
	string um_file_name = options.prefix + ".um.bam";
	string um_bai_name = options.prefix + ".um.bam.bai";
	string um_bni_name = options.prefix + ".um.bam.bni";
	string um_out_file_name = options.prefix + ".um.fq.gz";
	string um_1_out_file_name = options.prefix + ".um_1.fq.gz";
	string um_2_out_file_name = options.prefix + ".um_2.fq.gz";

	if (!options.working_dir.empty()) {
		um_file_name = options.working_prefix + ".um.bam";
		um_bai_name = options.working_prefix + ".um.bam.bai";
		um_bni_name = options.working_prefix + ".um.bam.bni";
		um_out_file_name = options.working_prefix + ".um.fq.gz";
		um_out_file_name = options.working_prefix + ".um.fq.gz";
		um_1_out_file_name = options.working_prefix + ".um_1.fq.gz";
		um_2_out_file_name = options.working_prefix + ".um_2.fq.gz";
	}
	if (boost::filesystem::exists(um_out_file_name) && boost::filesystem::exists(um_1_out_file_name) && boost::filesystem::exists(um_2_out_file_name)) {
		return;
	}
//	BAM_to_FASTQ_serial(um_file_name, um_out_file_name, um_1_out_file_name, um_2_out_file_name);
	BAM_to_FASTQ(um_file_name, um_bai_name, um_bni_name, um_out_file_name, um_1_out_file_name, um_2_out_file_name);
}

void TEA::generate_va_bams() {
	string um_1_file_name = options.prefix + ".um_1.fq.gz";
	string um_1_out_sam_name = options.prefix + ".um_1.fq.gz.sam";
	string um_1_out_tmp_bam_name = options.prefix + ".um_1.va.tmp.bam";
	string um_1_out_bam_name = options.prefix + ".um_1.va.bam";
	string um_2_file_name = options.prefix + ".um_2.fq.gz";
	string um_2_out_sam_name = options.prefix + ".um_2.fq.gz.sam";
	string um_2_out_tmp_bam_name = options.prefix + ".um_2.va.tmp.bam";
	string um_2_out_bam_name = options.prefix + ".um_2.va.bam";
	if (!options.working_dir.empty()) {
		um_1_file_name = options.working_prefix + ".um_1.fq.gz";
		um_1_out_sam_name = options.working_prefix + ".um_1.fq.gz.sam";
		um_1_out_tmp_bam_name = options.working_prefix + ".um_1.va.tmp.bam";
		um_1_out_bam_name = options.working_prefix + ".um_1.va.bam";
		um_2_file_name = options.working_prefix + ".um_2.fq.gz";
		um_2_out_sam_name = options.working_prefix + ".um_2.fq.gz.sam";
		um_2_out_tmp_bam_name = options.working_prefix + ".um_2.va.tmp.bam";
		um_2_out_bam_name = options.working_prefix + ".um_2.va.bam";
	}

	if (boost::filesystem::exists(um_1_out_bam_name) && boost::filesystem::exists(um_2_out_bam_name)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_va_bams");
	checker.start();

	vector<function<void()> > tasks;
	meerkat::BWACaller bc;
	bc.set_n_cores(n_cores);

	int64_t n_blocks_1 = bc.split_FASTQ_alt(um_1_file_name, um_2_file_name);
	bc.collect_single_align_tasks(tasks, n_blocks_1, um_1_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, um_1_file_name);
	bc.collect_single_align_tasks(tasks, n_blocks_1, um_2_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, um_2_file_name);

	if(static_cast<int64_t>(tasks.size()) < n_cores) {
		bc.align_single_reads(um_1_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, um_1_file_name);
		bc.align_single_reads(um_2_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, um_2_file_name);
		string um_1_sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % um_1_out_tmp_bam_name % um_1_out_sam_name).str();
		string um_2_sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % um_2_out_tmp_bam_name % um_2_out_sam_name).str();
		system(um_1_sam_to_bam_cmd.c_str());
		system(um_2_sam_to_bam_cmd.c_str());
	} else {
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		vector<string> um_1_bwa_errs;
		vector<string> um_2_bwa_errs;

		for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string um_1_bwa_err = um_1_file_name + ".bwa.err." + str_block_id;
			string um_2_bwa_err = um_2_file_name + ".bwa.err." + str_block_id;
			um_1_bwa_errs.push_back(um_1_bwa_err);
			um_2_bwa_errs.push_back(um_2_bwa_err);
		}
		bc.merge_SAM(um_1_file_name, um_1_out_tmp_bam_name, um_1_out_sam_name, n_blocks_1);
		bc.merge_SAM(um_2_file_name, um_2_out_tmp_bam_name, um_2_out_sam_name, n_blocks_1);
		castle::IOUtils::plain_file_merge(um_1_file_name + ".bwa.err", um_1_bwa_errs, n_cores, true);
		castle::IOUtils::plain_file_merge(um_2_file_name + ".bwa.err", um_2_bwa_errs, n_cores, true);
	}
	string um_1_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % um_1_out_bam_name % um_1_out_tmp_bam_name).str();
	string um_2_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % um_2_out_bam_name % um_2_out_tmp_bam_name).str();

	system(um_1_sort_cmd.c_str());
	system(um_2_sort_cmd.c_str());

	cout << checker;
}

void TEA::generate_vam_files() {
	// test.disc.bam test.disc.ram.bam test.disc.ram.bz2 test.disc_1.ra.bam test.disc_2.ra.bam
	// test.cl.disc.bam test.cl.disc.ram.bam test.cl.disc.ram.bz2 test.cl.disc_1.ra.bam test.cl.disc_2.ra.bam
	string refbam = options.prefix + ".um.bam";
	string vbamf = options.prefix + ".um.vam.bam";
	string vamf = options.prefix + ".um.vam";
	string um_1_va_bam = options.prefix + ".um_1.va.bam";
	string um_2_va_bam = options.prefix + ".um_2.va.bam";

	string the_merged_vbam = options.prefix + ".merged.vam.bam";
	string the_raw_vbam = options.prefix + ".vam.raw.bam";
	string the_vam = options.prefix + ".vam";
	string the_vbam = options.prefix + ".vam.bam";

	if (!options.working_dir.empty()) {
		refbam = options.working_prefix + ".um.bam";
		vbamf = options.working_prefix + ".um.vam.bam";
		vamf = options.working_prefix + ".um.vam";
		um_1_va_bam = options.working_prefix + ".um_1.ra.bam";
		um_2_va_bam = options.working_prefix + ".um_2.ra.bam";
		the_merged_vbam = options.working_prefix + ".merged.vam.bam";
		the_raw_vbam = options.working_prefix + ".vam.raw.bam";
		the_vam = options.working_prefix + ".vam";
		the_vbam = options.working_prefix + ".vam.bam";
	}
//	if (boost::filesystem::exists(the_vam) && boost::filesystem::exists(the_vbam)) {
//		return;
//	}
	generate_ram_file(refbam, vbamf, vamf, um_1_va_bam, um_2_va_bam, true);

	vector<string> vbam_files;
	vbam_files.push_back(vbamf);

	vector<string> vam_files;
	vam_files.push_back(vamf);

	castle::IOUtils::plain_file_merge_serial(the_merged_vbam, vbam_files, n_cores, false);
	castle::IOUtils::plain_file_merge_serial(the_vam, vam_files, n_cores, false);

	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % the_raw_vbam % the_merged_vbam).str();
	system(sam_to_bam_cmd.c_str());
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % the_vbam % the_raw_vbam).str();
	system(sambamba_sort_cmd.c_str());
}

void TEA::generate_um_bams(const string& a_path, const string& a_bai_path, const string& a_bni_path) {
	string bni_path;
	create_bni_even_index(a_path, a_bai_path, bni_path);
	int64_t size_block = 8192000;
//	vector<meerkat::BlockBoundary> fixed_size_blocks;
//	collect_boundaries_un(fixed_size_blocks, a_path, size_block);
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, a_path, a_bai_path, bni_path, size_block);

//	generate_um_raw_bam_serial();

	generate_um_raw_bam(local_unmapped_included_blocks);
	string in_um_raw_bam_name = options.prefix + ".um.raw.bam";
	string in_um_raw_bai_name = options.prefix + ".um.raw.bam.bai";
	string in_um_raw_bni_name = options.prefix + ".um.raw.bam.bni";
	string out_um_sam_name = options.prefix + ".um.sam";
	string out_um_bam_name = options.prefix + ".um.bam";
	string out_umm_name = options.prefix + ".umm";
	if (!options.working_dir.empty()) {
		in_um_raw_bam_name = options.working_prefix + ".um.raw.bam";
		in_um_raw_bai_name = options.working_prefix + ".um.raw.bam.bai";
		in_um_raw_bni_name = options.working_prefix + ".um.raw.bam.bni";
		out_um_sam_name = options.working_prefix + ".um.sam";
		out_um_bam_name = options.working_prefix + ".um.bam";
		out_umm_name = options.working_prefix + ".umm";
	}
	remove_duplicates(in_um_raw_bam_name, in_um_raw_bai_name, in_um_raw_bni_name, out_um_bam_name, out_um_sam_name);
//	write_duplicates(in_um_raw_bam_name, out_um_bam_name);
	generate_umm(out_um_bam_name, out_umm_name);
}

void TEA::generate_um_raw_bam_serial() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_um_raw_bam_serial");
	checker.start();
	string a_path = options.prefix + ".cl.sorted.bam";
	string out_um_raw_bam_name = options.prefix + ".um.raw.bam";
	string out_um_raw_sam_name = options.prefix + ".um.raw.sam";
	if (!options.working_dir.empty()) {
		out_um_raw_bam_name = options.working_prefix + ".um.raw.bam";
		out_um_raw_sam_name = options.working_prefix + ".um.raw.sam";
	}

	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		return;
	}

	const RefVector& refnames = local_reader.GetReferenceData();
	BamWriter um_raw_bam_file;
	if (!um_raw_bam_file.SAMOpen(out_um_raw_sam_name, local_reader.GetHeaderText(), refnames)) {
		cout << "ERROR: could not open output SAM file '" << out_um_raw_sam_name << "' for writing\n";
		exit(1);
	}
	BamAlignment local_alignment_entry;
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		if (!local_alignment_entry.IsMapped() || !local_alignment_entry.IsMateMapped()) {
			um_raw_bam_file.SaveSAMAlignment(local_alignment_entry);
		}
	}
	local_reader.Close();

	string sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % out_um_raw_bam_name % out_um_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	string bam_index_name = out_um_raw_bam_name + ".bai";
	string bam_index_cmd = (boost::format("sambamba index -t %d %s %s") % n_cores % out_um_raw_bam_name % bam_index_name).str();
	system(bam_index_cmd.c_str());
	cout << checker;
}

void TEA::generate_um_raw_bam(vector<meerkat::BlockBoundary>& actual_blocks) {
	string a_path = options.prefix + ".cl.sorted.bam";
	string out_um_raw_bam_name = options.prefix + ".um.raw.bam";
	if (!options.working_dir.empty()) {
		out_um_raw_bam_name = options.working_prefix + ".um.raw.bam";
	}

	if (boost::filesystem::exists(out_um_raw_bam_name)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_um_raw_bam");
	checker.start();
	vector<function<void()> > tasks;

	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	int64_t calculated_n_blocks = actual_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');
	vector<string> out_file_names(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, an_index_path)) {
				return;
			}
			int64_t num_total = 0;
//			int32_t the_current_ref_id = independent_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = independent_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = independent_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = independent_blocks[block_id + 1].pos;
//			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;

				string str_block_id = boost::lexical_cast<string>(block_id);

				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
				int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
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
//					cout << (boost::format("[TEA.generate_um_raw_bams] (start) Block-%d (%d/%d)-(%d/%d)\n")
//							% block_id % the_current_ref_id % the_current_ref_pos
//							% the_next_ref_id % the_next_ref_pos).str();
//				}
//				string rg_name;
				string out_um_raw_bam_name = options.prefix + ".um.raw.bam." + str_block_id;
				if (!options.working_dir.empty()) {
					out_um_raw_bam_name = options.working_prefix + ".um.raw.bam." + str_block_id;
				}
				out_file_names[block_id] = out_um_raw_bam_name;
				const RefVector& refnames = local_reader.GetReferenceData();
				BamWriter um_raw_bam_file;
				if(0 == block_id) {
					if (!um_raw_bam_file.SAMOpen(out_um_raw_bam_name,
									local_reader.GetHeaderText(), refnames)) {
						cout << "ERROR: could not open output SAM file '"
						<< out_um_raw_bam_name << "' for writing\n";
						exit(1);
					}
				} else {
					if (!um_raw_bam_file.SAMOpenNoHeader(out_um_raw_bam_name, refnames)) {
						cout << "ERROR: could not open output SAM file '"
						<< out_um_raw_bam_name << "' for writing\n";
						exit(1);
					}
				}
				BamAlignment local_alignment_entry;
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						cout << (boost::format("[TEA.generate_um_raw_bams] (first) Block-%d %s\n")
//								% block_id % a_block_boundary_str).str();
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
					if(!local_alignment_entry.IsMapped() || !local_alignment_entry.IsMateMapped()) {
						um_raw_bam_file.SaveSAMAlignment(local_alignment_entry);
					}
				}
				local_reader.Close();
				done_vector[block_id] = 'D';
//				if(prev_offset != the_next_ref_offset) {
//					cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] block-%d (last) offset: jump: prev(%d) cur(%d), calc: %d\n") % block_id % prev_offset % cur_offset % the_next_ref_offset).str();
//				}

//					if(verbose) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						cout << (boost::format("[ParallelDiscordExtractor.collect_clipped_discordants] (last) Block-%d %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//					else {
//						size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//						double processed = n/(double)done_vector.size() * 100.0;
//						cout << (boost::format("%.2f %%\n") % processed).str();
//					}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string out_um_raw_sam_name = options.prefix + ".um.raw.sam";
	if (!options.working_dir.empty()) {
		out_um_raw_sam_name = options.working_prefix + ".um.raw.sam";
	}

	castle::IOUtils::plain_file_merge(out_um_raw_sam_name, out_file_names, n_cores, true);

	string sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % out_um_raw_bam_name % out_um_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	string bam_index_name = out_um_raw_bam_name + ".bai";
	string bam_index_cmd = (boost::format("sambamba index -t %d %s %s") % n_cores % out_um_raw_bam_name % bam_index_name).str();
	system(bam_index_cmd.c_str());
	cout << checker;
}

void TEA::remove_duplicates_serial(const string& in_file_name, const string& out_file_name) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.remove_duplicates_serial");
	checker.start();
	string a_path(in_file_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	// <position, read_name, urec|mrec, alignment>
	boost::unordered_map<string, map<string, boost::unordered_map<string, BamAlignment>>> h;
	string pos;
	string oldpos;
	BamAlignment local_alignment_entry;
	BamWriter O;
	O.SAMOpen(out_file_name, reader.GetHeaderText(), reader.GetReferenceData());
	while (reader.LoadNextAlignmentCore(local_alignment_entry)) {
//        if (m/^@/) { print O; next }; # print header
//        chomp; @a=split(/\t/);
//        if ($a[1] =~ m/u/ && $a[1] =~ m/U/) { next }
		if (!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped()) {
			continue;
		}
		pos = (boost::format("%s:%s") % local_alignment_entry.RefID % local_alignment_entry.Position).str();
//		cout << (boost::format("[TEA.remove_duplicates] pos: %s\n") % pos).str();
//        $pos = "$a[2]:$a[3]";
//        my @x = keys %h;
//        #print "$a[0]:$pos:hash: ".scalar @x."\n";
		string xt_tag_str;
		local_alignment_entry.GetTag("XT", xt_tag_str);
//    	cout << (boost::format("[TEA.remove_duplicates] XT: %s\n") % xt_tag_str).str();
//    	string& the_read_name = local_alignment_entry.Name;
		auto an_itr = h.find(pos);
		if (h.end() != an_itr) {
//    		cout << (boost::format("[TEA.remove_duplicates] occurred before: %s\n") % pos).str();
			//        	# the same locus appeared before
			//            print "happened before\n" if $verbose;
			if (!local_alignment_entry.IsMapped()) {
				//# unmapped read
//                # check for the identical unmapped sequences
//                print "unmapped.. " if $verbose;
				bool dup = false;
//                auto& local_pos = an_itr->first;
				auto& read_name_map = an_itr->second;
				for (auto& read_name_entry : read_name_map) {
					auto& umrec_map = read_name_entry.second;
					for (auto& um_rec_entry : umrec_map) {
						auto& um_type = um_rec_entry.first;
						if ("urec" != um_type) {
							continue;
						}
						auto& b = um_rec_entry.second;
						if (local_alignment_entry.QueryBases == b.QueryBases) {
							dup = true;
							break;
						}
					}
				}
				if (!dup) {
					h[pos][local_alignment_entry.Name]["urec"] = local_alignment_entry;
//                    cout << (boost::format("[TEA.remove_duplicates] no duplicate. inserted into hash for locus %s:%s urec..\n") % pos % local_alignment_entry.Name).str();
				}
			} else if (!local_alignment_entry.IsMateMapped() && local_alignment_entry.MapQuality > 0 && (xt_tag_str.empty() || string::npos != xt_tag_str.find("U"))) {
				h[pos][local_alignment_entry.Name]["mrec"] = local_alignment_entry;
//            	  cout << (boost::format("[TEA.remove_duplicates] mapped(1) inserted into hash for %s:%s mrec..\n") % pos % local_alignment_entry.Name).str();
			}

		} else {
			//# the new locus
//            print "new locus.." if $verbose;
//            print "oldpos: $oldpos\n" if $verbose;
//        	cout << (boost::format("[TEA.remove_duplicates] new locus, oldpos: %s\n") % oldpos).str();
			auto an_old_itr = h.find(oldpos);
			if (h.end() != an_old_itr) {
//            if (exists $h{$oldpos}) {
//				cout << "[TEA.remove_duplicates] hash exists\n";
//                print "hash exists\n" if $verbose;
//                # print each pair if both mapped and unmapped read exist in the hash
//                auto& local_pos = an_old_itr->first;
				auto& read_name_map = an_old_itr->second;
				for (auto& read_name_entry : read_name_map) {
					auto& umrec_map = read_name_entry.second;
					auto mrec_itr = umrec_map.find("mrec");
					auto urec_itr = umrec_map.find("urec");
					if (umrec_map.end() != mrec_itr && umrec_map.end() != urec_itr) {
//						cout << (boost::format("[TEA.remove_duplicates] read name(0): %s\n") % mrec_itr->second.Name).str();
						if (mrec_itr->second.IsFirstMate()) {
							O.SaveSAMAlignment(mrec_itr->second);
							O.SaveSAMAlignment(urec_itr->second);
						} else {
							O.SaveSAMAlignment(urec_itr->second);
							O.SaveSAMAlignment(mrec_itr->second);
						}
					}
				}
//                # remove the items for the previous loci
//                delete $h{$oldpos};
				h.erase(oldpos);
//                print "deleted the hash item for locus $oldpos\n" if $verbose;
			}
			if (!local_alignment_entry.IsMapped()) {
//            if ($a[1] =~ m/u/) {
				//# unmapped read
				h[pos][local_alignment_entry.Name]["urec"] = local_alignment_entry;
				oldpos = pos;
//                cout << (boost::format("[TEA.remove_duplicates] unmapped(2) inserted into hash for %s:%s urec..\n") % pos % local_alignment_entry.Name).str();
			} else if (!local_alignment_entry.IsMateMapped() && local_alignment_entry.MapQuality > 0 && (xt_tag_str.empty() || string::npos != xt_tag_str.find("U"))) {
				// # unique map
				h[pos][local_alignment_entry.Name]["mrec"] = local_alignment_entry;
				oldpos = pos;
//                cout << (boost::format("[TEA.remove_duplicates] mapped(2) inserted into hash for %s:%s mrec..\n") % pos % local_alignment_entry.Name).str();
//                print "mapped.. inserted into hash for $pos:$a[0] mrec..\n" if $verbose;
			}
		}
	}
//# process the last loci
	for (auto& read_name_entry : h[oldpos]) {
//                	auto& read_name = read_name_entry.first;
		auto& umrec_map = read_name_entry.second;
		auto mrec_itr = umrec_map.find("mrec");
		auto urec_itr = umrec_map.find("urec");
		if (umrec_map.end() != mrec_itr && umrec_map.end() != urec_itr) {
//			cout << (boost::format("[TEA.remove_duplicates] read name(1): %s\n") % mrec_itr->second.Name).str();
			if (mrec_itr->second.IsFirstMate()) {
				O.SaveSAMAlignment(mrec_itr->second);
				O.SaveSAMAlignment(urec_itr->second);
			} else {
				O.SaveSAMAlignment(urec_itr->second);
				O.SaveSAMAlignment(mrec_itr->second);
			}
		}
	}
	O.Close();
	cout << checker;
}

void TEA::remove_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_bam_file_name, const string& out_sam_file_name) {
	if (boost::filesystem::exists(out_bam_file_name)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.remove_duplicates");
	checker.start();
	BamTools::BamReader reader;
	if (!reader.Open(a_path, a_bai_path)) {
		std::cout << "[TEA.remove_duplicates] ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	vector<function<void()> > tasks;
	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, a_path, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> independent_blocks = local_unmapped_included_blocks;
	int64_t calculated_n_blocks = independent_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');
	vector<string> out_file_names(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
				std::cout << "[TEA.remove_duplicates] ERROR: could not open BAM file '" << a_path << "'\n";
				exit(1);
			}

			string str_block_id = boost::lexical_cast<string>(block_id);

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

			// <position, read_name, urec|mrec, alignment>
				boost::unordered_map<string, map<string, boost::unordered_map<string, BamAlignment>>> h;
				string pos;
				string oldpos;

				BamAlignment local_alignment_entry;

				string local_out_file_name = out_sam_file_name + "." + str_block_id;
				out_file_names[block_id] = local_out_file_name;

				BamWriter O;
				if(0 == block_id) {
					O.SAMOpen(local_out_file_name, reader.GetHeaderText(), reader.GetReferenceData());
				} else {
					O.SAMOpenNoHeader(local_out_file_name, reader.GetReferenceData());
				}

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;

					if(!local_alignment_entry.IsMapped() && !local_alignment_entry.IsMateMapped()) {
						continue;
					}
					pos = (boost::format("%s:%s") % local_alignment_entry.RefID % local_alignment_entry.Position).str();
					string xt_tag_str;
					local_alignment_entry.GetTag("XT", xt_tag_str);
					auto an_itr = h.find(pos);
					if (h.end() != an_itr) {
						if (!local_alignment_entry.IsMapped()) {
							bool dup = false;
							auto& read_name_map = an_itr->second;
							for(auto& read_name_entry : read_name_map) {
								auto& umrec_map = read_name_entry.second;
								for(auto& um_rec_entry : umrec_map) {
									auto& um_type = um_rec_entry.first;
									if("urec" != um_type) {
										continue;
									}
									auto& b = um_rec_entry.second;
									if (local_alignment_entry.QueryBases == b.QueryBases) {
										dup = true;
										break;
									}
								}
							}
							if (!dup) {
								h[pos][local_alignment_entry.Name]["urec"] = local_alignment_entry;
							}
						} else if (!local_alignment_entry.IsMateMapped() && local_alignment_entry.MapQuality > 0 && (xt_tag_str.empty() || string::npos != xt_tag_str.find("U"))) {
							h[pos][local_alignment_entry.Name]["mrec"] = local_alignment_entry;
						}

					} else {
						auto an_old_itr = h.find(oldpos);
						if (h.end() != an_old_itr) {
							auto& read_name_map = an_old_itr->second;
							for(auto& read_name_entry: read_name_map) {
								auto& umrec_map = read_name_entry.second;
								auto mrec_itr = umrec_map.find("mrec");
								auto urec_itr = umrec_map.find("urec");
								if(umrec_map.end() != mrec_itr && umrec_map.end() != urec_itr) {
									if(mrec_itr->second.IsFirstMate()) {
										O.SaveSAMAlignment(mrec_itr->second);
										O.SaveSAMAlignment(urec_itr->second);
									} else {
										O.SaveSAMAlignment(urec_itr->second);
										O.SaveSAMAlignment(mrec_itr->second);
									}
								}
							}
							h.erase(oldpos);
						}
						if(!local_alignment_entry.IsMapped()) {
							h[pos][local_alignment_entry.Name]["urec"] = local_alignment_entry;
							oldpos = pos;
						} else if (!local_alignment_entry.IsMateMapped() && local_alignment_entry.MapQuality > 0 && (xt_tag_str.empty() || string::npos != xt_tag_str.find("U") )) {
							h[pos][local_alignment_entry.Name]["mrec"] = local_alignment_entry;
							oldpos = pos;
						}
					}
				}
				for (auto& read_name_entry : h[oldpos]) {
					auto& umrec_map = read_name_entry.second;
					auto mrec_itr = umrec_map.find("mrec");
					auto urec_itr = umrec_map.find("urec");
					if(umrec_map.end() != mrec_itr && umrec_map.end() != urec_itr) {
						if(mrec_itr->second.IsFirstMate()) {
							O.SaveSAMAlignment(mrec_itr->second);
							O.SaveSAMAlignment(urec_itr->second);
						} else {
							O.SaveSAMAlignment(urec_itr->second);
							O.SaveSAMAlignment(mrec_itr->second);
						}
					}
				}
				local_reader.Close();
				O.Close();
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge(out_sam_file_name, out_file_names, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % out_bam_file_name % out_sam_file_name).str();
	system(sam_to_bam_cmd.c_str());
	string bam_index_name = out_bam_file_name + ".bai";
	string bam_index_cmd = (boost::format("sambamba index -t %d %s %s") % n_cores % out_bam_file_name % bam_index_name).str();
	system(bam_index_cmd.c_str());
	cout << checker;
}

// After removing duplicate entries from a BAM in which entries with u or U flag presents
void TEA::generate_umm(const string& in_file_name, const string& out_file_name) {
	if (boost::filesystem::exists(out_file_name)) {
		return;
	}
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_umm");
	checker.start();
	string a_path(in_file_name);
	string an_index_path;
	get_bai_index_path(a_path, an_index_path);
	BamTools::BamReader reader;
	if (!reader.Open(a_path, an_index_path)) {
		std::cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	boost::unordered_set<string> h;
	BamAlignment local_alignment_entry;
	auto& ref_vec = reader.GetReferenceData();
	ofstream out(out_file_name, ios::binary);
	while (reader.LoadNextAlignmentCore(local_alignment_entry)) {
		string xt_tag_str;
		local_alignment_entry.GetTag("XT", xt_tag_str);
		if (!local_alignment_entry.IsMateMapped() && local_alignment_entry.IsMapped() && local_alignment_entry.MapQuality > 0 && (xt_tag_str.empty() || string::npos != xt_tag_str.find("U"))) {
			string r(local_alignment_entry.Name);
			boost::replace_all(r, "sc", "");
			boost::replace_all(r, "mu1", "");
			boost::replace_all(r, "mu2", "");
			if (h.end() == h.find(r)) {
				string ref_name = ref_vec[local_alignment_entry.RefID].RefName;
				boost::replace_all(ref_name, "chr", "");
				h.insert(r);
				// records coordinates of reads if their mates are unmapped
				if (local_alignment_entry.IsReverseStrand()) {
					out << (boost::format("%s\t%s\t-%s\tx\n") % local_alignment_entry.Name % ref_name % (local_alignment_entry.Position + 1)).str();
				} else {
					out << (boost::format("%s\t%s\t%s\tx\n") % local_alignment_entry.Name % ref_name % (local_alignment_entry.Position + 1)).str();
				}
			}
		}
	}
	cout << checker;
}

void TEA::write_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_file_name) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.write_duplicates");
	checker.start();
	BamTools::BamReader reader;
	if (!reader.Open(a_path, a_bai_path)) {
		std::cout << "[TEA.write_duplicates] ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	vector<function<void()> > tasks;
	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, a_path, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> independent_blocks = local_unmapped_included_blocks;
	int64_t calculated_n_blocks = independent_blocks.size();
	string done_vector(calculated_n_blocks - 1, 'U');
	vector<string> out_file_names(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
				std::cout << "[TEA.write_duplicates] ERROR: could not open BAM file '" << a_path << "'\n";
				exit(1);
			}

			string str_block_id = boost::lexical_cast<string>(block_id);

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

			BamAlignment local_alignment_entry;

			string local_out_file_name = out_file_name + "." + str_block_id;
			out_file_names[block_id] = local_out_file_name;
			BamWriter O;
			if(0 == block_id) {
				O.SAMOpen(local_out_file_name, reader.GetHeaderText(), reader.GetReferenceData());
			} else {
				O.SAMOpenNoHeader(local_out_file_name, reader.GetReferenceData());
			}

			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				cur_offset = m_bgzf.Tell();
				if(prev_offset >= the_next_ref_offset) {
					break;
				}
				prev_offset = cur_offset;
				O.SaveSAMAlignment(local_alignment_entry);
			}
			local_reader.Close();
			O.Close();
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void TEA::generate_ra_bams() {
	vector<function<void()> > tasks;
	string the_prefix = options.prefix;
	if(!options.working_dir.empty()) {
		the_prefix = options.working_prefix;
	}
	string disc_1_file_name = the_prefix + ".disc_1.fq.gz";
	string disc_1_out_sam_name = the_prefix + ".disc_1.fq.gz.sam";
	string disc_1_out_tmp_bam_name = the_prefix + ".disc_1.ra.tmp.bam";
	string disc_1_out_bam_name = the_prefix + ".disc_1.ra.bam";
	string disc_2_file_name = the_prefix + ".disc_2.fq.gz";
	string disc_2_out_sam_name = the_prefix + ".disc_2.fq.gz.sam";
	string disc_2_out_tmp_bam_name = the_prefix + ".disc_2.ra.tmp.bam";
	string disc_2_out_bam_name = the_prefix + ".disc_2.ra.bam";
	string cl_disc_1_file_name = the_prefix + ".cl.disc_1.fq.gz";
	string cl_disc_1_out_sam_name = the_prefix + ".cl.disc_1.fq.gz.sam";
	string cl_disc_1_out_tmp_bam_name = the_prefix + ".cl.disc_1.ra.tmp.bam";
	string cl_disc_1_out_bam_name = the_prefix + ".cl.disc_1.ra.bam";
	string cl_disc_2_file_name = the_prefix + ".cl.disc_2.fq.gz";
	string cl_disc_2_out_sam_name = the_prefix + ".cl.disc_2.fq.gz.sam";
	string cl_disc_2_out_tmp_bam_name = the_prefix + ".cl.disc_2.ra.tmp.bam";
	string cl_disc_2_out_bam_name = the_prefix + ".cl.disc_2.ra.bam";

	if (!options.is_force && boost::filesystem::exists(disc_1_out_bam_name) && boost::filesystem::exists(disc_2_out_bam_name) && boost::filesystem::exists(cl_disc_1_out_bam_name) && boost::filesystem::exists(cl_disc_2_out_bam_name)) {
		return;
	}

	meerkat::BWACaller bc;
	bc.set_n_cores(n_cores);
	bc.split_FASTQ_alt(disc_1_file_name, disc_2_file_name);
	int64_t n_blocks_disc_1 = bc.n_blocks_1;
	int64_t n_blocks_disc_2 = bc.n_blocks_2;
	int64_t n_blocks_cl_disc_1 = 0;
	int64_t n_blocks_cl_disc_2 = 0;

	bc.collect_single_align_tasks(tasks, n_blocks_disc_1, disc_1_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, disc_1_file_name);
	bc.collect_single_align_tasks(tasks, n_blocks_disc_2, disc_2_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, disc_2_file_name);
	if(boost::filesystem::exists(cl_disc_1_file_name) && boost::filesystem::exists(cl_disc_2_file_name)) {
		bc.split_FASTQ_alt(cl_disc_1_file_name, cl_disc_2_file_name);
		n_blocks_cl_disc_1 = bc.n_blocks_1;
		n_blocks_cl_disc_2 = bc.n_blocks_2;
		bc.collect_single_align_tasks(tasks, n_blocks_cl_disc_1, cl_disc_1_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, cl_disc_1_file_name);
		bc.collect_single_align_tasks(tasks, n_blocks_cl_disc_2, cl_disc_2_out_sam_name, options.repeat_reference, options.aln_param, options.samse_param, cl_disc_2_file_name);
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	vector<string> disc_1_bwa_errs;
	vector<string> disc_2_bwa_errs;

	for (int64_t block_id = 0; block_id < n_blocks_disc_1; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string disc_1_bwa_err = disc_1_file_name + ".bwa.err." + str_block_id;
		disc_1_bwa_errs.push_back(disc_1_bwa_err);
	}

	for (int64_t block_id = 0; block_id < n_blocks_disc_2; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string disc_2_bwa_err = disc_2_file_name + ".bwa.err." + str_block_id;
		disc_2_bwa_errs.push_back(disc_2_bwa_err);
	}

	bc.merge_SAM(disc_1_file_name, disc_1_out_tmp_bam_name, disc_1_out_sam_name, n_blocks_disc_1);
	bc.merge_SAM(disc_2_file_name, disc_2_out_tmp_bam_name, disc_2_out_sam_name, n_blocks_disc_2);
	string disc_1_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % disc_1_out_bam_name % disc_1_out_tmp_bam_name).str();
	string disc_2_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % disc_2_out_bam_name % disc_2_out_tmp_bam_name).str();

	if(boost::filesystem::exists(cl_disc_1_file_name) && boost::filesystem::exists(cl_disc_2_file_name)) {
		vector<string> cl_disc_1_bwa_errs;
		vector<string> cl_disc_2_bwa_errs;
		for (int64_t block_id = 0; block_id < n_blocks_cl_disc_1; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string cl_disc_1_bwa_err = cl_disc_1_file_name + ".bwa.err." + str_block_id;
			cl_disc_1_bwa_errs.push_back(cl_disc_1_bwa_err);
		}

		for (int64_t block_id = 0; block_id < n_blocks_cl_disc_2; ++block_id) {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string cl_disc_2_bwa_err = cl_disc_2_file_name + ".bwa.err." + str_block_id;
			cl_disc_2_bwa_errs.push_back(cl_disc_2_bwa_err);
		}
		bc.merge_SAM(cl_disc_1_file_name, cl_disc_1_out_tmp_bam_name, cl_disc_1_out_sam_name, n_blocks_cl_disc_1);
		bc.merge_SAM(cl_disc_2_file_name, cl_disc_2_out_tmp_bam_name, cl_disc_2_out_sam_name, n_blocks_cl_disc_2);
		string cl_disc_1_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % cl_disc_1_out_bam_name % cl_disc_1_out_tmp_bam_name).str();
		string cl_disc_2_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % cl_disc_2_out_bam_name % cl_disc_2_out_tmp_bam_name).str();
		system(cl_disc_1_sort_cmd.c_str());
		system(cl_disc_2_sort_cmd.c_str());
		castle::IOUtils::plain_file_merge(cl_disc_1_file_name + ".bwa.err", cl_disc_1_bwa_errs, n_cores, true);
		castle::IOUtils::plain_file_merge(cl_disc_2_file_name + ".bwa.err", cl_disc_2_bwa_errs, n_cores, true);
	}

	system(disc_1_sort_cmd.c_str());
	system(disc_2_sort_cmd.c_str());
	castle::IOUtils::plain_file_merge(disc_1_file_name + ".bwa.err", disc_1_bwa_errs, n_cores, true);
	castle::IOUtils::plain_file_merge(disc_2_file_name + ".bwa.err", disc_2_bwa_errs, n_cores, true);
}

void TEA::generate_ram_files() {
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_ram_files");
	checker.start();
//	  Usage:   bamram [options] <ref bam> <rambam> <ram bz2> <rbam1> <rbam2>
	// test.disc.bam test.disc.ram.bam test.disc.ram.bz2 test.disc_1.ra.bam test.disc_2.ra.bam
	// test.cl.disc.bam test.cl.disc.ram.bam test.cl.disc.ram.bz2 test.cl.disc_1.ra.bam test.cl.disc_2.ra.bam
	string the_prefix = options.prefix;
	if(!options.working_dir.empty()) {
		the_prefix = options.working_prefix;
	}

	string refbam = the_prefix + ".disc.num.bam";
	string rbamf = the_prefix + ".disc.ram.bam";
	string ramf = the_prefix + ".disc.ram";
	string disc_1_ra_bam = the_prefix + ".disc_1.ra.bam";
	string disc_2_ra_bam = the_prefix + ".disc_2.ra.bam";
	string cl_refbam = the_prefix + ".cl.sorted.disc.num.bam";
	string cl_rbamf = the_prefix + ".cl.disc.ram.bam";
	string cl_ramf = the_prefix + ".cl.disc.ram";
	string cl_disc_1_ra_bam = the_prefix + ".cl.disc_1.ra.bam";
	string cl_disc_2_ra_bam = the_prefix + ".cl.disc_2.ra.bam";

	string the_merged_rbam = the_prefix + ".merged.ram.bam";
	string the_raw_rbam = the_prefix + ".ram.raw.bam";
	string the_ram = the_prefix + ".ram";
	string the_rbam = the_prefix + ".ram.bam";

	if (!options.is_force && castle::IOUtils::get_file_size(the_ram) > 0 && castle::IOUtils::get_file_size(the_rbam) > 0) {
		cout << checker;
		return;
	}
	generate_ram_file(refbam, rbamf, ramf, disc_1_ra_bam, disc_2_ra_bam, options.exo, false);
	if(boost::filesystem::exists(cl_refbam) && boost::filesystem::exists(cl_disc_1_ra_bam) && boost::filesystem::exists(cl_disc_2_ra_bam)) {
		generate_ram_file(cl_refbam, cl_rbamf, cl_ramf, cl_disc_1_ra_bam, cl_disc_2_ra_bam, options.exo, true);
	}

	vector<string> rbam_files;
	rbam_files.push_back(rbamf);
	if(boost::filesystem::exists(cl_rbamf)) {
		rbam_files.push_back(cl_rbamf);
	}

	vector<string> ram_files;
	ram_files.push_back(ramf);
	if(boost::filesystem::exists(cl_ramf)) {
		ram_files.push_back(cl_ramf);
	}

	castle::IOUtils::plain_file_merge_serial(the_merged_rbam, rbam_files, n_cores, false);
	castle::IOUtils::plain_file_merge_serial(the_ram, ram_files, n_cores, false);

	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % the_raw_rbam % the_merged_rbam).str();
	system(sam_to_bam_cmd.c_str());
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % the_rbam % the_raw_rbam).str();
	system(sambamba_sort_cmd.c_str());
	cout << checker;
}

void TEA::generate_ram_file(const string& refbam, const string& rbamf, const string& ramf, const string& disc_1_ra_bam, const string& disc_2_ra_bam, bool exo, bool headless) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_ram_file");
	checker.start();
//	int64_t size_block = 8192000;
//	vector<meerkat::BlockBoundary> fixed_size_blocks;
//	collect_boundaries(fixed_size_blocks, refbam, size_block);

	string a_path(refbam);
	string an_index_path(a_path);
	an_index_path += ".bai";
	if (!boost::filesystem::exists(a_path) || !boost::filesystem::exists(an_index_path)) {
		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}

	vector<string> rabam_files;
	rabam_files.push_back(disc_1_ra_bam);
	rabam_files.push_back(disc_2_ra_bam);

//	# load repeat mappings into a hash
//	my %h = (); //# record end(1/2):READNAME\tREPAT_NAME1:fr1,REPEAT_NAME2:fr2, ...
	boost::unordered_map<string, string> h;
//	# save [strand:readname] -> [repeatname:seq] into a hash h
	load_repeat_mapping(h, exo, rabam_files);

//	# matching read ids from the reference bam and generate a ram file and a bam only with rams
	if(options.is_sampe) {
		write_ram_and_bam_serial(refbam, rbamf, ramf, h, exo, headless);
	} else {
		write_ram_and_bam_mem_serial(refbam, rbamf, ramf, h, exo, headless);
	}

	cout << checker;
}

void TEA::load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, vector<string>& rabam_files) {
	for (uint32_t file_id = 0; file_id < rabam_files.size(); ++file_id) {
		auto& rabam = rabam_files[file_id];
		int32_t end = file_id + 1;
		_load_repeat_mapping(h, exo, end, rabam);
	}
}

void TEA::_load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, const int32_t end, const string& rabam) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.load_repeat_mapping");
	checker.start();
	if (string::npos != rabam.find("./um")) {
		exo = true;
		cout << "[TEA.load_repeat_mapping] exogenous ..\n";
	}
	cout << (boost::format("[TEA.load_repeat_mapping] reading bamfile: %s, end: %d\n") % rabam % end).str();
	string a_bai_path = rabam + ".bai";
	string a_bni_path = rabam + ".bni";

	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;

	if(options.is_force) {
		boost::filesystem::remove_all(a_bni_path);
	}
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, rabam, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	collect_boundaries(actual_blocks, rabam, size_block);


	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

	// only for debugging
	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	if (verbose) {
//		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
//			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
//		}
//	}
	vector<boost::unordered_map<string, string>> h_lists(calculated_n_blocks - 1);
	vector<int64_t> mapcnt_lists(calculated_n_blocks - 1);
	vector<int64_t> dupmapcnt_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(rabam, a_bai_path)) {
				return;
			}
			int64_t num_total = 0;
			bool debug = false;
			vector<string> b;
			vector<string> rnames;
			const char* delim_semi_colon = ";";
			const char* delim_comma = ",";

			BamTools::BamAlignment local_alignment_entry;
//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.load_repeat_mapping] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.load_repeat_mapping] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			auto& ref_vec = local_reader.GetReferenceData();
			auto& local_h = h_lists[block_id];
			int64_t mapcnt = 0;
			int64_t dupmapcnt = 0;
			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				if(verbose && 0 == num_total) {
					string a_block_boundary_str = (boost::format("%s %d-%d %d")
							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
							% local_alignment_entry.AlignmentFlag).str();
					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
						cout << (boost::format("[TEA.load_repeat_mapping] Block-%d (first-wrong) %s\n")
								% block_id % a_block_boundary_str).str();
					} else {
						cout << (boost::format("[TEA.load_repeat_mapping] Block-%d (first) %s\n")
								% block_id % a_block_boundary_str).str();
					}
				}
//				if(local_alignment_entry.RefID == the_next_ref_id
//						&& local_alignment_entry.Position == the_next_ref_pos
//						&& local_alignment_entry.AlignmentFlag == the_next_aln_flag
//						&& local_alignment_entry.Name == the_next_block_read_name
//				) {
//					break;
//				}
//				debug = string::npos != local_alignment_entry.Name.find("HS2000-654_75:4:2206:13960:146767");
				if(debug) {
					cout << "[TEA.load_repeat_mapping] here-0\n";
				}
				cur_offset = m_bgzf.Tell();
				if(prev_offset >= the_next_ref_offset) {
					break;
				}

				prev_offset = cur_offset;
				++num_total;

				uint32_t nm = 0;
				local_alignment_entry.GetEditDistance(nm);
				//TODO: stringent mode
//				int64_t len = local_alignment_entry.Length;
				if(debug) {
					cout << "[TEA.load_repeat_mapping] here-1: " << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) <<"\n";
				}

				if (!local_alignment_entry.IsMapped()) {
				//TODO: stringent mode
//				if (!local_alignment_entry.IsMapped() || (len <= 50 && nm > 2)) {
					continue;
				}
				if(debug) {
					cout << "[TEA.load_repeat_mapping] here-2\n";
				}
				string key;
				string value;
				// check alternative mappings and list all TEs with similar mapping quality PolyA is reported with the top priority
				boost::unordered_map<string, bool> hte;
				string te = ref_vec[local_alignment_entry.RefID].RefName;
				hte[te] = true;//# a representative te type into a hash
				string xa_tag_str;
				if (local_alignment_entry.GetTag("XA", xa_tag_str)) {
					if (debug) {
						cout << xa_tag_str << "\n";
					}
					castle::StringUtils::c_string_multi_split(xa_tag_str, delim_semi_colon, b);
					// check alternative TE mappings
					for (uint64_t i = 0; i < b.size(); ++i) {
						castle::StringUtils::c_string_multi_split(b[i], delim_comma, rnames);
						if (rnames.size() < 4) {
							continue;
						}
						uint32_t local_nm = boost::lexical_cast<uint32_t>(rnames[3]);
						if (local_nm <= nm && hte.end() == hte.find(rnames[0])) {
							if ("PolyA" == rnames[0]) {
								te = "PolyA";
								break;
							}
							hte[rnames[0]] = true;
						}
					}
					if ("PolyA" != te) {
						vector<string> hte_keys;
						for (auto an_entry : hte) {
							hte_keys.push_back(an_entry.first);
						}
						sort(hte_keys.begin(), hte_keys.end());
						te = castle::StringUtils::join(hte_keys, ",");
					}
				}
				if(debug) {
					cout << "[TEA.load_repeat_mapping] here-3\n";
				}
				key = (boost::format("%s:%s") % end % local_alignment_entry.Name).str();
				value = (boost::format("%s:%s") % te % local_alignment_entry.QueryBases).str();
				local_h[key] = value;
				++mapcnt;
				if (debug) {
					cout << (boost::format("[TEA.load_repeat_mapping] te:%s\n%s\t%s\n") % te % key % value).str();
				}
				if (string::npos != te.find(",")) {
					++dupmapcnt;
				}
			}
			local_reader.Close();
			mapcnt_lists[block_id] = mapcnt;
			dupmapcnt_lists[block_id] = dupmapcnt;
			//				done_vector[block_id] = 'D';

			//				if(verbose) {
			//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
			//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
			//							% local_alignment_entry.AlignmentFlag).str();
			//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
			//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
			//								% block_id % a_block_boundary_str).str();
			//					} else {
			//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
			//								% block_id % a_block_boundary_str).str();
			//					}
			//				} else {
			//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
			//					double processed = n/(double)done_vector.size() * 100.0;
			//					cout << (boost::format("%.2f %%\n") % processed).str();
			//				}
		});
	}
	// sometimes the unaligned reads are in the last portion of BAM file, hence
	// changing the order.
	if (tasks.size() > 0) {
		swap(tasks[0], tasks.back());
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//	cout << "[TEA.load_repeat_mapping] gather scattered information\n";
	int64_t mapcnt = 0;
	int64_t dupmapcnt = 0;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		mapcnt += mapcnt_lists[block_id];
		dupmapcnt += dupmapcnt_lists[block_id];
		auto& local_h = h_lists[block_id];
		h.insert(local_h.begin(), local_h.end());
	}

	cout << (boost::format("[TEA.load_repeat_mapping] done processing %s: mapcnt: %s: dupmapcnt: %s\n") % rabam % mapcnt % dupmapcnt).str();
	cout << checker;
}

void TEA::write_ram_and_bam(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.write_ram_and_bam");
	checker.start();
	cout << (boost::format("[TEA.write_ram_and_bam] start matching the id to generate a BAM from %s ...\n") % refbam).str();
	string a_bai_path = refbam + ".bai";
	string a_bni_path = refbam + ".bni";

	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, refbam, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	vector<meerkat::BlockBoundary> actual_blocks;
//	collect_boundaries(actual_blocks, refbam, size_block);

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

	// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	if (verbose) {
//		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
//			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
//		}
//	}
	vector<boost::unordered_map<string, boost::unordered_map<string, RAMEntry>>> h_pos_lists(calculated_n_blocks - 1);
	vector<boost::unordered_map<string, string>> ram_lists(calculated_n_blocks - 1);
	vector<int64_t> cnt_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(refbam, a_bai_path)) {
				return;
			}

//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;
//			const bool debug = false;
			const bool debug = (0 == block_id);

			string previous_pos;
			string current_pos;
			//	h_pos[ram pos]->[TE sequence]->{bam_record}

				int64_t cnt = 0;
				string str;
				string map;
				string prname;
				//	# in the case of cl=1, generate unique map positions of rams among the potentially two read pairs
				//	# record read name => chr:pos:rname to generate unique pos
				//	# for the two pairs of reads originated from one pair
				//	# e.g. HWI-ST1001:7:2210:14748:43848mu1 and HWI-ST1001:7:2210:14748:43848mu2
				//	# with the exactly same mpos to the same repeat type

				vector<string> b;
				vector<string> c;
				const char* delim_colon = ":";
				const char* delim_tab = "\t";

				BamTools::BamAlignment local_alignment_entry;
				const auto& header = local_reader.GetHeaderText();
				string local_ram_filename = ramf + "." + str_block_id;
				string local_rbam_filename = rbamf + "." + str_block_id;

				ofstream O(local_ram_filename, ios::binary);
				BamTools::BamWriter O2;
				if(0 == block_id) {
					O2.SAMOpen(local_rbam_filename, header, local_reader.GetReferenceData());
				} else {
					O2.SAMOpenNoHeader(local_rbam_filename, local_reader.GetReferenceData());
				}
				auto& ref_vec = local_reader.GetReferenceData();
				auto& h_pos = h_pos_lists[block_id];
				auto& ram = ram_lists[block_id];
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first) %s\n")
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

					//		my @a = split(/\t/);
					//		$a[2] =~ s/chr//; # drop chr
					int64_t sam_flag = local_alignment_entry.AlignmentFlag;
					if ((sam_flag & 0x4) || 0 == local_alignment_entry.MapQuality) {
						continue;
					}

					string mate_key;
					//# end:1
					if (sam_flag & 0x40) {
						mate_key = (boost::format("2:%s") % local_alignment_entry.Name).str();
					} else if (sam_flag & 0x80) {
						mate_key = (boost::format("1:%s") % local_alignment_entry.Name).str();
					}
					if(debug) {
						cout << mate_key << "\n";
					}

					auto h_itr = h.find(mate_key);
					if (h.end() != h_itr) {
						string xt_tag_str;
						if (local_alignment_entry.GetTag("XT", xt_tag_str)) {
							//    		cout << xt_tag_str << "\n";
							if (string::npos == xt_tag_str.find("U")) {
								continue;
							}
							//    		(m/XT:A:/ && !m/XT:A:U/)
						}
						//	# check if its mate is mapped to TE
						auto mate = h_itr->second;

						//the mate variable: Alu:seq,...
						castle::StringUtils::c_string_multi_split(mate, delim_colon, b);
						current_pos = (boost::format("%s\t%s") % ref_vec[local_alignment_entry.RefID].RefName % local_alignment_entry.Position).str();
						map = (boost::format("%s\t%s\t%s") % local_alignment_entry.Name % current_pos % b[0]).str();//# append the read and  TE name

						if (debug) {
							cout << "mate: " << mate << "\n";
							cout << "current_pos: " << current_pos << "\n";
							cout << "previous_pos: " << previous_pos << "\n";
							cout << "map: " << map << "\n";
						}

						if (current_pos == previous_pos) {
							if (debug) {
								cout << "##the same pos\n";
							}
							// check the mate TE sequence is already in the hash, i.e., the same DNA fragment appeared before not
							auto h_pos_itr = h_pos.find(current_pos);
							if (h_pos.end() != h_pos_itr) {
								for (auto an_entry : h_pos_itr->second) {
									// s represents a TE sequence that appeared before
									auto& s = an_entry.first;
									//					}
									if (debug) {
										cout << "##the same pos has rams\n";
									}
									bool dup = false;
									if (debug) {
										cout << s << "\n";
									}
									if (s == b[1]) {
										if (debug) {
											cout << "matches with " << b[1] << "\n";
										}

										dup = true;
										auto& record = an_entry.second.bam;
										//							# compare the mapping quality of the current ram and the rams that appeared before
										if (local_alignment_entry.MapQuality > record.MapQuality) {
											if (debug) {
												cout << "choose the best quality\n" << map << "\n";
											}
											record = local_alignment_entry;
											an_entry.second.map = map;
										}
										break;
									}
									if (!dup) {
										// no match, thus add this record
										if (debug) {
											cout << "no dup\n";
										}

										h_pos[current_pos][b[1]].bam = local_alignment_entry;
										h_pos[current_pos][b[1]].map = map;
									}
								}
							} else {
								// if there is no hash, create a hash entry
								if (debug) {
									cout << "no ram at the same position before. create a new record\n" << b[1] << "\n" << local_alignment_entry.Name << "\n" << map << "\n";
								}
								h_pos[current_pos][b[1]].bam = local_alignment_entry;
								h_pos[current_pos][b[1]].map = map;
							}
						} else if (!previous_pos.empty()) {
							// print out rams for the previous_pos
							if (debug) {
								cout << "##new pos: print out rams for " << previous_pos << "\n";
							}

							auto seq_itr = h_pos.find(previous_pos);
							if (h_pos.end() != seq_itr) {
								// seq represents the TE sequences
								auto& seq = seq_itr->second;
								for (auto& s : seq) {
									auto& the_RAM_entry = s.second;
									auto& record1 = the_RAM_entry.map;
									auto& record2 = the_RAM_entry.bam;
									castle::StringUtils::c_string_multi_split(record1, delim_tab, c);
									// check whether the record with the same read name and ram pos exists
									// (only differ in the last suffix 1 or 2 for mu & sc) same read name exists
									prname = record2.Name;
									auto the_mu_pos = prname.rfind("mu");
									if (string::npos != the_mu_pos) {
										int64_t the_pair_id_pos = the_mu_pos + 2;
										char the_pair_id = prname[the_pair_id_pos];
										if ('1' == the_pair_id) {
											prname[the_pair_id_pos] = '2';
										} else if ('2' == the_pair_id) {
											prname[the_pair_id_pos] = '1';
										}
									}

									auto ram_itr = ram.find(prname);
									if (!exo || (exo && ram.end() == ram_itr) || (exo && ram_itr->second != record1)) {
										if (exo) {
											ram[prname] = map;
										}
										str = (boost::format("%s:\"%s\"") % record2.Name % c[3]).str();
										if (str.size() > 250) {
											str = str.substr(0, 250);
										}
										record2.Name = str;
										O << record1 << "\n";
										O2.SaveSAMAlignment(record2);
										++cnt;
									}
								}
							}
							h_pos[current_pos][b[1]].bam = local_alignment_entry;
							h_pos[current_pos][b[1]].map = map;
						} else {
							//				#print "adding $current_pos: $b[1]; $_; $map\n";
							h_pos[current_pos][b[1]].bam = local_alignment_entry;
							h_pos[current_pos][b[1]].map = map;
						}
					}
					previous_pos = current_pos;
				}
				local_reader.Close();
				O2.Close();
				cnt_lists[block_id] = cnt;
				//				done_vector[block_id] = 'D';

				//				if(verbose) {
				//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
				//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
				//							% local_alignment_entry.AlignmentFlag).str();
				//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					} else {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					}
				//				} else {
				//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
				//					double processed = n/(double)done_vector.size() * 100.0;
				//					cout << (boost::format("%.2f %%\n") % processed).str();
				//				}
			});
	}
	// sometimes the unaligned reads are in the last portion of BAM file, hence
	// changing the order.
	if (tasks.size() > 0) {
		swap(tasks[0], tasks.back());
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

//	cout << "[TEA.load_repeat_mapping] gather scattered information\n";
	int64_t cnt = 0;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		cnt += cnt_lists[block_id];
	}

	cout << (boost::format("[TEA.write_ram_and_bam] done generating ram and ram.bam with %d rams.\n") % cnt).str();
	cout << checker;
}
void TEA::write_ram_and_bam_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.write_ram_and_bam_serial");
	checker.start();
	cout << (boost::format("[TEA.write_ram_and_bam_serial] start matching the id to generate a BAM from %s ...\n") % refbam).str();
	string a_path(refbam);
	string an_index_path(a_path);
	an_index_path += ".bai";

	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	boost::unordered_map<string, string> ram;
	int64_t cnt = 0;
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}

	const bool debug = false;
//	const bool debug = true;

	string previous_pos;
	string current_pos;
	//	h_pos[ram pos]->[TE sequence]->{bam_record}

	string str;
	string map;
	string prname;
	//	# in the case of cl=1, generate unique map positions of rams among the potentially two read pairs
	//	# record read name => chr:pos:rname to generate unique pos
	//	# for the two pairs of reads originated from one pair
	//	# e.g. HWI-ST1001:7:2210:14748:43848mu1 and HWI-ST1001:7:2210:14748:43848mu2
	//	# with the exactly same mpos to the same repeat type

	vector<string> b;
	vector<string> c;
	const char* delim_colon = ":";
	const char* delim_tab = "\t";

	const auto& header = local_reader.GetHeaderText();
	string local_ram_filename = ramf;
	string local_rbam_filename = rbamf;

	ofstream out_ram(local_ram_filename, ios::binary);
	BamTools::BamWriter out_rbam;
	if (headless) {
		out_rbam.SAMOpenNoHeader(local_rbam_filename, local_reader.GetReferenceData());
	} else {
		out_rbam.SAMOpen(local_rbam_filename, header, local_reader.GetReferenceData());
	}

	auto& ref_vec = local_reader.GetReferenceData();
	BamTools::BamAlignment local_alignment_entry;
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		//		my @a = split(/\t/);
		//		$a[2] =~ s/chr//; # drop chr
//		debug = string::npos != local_alignment_entry.Name.find("HS2000-654_75:4:2206:13960:146767");
		if(debug) {
			cout << "[TEA.write_ram_and_bam_serial] " << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
		}
		int64_t sam_flag = local_alignment_entry.AlignmentFlag;
		if ((sam_flag & 0x4) || 0 == local_alignment_entry.MapQuality) {
			continue;
		}

		string mate_key;
		//# end:1
		if (local_alignment_entry.IsFirstMate()) {
			mate_key = (boost::format("2:%s") % local_alignment_entry.Name).str();
		} else if (local_alignment_entry.IsSecondMate()) {
			mate_key = (boost::format("1:%s") % local_alignment_entry.Name).str();
		}

		if(debug) {
			cout << "[TEA.write_ram_and_bam_serial] " << mate_key << "\n";
		}

		auto h_itr = h.find(mate_key);
		if (h.end() != h_itr) {
			string xt_tag_str;
			if (local_alignment_entry.GetTag("XT", xt_tag_str)) {
				if(debug) {
					cout << "[TEA.write_ram_and_bam_serial] " << xt_tag_str << "\n";
				}
				if (string::npos == xt_tag_str.find("U")) {
					continue;
				}
				//    		(m/XT:A:/ && !m/XT:A:U/)
			}
			//	# check if its mate is mapped to TE
			auto mate = h_itr->second;


			//the mate variable: Alu:seq,...
			castle::StringUtils::c_string_multi_split(mate, delim_colon, b);
			// the BAM is 0-index system, hence the position should be increased by 1
			if (local_alignment_entry.IsReverseStrand()) {
				current_pos = (boost::format("%s\t-%s") % ref_vec[local_alignment_entry.RefID].RefName % (local_alignment_entry.Position + 1)).str();
			} else {
				current_pos = (boost::format("%s\t%s") % ref_vec[local_alignment_entry.RefID].RefName % (local_alignment_entry.Position + 1)).str();
			}
			map = (boost::format("%s\t%s\t%s") % local_alignment_entry.Name % current_pos % b[0]).str();					//# append the read and  TE name

			if (debug) {
				cout << "[TEA.write_ram_and_bam_serial] mate: " << mate << "\n";
				cout << "[TEA.write_ram_and_bam_serial] current_pos: " << current_pos << "\n";
				cout << "[TEA.write_ram_and_bam_serial] previous_pos: " << previous_pos << "\n";
				cout << "[TEA.write_ram_and_bam_serial] map: " << map << "\n";
			}

			if (current_pos == previous_pos) {
				if (debug) {
					cout << "[TEA.write_ram_and_bam_serial] ##the same pos\n";
				}
				// check the mate TE sequence is already in the hash, i.e., the same DNA fragment appeared before not
				auto h_pos_itr = h_pos.find(current_pos);
				if (h_pos.end() != h_pos_itr) {
					for (auto an_entry : h_pos_itr->second) {
						// s represents a TE sequence that appeared before
						auto& s = an_entry.first;
						if (debug) {
							cout << "[TEA.write_ram_and_bam_serial] ##the same pos has rams\n";
						}
						bool dup = false;
						if (debug) {
							cout << "[TEA.write_ram_and_bam_serial] " << s << "/" << b[1] << "\n";
						}
						if (s == b[1]) {
							if (debug) {
								cout << "[TEA.write_ram_and_bam_serial] matches with " << b[1] << "\n";
							}

							dup = true;
							auto& record = an_entry.second.bam;
//							# compare the mapping quality of the current ram and the rams that appeared before
							if (local_alignment_entry.MapQuality > record.MapQuality) {
								if (debug) {
									cout << "[TEA.write_ram_and_bam_serial] choose the best quality\n" << map << "\n";
								}
								record = local_alignment_entry;
								an_entry.second.map = map;
							}
							break;
						}
						if (!dup) {
							// no match, thus add this record
							if (debug) {
								cout << "[TEA.write_ram_and_bam_serial] no dup\n";
							}

							h_pos[current_pos][b[1]].bam = local_alignment_entry;
							h_pos[current_pos][b[1]].map = map;
						}
					}
				} else {
					// if there is no hash, create a hash entry
					if (debug) {
						cout << "[TEA.write_ram_and_bam_serial] no ram at the same position before. create a new record\n" << b[1] << "\n" << local_alignment_entry.Name << "\n" << map << "\n";
					}
					h_pos[current_pos][b[1]].bam = local_alignment_entry;
					h_pos[current_pos][b[1]].map = map;
				}
			} else if (!previous_pos.empty()) {
				// print out rams for the previous_pos
				if (debug) {
					cout << "[TEA.write_ram_and_bam_serial] ##new pos: print out rams for " << previous_pos << "\n";
				}

				auto seq_itr = h_pos.find(previous_pos);
				if (h_pos.end() != seq_itr) {
					// seq represents the TE sequences
					auto& seq = seq_itr->second;
					for (auto& s : seq) {
						auto& the_RAM_entry = s.second;
						auto& record1 = the_RAM_entry.map;
						auto& record2 = the_RAM_entry.bam;
						castle::StringUtils::c_string_multi_split(record1, delim_tab, c);
						// check whether the record with the same read name and ram pos exists
						// (only differ in the last suffix 1 or 2 for mu & sc) same read name exists
						prname = record2.Name;
						auto the_mu_pos = prname.rfind("mu");
						if (string::npos != the_mu_pos) {
							int64_t the_pair_id_pos = the_mu_pos + 2;
							char the_pair_id = prname[the_pair_id_pos];
							if ('1' == the_pair_id) {
								prname[the_pair_id_pos] = '2';
							} else if ('2' == the_pair_id) {
								prname[the_pair_id_pos] = '1';
							}
						}

						auto ram_itr = ram.find(prname);
						if (!exo || (exo && ram.end() == ram_itr) || (exo && ram_itr->second != record1)) {
							if (exo) {
//								if(ram.end() != ram_itr) {
//									continue;
//								}
								ram[prname] = map;
							}
							str = (boost::format("%s:\"%s\"") % record2.Name % c[3]).str();
							if (str.size() > 250) {
								str = str.substr(0, 250);
							}
							record2.Name = str;
							if(debug) {
								cout << "[TEA.write_ram_and_bam_serial] ram: " << record1 << "\n";
							}
							out_ram << record1 << "\n";
							out_rbam.SaveSAMAlignment(record2);
							++cnt;
						}
					}
				}
				h_pos[current_pos][b[1]].bam = local_alignment_entry;
				h_pos[current_pos][b[1]].map = map;
			} else {
				//				#print "adding $current_pos: $b[1]; $_; $map\n";
				h_pos[current_pos][b[1]].bam = local_alignment_entry;
				h_pos[current_pos][b[1]].map = map;
			}
		}
		previous_pos = current_pos;
	}
	local_reader.Close();
	out_rbam.Close();

	cout << (boost::format("[TEA.write_ram_and_bam_serial] done generating ram and ram.bam with %d rams.\n") % cnt).str();
	cout << checker;
}

void TEA::write_ram_and_bam_mem_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.write_ram_and_bam_mem_serial");
	checker.start();
	cout << (boost::format("[TEA.write_ram_and_bam_mem_serial] start matching the id to generate a BAM from %s ...\n") % refbam).str();
	string a_path(refbam);
	string an_index_path(a_path);
	an_index_path += ".bai";

	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	boost::unordered_map<string, string> ram;
	int64_t cnt = 0;
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}

//	const bool debug = false;
//	const bool debug = true;

	string previous_pos;
	string current_pos;
	//	h_pos[ram pos]->[TE sequence]->{bam_record}

	string str;
	string map;
	string prname;
	//	# in the case of cl=1, generate unique map positions of rams among the potentially two read pairs
	//	# record read name => chr:pos:rname to generate unique pos
	//	# for the two pairs of reads originated from one pair
	//	# e.g. HWI-ST1001:7:2210:14748:43848mu1 and HWI-ST1001:7:2210:14748:43848mu2
	//	# with the exactly same mpos to the same repeat type

	vector<string> b;
	vector<string> c;
	const char* delim_colon = ":";
	const char* delim_tab = "\t";

	const auto& header = local_reader.GetHeaderText();
	string local_ram_filename = ramf;
	string local_rbam_filename = rbamf;

	ofstream out_ram(local_ram_filename, ios::binary);
	BamTools::BamWriter out_rbam;
	if (headless) {
		out_rbam.SAMOpenNoHeader(local_rbam_filename, local_reader.GetReferenceData());
	} else {
		out_rbam.SAMOpen(local_rbam_filename, header, local_reader.GetReferenceData());
	}

	auto& ref_vec = local_reader.GetReferenceData();
	BamTools::BamAlignment local_alignment_entry;
	set<string> written_set;
	const bool debug = false;
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {

		//		my @a = split(/\t/);
		//		$a[2] =~ s/chr//; # drop chr
//		const bool debug = string::npos != local_alignment_entry.Name.find("HSQ700642:191:D18JJACXX:1:2312:4453:42941");
		if(debug) {
			cout << "[TEA.write_ram_and_bam_mem_serial] " << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
		}
		int64_t sam_flag = local_alignment_entry.AlignmentFlag;
		if ((sam_flag & 0x4) || 0 == local_alignment_entry.MapQuality) {
			continue;
		}

		string mate_key_prefix;
		//# end:1
		if (local_alignment_entry.IsFirstMate()) {
			mate_key_prefix = (boost::format("2:%s") % local_alignment_entry.Name).str();
		} else if (local_alignment_entry.IsSecondMate()) {
			mate_key_prefix = (boost::format("1:%s") % local_alignment_entry.Name).str();
		}

		for(int64_t seq_id = 1; seq_id < 20; ++seq_id) {
			string mate_key = mate_key_prefix + ":m" + boost::lexical_cast<string>(seq_id);
			auto h_itr = h.find(mate_key);
			if (h.end() != h_itr) {
				if(debug) {
					cout << "[TEA.write_ram_and_bam_mem_serial] " << mate_key << "\n";
				}
				string xt_tag_str;
				if (local_alignment_entry.GetTag("XT", xt_tag_str)) {
					if(debug) {
						cout << "[TEA.write_ram_and_bam_mem_serial] " << xt_tag_str << "\n";
					}
					if (string::npos == xt_tag_str.find("U")) {
						continue;
					}
					//    		(m/XT:A:/ && !m/XT:A:U/)
				}
				//	# check if its mate is mapped to TE
				auto mate = h_itr->second;

				//the mate variable: Alu:seq,...
				castle::StringUtils::c_string_multi_split(mate, delim_colon, b);
				// the BAM is 0-index system, hence the position should be increased by 1
				if (local_alignment_entry.IsReverseStrand()) {
					current_pos = (boost::format("%s\t-%s\t%s") % ref_vec[local_alignment_entry.RefID].RefName % (local_alignment_entry.Position + 1) % seq_id).str();
				} else {
					current_pos = (boost::format("%s\t%s\t%s") % ref_vec[local_alignment_entry.RefID].RefName % (local_alignment_entry.Position + 1) % seq_id).str();
				}
				map = (boost::format("%s\t%s\t%s") % local_alignment_entry.Name % current_pos % b[0]).str();					//# append the read and  TE name

				if (debug) {
					cout << "[TEA.write_ram_and_bam_mem_serial] mate: " << mate << "\n";
					cout << "[TEA.write_ram_and_bam_mem_serial] current_pos: " << current_pos << "\n";
					cout << "[TEA.write_ram_and_bam_mem_serial] previous_pos: " << previous_pos << "\n";
					cout << "[TEA.write_ram_and_bam_mem_serial] map: " << map << "\n";
				}

				if (current_pos == previous_pos) {
					if (debug) {
						cout << "[TEA.write_ram_and_bam_mem_serial] ##the same pos\n";
					}
					// check the mate TE sequence is already in the hash, i.e., the same DNA fragment appeared before not
					auto h_pos_itr = h_pos.find(current_pos);
					if (h_pos.end() != h_pos_itr) {
						for (auto an_entry : h_pos_itr->second) {
							// s represents a TE sequence that appeared before
							auto& s = an_entry.first;
							if (debug) {
								cout << "[TEA.write_ram_and_bam_mem_serial] ##the same pos has rams\n";
							}
							bool dup = false;
							if (debug) {
								cout << "[TEA.write_ram_and_bam_mem_serial] " << s << "/" << b[1] << "\n";
							}
							if (s == b[1]) {
								if (debug) {
									cout << "[TEA.write_ram_and_bam_mem_serial] matches with " << b[1] << "\n";
								}

								dup = true;
								auto& record = an_entry.second.bam;
	//							# compare the mapping quality of the current ram and the rams that appeared before
								if (local_alignment_entry.MapQuality > record.MapQuality) {
									if (debug) {
										cout << "[TEA.write_ram_and_bam_mem_serial] choose the best quality\n" << map << "\n";
									}
									record = local_alignment_entry;
									an_entry.second.map = map;
								}
								break;
							}
							if (!dup) {
								// no match, thus add this record
								if (debug) {
									cout << "[TEA.write_ram_and_bam_mem_serial] no dup\n";
								}

								h_pos[current_pos][b[1]].bam = local_alignment_entry;
								h_pos[current_pos][b[1]].map = map;
							}
						}
					} else {
						// if there is no hash, create a hash entry
						if (debug) {
							cout << "[TEA.write_ram_and_bam_mem_serial] no ram at the same position before. create a new record\n" << b[1] << "\n" << local_alignment_entry.Name << "\n" << map << "\n";
						}
						h_pos[current_pos][b[1]].bam = local_alignment_entry;
						h_pos[current_pos][b[1]].map = map;
					}
				} else if (!previous_pos.empty()) {
					auto seq_itr = h_pos.find(previous_pos);
					if (h_pos.end() != seq_itr) {
						// print out rams for the previous_pos
						if (debug) {
							cout << "[TEA.write_ram_and_bam_mem_serial] ##new pos: print out rams for " << previous_pos << "\n";
						}
						// seq represents the TE sequences
						auto& seq = seq_itr->second;
						for (auto& s : seq) {
							auto& the_RAM_entry = s.second;
							auto& record1 = the_RAM_entry.map;
							auto& record2 = the_RAM_entry.bam;
							castle::StringUtils::c_string_multi_split(record1, delim_tab, c);
							// check whether the record with the same read name and ram pos exists
							// (only differ in the last suffix 1 or 2 for mu & sc) same read name exists
							prname = record2.Name;
							auto the_mu_pos = prname.rfind("mu");
							if (string::npos != the_mu_pos) {
								int64_t the_pair_id_pos = the_mu_pos + 2;
								char the_pair_id = prname[the_pair_id_pos];
								if ('1' == the_pair_id) {
									prname[the_pair_id_pos] = '2';
								} else if ('2' == the_pair_id) {
									prname[the_pair_id_pos] = '1';
								}
							}

							auto ram_itr = ram.find(prname);
							if (!exo || (exo && ram.end() == ram_itr) || (exo && ram_itr->second != record1)) {
								if (exo) {
	//								if(ram.end() != ram_itr) {
	//									continue;
	//								}
									ram[prname] = map;
								}
								str = (boost::format("%s:\"%s\"") % record2.Name % c[c.size() - 1]).str();
								if (str.size() > 250) {
									str = str.substr(0, 250);
								}
								record2.Name = str;
								string key = c[0] + c[1] + c[2] + c[4];
								if(written_set.end() == written_set.find(key)) {
									if(debug) {
										cout << "[TEA.write_ram_and_bam_mem_serial] key: " << key << "\n";
										cout << "[TEA.write_ram_and_bam_mem_serial] ram: " << record1 << "\n";
									}
									out_ram << record1 << "\n";
									out_rbam.SaveSAMAlignment(record2);
									written_set.insert(key);
									++cnt;
								}
							}
						}
					}
					h_pos[current_pos][b[1]].bam = local_alignment_entry;
					h_pos[current_pos][b[1]].map = map;
				} else {
					//				#print "adding $current_pos: $b[1]; $_; $map\n";
					h_pos[current_pos][b[1]].bam = local_alignment_entry;
					h_pos[current_pos][b[1]].map = map;
				}
			} else {
				if(debug) {
					cout << "[TEA.write_ram_and_bam_mem_serial] =================== BREAK =================== " << h_pos.size() << "\n";
				}
//				h_pos.clear();
				break;
			}
			previous_pos = current_pos;
		}

	}
	local_reader.Close();
	out_rbam.Close();

	cout << (boost::format("[TEA.write_ram_and_bam_mem_serial] done generating ram and ram.bam with %d rams.\n") % cnt).str();
	cout << checker;
}

void TEA::generate_cbam_files() {
	if(options.is_sampe) {
		_generate_cbam_files_sampe();
	} else {
		_generate_cbam_files_mem();
//		_generate_cbam_files_mem_alt();
//		_generate_cbam_files_mem_alt2();
//		_generate_cbam_files_mem_org();
	}
}

void TEA::_generate_cbam_files_mem() {
	string original_input_bam_name = options.prefix + ".bam";
	string firststat_name = options.prefix + ".firststat";
	string input_bam_name = options.prefix + ".disc.sorted.bam";
	string input_bam_bfi_name = options.prefix + ".disc.sorted.bam.bfi";
	string consd_raw_sam_name = options.prefix + ".softclips.consd.raw.sam";
	string consd_raw_bam_name = options.prefix + ".softclips.consd.raw.bam";
	string consd_sorted_bam_name = options.prefix + ".softclips.consd.sorted.bam";
	string consd_sorted_sam_name = options.prefix + ".softclips.consd.sorted.sam";
	string consd_sorted_bam_bni_name = options.prefix + ".softclips.consd.sorted.bam.bni";
	string consd_sorted_raw_sam_name = options.prefix + ".softclips.consd.sorted.raw.sam";
	string consd_sorted_raw_bam_name = options.prefix + ".softclips.consd.sorted.raw.bam";
	string out_softclips_consd_bam_name = options.prefix + ".softclips.consd.bam";
	string out_softclips_consd_cpos_name = options.prefix + ".softclips.consd.cpos";
	if (!options.working_dir.empty()) {
		firststat_name = options.working_prefix + ".firststat";
		input_bam_name = options.working_prefix + ".disc.sorted.bam";
		input_bam_bfi_name = options.working_prefix + ".disc.sorted.bam.bfi";
		consd_raw_sam_name = options.working_prefix + ".softclips.consd.raw.sam";
		consd_raw_bam_name = options.working_prefix + ".softclips.consd.raw.bam";
		consd_sorted_bam_name = options.working_prefix + ".softclips.consd.sorted.bam";
		consd_sorted_sam_name = options.working_prefix + ".softclips.consd.sorted.sam";
		consd_sorted_bam_bni_name = options.working_prefix + ".softclips.consd.sorted.bam.bni";
		consd_sorted_raw_sam_name = options.working_prefix + ".softclips.consd.sorted.raw.sam";
		consd_sorted_raw_bam_name = options.working_prefix + ".softclips.consd.sorted.raw.bam";
		out_softclips_consd_bam_name = options.working_prefix + ".softclips.consd.bam";
		out_softclips_consd_cpos_name = options.working_prefix + ".softclips.consd.cpos";
	}

	if (!options.is_force && castle::IOUtils::get_file_size(out_softclips_consd_bam_name) > 0 && castle::IOUtils::get_file_size(out_softclips_consd_cpos_name) > 0) {
		return;
	}

	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_cbam_files_mem");
	checker.start();

	vector<int64_t> block_boundary;
	string an_original_index_path = original_input_bam_name + ".bai";

	if(!boost::filesystem::exists(original_input_bam_name)) {
		cout << (boost::format("[TEA.generate_cbam_files_mem] cannot read: %s\n") % original_input_bam_name).str();
		cout << checker;
		return;
	}

	// create a .softclips.consd.raw.bam file, which contains a BAM having converted of CIGAR string 'H' -> 'S'.
	string an_index_path = input_bam_name + ".bfi";

	if (boost::filesystem::exists(an_index_path)) {
		castle::ParallelRunner::load_bfi_index(block_boundary, an_index_path);
	} else {
		const int64_t N_BLOCK_ENTRIES = 262144;
		castle::ParallelRunner::create_bfi_index(block_boundary, input_bam_name, input_bam_bfi_name, N_BLOCK_ENTRIES);
	}
	{
		string line;
		ifstream in(firststat_name, ios::binary);
		for(int64_t e_id = 0; e_id < 5; ++e_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}
	vector<function<void()> > tasks;

	int64_t n_block_boundaries = block_boundary.size();
	cout << (boost::format("[TEA.generate_cbam_files_mem] # blocks: %d of .disc.sorted.bam\n") % (n_block_boundaries - 1)).str();
	vector<string> raw_bam_name_lists(n_block_boundaries - 1);
	vector<string> softclips_consd_cpos_name_lists(n_block_boundaries - 1);

	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(input_bam_name)) {
				return;
			}
			const RefVector& m_ref = local_reader.GetReferenceData();
			boost::unordered_map<string, int64_t> ref_reverse_index;
			for (uint64_t ref_id = 0; ref_id < m_ref.size(); ++ref_id) {
				auto& a_ref = m_ref[ref_id];
				ref_reverse_index[a_ref.RefName] = ref_id;
			}
			auto& local_data = local_reader.GetBGZF();
			if(0 != block_boundary[block_id]) {
				local_data.Seek(block_boundary[block_id]);
			}
			int64_t next_block_pos = block_boundary[block_id + 1];
			int64_t cur_pos = -1;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();

			string local_raw_sam_name = consd_raw_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;
//			bool debug = false;
			BamAlignment local_alignment_entry;
			BamWriter out_raw_sam;
			if(0 == block_id) {
				out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
			} else {
				out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
			}

//			ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);
			vector<BamAlignment> pair1_alns;
			vector<BamAlignment> pair2_alns;
			string prev_name;

			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				auto& cur_name = local_alignment_entry.Name;
				if (prev_name != cur_name) {
					if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
						remove_entry_enclosed_with_large_H(pair1_alns);
						remove_entry_enclosed_with_large_H(pair2_alns);
						convert_H_to_S(pair1_alns);
						convert_H_to_S(pair2_alns);
						for(auto& an_aln : pair1_alns) {
							if(has_S_or_H(an_aln)) {
								out_raw_sam.SaveSAMAlignment(an_aln);
							}
						}
						for(auto& an_aln : pair2_alns) {
							if(has_S_or_H(an_aln)) {
								out_raw_sam.SaveSAMAlignment(an_aln);
							}
						}
						pair1_alns.clear();
						pair2_alns.clear();
					}
					prev_name = cur_name;
				}
				if(local_alignment_entry.IsFirstMate()) {
					pair1_alns.push_back(local_alignment_entry);
				} else {
					pair2_alns.push_back(local_alignment_entry);
				}
				cur_pos = local_data.Tell();
				if(cur_pos >= next_block_pos) {
					break;
				}
			}
			if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
				remove_entry_enclosed_with_large_H(pair1_alns);
				remove_entry_enclosed_with_large_H(pair2_alns);
				convert_H_to_S(pair1_alns);
				convert_H_to_S(pair2_alns);
				for(auto& an_aln : pair1_alns) {
					if(has_S_or_H(an_aln)) {
						out_raw_sam.SaveSAMAlignment(an_aln);
					}
				}
				for(auto& an_aln : pair2_alns) {
					if(has_S_or_H(an_aln)) {
						out_raw_sam.SaveSAMAlignment(an_aln);
					}
				}
				pair1_alns.clear();
				pair2_alns.clear();
			}
			local_reader.Close();
			out_raw_sam.Close();
		});
	}


	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << (boost::format("[TEA.generate_cbam_files_mem] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + options.qcutoff)).str();
	castle::IOUtils::plain_file_merge(consd_raw_sam_name, raw_bam_name_lists, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_raw_bam_name % consd_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem] done generating consd raw bam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem] start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % consd_sorted_bam_name % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());

	// create .cpos .softclips.consd.bam files
	if(options.is_force) {
		boost::filesystem::remove_all(consd_sorted_bam_bni_name);
	}
	int64_t size_block = 8192000;
	vector<int64_t> unmapped_block_boundaries;

	castle::ParallelRunner::load_bai_bni_index(unmapped_block_boundaries, consd_sorted_bam_name, consd_sorted_bam_bni_name, size_block, n_cores);

	int64_t calculated_n_blocks = unmapped_block_boundaries.size();
	cout << (boost::format("[TEA.generate_cbam_files_mem] # blocks: %d of .consd.sorted.bam\n") % (calculated_n_blocks - 1)).str();
	raw_bam_name_lists.clear();
	softclips_consd_cpos_name_lists.clear();
	raw_bam_name_lists.resize(calculated_n_blocks - 1);
	softclips_consd_cpos_name_lists.resize(calculated_n_blocks - 1);
	boost::mutex print_mutex;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(consd_sorted_bam_name, consd_sorted_bam_name + ".bai")) {
				return;
			}
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t the_current_ref_offset = unmapped_block_boundaries[block_id];
			int64_t the_next_ref_offset = unmapped_block_boundaries[block_id + 1];
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;

			uint32_t qcutoff = options.qcutoff;
			uint32_t max_mismatches = options.max_mismatches;
			int32_t min_matches = options.min_matches;
			int32_t min_polyAT = options.min_polyAT;

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();
			string cseq1;
			string cseq2;
			uint32_t clen = 0;
			bool polyAT = false;
			uint32_t nm = 0;

			string qual1;
			string qual2;

			bool selected = false;
			bool selected2 = false;
			int32_t cpos;

			string poly_a_str(min_polyAT, 'A');
			string poly_t_str(min_polyAT, 'T');
			string a_delim_str("^");
			for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			vector<string> b;
			const char* delim_all_stopwords = a_delim_str.c_str();

			string local_raw_sam_name = consd_sorted_raw_sam_name + "." + str_block_id;
			string local_sorted_sam_name = consd_sorted_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;

			BamTools::BamAlignment local_alignment_entry;
			BamTools::BamWriter out_raw_sam;
//			BamTools::BamWriter out_sorted_sam;
			if(0 == block_id) {
				out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
//				out_sorted_sam.SAMOpen(local_sorted_sam_name, header, ref_vec);
			} else {
				out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
//				out_sorted_sam.SAMOpenNoHeader(local_sorted_sam_name, ref_vec);
			}
//			bool debug = (0 == block_id);
			const bool debug = false;
//			const bool debug = true;
			ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);

			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				cur_offset = m_bgzf.Tell();
				if(prev_offset >= the_next_ref_offset) {
					break;
				}
				prev_offset = cur_offset;
				++num_total;

//				debug = "DHFC08P1:370:C1784ACXX:2:1214:21254:59267" == local_alignment_entry.Name;
//				if(!debug) {
//					continue;
//				}

//				out_sorted_sam.SaveSAMAlignment(local_alignment_entry);
				//	while(<F>) {
				//		chomp; my @a=split(/\t/);
				//		if ($a[1] & 0x400 | $a[1] & 0x4) { next }; # skip PCR duplicate or unmapped reads
				uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
				if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
					continue;
				}

				auto& seq = local_alignment_entry.QueryBases;
				auto& qual = local_alignment_entry.Qualities;
				//		$s1 = $cigar;
				//		$s2 = $cigar;
				polyAT = selected = selected2 = false;
				auto& cigar = local_alignment_entry.CigarData;
				// # clipped in the beginning
				try {
				if (!cigar.empty() && ('S' == cigar.front().Type)) {
					if(debug) {
						boost::lock_guard<boost::mutex> a_lock(print_mutex);
						cout << (boost::format("[TEA.generate_cbam_files_mem] initial clipped: %s\n") % BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec)).str();
					}
					clen = cigar.front().Length;
					if (clen >= 5) {
						qual1 = qual.substr(0, clen);
						if(debug) {
							boost::lock_guard<boost::mutex> a_lock(print_mutex);
							cout << "[TEA.generate_cbam_files_mem] qual1: " << qual1 << "\n";
						}
						//# more than 5 good quality bases
						int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
						if (n_good_quals >= 5) {
							if(debug) {
								boost::lock_guard<boost::mutex> a_lock(print_mutex);
								cout << "[TEA.generate_cbam_files_mem] selected 1-1\n";
							}
							selected = true;
//							if (!(sam_flag & 0x0010)) { // # positive strand mapping
							if (!local_alignment_entry.IsReverseStrand()) { // # positive strand mapping
								if(debug) {
									boost::lock_guard<boost::mutex> a_lock(print_mutex);
									cout << "[TEA.generate_cbam_files_mem] selected 2-1\n";
								}
								selected2 = true;
								cpos = local_alignment_entry.Position;
								//								# check whether the clipped seq has >=10
								cseq1 = seq.substr(0, clen);
								if (string::npos != cseq1.find(poly_a_str)) {
									if(debug) {
										boost::lock_guard<boost::mutex> a_lock(print_mutex);
										cout << "[TEA.generate_cbam_files_mem] polyAT 1\n";
									}
									polyAT = true;
								}
							}
						}
					}
				}
				// # clipping in the end
				clen = 0;
//				for (auto& a_cigar : cigar) {
//					if ('D' == a_cigar.Type) {
//						clen = a_cigar.Length;
//						break;
//					}
//				}
				if(!cigar.empty() && ('S' == cigar.back().Type)) {
					clen = cigar.back().Length;
				}
				if (clen >= 5) {
					int64_t the_qual_sub_str_pos = qual.size() - clen;
					if (the_qual_sub_str_pos > 0) {
						qual2 = qual.substr(the_qual_sub_str_pos);
					} else {
						qual2 = qual;
					}
					if(debug) {
						boost::lock_guard<boost::mutex> a_lock(print_mutex);
						cout << "[TEA.generate_cbam_files_mem] qual2: " << qual2 << "\n";
					}
					int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
					//		if ($qual2 =~ m/([^#]+)/) {
					//				if (length($1) >= 5 ) {
					if (n_good_quals >= 5) {
						if(debug) {
							boost::lock_guard<boost::mutex> a_lock(print_mutex);
							cout << "[TEA.generate_cbam_files_mem] selected 1-2\n";
						}
						selected = true;
//						if (sam_flag & 0x0010) { // # negative strand mapping
						if (local_alignment_entry.IsReverseStrand()) {
							if(debug) {
								boost::lock_guard<boost::mutex> a_lock(print_mutex);
								cout << "[TEA.generate_cbam_files_mem] selected 2-2\n";
							}
							selected2 = true;
							cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
							//							# check whether the clipped seq has >=10
							int64_t the_seq_sub_str_pos = seq.size() - clen;
							cseq2 = seq.substr(the_seq_sub_str_pos);
							//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
							if (string::npos != cseq2.find(poly_t_str)) {
								if(debug) {
									boost::lock_guard<boost::mutex> a_lock(print_mutex);
									cout << "[TEA.generate_cbam_files_mem] polyAT 2\n";
								}
								polyAT = true;
							}
						}
					}
				}

				//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
				if (!polyAT) {
					//        # check max_mismatches
					local_alignment_entry.GetEditDistance(nm);
					if (nm > max_mismatches) {
						continue;
					};

					string md_str;
					local_alignment_entry.GetTag("MD", md_str);
					//        $md = $_;

					//        # check min_matches is satisfied

					//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
					//        my @b=split(/[A-Za-z^]/, $md);
					//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
					castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
					int64_t summed_n_bases = 0;
					for (auto a_value_str : b) {
						summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
					}
					if (summed_n_bases < min_matches) {
						continue;
					}
				}

				if (selected2) {
					out_raw_sam.SaveSAMAlignment(local_alignment_entry);
					out_consd_cpos << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
				}
				} catch(exception& ex) {
					cout << ex.what() << "\n";
					cout << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
				}

			}
			local_reader.Close();
			out_raw_sam.Close();
		});
	}

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	castle::IOUtils::plain_file_merge(consd_sorted_raw_sam_name, raw_bam_name_lists, n_cores, true);
	string sam_to_sorted_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_sorted_raw_bam_name % consd_sorted_raw_sam_name).str();
	system(sam_to_sorted_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem] done generating consd sorted raw bam: %s\n") % consd_sorted_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem] start sorting and generating the index for %s ...\n") % consd_sorted_raw_bam_name).str();
	string sambamba_last_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_softclips_consd_bam_name % consd_sorted_raw_bam_name).str();
	system(sambamba_last_sort_cmd.c_str());
	castle::IOUtils::plain_file_merge(out_softclips_consd_cpos_name, softclips_consd_cpos_name_lists, n_cores, true);
	cout << "[TEA.generate_cbam_files_mem] done generating cbam and its index.\n";

	cout << checker;
}

void TEA::_generate_cbam_files_mem_alt() {
	string the_prefix = options.prefix;
	string input_bam_name = the_prefix + ".bam";
	if (!options.working_dir.empty()) {
		the_prefix = options.working_prefix;
	}
	string firststat_name = the_prefix + ".firststat";
	string consd_raw_sam_name = the_prefix + ".softclips.consd.raw.sam";
	string consd_raw_bam_name = the_prefix + ".softclips.consd.raw.bam";
	string out_softclips_consd_bam_name = the_prefix + ".softclips.consd.bam";
	string out_softclips_consd_cpos_name = the_prefix + ".softclips.consd.cpos";

//	if (!options.is_force && boost::filesystem::exists(out_softclips_consd_bam_name) && boost::filesystem::exists(out_softclips_consd_cpos_name)) {
//		return;
//	}

	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_cbam_files_mem_alt");
	checker.start();

	string a_bai_path;
	get_bai_index_path(input_bam_name, a_bai_path);
	string a_bni_path = options.prefix + ".bam.bni";
	if(!options.working_dir.empty()) {
		a_bni_path = options.working_prefix + ".bam.bni";
	}
	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, input_bam_name, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;


//	vector<meerkat::BlockBoundary> actual_blocks;
//	collect_boundaries(actual_blocks, input_bam_name, size_block);

//	BamTools::BamReader local_reader;
//	if (!local_reader.Open(a_path, an_index_path)) {
//		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
//		exit(1);
//	}
//	find_quality_standard(local_reader);
	{
		string line;
		ifstream in(firststat_name, ios::binary);
		for(int64_t e_id = 0; e_id < 5; ++e_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + options.qcutoff)).str();

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;
	boost::unordered_map<string, AlnSeqQualEntry> aln_seq_map;
	{
		castle::TimeChecker sub_checker_set;
		sub_checker_set.setTarget("TEA.generate_cbam_files_mem_alt - collect 'H' set");
		sub_checker_set.start();

		//collect 'H' containing seq ids
		vector<boost::unordered_set<string>> aln_seq_qual_set_lists(calculated_n_blocks - 1);

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				BamTools::BamReader local_reader;
				if (!local_reader.Open(input_bam_name, a_bai_path)) {
					return;
				}
				string str_block_id = boost::lexical_cast<string>(block_id);
				string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
				int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;

				auto& aln_seq_qual_set = aln_seq_qual_set_lists[block_id];
				BamTools::BamAlignment local_alignment_entry;

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;
					if(!local_alignment_entry.IsMapped()) {
						continue;
					}
	//				if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
	//					continue;
	//				}

					auto& seq_id = local_alignment_entry.Name;
					auto& cigar = local_alignment_entry.CigarData;
					// # clipped in the beginning
					if (cigar.size() > 1 && ('H' == cigar.front().Type || 'H' == cigar.back().Type)) {
						aln_seq_qual_set.insert(seq_id);
//						if(local_alignment_entry.IsFirstMate()) {
//							aln_seq_qual_map[seq_id].seq1 = seq;
//							aln_seq_qual_map[seq_id].qual1 = qual;
//						} else if(local_alignment_entry.IsSecondMate()) {
//							aln_seq_qual_map[seq_id].seq2 = seq;
//							aln_seq_qual_map[seq_id].qual2 = qual;
//						}
					}
				}
				local_reader.Close();
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

		cout << sub_checker_set;

		castle::TimeChecker sub_checker_map;
		sub_checker_map.setTarget("TEA.generate_cbam_files_mem_alt - collect 'H' map");
		sub_checker_map.start();

		//collect 'H' containing seq ids
		vector<boost::unordered_map<string, AlnSeqQualEntry>> aln_seq_qual_map_lists(calculated_n_blocks - 1);

		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				BamTools::BamReader local_reader;
				if (!local_reader.Open(input_bam_name, a_bai_path)) {
					return;
				}
				string str_block_id = boost::lexical_cast<string>(block_id);
				string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
				int64_t the_current_ref_offset = actual_blocks[block_id].offset;
				int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
				auto& m_bgzf = local_reader.GetBGZF();
				if(0 != block_id) {
					if(!m_bgzf.Seek(the_current_ref_offset)) {
						local_reader.Close();
						return;
					}
				}
				int64_t cur_offset = m_bgzf.Tell();
				int64_t prev_offset = cur_offset;


				auto& aln_seq_qual_map = aln_seq_qual_map_lists[block_id];
				BamTools::BamAlignment local_alignment_entry;

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
					cur_offset = m_bgzf.Tell();
					if(prev_offset >= the_next_ref_offset) {
						break;
					}
					prev_offset = cur_offset;
					if(!local_alignment_entry.IsMapped()) {
						continue;
					}
	//				if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
	//					continue;
	//				}

					auto& seq_id = local_alignment_entry.Name;
					auto& seq = local_alignment_entry.QueryBases;
					auto& qual = local_alignment_entry.Qualities;
					auto& cigar = local_alignment_entry.CigarData;
					// # clipped in the beginning
					if (cigar.size() > 1 && ('H' == cigar.front().Type || 'H' == cigar.back().Type)) {
						if(local_alignment_entry.IsFirstMate()) {
							if(aln_seq_qual_map[seq_id].seq1.size() < seq.size()) {
								aln_seq_qual_map[seq_id].seq1 = seq;
								aln_seq_qual_map[seq_id].qual1 = qual;
							}
						} else if(local_alignment_entry.IsSecondMate()) {
							if(aln_seq_qual_map[seq_id].seq2.size() < seq.size()) {
								aln_seq_qual_map[seq_id].seq2 = seq;
								aln_seq_qual_map[seq_id].qual2 = qual;
							}
						}
					}
				}
				local_reader.Close();
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& aln_seq_qual_map = aln_seq_qual_map_lists[block_id];
			for(auto& an_entry : aln_seq_qual_map) {
				auto& seq_id = an_entry.first;
				auto& seq_qual_entry = an_entry.second;
				if(aln_seq_map[seq_id].seq1.size() < seq_qual_entry.seq1.size()) {
					aln_seq_map[seq_id].seq1 = seq_qual_entry.seq1;
					aln_seq_map[seq_id].qual1 = seq_qual_entry.qual1;
				}
				if(aln_seq_map[seq_id].seq2.size() < seq_qual_entry.seq2.size()) {
					aln_seq_map[seq_id].seq2 = seq_qual_entry.seq2;
					aln_seq_map[seq_id].qual2 = seq_qual_entry.qual2;
				}
			}
		}

		cout << sub_checker_map;
	}
	{
		castle::TimeChecker sub_checker_mem;
		sub_checker_mem.setTarget("TEA.generate_cbam_files_mem_alt - check memory consumption");
		sub_checker_mem.start();
		cout << (boost::format("[TEA.generate_cbam_files_mem_alt] # stored [seq, qual] entries: %d\n") % aln_seq_map.size()).str();
		cout << sub_checker_mem;
	}
	// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	vector<boost::unordered_map<string, boost::unordered_map<string, RAMEntry>>>h_pos_lists(calculated_n_blocks - 1);
//	vector<boost::unordered_map<string, string>> ram_lists(calculated_n_blocks - 1);
	//	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	//	boost::unordered_map<string, string> ram;

	vector<string> raw_bam_name_lists(calculated_n_blocks - 1);
	vector<string> softclips_consd_cpos_name_lists(calculated_n_blocks - 1);
//	vector<int64_t> cnt_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(input_bam_name, a_bai_path)) {
				return;
			}

//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;

			uint32_t qcutoff = options.qcutoff;
			uint32_t max_mismatches = options.max_mismatches;
			int32_t min_matches = options.min_matches;
			int32_t min_polyAT = options.min_polyAT;

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();
			string cseq1;
			string cseq2;
			uint32_t clen = 0;
			bool polyAT = false;
			uint32_t nm = 0;

			string qual1;
			string qual2;

			bool selected = false;
			bool selected2 = false;
			int32_t cpos;

			string poly_a_str(min_polyAT, 'A');
			string poly_t_str(min_polyAT, 'T');
			string a_delim_str("^");
			for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			vector<string> b;
			const char* delim_all_stopwords = a_delim_str.c_str();

			string local_raw_sam_name = consd_raw_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;
//			bool debug = false;
				BamTools::BamAlignment local_alignment_entry;
				BamTools::BamWriter out_raw_sam;
				if(0 == block_id) {
					out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
				} else {
					out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
				}
				ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first) %s\n")
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
//				debug = "DHFC08P1:370:C1784ACXX:2:1214:21254:59267" == local_alignment_entry.Name;
//				if(!debug) {
//					continue;
//				}

					//	while(<F>) {
					//		chomp; my @a=split(/\t/);
					//		if ($a[1] & 0x400 | $a[1] & 0x4) { next }; # skip PCR duplicate or unmapped reads
//					uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
//					if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
//						continue;
//					}
					if(!local_alignment_entry.IsMapped()) {
						continue;
					}

					try {
						auto& seq_id = local_alignment_entry.Name;
						auto& seq = local_alignment_entry.QueryBases;
						auto& qual = local_alignment_entry.Qualities;
						//		$s1 = $cigar;
						//		$s2 = $cigar;
						polyAT = selected = selected2 = false;
						auto& cigar = local_alignment_entry.CigarData;
						if (cigar.size() > 0 && ('H' == cigar.front().Type || 'H' == cigar.back().Type)) {
							auto a_seq_qual_entry = aln_seq_map.find(seq_id);
							if(aln_seq_map.end() != a_seq_qual_entry) {
								fill_H_to_S(local_alignment_entry, a_seq_qual_entry->second);
							}
						}
						// # clipped in the beginning
						if (!cigar.empty() && ('S' == cigar.front().Type)) {
	//					if(debug) {
	//						cout << "initial clipped\n";
	//					}
							clen = cigar.front().Length;
							if (clen >= 5) {
								qual1 = qual.substr(0, clen);
	//						if(debug) {
	//							cout << "qual1: " << qual1 << "\n";
	//						}
								//# more than 5 good quality bases
								int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
								if (n_good_quals >= 5) {
	//							if(debug) {
	//								cout << "selected 1\n";
	//								exit(0);
	//							}
									selected = true;
	//							if (!(sam_flag & 0x0010)) { // # positive strand mapping
									if (!local_alignment_entry.IsReverseStrand()) { // # positive strand mapping
										selected2 = true;
										cpos = local_alignment_entry.Position;
										//								# check whether the clipped seq has >=10
										cseq1 = seq.substr(0, clen);
										if (string::npos != cseq1.find(poly_a_str)) {
											polyAT = true;
										}
									}
								}
							}
						}
						// # clipping in the end
						clen = 0;
						//				for (auto& a_cigar : cigar) {
						//					if ('D' == a_cigar.Type) {
						//						clen = a_cigar.Length;
						//						break;
						//					}
						//				}
	//				for (auto& a_cigar : cigar) {
	//					if ('D' == a_cigar.Type) {
	//						clen = a_cigar.Length;
	//						break;
	//					}
	//				}
					if(!cigar.empty() && ('S' == cigar.back().Type)) {
						clen = cigar.back().Length;
					}
					if (clen >= 5) {
						int64_t the_qual_sub_str_pos = qual.size() - clen;
						if (the_qual_sub_str_pos > 0) {
							qual2 = qual.substr(the_qual_sub_str_pos);
						} else {
							qual2 = qual;
						}
	//					if(debug) {
	//						cout << "qual2: " << qual2 << "\n";
	//					}
						int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
						//		if ($qual2 =~ m/([^#]+)/) {
						//				if (length($1) >= 5 ) {
						if (n_good_quals >= 5) {
	//						if(debug) {
	//							cout << "selected: true\n";
	//						}
							selected = true;
	//						if (sam_flag & 0x0010) { // # negative strand mapping
							if (local_alignment_entry.IsReverseStrand()) {
								selected2 = true;
								cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
								//							# check whether the clipped seq has >=10
								int64_t the_seq_sub_str_pos = seq.size() - clen;
								cseq2 = seq.substr(the_seq_sub_str_pos);
								//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
								if (string::npos != cseq2.find(poly_t_str)) {
									polyAT = true;
								}
							}
						}
					}

					//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
					if (!polyAT) {
						//        # check max_mismatches
						local_alignment_entry.GetEditDistance(nm);
						if (nm > max_mismatches) {
							continue;
						}

						string md_str;
						local_alignment_entry.GetTag("MD", md_str);
						//        $md = $_;

						//        # check min_matches is satisfied

						//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
						//        my @b=split(/[A-Za-z^]/, $md);
						//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
						castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
						int64_t summed_n_bases = 0;
						for (auto a_value_str : b) {
							summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
						}
						if (summed_n_bases < min_matches) {
							continue;
						}
					}

					if (selected2) {
						out_raw_sam.SaveSAMAlignment(local_alignment_entry);
						out_consd_cpos << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
					}
				} catch(exception& ex) {
					cout << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
					cout << ex.what() << "\n";
				}
			}
			local_reader.Close();
			out_raw_sam.Close();
			//				done_vector[block_id] = 'D';

			//				if(verbose) {
			//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
			//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
			//							% local_alignment_entry.AlignmentFlag).str();
			//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
			//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
			//								% block_id % a_block_boundary_str).str();
			//					} else {
			//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
			//								% block_id % a_block_boundary_str).str();
			//					}
			//				} else {
			//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
			//					double processed = n/(double)done_vector.size() * 100.0;
			//					cout << (boost::format("%.2f %%\n") % processed).str();
			//				}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	//	cout << "[TEA.load_repeat_mapping] gather scattered information\n";

	castle::IOUtils::plain_file_merge(consd_raw_sam_name, raw_bam_name_lists, n_cores, true);
	castle::IOUtils::plain_file_merge(out_softclips_consd_cpos_name, softclips_consd_cpos_name_lists, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_raw_bam_name % consd_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt] done generating cbam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt] start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();

	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_softclips_consd_bam_name % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());
	cout << "[TEA.generate_cbam_files_mem_alt] done generating cbam and its index.\n";
	cout << checker;
}

void TEA::_generate_cbam_files_mem_alt2() {
	string original_input_bam_name = options.prefix + ".bam";
	string the_prefix = options.prefix;
	if (!options.working_dir.empty()) {
		the_prefix = options.working_prefix;
	}
	string firststat_name = the_prefix + ".firststat";
	string input_bam_name = the_prefix + ".disc.sorted.bam";
	string input_bam_bfi_name = the_prefix + ".disc.sorted.bam.bfi";
	string consd_raw_sam_name = the_prefix + ".softclips.consd.raw.sam";
	string consd_raw_bam_name = the_prefix + ".softclips.consd.raw.bam";
	string consd_sorted_bam_name = the_prefix + ".softclips.consd.sorted.bam";
	string consd_sorted_sam_name = the_prefix + ".softclips.consd.sorted.sam";
	string consd_sorted_bam_bni_name = the_prefix + ".softclips.consd.sorted.bam.bni";
	string consd_sorted_raw_sam_name = the_prefix + ".softclips.consd.sorted.raw.sam";
	string consd_sorted_raw_bam_name = the_prefix + ".softclips.consd.sorted.raw.bam";
	string out_softclips_consd_bam_name = the_prefix + ".softclips.consd.bam";
	string out_softclips_consd_cpos_name = the_prefix + ".softclips.consd.cpos";

//	if (castle::IOUtils::get_file_size(out_softclips_consd_bam_name) > 0 && castle::IOUtils::get_file_size(out_softclips_consd_cpos_name) > 0) {
//		return;
//	}

	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_cbam_files_mem_alt2");
	checker.start();

	vector<int64_t> block_boundary;
	string an_original_index_path;
	get_bai_index_path(original_input_bam_name, an_original_index_path);

	if(!boost::filesystem::exists(original_input_bam_name)) {
		cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] cannot read: %s\n") % original_input_bam_name).str();
		cout << checker;
		return;
	}

	{
		string line;
		ifstream in(firststat_name, ios::binary);
		for(int64_t e_id = 0; e_id < 5; ++e_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}

	// create a .softclips.consd.raw.bam file, which contains a BAM having converted of CIGAR string 'H' -> 'S'.

	if (boost::filesystem::exists(input_bam_bfi_name)) {
		castle::ParallelRunner::load_bfi_index(block_boundary, input_bam_bfi_name);
	} else {
		const int64_t N_BLOCK_ENTRIES = 262144;
		castle::ParallelRunner::create_bfi_index(block_boundary, input_bam_name, input_bam_bfi_name, N_BLOCK_ENTRIES);
	}

	vector<function<void()> > tasks;

	int64_t n_block_boundaries = block_boundary.size();
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] # blocks: %d of .disc.sorted.bam\n") % (n_block_boundaries - 1)).str();
	vector<string> raw_bam_name_lists(n_block_boundaries - 1);
	vector<string> softclips_consd_cpos_name_lists(n_block_boundaries - 1);

	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(input_bam_name)) {
				return;
			}
			const RefVector& m_ref = local_reader.GetReferenceData();
			boost::unordered_map<string, int64_t> ref_reverse_index;
			for (uint64_t ref_id = 0; ref_id < m_ref.size(); ++ref_id) {
				auto& a_ref = m_ref[ref_id];
				ref_reverse_index[a_ref.RefName] = ref_id;
			}
			auto& local_data = local_reader.GetBGZF();
			if(0 != block_boundary[block_id]) {
				local_data.Seek(block_boundary[block_id]);
			}
			int64_t next_block_pos = block_boundary[block_id + 1];
			int64_t cur_pos = -1;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();

			string local_raw_sam_name = consd_raw_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;
//			bool debug = false;
			BamAlignment local_alignment_entry;
			BamWriter out_raw_sam;
			if(0 == block_id) {
				out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
			} else {
				out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
			}

//			ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);
			vector<BamAlignment> pair1_alns;
			vector<BamAlignment> pair2_alns;
			string prev_name;

			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				auto& cur_name = local_alignment_entry.Name;
				if (prev_name != cur_name) {
					if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
						remove_entry_enclosed_with_large_H(pair1_alns);
						remove_entry_enclosed_with_large_H(pair2_alns);
						convert_H_to_S(pair1_alns);
						convert_H_to_S(pair2_alns);
						for(auto& an_aln : pair1_alns) {
							if(has_S_or_H(an_aln)) {
								out_raw_sam.SaveSAMAlignment(an_aln);
							}
						}
						for(auto& an_aln : pair2_alns) {
							if(has_S_or_H(an_aln)) {
								out_raw_sam.SaveSAMAlignment(an_aln);
							}
						}
						pair1_alns.clear();
						pair2_alns.clear();
					}
					prev_name = cur_name;
				}
				if(local_alignment_entry.IsFirstMate()) {
					pair1_alns.push_back(local_alignment_entry);
				} else {
					pair2_alns.push_back(local_alignment_entry);
				}
				cur_pos = local_data.Tell();
				if(cur_pos >= next_block_pos) {
					break;
				}
			}
			if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
				remove_entry_enclosed_with_large_H(pair1_alns);
				remove_entry_enclosed_with_large_H(pair2_alns);
				convert_H_to_S(pair1_alns);
				convert_H_to_S(pair2_alns);
				for(auto& an_aln : pair1_alns) {
					if(has_S_or_H(an_aln)) {
						out_raw_sam.SaveSAMAlignment(an_aln);
					}
				}
				for(auto& an_aln : pair2_alns) {
					if(has_S_or_H(an_aln)) {
						out_raw_sam.SaveSAMAlignment(an_aln);
					}
				}
				pair1_alns.clear();
				pair2_alns.clear();
			}
			local_reader.Close();
			out_raw_sam.Close();
		});
	}

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + options.qcutoff)).str();
	castle::IOUtils::plain_file_merge(consd_raw_sam_name, raw_bam_name_lists, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_raw_bam_name % consd_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] done generating consd raw bam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % consd_sorted_bam_name % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());

	// create .cpos .softclips.consd.bam files

	int64_t size_block = 8192000;
	vector<int64_t> unmapped_block_boundaries;
	castle::ParallelRunner::load_bai_bni_index(unmapped_block_boundaries, consd_sorted_bam_name, consd_sorted_bam_bni_name, size_block, n_cores);

	int64_t calculated_n_blocks = unmapped_block_boundaries.size();
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] # blocks: %d of .consd.sorted.bam\n") % (calculated_n_blocks - 1)).str();
	raw_bam_name_lists.clear();
	softclips_consd_cpos_name_lists.clear();
	raw_bam_name_lists.resize(calculated_n_blocks - 1);
	softclips_consd_cpos_name_lists.resize(calculated_n_blocks - 1);
	boost::mutex print_mutex;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamReader local_reader;
			if (!local_reader.Open(consd_sorted_bam_name, consd_sorted_bam_name + ".bai")) {
				return;
			}
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t the_current_ref_offset = unmapped_block_boundaries[block_id];
			int64_t the_next_ref_offset = unmapped_block_boundaries[block_id + 1];
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;

			uint32_t qcutoff = options.qcutoff;
			uint32_t max_mismatches = options.max_mismatches;
			int32_t min_matches = options.min_matches;
			int32_t min_polyAT = options.min_polyAT;

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();
			string cseq1;
			string cseq2;
			uint32_t clen = 0;
			bool polyAT = false;
			uint32_t nm = 0;

			string qual1;
			string qual2;

			bool selected = false;
			bool selected2 = false;
			int32_t cpos;

			string poly_a_str(min_polyAT, 'A');
			string poly_t_str(min_polyAT, 'T');
			string a_delim_str("^");
			for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			vector<string> b;
			const char* delim_all_stopwords = a_delim_str.c_str();

			string local_raw_sam_name = consd_sorted_raw_sam_name + "." + str_block_id;
			string local_sorted_sam_name = consd_sorted_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;

			BamTools::BamAlignment local_alignment_entry;
			BamTools::BamWriter out_raw_sam;
//			BamTools::BamWriter out_sorted_sam;
			if(0 == block_id) {
				out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
//				out_sorted_sam.SAMOpen(local_sorted_sam_name, header, ref_vec);
			} else {
				out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
//				out_sorted_sam.SAMOpenNoHeader(local_sorted_sam_name, ref_vec);
			}
//			bool debug = (0 == block_id);
			const bool debug = false;
//			const bool debug = true;
			ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);

			while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
				cur_offset = m_bgzf.Tell();
				if(prev_offset >= the_next_ref_offset) {
					break;
				}
				prev_offset = cur_offset;
				++num_total;
				uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
				if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
					continue;
				}

				auto& seq = local_alignment_entry.QueryBases;
				auto& qual = local_alignment_entry.Qualities;
				//		$s1 = $cigar;
				//		$s2 = $cigar;
				polyAT = selected = selected2 = false;
				auto& cigar = local_alignment_entry.CigarData;
				// # clipped in the beginning
				try {
				if (!cigar.empty() && ('S' == cigar.front().Type)) {
					if(debug) {
						boost::lock_guard<boost::mutex> a_lock(print_mutex);
						cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] initial clipped: %s\n") % BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec)).str();
					}
					clen = cigar.front().Length;
					if (clen >= 5) {
						qual1 = qual.substr(0, clen);
						if(debug) {
							boost::lock_guard<boost::mutex> a_lock(print_mutex);
							cout << "[TEA.generate_cbam_files_mem_alt2] qual1: " << qual1 << "\n";
						}
						//# more than 5 good quality bases
						int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
						if (n_good_quals >= 5) {
							if(debug) {
								boost::lock_guard<boost::mutex> a_lock(print_mutex);
								cout << "[TEA.generate_cbam_files_mem_alt2] selected 1-1\n";
							}
							selected = true;
//							if (!(sam_flag & 0x0010)) { // # positive strand mapping
							if (!local_alignment_entry.IsReverseStrand()) { // # positive strand mapping
								if(debug) {
									boost::lock_guard<boost::mutex> a_lock(print_mutex);
									cout << "[TEA.generate_cbam_files_mem_alt2] selected 2-1\n";
								}
								selected2 = true;
								cpos = local_alignment_entry.Position;
								//								# check whether the clipped seq has >=10
								cseq1 = seq.substr(0, clen);
								if (string::npos != cseq1.find(poly_a_str)) {
									if(debug) {
										boost::lock_guard<boost::mutex> a_lock(print_mutex);
										cout << "[TEA.generate_cbam_files_mem_alt2] polyAT 1\n";
									}
									polyAT = true;
								}
							}
						}
					}
				}
				// # clipping in the end
				clen = 0;
//				for (auto& a_cigar : cigar) {
//					if ('D' == a_cigar.Type) {
//						clen = a_cigar.Length;
//						break;
//					}
//				}
				if(!cigar.empty() && ('S' == cigar.back().Type)) {
					clen = cigar.back().Length;
				}
				if (clen >= 5) {
					int64_t the_qual_sub_str_pos = qual.size() - clen;
					if (the_qual_sub_str_pos > 0) {
						qual2 = qual.substr(the_qual_sub_str_pos);
					} else {
						qual2 = qual;
					}
					if(debug) {
						boost::lock_guard<boost::mutex> a_lock(print_mutex);
						cout << "[TEA.generate_cbam_files_mem_alt2] qual2: " << qual2 << "\n";
					}
					int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
					//		if ($qual2 =~ m/([^#]+)/) {
					//				if (length($1) >= 5 ) {
					if (n_good_quals >= 5) {
						if(debug) {
							boost::lock_guard<boost::mutex> a_lock(print_mutex);
							cout << "[TEA.generate_cbam_files_mem_alt2] selected 1-2\n";
						}
						selected = true;
//						if (sam_flag & 0x0010) { // # negative strand mapping
						if (local_alignment_entry.IsReverseStrand()) {
							if(debug) {
								boost::lock_guard<boost::mutex> a_lock(print_mutex);
								cout << "[TEA.generate_cbam_files_mem_alt2] selected 2-2\n";
							}
							selected2 = true;
							cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
							//							# check whether the clipped seq has >=10
							int64_t the_seq_sub_str_pos = seq.size() - clen;
							cseq2 = seq.substr(the_seq_sub_str_pos);
							//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
							if (string::npos != cseq2.find(poly_t_str)) {
								if(debug) {
									boost::lock_guard<boost::mutex> a_lock(print_mutex);
									cout << "[TEA.generate_cbam_files_mem_alt2] polyAT 2\n";
								}
								polyAT = true;
							}
						}
					}
				}

				//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
				if (!polyAT) {
					//        # check max_mismatches
					local_alignment_entry.GetEditDistance(nm);
					if (nm > max_mismatches) {
						continue;
					};

					string md_str;
					local_alignment_entry.GetTag("MD", md_str);
					//        $md = $_;

					//        # check min_matches is satisfied

					//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
					//        my @b=split(/[A-Za-z^]/, $md);
					//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
					castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
					int64_t summed_n_bases = 0;
					for (auto a_value_str : b) {
						summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
					}
					if (summed_n_bases < min_matches) {
						continue;
					}
				}

				if (selected2) {
					out_raw_sam.SaveSAMAlignment(local_alignment_entry);
					out_consd_cpos << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
				}
				} catch(exception& ex) {
					cout << ex.what() << "\n";
					cout << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
				}

			}
			local_reader.Close();
			out_raw_sam.Close();
		});
	}

	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	castle::IOUtils::plain_file_merge(consd_sorted_raw_sam_name, raw_bam_name_lists, n_cores, true);
	string sam_to_sorted_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_sorted_raw_bam_name % consd_sorted_raw_sam_name).str();
	system(sam_to_sorted_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] done generating consd sorted raw bam: %s\n") % consd_sorted_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem_alt2] start sorting and generating the index for %s ...\n") % consd_sorted_raw_bam_name).str();
	string sambamba_last_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_softclips_consd_bam_name % consd_sorted_raw_bam_name).str();
	system(sambamba_last_sort_cmd.c_str());
	castle::IOUtils::plain_file_merge(out_softclips_consd_cpos_name, softclips_consd_cpos_name_lists, n_cores, true);
	cout << "[TEA.generate_cbam_files_mem_alt2] done generating cbam and its index.\n";

	cout << checker;
}

void TEA::_generate_cbam_files_mem_org() {
	string input_bam_name = options.prefix + ".bam";
	string firststat_name = options.prefix + ".firststat";
	string consd_raw_sam_name = options.prefix + ".softclips.consd.raw.sam";
	string consd_raw_bam_name = options.prefix + ".softclips.consd.raw.bam";
	string out_softclips_consd_bam_name = options.prefix + ".softclips.consd.bam";
	string out_softclips_consd_cpos_name = options.prefix + ".softclips.consd.cpos";
	if (!options.working_dir.empty()) {
		firststat_name = options.working_prefix + ".firststat";
		consd_raw_sam_name = options.working_prefix + ".softclips.consd.raw.sam";
		consd_raw_bam_name = options.working_prefix + ".softclips.consd.raw.bam";
		out_softclips_consd_bam_name = options.working_prefix + ".softclips.consd.bam";
		out_softclips_consd_cpos_name = options.working_prefix + ".softclips.consd.cpos";
	}

//	if (!options.is_force && boost::filesystem::exists(out_softclips_consd_bam_name) && boost::filesystem::exists(out_softclips_consd_cpos_name)) {
//		return;
//	}

	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_cbam_files_mem_org");
	checker.start();

	string a_bai_path;
	get_bai_index_path(input_bam_name, a_bai_path);
	string a_bni_path = options.prefix + ".bam.bni";
	if(!options.working_dir.empty()) {
		a_bni_path = options.working_prefix + ".bam.bni";
	}
	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, input_bam_name, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	vector<meerkat::BlockBoundary> actual_blocks;
//	collect_boundaries(actual_blocks, input_bam_name, size_block);

//	BamTools::BamReader local_reader;
//	if (!local_reader.Open(a_path, an_index_path)) {
//		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
//		exit(1);
//	}
//	find_quality_standard(local_reader);
	{
		string line;
		ifstream in(firststat_name, ios::binary);
		for(int64_t e_id = 0; e_id < 5; ++e_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}
	cout << (boost::format("[TEA._generate_cbam_files_mem_org] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + options.qcutoff)).str();

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

	// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	vector<boost::unordered_map<string, boost::unordered_map<string, RAMEntry>>>h_pos_lists(calculated_n_blocks - 1);
//	vector<boost::unordered_map<string, string>> ram_lists(calculated_n_blocks - 1);
	//	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	//	boost::unordered_map<string, string> ram;

	vector<string> raw_bam_name_lists(calculated_n_blocks - 1);
	vector<string> softclips_consd_cpos_name_lists(calculated_n_blocks - 1);
//	vector<int64_t> cnt_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(input_bam_name, a_bai_path)) {
				return;
			}

//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;

			uint32_t qcutoff = options.qcutoff;
			uint32_t max_mismatches = options.max_mismatches;
			int32_t min_matches = options.min_matches;
			int32_t min_polyAT = options.min_polyAT;

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();
			string cseq1;
			string cseq2;
			uint32_t clen = 0;
			bool polyAT = false;
			uint32_t nm = 0;

			string qual1;
			string qual2;

			bool selected = false;
			bool selected2 = false;
			int32_t cpos;

			string poly_a_str(min_polyAT, 'A');
			string poly_t_str(min_polyAT, 'T');
			string a_delim_str("^");
			for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			vector<string> b;
			const char* delim_all_stopwords = a_delim_str.c_str();

			string local_raw_sam_name = consd_raw_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;
//			bool debug = false;
				BamTools::BamAlignment local_alignment_entry;
				BamTools::BamWriter out_raw_sam;
				if(0 == block_id) {
					out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
				} else {
					out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
				}
				ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first) %s\n")
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
//				debug = "DHFC08P1:370:C1784ACXX:2:1214:21254:59267" == local_alignment_entry.Name;
//				if(!debug) {
//					continue;
//				}

					//	while(<F>) {
					//		chomp; my @a=split(/\t/);
					//		if ($a[1] & 0x400 | $a[1] & 0x4) { next }; # skip PCR duplicate or unmapped reads
					uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
					if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
						continue;
					}

					auto& seq = local_alignment_entry.QueryBases;
					auto& qual = local_alignment_entry.Qualities;
					//		$s1 = $cigar;
					//		$s2 = $cigar;
					polyAT = selected = selected2 = false;
					auto& cigar = local_alignment_entry.CigarData;
					// # clipped in the beginning
					if (!cigar.empty() && ('S' == cigar.front().Type)) {
//					if(debug) {
//						cout << "initial clipped\n";
//					}
						clen = cigar.front().Length;
						if (clen >= 5) {
							qual1 = qual.substr(0, clen);
//						if(debug) {
//							cout << "qual1: " << qual1 << "\n";
//						}
							//# more than 5 good quality bases
							int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
							if (n_good_quals >= 5) {
//							if(debug) {
//								cout << "selected 1\n";
//								exit(0);
//							}
								selected = true;
//							if (!(sam_flag & 0x0010)) { // # positive strand mapping
								if (!local_alignment_entry.IsReverseStrand()) { // # positive strand mapping
									selected2 = true;
									cpos = local_alignment_entry.Position;
									//								# check whether the clipped seq has >=10
									cseq1 = seq.substr(0, clen);
									if (string::npos != cseq1.find(poly_a_str)) {
										polyAT = true;
									}
								}
							}
						}
					}
					// # clipping in the end
					clen = 0;
//				for (auto& a_cigar : cigar) {
//					if ('D' == a_cigar.Type) {
//						clen = a_cigar.Length;
//						break;
//					}
//				}
					if(!cigar.empty() && ('S' == cigar.back().Type)) {
						clen = cigar.back().Length;
					}
					if (clen >= 5) {
						int64_t the_qual_sub_str_pos = qual.size() - clen;
						if (the_qual_sub_str_pos > 0) {
							qual2 = qual.substr(the_qual_sub_str_pos);
						} else {
							qual2 = qual;
						}
//					if(debug) {
//						cout << "qual2: " << qual2 << "\n";
//					}
						int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
						//		if ($qual2 =~ m/([^#]+)/) {
						//				if (length($1) >= 5 ) {
						if (n_good_quals >= 5) {
//						if(debug) {
//							cout << "selected: true\n";
//						}
							selected = true;
//						if (sam_flag & 0x0010) { // # negative strand mapping
							if (local_alignment_entry.IsReverseStrand()) {
								selected2 = true;
								cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
								//							# check whether the clipped seq has >=10
								int64_t the_seq_sub_str_pos = seq.size() - clen;
								cseq2 = seq.substr(the_seq_sub_str_pos);
								//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
								if (string::npos != cseq2.find(poly_t_str)) {
									polyAT = true;
								}
							}
						}
					}

					//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
					if (!polyAT) {
						//        # check max_mismatches
						local_alignment_entry.GetEditDistance(nm);
						if (nm > max_mismatches) {
							continue;
						};

						string md_str;
						local_alignment_entry.GetTag("MD", md_str);
						//        $md = $_;

						//        # check min_matches is satisfied

						//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
						//        my @b=split(/[A-Za-z^]/, $md);
						//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
						castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
						int64_t summed_n_bases = 0;
						for (auto a_value_str : b) {
							summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
						}
						if (summed_n_bases < min_matches) {
							continue;
						}
					}

					if (selected2) {
						out_raw_sam.SaveSAMAlignment(local_alignment_entry);
						out_consd_cpos << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
					}

				}
				local_reader.Close();
				out_raw_sam.Close();
				//				done_vector[block_id] = 'D';

				//				if(verbose) {
				//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
				//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
				//							% local_alignment_entry.AlignmentFlag).str();
				//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					} else {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					}
				//				} else {
				//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
				//					double processed = n/(double)done_vector.size() * 100.0;
				//					cout << (boost::format("%.2f %%\n") % processed).str();
				//				}
			});
	}
	// sometimes the unaligned reads are in the last portion of BAM file, hence
	// changing the order.
//	if (tasks.size() > 0) {
//		swap(tasks[0], tasks.back());
//	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	//	cout << "[TEA.load_repeat_mapping] gather scattered information\n";

	castle::IOUtils::plain_file_merge(consd_raw_sam_name, raw_bam_name_lists, n_cores, true);
	castle::IOUtils::plain_file_merge(out_softclips_consd_cpos_name, softclips_consd_cpos_name_lists, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_raw_bam_name % consd_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_mem_org] done generating cbam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_mem_org] start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();

	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_softclips_consd_bam_name % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());
	cout << "[TEA.generate_cbam_files_mem_org] done generating cbam and its index.\n";
	cout << checker;

}

void TEA::_generate_cbam_files_sampe() {
	string input_bam_name = options.prefix + ".bam";
	string firststat_name = options.prefix + ".firststat";
	string consd_raw_sam_name = options.prefix + ".softclips.consd.raw.sam";
	string consd_raw_bam_name = options.prefix + ".softclips.consd.raw.bam";
	string out_softclips_consd_bam_name = options.prefix + ".softclips.consd.bam";
	string out_softclips_consd_cpos_name = options.prefix + ".softclips.consd.cpos";
	if (!options.working_dir.empty()) {
		firststat_name = options.working_prefix + ".firststat";
		consd_raw_sam_name = options.working_prefix + ".softclips.consd.raw.sam";
		consd_raw_bam_name = options.working_prefix + ".softclips.consd.raw.bam";
		out_softclips_consd_bam_name = options.working_prefix + ".softclips.consd.bam";
		out_softclips_consd_cpos_name = options.working_prefix + ".softclips.consd.cpos";
	}

	if (!options.is_force && boost::filesystem::exists(out_softclips_consd_bam_name) && boost::filesystem::exists(out_softclips_consd_cpos_name)) {
		return;
	}

	castle::TimeChecker checker;
	checker.setTarget("TEA.generate_cbam_files_sampe");
	checker.start();

	string a_bai_path;
	get_bai_index_path(input_bam_name, a_bai_path);
	string a_bni_path = options.prefix + ".bam.bni";
	if(!options.working_dir.empty()) {
		a_bni_path = options.working_prefix + ".bam.bni";
	}
	int64_t size_block = 8192000;
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, input_bam_name, a_bai_path, a_bni_path, size_block);
	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	vector<meerkat::BlockBoundary> actual_blocks;
//	collect_boundaries(actual_blocks, input_bam_name, size_block);

//	BamTools::BamReader local_reader;
//	if (!local_reader.Open(a_path, an_index_path)) {
//		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
//		exit(1);
//	}
//	find_quality_standard(local_reader);
	{
		string line;
		ifstream in(firststat_name, ios::binary);
		for(int64_t e_id = 0; e_id < 5; ++e_id) {
			getline(in, line, '\n');
		}
		qenc = boost::lexical_cast<int>(line);
		min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	}
	cout << (boost::format("[TEA.generate_cbam_files_sampe] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + options.qcutoff)).str();

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

	// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	vector<boost::unordered_map<string, boost::unordered_map<string, RAMEntry>>>h_pos_lists(calculated_n_blocks - 1);
//	vector<boost::unordered_map<string, string>> ram_lists(calculated_n_blocks - 1);
	//	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	//	boost::unordered_map<string, string> ram;

	vector<string> raw_bam_name_lists(calculated_n_blocks - 1);
	vector<string> softclips_consd_cpos_name_lists(calculated_n_blocks - 1);
//	vector<int64_t> cnt_lists(calculated_n_blocks - 1);
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(input_bam_name, a_bai_path)) {
				return;
			}

//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA.write_ram_and_bam] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;
			int64_t num_total = 0;

			uint32_t qcutoff = options.qcutoff;
			uint32_t max_mismatches = options.max_mismatches;
			int32_t min_matches = options.min_matches;
			int32_t min_polyAT = options.min_polyAT;

			string previous_pos;
			string current_pos;

			const auto& header = local_reader.GetHeaderText();
			auto& ref_vec = local_reader.GetReferenceData();
			string cseq1;
			string cseq2;
			uint32_t clen = 0;
			bool polyAT = false;
			uint32_t nm = 0;

			string qual1;
			string qual2;

			bool selected = false;
			bool selected2 = false;
			int32_t cpos;

			string poly_a_str(min_polyAT, 'A');
			string poly_t_str(min_polyAT, 'T');
			string a_delim_str("^");
			for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
				a_delim_str += static_cast<char>(c);
			}
			vector<string> b;
			const char* delim_all_stopwords = a_delim_str.c_str();

			string local_raw_sam_name = consd_raw_sam_name + "." + str_block_id;
			string local_softclips_consd_cpos_name = out_softclips_consd_cpos_name + "." + str_block_id;

			raw_bam_name_lists[block_id] = local_raw_sam_name;
			softclips_consd_cpos_name_lists[block_id] = local_softclips_consd_cpos_name;
			const bool debug = false;
//			bool debug = false;
				BamTools::BamAlignment local_alignment_entry;
				BamTools::BamWriter out_raw_sam;
				if(0 == block_id) {
					out_raw_sam.SAMOpen(local_raw_sam_name, header, ref_vec);
				} else {
					out_raw_sam.SAMOpenNoHeader(local_raw_sam_name, ref_vec);
				}
				ofstream out_consd_cpos(local_softclips_consd_cpos_name, ios::binary);

				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA.write_ram_and_bam] Block-%d (first) %s\n")
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
//				debug = ("HSQ700642:192:C13FVACXX:1:2101:19274:55222" == local_alignment_entry.Name || "HSQ700642:188:C0U39ACXX:8:1101:7056:89344" == local_alignment_entry.Name);
//					debug = ("HSQ700642:192:C13FVACXX:1:2301:15557:8813" == local_alignment_entry.Name);
//				if(!debug) {
//					continue;
//				}

					//	while(<F>) {
					//		chomp; my @a=split(/\t/);
					//		if ($a[1] & 0x400 | $a[1] & 0x4) { next }; # skip PCR duplicate or unmapped reads
					uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
					if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
						continue;
					}

					auto& seq = local_alignment_entry.QueryBases;
					auto& qual = local_alignment_entry.Qualities;
					//		$s1 = $cigar;
					//		$s2 = $cigar;
					polyAT = selected = selected2 = false;
					auto& cigar = local_alignment_entry.CigarData;
					// # clipped in the beginning
					if (!cigar.empty() && ('S' == cigar.front().Type)) {
					if(debug) {
						cout << "[TEA.generate_cbam_files_sampe] initial clipped: " << block_id << "\n";
					}
						clen = cigar.front().Length;
						if (clen >= 5) {
							qual1 = qual.substr(0, clen);
						if(debug) {
							cout << "[TEA.generate_cbam_files_sampe] qual1: " << qual1 << "\n";
						}
							//# more than 5 good quality bases
							int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
							if (n_good_quals >= 5) {
							if(debug) {
								cout << "[TEA.generate_cbam_files_sampe] selected 1-1: qual > 5\n";
							}
								selected = true;
//							if (!(sam_flag & 0x0010)) { // # positive strand mapping
								if (!local_alignment_entry.IsReverseStrand()) { // # positive strand mapping
									selected2 = true;
									cpos = local_alignment_entry.Position;
									//								# check whether the clipped seq has >=10
									cseq1 = seq.substr(0, clen);
									if (string::npos != cseq1.find(poly_a_str)) {
										if(debug) {
											cout << "[TEA.generate_cbam_files_sampe] selected 1-2: poly a\n";
										}
										polyAT = true;
									}
								}
							}
						}
					}
					// # clipping in the end
					clen = 0;
//				for (auto& a_cigar : cigar) {
//					if ('D' == a_cigar.Type) {
//						clen = a_cigar.Length;
//						break;
//					}
//				}
					if(!cigar.empty() && ('S' == cigar.back().Type)) {
						clen = cigar.back().Length;
					}
					if (clen >= 5) {
						int64_t the_qual_sub_str_pos = qual.size() - clen;
						if (the_qual_sub_str_pos > 0) {
							qual2 = qual.substr(the_qual_sub_str_pos);
						} else {
							qual2 = qual;
						}
					if(debug) {
						cout << "[TEA.generate_cbam_files_sampe] qual2: " << qual2 << "\n";
					}
						int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
						//		if ($qual2 =~ m/([^#]+)/) {
						//				if (length($1) >= 5 ) {
						if (n_good_quals >= 5) {
						if(debug) {
							cout << "[TEA.generate_cbam_files_sampe] selected 2-1: qual\n";
						}
							selected = true;
//						if (sam_flag & 0x0010) { // # negative strand mapping
							if (local_alignment_entry.IsReverseStrand()) {
								if(debug) {
									cout << "[TEA.generate_cbam_files_sampe] selected 2-2: qual\n";
								}
								selected2 = true;
								cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
								//							# check whether the clipped seq has >=10
								int64_t the_seq_sub_str_pos = seq.size() - clen;
								cseq2 = seq.substr(the_seq_sub_str_pos);
								//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
								if (string::npos != cseq2.find(poly_t_str)) {
									if(debug) {
										cout << "[TEA.generate_cbam_files_sampe] selected 2-3: poly a\n";
									}
									polyAT = true;
								}
							}
						}
					}

					//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
					if (!polyAT) {
						//        # check max_mismatches
						local_alignment_entry.GetEditDistance(nm);
						if (nm > max_mismatches) {
							if(debug) {
								cout << "[TEA.generate_cbam_files_sampe] wrong mismatches\n";
							}
							continue;
						};

						string md_str;
						local_alignment_entry.GetTag("MD", md_str);
						//        $md = $_;

						//        # check min_matches is satisfied

						//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
						//        my @b=split(/[A-Za-z^]/, $md);
						//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
						castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
						int64_t summed_n_bases = 0;
						for (auto a_value_str : b) {
							summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
						}
						if (summed_n_bases < min_matches) {
							if(debug) {
								cout << "[TEA.generate_cbam_files_sampe] MD: " << summed_n_bases << "\n";
							}
							continue;
						}
					}

					if (selected2) {
						if(debug) {
							cout << "[TEA.generate_cbam_files_sampe] selected 3: write\n";
						}
						out_raw_sam.SaveSAMAlignment(local_alignment_entry);
						out_consd_cpos << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
					}

				}
				local_reader.Close();
				out_raw_sam.Close();
				//				done_vector[block_id] = 'D';

				//				if(verbose) {
				//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
				//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
				//							% local_alignment_entry.AlignmentFlag).str();
				//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					} else {
				//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
				//								% block_id % a_block_boundary_str).str();
				//					}
				//				} else {
				//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
				//					double processed = n/(double)done_vector.size() * 100.0;
				//					cout << (boost::format("%.2f %%\n") % processed).str();
				//				}
			});
	}
	// sometimes the unaligned reads are in the last portion of BAM file, hence
	// changing the order.
	if (tasks.size() > 0) {
		swap(tasks[0], tasks.back());
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

	//	cout << "[TEA.load_repeat_mapping] gather scattered information\n";

	castle::IOUtils::plain_file_merge(consd_raw_sam_name, raw_bam_name_lists, n_cores, true);
	castle::IOUtils::plain_file_merge(out_softclips_consd_cpos_name, softclips_consd_cpos_name_lists, n_cores, true);
	string sam_to_bam_cmd = (boost::format("samtools view -@ %d -o %s -Sb %s") % n_cores % consd_raw_bam_name % consd_raw_sam_name).str();
	system(sam_to_bam_cmd.c_str());
	cout << (boost::format("[TEA.generate_cbam_files_sampe] done generating cbam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("[TEA.generate_cbam_files_sampe] start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();

	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_softclips_consd_bam_name % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());
	cout << "[TEA.generate_cbam_files_sampe] done generating cbam and its index.\n";
	cout << checker;
}

void TEA::generate_cbam_files_serial() {
	uint32_t qcutoff = options.qcutoff;
	uint32_t max_mismatches = options.max_mismatches;
	int32_t min_matches = options.min_matches;
	int32_t min_polyAT = options.min_polyAT;
	string bamf = options.prefix + ".bam";
	string consd_raw_bam_name = options.prefix + ".softclips.consd.raw.bam";
	string out_bam = options.prefix + ".softclips.consd.bam";
	string outf3 = options.prefix + ".softclips.consd.cpos.bz2";
	if (!options.working_dir.empty()) {
		bamf = options.working_prefix + ".bam";
		consd_raw_bam_name = options.working_prefix + ".softclips.consd.raw.bam";
		out_bam = options.working_prefix + ".softclips.consd.bam";
		outf3 = options.working_prefix + ".softclips.consd.cpos.bz2";
	}

	string a_path(bamf);
	string an_index_path(a_path);
	an_index_path += ".bai";

	boost::unordered_map<string, boost::unordered_map<string, RAMEntry>> h_pos;
	boost::unordered_map<string, string> ram;
//	int64_t cnt = 0;
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		cout << "ERROR: could not open BAM file '" << a_path << "'\n";
		exit(1);
	}
	find_quality_standard(local_reader);

	cout << (boost::format("[TEA.generate_cbam_files] Actual Qual Cutoff: %c\n") % static_cast<char>(min_qual + qcutoff)).str();

	BamTools::BamWriter O2;
	O2.SAMOpen(consd_raw_bam_name, local_reader.GetHeaderText(), local_reader.GetReferenceData());
	ofstream O3(outf3, ios::binary);
//	open(O2, "| $samtools view -bS - > $outf2") || die "Can't create $outf2";
//	open(O3, "| bzip2 - > $outf3") || die "Can't create $outf3";

//	# print the bam header
//  open(F, "$samtools view -H $bamf |");
//	my @header = <F>;
//	print O2 @header;
//  close F;

//	open(F, "$samtools view $bamf|") || die "Can't open $bamf";
//	my ($seq, $cseq1, $cseq2, $clen, $polyAT, $nm, $md, $cigar, $qual, $qual1, $qual2, $s1, $s2, $selected, $selected2, $cpos);
	string cseq1;
	string cseq2;
	uint32_t clen = 0;
	bool polyAT = false;
	uint32_t nm = 0;
//	int64_t md;

	string qual1;
	string qual2;

	bool selected = false;
	bool selected2 = false;
	int32_t cpos;

	string poly_a_str(min_polyAT, 'A');
	string poly_t_str(min_polyAT, 'T');
	string a_delim_str("^");
	for (uint8_t c = 'A'; c < (static_cast<uint8_t>('Z') + 1); ++c) {
		a_delim_str += static_cast<char>(c);
	}
	for (uint8_t c = 'a'; c < (static_cast<uint8_t>('z') + 1); ++c) {
		a_delim_str += static_cast<char>(c);
	}
	int64_t n_entries = 0;
	vector<string> b;
	const char* delim_all_stopwords = a_delim_str.c_str();
	auto& ref_vec = local_reader.GetReferenceData();
	BamTools::BamAlignment local_alignment_entry;
	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		++n_entries;
		if (0 == (n_entries & 1048575)) {
			cout << (boost::format("%d processed\n") % n_entries).str();
		}
//	while(<F>) {
//		chomp; my @a=split(/\t/);
//		if ($a[1] & 0x400 | $a[1] & 0x4) { next }; # skip PCR duplicate or unmapped reads
		uint32_t sam_flag = local_alignment_entry.AlignmentFlag;
		if ((sam_flag & 0x400) || (sam_flag & 0x4)) {
			continue;
		}

//		$seq = $a[9]; $cigar = $a[5]; $qual = $a[10];
		auto& seq = local_alignment_entry.QueryBases;
		auto& qual = local_alignment_entry.Qualities;
//		$s1 = $cigar;
//		$s2 = $cigar;
		polyAT = selected = selected2 = false;
		auto& cigar = local_alignment_entry.CigarData;
		// # clipped in the beginning
		if (!cigar.empty() && 'S' == cigar.front().Type) {
//			$s1 =~ /^(\d+)S\S+/) {
			clen = cigar.front().Length;
//		$clen = $1;
			if (clen >= 5) {
				qual1 = qual.substr(0, clen);
				int64_t n_good_quals = get_number_of_good_qualities(qual1, qcutoff);
//        	if (qual1 =~ m/([^#]+)/) {
//						 if (length($1) >= 5) {  //# more than 5 good quality bases
				if (n_good_quals > 5) {
					selected = true;
					if (!(sam_flag & 0x0010)) { // # positive strand mapping
						selected2 = true;
						cpos = local_alignment_entry.Position;
//								# check whether the clipped seq has >=10
						cseq1 = seq.substr(0, clen);
						if (string::npos != cseq1.find(poly_a_str)) {
							polyAT = true;
						}
					}
				}
			}
//				}
		}
		// # clipping in the end
		clen = 0;
//		for (auto& a_cigar : cigar) {
//			if ('D' == a_cigar.Type) {
//				clen = a_cigar.Length;
//				break;
//			}
//		}
		if ('S' == cigar.back().Type) {
			clen = cigar.back().Length;
		}
		if (clen >= 5) {
			int64_t the_qual_sub_str_pos = qual.size() - clen;
			if (the_qual_sub_str_pos > 0) {
				qual2 = qual.substr(the_qual_sub_str_pos);
			} else {
				qual2 = qual;
			}
			int64_t n_good_quals = get_number_of_good_qualities(qual2, qcutoff);
//		if ($qual2 =~ m/([^#]+)/) {
//				if (length($1) >= 5 ) {
			if (n_good_quals > 5) {
				selected = true;
				if (sam_flag & 0x0010) { // # negative strand mapping
					selected2 = true;
					cpos = get_cpos(local_alignment_entry.Position, cigar, qual2, -1);
//							# check whether the clipped seq has >=10
					int64_t the_seq_sub_str_pos = seq.size() - clen;
					cseq2 = seq.substr(the_seq_sub_str_pos);
//							if ($cseq2 =~ m/T{$min_polyAT,}/) { $polyAT = 1 }
					if (string::npos != cseq2.find(poly_t_str)) {
						polyAT = true;
					}
				}
			}
		}

//	# if there is not >=10 A/T bases in the clipped sequence, apply NM and MD filters
		if (!polyAT) {
//        # check max_mismatches
			local_alignment_entry.GetEditDistance(nm);
			if (nm > max_mismatches) {
				continue;
			};

			string md_str;
			local_alignment_entry.GetTag("MD", md_str);
//        $md = $_;

//        # check min_matches is satisfied

//        $md =~  s/.*\sMD:Z:([\w^]+)($|\s.*)/$1/;
//        my @b=split(/[A-Za-z^]/, $md);
//				cout << local_alignment_entry.Name << "/" << md_str << "\n";
			castle::StringUtils::c_string_multi_split(md_str, delim_all_stopwords, b);
			int64_t summed_n_bases = 0;
			for (auto a_value_str : b) {
				summed_n_bases += boost::lexical_cast<int64_t>(a_value_str);
			}
			if (summed_n_bases < min_matches) {
				continue;
			}
		}

		if (selected2) {
			O2.SaveSAMAlignment(local_alignment_entry);
			O3 << (boost::format("%s\t%s\n") % ref_vec[local_alignment_entry.RefID].RefName % cpos).str();
		}
	}
	cout << (boost::format("[TEA.generate_cbam_files] %d processed\n") % n_entries).str();
//	close F;
	O2.Close();

	cout << (boost::format("done generating cbam: %s\n") % consd_raw_bam_name).str();
	cout << (boost::format("start sorting and generating the index for %s ...\n") % consd_raw_bam_name).str();

//	string prefix = outf2;
//	$prefix =~ s/\.raw\.bam$ //;
//	system("samtools sort $outf2 $prefix; $samtools index $prefix.bam; rm $outf2");
	string sambamba_sort_cmd = (boost::format("sambamba sort -l 1 -t %d -o %s %s") % n_cores % out_bam % consd_raw_bam_name).str();
	system(sambamba_sort_cmd.c_str());
	cout << "done generating cbam and its index.\n";
}

void TEA::is_containing_S_or_H(bool& has_S, bool& has_H, vector<BamAlignment>& alns) {
	for(auto& an_aln : alns) {
		auto& cigars = an_aln.CigarData;
		if(cigars.empty()) {
			continue;
		}
		char front_type = cigars.front().Type;
		char back_type = cigars.back().Type;
		if('S' == front_type || 'S' == back_type) {
			has_S = true;
		} else if('S' == back_type || 'H' == back_type) {
			has_H = true;
		}
		if(has_S && has_H) {
			break;
		}
	}
}
bool TEA::has_S_or_H(const BamAlignment& aln) {
	auto& cigars = aln.CigarData;
	if(cigars.empty()) {
		return false;
	}
	char front_type = cigars.front().Type;
	char back_type = cigars.back().Type;
	return 'S' == front_type || 'S' == back_type || 'H' == front_type || 'H' == back_type;
}
bool TEA::has_tag(const BamAlignment& aln, const char tag) {
	auto& cigars = aln.CigarData;
	if(cigars.empty()) {
		return false;
	}
	char front_type = cigars.front().Type;
	char back_type = cigars.back().Type;
	return tag == front_type || tag == back_type;
}
void TEA::find_quality_standard(BamTools::BamReader& a_reader) {
	a_reader.Rewind();
	meerkat::ReadGroup rg;
	vector<uint64_t> quality_sample_space(255);
	BamTools::BamAlignment local_alignment_entry;
	int64_t n_samples = 0;
	while (a_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		auto& qual = local_alignment_entry.Qualities;
		for (uint64_t qual_id = 0; qual_id < qual.size(); ++qual_id) {
			int32_t the_chr_id = qual[qual_id];
			++quality_sample_space[the_chr_id];
		}
		++n_samples;
		if (n_samples > 600000) {
			break;
		}
	}
	int a = 0;
	int b = 0;
	for (int i = 0; i < 255; ++i) {
		if (quality_sample_space[i] != 0) {
			a = i;
			break;
		}
	}
	for (int i = 255 - 1; i >= 0; --i) {
		if (quality_sample_space[i] != 0) {
			b = i;
			break;
		}
	}

	/* Simple guess at the encoding type */
	int best = 0, best_diff = INT_MAX;
	for (unsigned int i = 0; i < sizeof(meerkat::ReadGroup::known_qencodings) / sizeof(meerkat::QENCODING); ++i) {
		int d = abs(meerkat::ReadGroup::known_qencodings[i].min - a) + abs(meerkat::ReadGroup::known_qencodings[i].max - b);
		if (d < best_diff) {
			best = i;
			best_diff = d;
		}
	}
	qenc = best;
	for (int i = 0; i < static_cast<int>(sizeof(meerkat::ReadGroup::known_qencodings) / sizeof(meerkat::QENCODING)); ++i) {
		if (qenc == i)
			cout << "    * ";
		else
			cout << "      ";

		cout.width(25);
		cout << left << meerkat::ReadGroup::known_qencodings[i].name << "[" << meerkat::ReadGroup::known_qencodings[i].min << ", " << meerkat::ReadGroup::known_qencodings[i].max << "]\n";
	}
	min_qual = meerkat::ReadGroup::known_qencodings[qenc].min;
	a_reader.Rewind();
}

void TEA::load_repeat_annotation(boost::unordered_map<string, pair<string, string>>& rannot) {
	string rannot_file = options.rannot_file;
	if (rannot_file.empty()) {
		rannot_file = (boost::format("%s/lib/repeats.txt") % string(getenv("tea_base"))).str();
	}
	cout << "[TEA.load_repeat_annotation] reads repeat annotation from " << rannot_file << "\n";
	const char* delim_tab = "\t";
	const char* delim_slash = "/";
	vector<string> data;
	vector<string> cl_data;

	string line;
	ifstream in_rannot(rannot_file, ios::binary);
	while (getline(in_rannot, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		if (data.size() < 2) {
			continue;
		}
		auto& the_key = data[0];
		castle::StringUtils::c_string_multi_split(data[1], delim_slash, cl_data);
		if (cl_data.size() > 1) {
			rannot[the_key].first = cl_data[0];
			rannot[the_key].second = cl_data[1];
		} else {
			rannot[the_key].first = "NA";
		}
	}
}

void TEA::load_virus_annotation(map<int64_t, string>& vannot) {
	string vannot_file = options.rannot_file;
	if (vannot_file.empty()) {
		vannot_file = (boost::format("%s/lib/viruses.txt") % string(getenv("tea_base"))).str();
	}
	cout << "[TEA.load_virus_annotation] reads virus annotation from " << vannot_file << "\n";
	const char* delim_ncbi = "|,";
	vector<string> data;
	vector<string> cl_data;

	string line;
	ifstream in_rannot(vannot_file, ios::binary);
	while (getline(in_rannot, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_ncbi, data);
		if (data.size() < 2) {
			continue;
		}
		auto the_key = boost::lexical_cast<int64_t>(data[1]);
		if (data.size() > 4) {
			vannot[the_key] = data[4];
		}
	}
}

void TEA::load_ref_annotation(set<string>& chrl, boost::unordered_map<string, pair<string, string>>& rannot,
//		boost::unordered_map<string, boost::unordered_map<string, vector<pair<int64_t, int64_t> > > >& ril_annot,
		boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt, boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot,
boost::unordered_map<string, GeneIntervalVector>& gene_annot, bool out_chrl, bool out_gap) {
	//ref, rmasker.rfile, out.chrl=T, out.ril=!exo&annot.oi, out.gene=annot.gene, out.gap=F, verbose
	castle::TimeChecker checker;
	checker.setTarget("TEA.load_ref_annotation");
	checker.start();
	string& ref = options.ref;
	string& rmasker_rfile = options.rmasker_rfile;
//	const bool out_ril = !options.exo && options.annot_oi;
//	const bool out_ril = !rmasker_rfile.empty();
	const bool out_gene = options.annot_gene;

	if (out_chrl) {
		load_chr(chrl);
	}
//	if(out_ril) {
		if (rmasker_rfile.empty()) {
			rmasker_rfile = (boost::format("%s/lib/rmasker/rmasker.%s.merged.RData.txt") % string(getenv("tea_base")) % ref).str();
		}
		if (boost::filesystem::exists(rmasker_rfile)) {
			cout << (boost::format("[TEA.load_ref_annotation] reads rmasker te annotation %s\n") % rmasker_rfile).str();
			const char* delim_tab = " \t";
			vector<string> data;
			string current_repeat_name;
			string current_repeat_class;
			string current_repeat_family;
			string current_chromosome;
			string line;
			ifstream in_rmasker(rmasker_rfile, ios::binary);
			while(getline(in_rmasker, line, '\n')) {
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
				if(string::npos != line.find("key:")) {
					current_repeat_name = data[1];
					for(auto& c : data[1]) {
						c = toupper(c);
					}
					auto a_repeat_entry = rannot.find(data[1]);
					if(rannot.end() != a_repeat_entry) {
						current_repeat_class = a_repeat_entry->second.first;
						current_repeat_family = a_repeat_entry->second.second;
					}

					current_chromosome = data[2];
					continue;
				}
				if(data.size() < 2) {
					continue;
				}
				int64_t start_pos = boost::lexical_cast<int64_t>(data[0]);
				int64_t end_pos = boost::lexical_cast<int64_t>(data[1]);
//				ril_annot[current_repeat_class][current_chromosome].push_back(make_pair(start_pos, end_pos));
				RefRepeatEntry an_entry;
				an_entry.chromosome = current_chromosome;
				an_entry.repeat_name = current_repeat_name;
				an_entry.repeat_class = current_repeat_class;
				an_entry.repeat_family = current_repeat_family;
				RefRepeatIntervalEntry an_interval_entry(start_pos, end_pos, an_entry);
				ril_annot_alt[current_chromosome].push_back(an_interval_entry);
			}
//				ril = readRDS(rmasker.rfile)
//				#outfile = sprintf("%s.txt", rmasker.rfile)
			cout << (boost::format("[TEA.load_ref_annotation] rmasker annotation loaded for %d chromosomes\n") % ril_annot_alt.size()).str();
		} else {
			cout << (boost::format("[TEA.load_ref_annotation] no rmasker annotation rfile: %s\n") % rmasker_rfile).str();
			exit(1);
		}
//	}
	// # reference gap annotation
	if (out_gap) {
		string gap_fname = (boost::format("%s/lib/gap/gap.%s.gz") % string(getenv("tea_base")) % ref).str();
		cout << (boost::format("[TEA.load_ref_annotation] reads %s\n") % gap_fname).str();
		read_gap_rfile(gap_annot, gap_fname, chrl);
		cout << (boost::format("[TEA.load_ref_annotation] # gapped regions: %d\n") % gap_annot.size()).str();
	}
	// # gene annotation
	if (out_gene) {
		string gene_rfile = (boost::format("%s/lib/gene/%s.genes.RData.txt") % string(getenv("tea_base")) % ref).str();
		if (boost::filesystem::exists(gene_rfile)) {
			cout << (boost::format("[TEA.load_ref_annotation] reads gene annotations from %s\n") % gene_rfile).str();
			read_gene_rfile(gene_annot, gene_rfile, chrl);
//			genes = readRDS(gene.rfile)
		} else {
			cout << (boost::format("[TEA.load_ref_annotation] no gene annotation file: %s\n") % gene_rfile).str();
			exit(1);
		}
	}
	cout << checker;
}

void TEA::load_chr(set<string>& chrl) const {
	const string& ref = options.ref;
	if (ref == "hg18" || ref == "hg19") {
		for(int32_t chr_id = 1; chr_id < 23; ++chr_id) {
			string str_chr_id = boost::lexical_cast<string>(chr_id);
			string a_chr_id = "chr" + str_chr_id;
			chrl.insert(a_chr_id);
		}
		chrl.insert("chrX");
		chrl.insert("chrY");
	} else if (ref == "ponAbe2") {
		for(int32_t chr_id = 3; chr_id < 23; ++chr_id) {
			string str_chr_id = boost::lexical_cast<string>(chr_id);
			string a_chr_id = "chr" + str_chr_id;
			chrl.insert(a_chr_id);
		}
		chrl.insert("chr1");
		chrl.insert("chr2a");
		chrl.insert("chr2b");
		chrl.insert("chrX");
	} else if (ref == "panTro3") {
		for(int32_t chr_id = 3; chr_id < 23; ++chr_id) {
			string str_chr_id = boost::lexical_cast<string>(chr_id);
			string a_chr_id = "chr" + str_chr_id;
			chrl.insert(a_chr_id);
		}
		chrl.insert("chr1");
		chrl.insert("chr2A");
		chrl.insert("chr2B");
		chrl.insert("chrX");
		chrl.insert("chrY");
	} else if (ref == "rheMac2") {
		for(int32_t chr_id = 1; chr_id < 21; ++chr_id) {
			string str_chr_id = boost::lexical_cast<string>(chr_id);
			string a_chr_id = "chr" + str_chr_id;
			chrl.insert(a_chr_id);
		}
		chrl.insert("chrX");
	}
}

void TEA::read_gap_rfile(boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot, const string& file_name, const set<string>& chrl) {
	string line;

	int32_t chr_idx = 0;
	int32_t chr_start_idx = 1;
	int32_t chr_end_idx = 2;

	const char* delim_tab = "\t";
	vector<string> data;
	ifstream in_gap(file_name, ios::binary);
	while(getline(in_gap, line, '\n')) {
		if(line.empty()) {
			continue;
		}
		if('#' == line[0]) {
			// remove the initial #
			line = line.substr(1);
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			for(uint64_t d_id = 0; d_id < data.size(); ++d_id) {
				if("chrom" == data[d_id]) {
					chr_idx = d_id;
				} else if("chromStart" == data[d_id]) {
					chr_start_idx = d_id;
				} else if("chromEnd" == data[d_id]) {
					chr_end_idx = d_id;
				}
			}
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		string& chr_name = data[chr_idx];
		if(!chrl.empty()) {
			if(chrl.end() == chrl.find(chr_name)) {
				continue;
			}
		}
		// note that the start coordinates are 0-based
		int64_t start_pos = boost::lexical_cast<int64_t>(data[chr_start_idx]) + 1;
		int64_t end_pos = boost::lexical_cast<int64_t>(data[chr_end_idx]) + 1;
		gap_annot[chr_name].push_back(make_pair(start_pos, end_pos));
	}
}

void TEA::read_gene_rfile(boost::unordered_map<string, GeneIntervalVector>& gene_annot, const string& file_name, const set<string>& chrl) {
	string line;
	const char* delim_tab = "\t";
	vector<string> data;
	ifstream in_gene(file_name, ios::binary);
	cout << "[TEA.read_gene_rfile] reading " << file_name << "\n";
	// ignore the first line
	getline(in_gene, line, '\n');
	while (getline(in_gene, line, '\n')) {
		if (line.empty()) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		if (data.size() < 5) {
			continue;
		}
		string& chr_name = data[0];
		GeneEntry an_entry;
		an_entry.strand = data[1][0];
		an_entry.s = boost::lexical_cast<int64_t>(data[2]);
		an_entry.e = boost::lexical_cast<int64_t>(data[3]);
		an_entry.name = data[4];
		an_entry.type = data[5];
		GeneIntervalEntry an_interval_entry(an_entry.s, an_entry.e, an_entry);
		gene_annot[chr_name].push_back(an_interval_entry);
	}
}

void TEA::load_read_length(boost::unordered_map<string, int32_t>& rl) {
	string rl_file = options.prefix + ".rl";
	if (!options.working_dir.empty()) {
		rl_file = options.working_prefix + ".rl";
	}
	cout << (boost::format("[TEA.load_read_length] reads rl from %s\n") % rl_file);
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in_rl(rl_file, ios::binary);
	while (getline(in_rl, line, '\n')) {
		if (line.empty()) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		if (data.size() < 2) {
			continue;
		}
		auto& rg_name = data[0];
		int32_t read_length = boost::lexical_cast<int32_t>(data[1]);
		rl[rg_name] = read_length;
	}
}

void TEA::load_insert_size(boost::unordered_map<string, boost::unordered_map<string, double>>& is, boost::unordered_map<string, int32_t>& rl) {
	string isize_file = options.prefix + ".isize";
	if (!options.working_dir.empty()) {
		isize_file = options.working_prefix + ".isize";
	}
	cout << (boost::format("[TEA.load_insert_size] reads isize from %s\n") % isize_file);
	if (!boost::filesystem::exists(isize_file)) {
		return;
	}
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in_is(isize_file, ios::binary);
	while (getline(in_is, line, '\n')) {
		if (line.empty()) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
		if (data.size() < 3) {
			continue;
		}
		auto& rg_name = data[0];
		int32_t mean = boost::lexical_cast<int32_t>(data[1]);
		int32_t stdev = boost::lexical_cast<int32_t>(data[2]);
		is[rg_name]["mu"] = mean;
		is[rg_name]["sd"] = stdev;
	}

	double rl_val = rl["all"];
	double mu = floor(is["all"]["mu"]);
	double sd = is["all"]["sd"];
	double fr = floor(mu + 3 * sd);
	double candidate_intra_gap = floor(static_cast<double>((fr - rl_val) * 2) / (double) 3);
	double intra_gap = min(static_cast<double>(250.0), candidate_intra_gap);
	double inter_gap = intra_gap * 2;
	double gap_size = floor(static_cast<double>(mu * 2) / (double) 3 + rl_val);
	double ins_margin = floor(1.5 * fr);
	auto& the_represent = is["all"];
	the_represent["mu"] = mu;
	the_represent["sd"] = sd;
	the_represent["fr"] = fr;
	the_represent["intra_gap"] = intra_gap;
	the_represent["inter_gap"] = inter_gap;
	the_represent["gap_size"] = gap_size;
	the_represent["ins_margin"] = ins_margin;
}
void TEA::load_ram(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram, boost::unordered_map<string, pair<string, string>>& rannot, const bool rm_dup) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.load_ram");
	checker.start();
	string ram_file = options.prefix + "." + options.rasym + "m";

	if(!options.working_dir.empty()) {
		ram_file = options.working_prefix + "." + options.rasym + "m";
	}
	cout << (boost::format("[TEA.load_ram] reads rams from %s\n") % ram_file);
	const char* delim_tab = "\t";
	const char* delim_comma = ",";
	vector<string> data;
	vector<string> ram_data;
	vector<string> repeat_class_data;

	string line;
	ifstream in_ram(ram_file, ios::binary);
	if(options.is_sampe) {
		while(getline(in_ram, line, '\n')) {
			if(line.empty()) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			if(3 == data.size()) {
				auto& chr = data[1];
				RAMRepeatEntry an_entry;
				an_entry.read_name = data[0];
				an_entry.pos = boost::lexical_cast<int64_t>(data[2]);
				if(an_entry.pos > 0) {
					ram[chr][1].push_back(an_entry);
				} else if(an_entry.pos < 0) {
					an_entry.pos = -an_entry.pos;
					ram[chr][-1].push_back(an_entry);
				}
				// read name, chromosome, position, repeat class
			} else if(4 == data.size()) {
				RAMRepeatEntry an_entry;
				an_entry.read_name = data[0];
				auto& chr = data[1];
				an_entry.pos = boost::lexical_cast<int64_t>(data[2]);
				auto& repeat_name = data[3];

				boost::replace_all(repeat_name, "/", "_");
	//			if(string::npos != repeat_name.find(",")) {
				castle::StringUtils::c_string_multi_split(repeat_name, delim_comma, repeat_class_data);
				for(auto& a_repeat_name: repeat_class_data) {
	//					repeat_name = repeat_class_data[0];
					string capitalized_repeat_name(a_repeat_name);
					for (auto &c: capitalized_repeat_name) {
						c = toupper(c);
					}
					an_entry.repeat_name = a_repeat_name;

					auto a_repeat_entry = rannot.find(capitalized_repeat_name);
					if(rannot.end() != a_repeat_entry) {
						auto& the_class = a_repeat_entry->second.first;
						auto& the_family = a_repeat_entry->second.second;
						an_entry.repeat_class = the_class;
						an_entry.repeat_family = the_family;
					} else {
						an_entry.repeat_class = "-";
						an_entry.repeat_family = "*";
					}

					if(an_entry.pos > 0) {
						ram[chr][1].push_back(an_entry);
					} else if(an_entry.pos < 0) {
						an_entry.pos = -an_entry.pos;
						ram[chr][-1].push_back(an_entry);
					}
				}
			}
		}
	} else {
		map<string, string> prev_entries;
		string prev_name;
		while(getline(in_ram, line, '\n')) {
			if(line.empty()) {
				continue;
			}
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			if(5 == data.size()) {
				auto& cur_name = data[0];
				if(prev_name != cur_name) {
					for(auto& ram_entry : prev_entries) {
						auto& a_line = ram_entry.second;
						castle::StringUtils::c_string_multi_split(a_line, delim_tab, ram_data);
						RAMRepeatEntry an_entry;
						an_entry.read_name = ram_data[0];
						auto& chr = ram_data[1];
						an_entry.pos = boost::lexical_cast<int64_t>(ram_data[2]);
						auto& repeat_name = ram_data[4];
						boost::replace_all(repeat_name, "/", "_");
						castle::StringUtils::c_string_multi_split(repeat_name, delim_comma, repeat_class_data);

						for(auto& a_repeat_name: repeat_class_data) {
							string capitalized_repeat_name(a_repeat_name);
							for (auto &c: capitalized_repeat_name) {
								c = toupper(c);
							}
							an_entry.repeat_name = a_repeat_name;

							auto a_repeat_entry = rannot.find(capitalized_repeat_name);
							if(rannot.end() != a_repeat_entry) {
								auto& the_class = a_repeat_entry->second.first;
								auto& the_family = a_repeat_entry->second.second;
								an_entry.repeat_class = the_class;
								an_entry.repeat_family = the_family;
							} else {
								an_entry.repeat_class = "-";
								an_entry.repeat_family = "-";
							}

							if(an_entry.pos > 0) {
								ram[chr][1].push_back(an_entry);
							} else if(an_entry.pos < 0) {
								an_entry.pos = -an_entry.pos;
								ram[chr][-1].push_back(an_entry);
							}
						}
					}
					prev_entries.clear();
					prev_name = cur_name;
				}
			prev_entries.insert(make_pair(data[0] + data[4], line));
			}
		}
	}
	int64_t n_ram = 0;
	for(auto& chr_entry: ram) {
		for(auto& sign_entry : chr_entry.second) {
			auto& the_pos_vec = sign_entry.second;
			sort(the_pos_vec.begin(), the_pos_vec.end());
			if(rm_dup) {
				the_pos_vec.erase(unique(the_pos_vec.begin(), the_pos_vec.end()), the_pos_vec.end());
			}
			n_ram += the_pos_vec.size();
		}
	}
	cout << "[TEA.load_ram] # rams: " << n_ram << "\n";
	cout << checker;
}

void TEA::get_cluster_alt(const string& chr, RAMIntervalVector& cl, vector<RAMRepeatEntry>& sram, boost::unordered_map<string, pair<string, string>>& rannot, const int32_t strand, int64_t gap_cutoff) {
//	castle::TimeChecker checker;
//	checker.setTarget("TEA.get_cluster_alt");
//	checker.start();
	int64_t max_id = sram.size();
	if (max_id < 2) {
		if (1 == max_id) {
			RepeatClusterEntry an_entry;
			an_entry.ram = 1;
			//		cout << "[TEA.get_cluster] ram size: " << an_entry.ram << "\n";
			//		cout << "[TEA.get_cluster] the parent family: " << sram[si].repeat_family << "\n";
			an_entry.s = sram[0].pos;
			an_entry.e = sram[0].pos;
			auto& r = sram[0];
			an_entry.rep_repeat.insert(r.repeat_name);
//			cout << "[TEA.get_cluster] the current family: " << r.repeat_family << "\n";
			an_entry.family.push_back(r.repeat_family);
			an_entry.repeat_class.push_back(r.repeat_class);
			an_entry.repeat_name.push_back(r.repeat_name);
			an_entry.pos.push_back(r.pos);
			an_entry.rname.push_back(r.read_name);
			RAMIntervalEntry an_interval_entry(an_entry.s, an_entry.e, an_entry);
			cl.push_back(an_interval_entry);
		}
		return;
	}
	// declares each break point if the gap between rams is larger than gap_cutoff or the repeat family name is not the same except for PolyA family
	// since Poly-A family can be attached to any family.
	boost::unordered_map<string, vector<vector<int64_t>>> breakpoint_ids;
	for (int64_t pos_id = 0; pos_id < max_id; ++pos_id) {
		auto& a_ram = sram[pos_id];
		auto& the_family = a_ram.repeat_family;
		if ("PolyA" == the_family) {
			continue;
		}
		breakpoint_ids[the_family].reserve(10000);
	}
//	const bool debug = false;
	for (int64_t pos_id = 0; pos_id < max_id; ++pos_id) {
		auto& a_ram = sram[pos_id];
		const bool debug = chr == "1" && a_ram.pos >= 93518457 && a_ram.pos <= 93519148;
		if(debug) {
			cout << a_ram.str() << "\n";
		}

		auto& the_family = a_ram.repeat_family;

		if ("PolyA" == the_family) {
			for (auto& a_cluster_parent_vec_entry : breakpoint_ids) {
//				auto& local_family = a_cluster_parent_vec_entry.first;
				auto& a_cluster_parent_vec = a_cluster_parent_vec_entry.second;
				if (!a_cluster_parent_vec.empty()) {
					auto& a_cluster_vec = a_cluster_parent_vec.back();
					if (a_cluster_vec.empty()) {
						continue;
					}
					int64_t the_last_pos_id = a_cluster_vec.back();
					int64_t delta = sram[pos_id].pos - sram[the_last_pos_id].pos;
					if (delta <= gap_cutoff) {
//						cout << a_ram.str() << "\n";
//						cout << "Insert(Poly) : " << local_family << ":" << a_ram.pos << "\n";
						a_cluster_vec.push_back(pos_id);
					}
				}
			}
		} else {
			if (breakpoint_ids[the_family].empty()) {
				vector<int64_t> the_next_cluster_vec;
				breakpoint_ids[the_family].push_back(the_next_cluster_vec);
			}
			if (breakpoint_ids[the_family].back().empty()) {
				if(debug) {
					cout << "Insert(New) : " << a_ram.pos << "\n";
				}
				breakpoint_ids[the_family].back().push_back(pos_id);
			} else {
				int64_t the_last_pos_id = breakpoint_ids[the_family].back().back();
				auto& prev_ram = sram[the_last_pos_id];
				int64_t delta = a_ram.pos - prev_ram.pos;
				if(debug) {
					cout << "prev: " << prev_ram.str() << ", current: " << a_ram.str() << ", delta: " << delta << "\n";
				}

				if (delta > gap_cutoff) {
					if(debug) {
						cout << "Insert(New Current BreakPoint) : " << a_ram.pos << ", gap_cutoff: " << gap_cutoff << ", delta: " << delta << "\n";
					}
					vector<int64_t> the_next_cluster_vec;
					breakpoint_ids[the_family].push_back(the_next_cluster_vec);
				}
//				else {
//					cout << "Insert(Normal) : " << a_ram.pos << "\n";
//				}
				breakpoint_ids[the_family].back().push_back(pos_id);
			}
		}
	}
//	int64_t max_sram_id = sram.size();
	vector<RepeatClusterEntry> cluster_entries;
	for (auto& a_family_entry : breakpoint_ids) {
//		auto& the_family = a_family_entry.first;
//		cout << "[TEA.get_cluster] the parent family: " << the_family << "\n";
		auto& the_parent_cluster_vec = a_family_entry.second;
		for (auto& a_cluster : the_parent_cluster_vec) {
			if (a_cluster.empty()) {
				continue;
			}
			int64_t si = a_cluster[0];
			int64_t ei = a_cluster.back();
//			for(uint64_t cl_id = 0; cl_id < a_cluster.size(); ++cl_id) {
//				int64_t point_id = a_cluster[cl_id];
//				cout << sram[point_id].str() << "\n";
//			}
//			cout << "\n";
			RepeatClusterEntry an_entry;
			an_entry.ram = a_cluster.size();
//		cout << "[TEA.get_cluster] ram size: " << an_entry.ram << "\n";
//		cout << "[TEA.get_cluster] the parent family: " << sram[si].repeat_family << "\n";
			an_entry.s = sram[si].pos;
			an_entry.e = sram[ei].pos;

			boost::unordered_set<string> name_checker;
			for (uint64_t cl_id = 0; cl_id < a_cluster.size(); ++cl_id) {
				int64_t point_id = a_cluster[cl_id];
				auto& r = sram[point_id];
				if(name_checker.end() != name_checker.find(r.read_name)) {
					continue;
				}
				an_entry.rep_repeat.insert(r.repeat_name);
				//			cout << "[TEA.get_cluster] the current family: " << r.repeat_family << "\n";
				an_entry.family.push_back(r.repeat_family);
				an_entry.repeat_class.push_back(r.repeat_class);
				an_entry.repeat_name.push_back(r.repeat_name);
				an_entry.pos.push_back(r.pos);
				an_entry.rname.push_back(r.read_name);
				name_checker.insert(r.read_name);
				string sc_read_name = r.read_name + "sc";
				name_checker.insert(sc_read_name);
			}
			cluster_entries.push_back(an_entry);
		}
	}
	sort(cluster_entries.begin(), cluster_entries.end());
	for (auto& an_entry : cluster_entries) {
		RAMIntervalEntry an_interval_entry(an_entry.s, an_entry.e, an_entry);
		cl.push_back(an_interval_entry);
	}

//	cout << checker;
}

void TEA::pair_cluster_alt(map<int64_t, int64_t>& pm_cl, RAMIntervalVector& p_cl, RAMIntervalVector& n_cl, const int64_t gap_cutoff, const int64_t read_length, bool stringent_pair) {
//	castle::TimeChecker checker;
//	checker.setTarget("TEA.pair_cluster_alt");
//	checker.start();
//	cout << (boost::format("[TEA.pair_cluster_alt] gap cutoff: %d\n") % gap_cutoff).str();
	RAMIntervalTree negative_strand_interval_tree(n_cl);
	RAMIntervalVector results;
	results.reserve(10000);
	const int64_t ram_boundary_limit = 2 * read_length + 10;
	for (uint64_t r_id = 0; r_id < p_cl.size(); ++r_id) {
		auto& a_ram_entry = p_cl[r_id];
		RAMIntervalVector result_selected;
		int64_t n_ram_supports = 0;
//		int64_t ram_boundary_size = (a_ram_entry.stop - a_ram_entry.start + 1);
//		if (ram_boundary_size < ram_boundary_limit) {
//			continue;
//		}
		n_ram_supports += a_ram_entry.value.ram;
		results.clear();
		boost::unordered_set<string> possible_family;
		possible_family.insert("PolyA");
		possible_family.insert(a_ram_entry.value.family.begin(), a_ram_entry.value.family.end());
//		const bool debug  = (37925043 == a_ram_entry.start);
		const bool debug = false;
		negative_strand_interval_tree.find_overlap(a_ram_entry.start, a_ram_entry.stop + gap_cutoff, results);
		if (debug) {
			cout << "[TEA.pair_cluster_alt] results: " << results.size() << "\n";
		}
		for (auto& a_negative_ram_entry : results) {

			bool found_incompatible = false;
			for (auto& a_family : a_negative_ram_entry.value.family) {
				if (possible_family.end() == possible_family.find(a_family)) {
					found_incompatible = true;
					break;
				}
			}
			if (debug) {
				cout << "[TEA.pair_cluster_alt] found_incompatible: " << found_incompatible << "\n";
			}
			if (found_incompatible) {
				continue;
			}

			const int64_t ram_boundary_size = a_negative_ram_entry.stop - a_ram_entry.start + read_length;
			if (debug) {
				cout << "[TEA.pair_cluster_alt] ram_boundary_size: " << ram_boundary_size << "/" << ram_boundary_limit << "\n";
			}
			// the ram boundary is smaller than the minimum ram boundary size
			if (ram_boundary_size > ram_boundary_limit) {
				if (stringent_pair) {
					if (a_ram_entry.start <= a_negative_ram_entry.start && a_ram_entry.stop <= a_negative_ram_entry.stop) {
						result_selected.push_back(a_negative_ram_entry);
					}
				} else {
					if (a_ram_entry.start < a_negative_ram_entry.stop) {
						result_selected.push_back(a_negative_ram_entry);
					}
				}
			}
		}

		if (result_selected.size() > 0) {
//			cout << "Positive: " << n_ram_supports << "/Negative: " << result_selected.size() << "\n";
			for (auto& a_negative_entry : result_selected) {
//				pm_cl[a_negative_entry.value.global_cluster_id] = result_selected;
				pm_cl[r_id] = a_negative_entry.value.global_cluster_id;
//				cout << r_id << "<->" << a_negative_entry.value.global_cluster_id << "\n";
			}
//			pm_cl[r_id] = result_selected;
		}
	}
//	cout << checker;
}

void TEA::count_clipped(boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt, boost::unordered_map<string, GeneIntervalVector>& gene_annot, const string& chr, const string& cl_prefix, const string& contig_dir, const map<int64_t, int64_t>& pm_cl, const RAMIntervalVector& p_cl,
		const RAMIntervalVector& n_cl, boost::unordered_set<int64_t>& positive_only, boost::unordered_set<int64_t>& negative_only, const int64_t read_length, const int64_t fragment_size, const int64_t rmasker_filter_margin, const int64_t gene_margin, const bool headless) {
//	castle::TimeChecker checker;
//	checker.setTarget("TEA.count_clipped");
//	checker.start();
	string tmp_chr_name(chr);
	if (string::npos == tmp_chr_name.find("chr")) {
		tmp_chr_name = "chr" + chr;
	}
	const string prefixed_chr = "chr" + chr;
	string input_softclips_consd_bam_name = options.prefix + ".softclips.consd.bam";
	string input_softclips_consd_cpos_name = options.prefix + ".softclips.consd.cpos";
	if (!options.working_dir.empty()) {
		input_softclips_consd_bam_name = options.working_prefix + ".softclips.consd.bam";
		input_softclips_consd_cpos_name = options.working_prefix + ".softclips.consd.cpos";
	}

	auto repeat_annot_itr = ril_annot_alt.find(chr);
	if (ril_annot_alt.end() == repeat_annot_itr) {
		repeat_annot_itr = ril_annot_alt.find(prefixed_chr);
	}
	if (ril_annot_alt.end() == repeat_annot_itr) {
//		cout << (boost::format("[TEA.count_clipped] locating repeat annotation failed: %s\n") % prefixed_chr).str();
		return;
	}

	auto gene_itr = gene_annot.find(chr);
	if (gene_annot.end() == gene_itr) {
		gene_itr = gene_annot.find(prefixed_chr);
	}
	if (gene_annot.end() == gene_itr) {
//		cout << (boost::format("[TEA.count_clipped] locating gene annotation failed: %s\n") % chr).str();
		return;
	}
	auto& the_target_repeat_annot = ril_annot_alt.end() == ril_annot_alt.find(chr) ? ril_annot_alt[prefixed_chr] : ril_annot_alt[chr];
	auto& the_target_gene_annot = gene_annot.end() == gene_annot.find(chr) ? gene_annot[prefixed_chr] : gene_annot[chr];
	string a_bai_path;
	get_bai_index_path(input_softclips_consd_bam_name, a_bai_path);

	BamTools::BamReader local_reader;
	if (!local_reader.Open(input_softclips_consd_bam_name, a_bai_path)) {
		cout << "[TEA.count_clipped] ERROR: could not open BAM file '" << input_softclips_consd_bam_name << "'\n";
		exit(1);
	}
	map<string, int64_t> ref_reverse_index;
	const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
	for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
		auto& a_ref = a_ref_vector[ref_id];
		ref_reverse_index[a_ref.RefName] = ref_id;
	}
	auto rev_itr = ref_reverse_index.find(chr);
	if (ref_reverse_index.end() == rev_itr) {
		cout << "the chromosome name is not the same\n";
		exit(0);
	}

	string header_cl = "chr\ts\te\tsize\ttsd\tpbp\tnbp\trep.repeat\tfamily\tclass\tram\tpram\tnram\tcr\tpcr\tncr\t"
			"acr\tpacr\tnacr\tacrr\tpram_start\tpram_end\tnram_start\tnram_end\tpgene\tngene\tscore\ts2n\toi\tdesc\tconf\n";
	string header_germline = "sample\t" + header_cl;

	string header_clipped = "chr\ts\te\trep.repeat\tcpos\taligned\tcigar\ttname\tref.seq\tclipped.seq\tclipped.qual\n";

	RefRepeatIntervalTree ref_repeat_interval_tree(the_target_repeat_annot);
	RefRepeatIntervalVector stat_results;
	GeneIntervalTree gene_interval_tree(the_target_gene_annot);
	GeneIntervalVector gene_results;
	const int64_t chr_ref_id = rev_itr->second;
//	cout << "[TEA.count_clipped] process paired\n";
	string cl_file = cl_prefix + "." + tmp_chr_name + ".cluster";
	string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";
	string p_mate_readname_file = cl_prefix + "." + tmp_chr_name + ".p.mate.rname";
	string n_mate_readname_file = cl_prefix + "." + tmp_chr_name + ".n.mate.rname";
	string p_clipped_filename_file = cl_prefix + "." + tmp_chr_name + ".p.clipped.fname";
	string n_clipped_filename_file = cl_prefix + "." + tmp_chr_name + ".n.clipped.fname";

	string clipped_file = cl_prefix + "." + tmp_chr_name + ".clipped";

	ofstream out_cl(cl_file, ios::binary);
	ofstream out_germline(germline_file, ios::binary);
	ofstream out_clipped(clipped_file, ios::binary);
	ofstream out_p_mate_rname(p_mate_readname_file, ios::binary);
	ofstream out_n_mate_rname(n_mate_readname_file, ios::binary);
	ofstream out_p_clipped_filename(p_clipped_filename_file, ios::binary);
	ofstream out_n_clipped_filename(n_clipped_filename_file, ios::binary);

	if(!headless) {
		out_cl << header_cl;
		out_germline << header_germline;
		out_clipped << header_clipped;
	}
	{
		for (auto& a_pair : pm_cl) {
			int64_t p_idx = a_pair.first;
			int64_t n_idx = a_pair.second;
			auto& positive_entry = p_cl[p_idx];
			auto& negative_entry = n_cl[n_idx];
			int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
			if (n_ram < options.min_ram) {
				continue;
			}
			// s column
			int64_t the_ram_boundary_start = positive_entry.start;
			// e column
			int64_t the_ram_boundary_end = negative_entry.stop + read_length;
			int64_t start_pos = the_ram_boundary_start + read_length;
			if (!positive_entry.value.pos.empty()) {
//				start_pos = min(start_pos, positive_entry.value.pos.back());
				start_pos = min(start_pos, positive_entry.value.pos.back() - read_length);
			}
//			int64_t end_pos = the_ram_boundary_end - read_length;
			int64_t end_pos = negative_entry.value.pos[0] + read_length;
			if (!negative_entry.value.pos.empty()) {
				end_pos = max(end_pos, the_ram_boundary_end - read_length);
			}
			const bool debug = (ref_reverse_index["chr2"] == chr_ref_id && the_ram_boundary_start == 31058709 && 31058709 == the_ram_boundary_end);
			if(debug) {
				cout << (boost::format("[TEA.count_clipped] pram start: %d, pram end: %d, nram start: %d, nram end: %d\n") % positive_entry.start % positive_entry.stop % negative_entry.start % negative_entry.stop).str();
				cout << (boost::format("[TEA.count_clipped] pram start-2: %d, pram end: %d, nram start: %d, nram end: %d\n") % positive_entry.value.pos[0] % positive_entry.value.pos.back() % negative_entry.value.pos[0] % negative_entry.value.pos.back()).str();
				cout << (boost::format("[TEA.count_clipped] BAM start: %d, end: %d, read len.: %d\n") % start_pos % end_pos % read_length).str();
				for(auto pos: positive_entry.value.pos) {
					cout << (boost::format("[TEA.count_clipped] pos position: %d\n") % pos).str();
				}
				for(auto pos: negative_entry.value.pos) {
					cout << (boost::format("[TEA.count_clipped] neg position: %d\n") % pos).str();
				}
			}

			int64_t mid_point = (positive_entry.stop + negative_entry.start + read_length) / (double)2;
			local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
			output_clipped_stat(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, gene_interval_tree, gene_results, local_reader, the_ram_boundary_start, the_ram_boundary_end,
					positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin, mid_point);
		}
	}
	RepeatClusterEntry an_empty_entry;
	RAMIntervalEntry empty_interval_entry(0, 0, an_empty_entry);
//	cout << "[TEA.count_clipped] process positive only\n";
	{
//		string cl_file = cl_prefix + "." + tmp_chr_name + ".p.cluster";
//		string clipped_file = cl_prefix + "." + tmp_chr_name + ".p.clipped";

//		ofstream out_cl(cl_file, ios::binary);
//		ofstream out_clipped(clipped_file, ios::binary);
//		out_cl << header_cl;
//		out_clipped << header_clipped;
		for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
			if (positive_only.end() == positive_only.find(r_id)) {
				continue;
			}
			auto& positive_entry = p_cl[r_id];
			auto& negative_entry = empty_interval_entry;
			int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
			if (n_ram < options.min_ram) {
				continue;
			}
			// s column
			int64_t the_ram_boundary_start = positive_entry.start;
			// e column
			int64_t the_ram_boundary_end = positive_entry.stop + read_length + fragment_size;

			int64_t start_pos = the_ram_boundary_start + read_length;
			if (!positive_entry.value.pos.empty()) {
				start_pos = min(start_pos, positive_entry.value.pos.back());
			}
			int64_t end_pos = the_ram_boundary_end - read_length;
//			if(!negative_entry.value.pos.empty()) {
//				end_pos = max(end_pos, negative_entry.value.pos[0]);
//			}

			local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
			int64_t mid_point = positive_entry.stop + read_length;
			output_clipped_stat(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, gene_interval_tree, gene_results, local_reader, the_ram_boundary_start, the_ram_boundary_end,
					positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin, mid_point);
		}
	}

//	cout << "[TEA.count_clipped] process negative only\n";
	{
//		string cl_file = cl_prefix + "." + tmp_chr_name + ".n.cluster";
//		string clipped_file = cl_prefix + "." + tmp_chr_name + ".n.clipped";
//		ofstream out_cl(cl_file, ios::binary);
//		ofstream out_clipped(clipped_file, ios::binary);
//		out_cl << header_cl;
//		out_clipped << header_clipped;
		for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
			if (negative_only.end() == negative_only.find(r_id)) {
				continue;
			}
			auto& positive_entry = empty_interval_entry;
			auto& negative_entry = n_cl[r_id];
			int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
			if (n_ram < options.min_ram) {
				continue;
			}
			// s column
			int64_t the_ram_boundary_start = negative_entry.start - fragment_size;
			// e column
			int64_t the_ram_boundary_end = negative_entry.stop + read_length;

			int64_t start_pos = the_ram_boundary_start + read_length;
//			if(!positive_entry.value.pos.empty()) {
//				start_pos = min(start_pos, positive_entry.value.pos.back());
//			}
			int64_t end_pos = the_ram_boundary_end - read_length;
			if (!negative_entry.value.pos.empty()) {
				end_pos = max(end_pos, negative_entry.value.pos[0]);
			}
			int64_t mid_point = negative_entry.start;
			local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
			output_clipped_stat(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, gene_interval_tree, gene_results, local_reader, the_ram_boundary_start, the_ram_boundary_end,
					positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin, mid_point);
		}
	}
	local_reader.Close();
//	cout << checker;
}

void TEA::count_clipped_v(
				boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
				map<int64_t, string>& vannot,
				const string& chr, const string& cl_prefix, const string& contig_dir,
				const map<int64_t, int64_t>& pm_cl, const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl,
				boost::unordered_set<int64_t>& positive_only, boost::unordered_set<int64_t>& negative_only,
				const int64_t read_length, const int64_t fragment_size,
				const int64_t rmasker_filter_margin, const int64_t gene_margin) {
//		castle::TimeChecker checker;
//		checker.setTarget("TEA.count_clipped");
//		checker.start();
		string tmp_chr_name(chr);
		if (string::npos == tmp_chr_name.find("chr")) {
			tmp_chr_name = "chr" + chr;
		}
		const string prefixed_chr = "chr" + chr;
		string input_softclips_consd_bam_name = options.prefix + ".softclips.consd.bam";
		string input_softclips_consd_cpos_name = options.prefix + ".softclips.consd.cpos";
		if (!options.working_dir.empty()) {
			input_softclips_consd_bam_name = options.working_prefix + ".softclips.consd.bam";
			input_softclips_consd_cpos_name = options.working_prefix + ".softclips.consd.cpos";
		}

		auto repeat_annot_itr = ril_annot_alt.find(chr);
		if (ril_annot_alt.end() == repeat_annot_itr) {
			repeat_annot_itr = ril_annot_alt.find(prefixed_chr);
		}
		if (ril_annot_alt.end() == repeat_annot_itr) {
			return;
		}

		auto& the_target_repeat_annot = ril_annot_alt.end() == ril_annot_alt.find(chr) ? ril_annot_alt[prefixed_chr] : ril_annot_alt[chr];
		string a_bai_path;
		get_bai_index_path(input_softclips_consd_bam_name, a_bai_path);

		BamTools::BamReader local_reader;
		if (!local_reader.Open(input_softclips_consd_bam_name, a_bai_path)) {
			cout << "ERROR: could not open BAM file '" << input_softclips_consd_bam_name << "'\n";
			exit(1);
		}
		map<string, int64_t> ref_reverse_index;
		const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
		for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
			auto& a_ref = a_ref_vector[ref_id];
			ref_reverse_index[a_ref.RefName] = ref_id;
		}
		auto rev_itr = ref_reverse_index.find(chr);
		if (ref_reverse_index.end() == rev_itr) {
			cout << "the chromosome name is not the same\n";
			exit(0);
		}

		string header_cl = "chr\ts\te\tsize\ttsd\tpbp\tnbp\trep.repeat\tfamily\tclass\tram\tpram\tnram\tcr\tpcr\tncr\t"
				"acr\tpacr\tnacr\tacrr\tpram_start\tpram_end\tnram_start\tnram_end\tpgene\tngene\tscore\ts2n\toi\tdesc\tconf\n";
		string header_germline = "sample\t" + header_cl;

		string header_clipped = "chr\ts\te\trep.repeat\tcpos\taligned\tcigar\ttname\tref.seq\tclipped.seq\tclipped.qual\n";

		RefRepeatIntervalTree ref_repeat_interval_tree(the_target_repeat_annot);
		RefRepeatIntervalVector stat_results;
		const int64_t chr_ref_id = rev_itr->second;
	//	cout << "[TEA.count_clipped] process paired\n";
		string cl_file = cl_prefix + "." + tmp_chr_name + ".cluster";
		string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";
		string p_mate_readname_file = cl_prefix + "." + tmp_chr_name + ".p.mate.rname";
		string n_mate_readname_file = cl_prefix + "." + tmp_chr_name + ".n.mate.rname";
		string p_clipped_filename_file = cl_prefix + "." + tmp_chr_name + ".p.clipped.fname";
		string n_clipped_filename_file = cl_prefix + "." + tmp_chr_name + ".n.clipped.fname";

		string clipped_file = cl_prefix + "." + tmp_chr_name + ".clipped";

		ofstream out_cl(cl_file, ios::binary);
		ofstream out_germline(germline_file, ios::binary);

		ofstream out_clipped(clipped_file, ios::binary);
		ofstream out_p_mate_rname(p_mate_readname_file, ios::binary);
		ofstream out_n_mate_rname(n_mate_readname_file, ios::binary);
		ofstream out_p_clipped_filename(p_clipped_filename_file, ios::binary);
		ofstream out_n_clipped_filename(n_clipped_filename_file, ios::binary);

		out_cl << header_cl;
		out_germline << header_germline;
		out_clipped << header_clipped;
		{
			for (auto& a_pair : pm_cl) {
				int64_t p_idx = a_pair.first;
				int64_t n_idx = a_pair.second;
				auto& positive_entry = p_cl[p_idx];
				auto& negative_entry = n_cl[n_idx];
				int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
				if (n_ram < options.min_ram) {
					continue;
				}
				// s column
				int64_t the_ram_boundary_start = positive_entry.start;
				// e column
				int64_t the_ram_boundary_end = negative_entry.stop + read_length;
				int64_t start_pos = the_ram_boundary_start + read_length;
				if (!positive_entry.value.pos.empty()) {
					start_pos = min(start_pos, positive_entry.value.pos.back());
				}
				int64_t end_pos = the_ram_boundary_end - read_length;
				if (!negative_entry.value.pos.empty()) {
					end_pos = max(end_pos, negative_entry.value.pos[0]);
				}
				local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
				output_clipped_stat_v(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, vannot, local_reader, the_ram_boundary_start, the_ram_boundary_end,
						positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin);
			}
		}
		RepeatClusterEntry an_empty_entry;
		RAMIntervalEntry empty_interval_entry(0, 0, an_empty_entry);
	//	cout << "[TEA.count_clipped] process positive only\n";
		{
	//		string cl_file = cl_prefix + "." + tmp_chr_name + ".p.cluster";
	//		string clipped_file = cl_prefix + "." + tmp_chr_name + ".p.clipped";

	//		ofstream out_cl(cl_file, ios::binary);
	//		ofstream out_clipped(clipped_file, ios::binary);
	//		out_cl << header_cl;
	//		out_clipped << header_clipped;
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(p_cl.size()); ++r_id) {
				if (positive_only.end() == positive_only.find(r_id)) {
					continue;
				}
				auto& positive_entry = p_cl[r_id];
				auto& negative_entry = empty_interval_entry;
				int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
				if (n_ram < options.min_ram) {
					continue;
				}
				// s column
				int64_t the_ram_boundary_start = positive_entry.start;
				// e column
				int64_t the_ram_boundary_end = positive_entry.stop + read_length + fragment_size;

				int64_t start_pos = the_ram_boundary_start + read_length;
				if (!positive_entry.value.pos.empty()) {
					start_pos = min(start_pos, positive_entry.value.pos.back());
				}
				int64_t end_pos = the_ram_boundary_end - read_length;
	//			if(!negative_entry.value.pos.empty()) {
	//				end_pos = max(end_pos, negative_entry.value.pos[0]);
	//			}

				local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
				output_clipped_stat_v(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, vannot, local_reader, the_ram_boundary_start, the_ram_boundary_end,
						positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin);
			}
		}

	//	cout << "[TEA.count_clipped] process negative only\n";
		{
	//		string cl_file = cl_prefix + "." + tmp_chr_name + ".n.cluster";
	//		string clipped_file = cl_prefix + "." + tmp_chr_name + ".n.clipped";
	//		ofstream out_cl(cl_file, ios::binary);
	//		ofstream out_clipped(clipped_file, ios::binary);
	//		out_cl << header_cl;
	//		out_clipped << header_clipped;
			for (int64_t r_id = 0; r_id < static_cast<int64_t>(n_cl.size()); ++r_id) {
				if (negative_only.end() == negative_only.find(r_id)) {
					continue;
				}
				auto& positive_entry = empty_interval_entry;
				auto& negative_entry = n_cl[r_id];
				int64_t n_ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
				if (n_ram < options.min_ram) {
					continue;
				}
				// s column
				int64_t the_ram_boundary_start = negative_entry.start - fragment_size;
				// e column
				int64_t the_ram_boundary_end = negative_entry.stop + read_length;

				int64_t start_pos = the_ram_boundary_start + read_length;
	//			if(!positive_entry.value.pos.empty()) {
	//				start_pos = min(start_pos, positive_entry.value.pos.back());
	//			}
				int64_t end_pos = the_ram_boundary_end - read_length;
				if (!negative_entry.value.pos.empty()) {
					end_pos = max(end_pos, negative_entry.value.pos[0]);
				}
				local_reader.SetRegion(chr_ref_id, start_pos, chr_ref_id, end_pos);
				output_clipped_stat_v(out_p_clipped_filename, out_n_clipped_filename, out_p_mate_rname, out_n_mate_rname, out_cl, out_germline, out_clipped, contig_dir, ref_repeat_interval_tree, stat_results, vannot, local_reader, the_ram_boundary_start, the_ram_boundary_end,
						positive_entry, negative_entry, chr, prefixed_chr, read_length, rmasker_filter_margin, gene_margin);
			}
		}
		local_reader.Close();
//		cout << checker;
}

void TEA::get_clipped_entries(vector<ClippedEntry>& clipped_entries, int64_t& max_pos_positive, int64_t& max_pos_negative, int64_t& n_positive_clipped_reads, int64_t& n_negative_clipped_reads, int64_t& n_aligned_clipped_positive, int64_t& n_aligned_clipped_negative,
		BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const int64_t read_length, const int64_t mid_point) {
	string tmp_chr_name(chr);
	if (string::npos == tmp_chr_name.find("chr")) {
		tmp_chr_name = "chr" + chr;
	}
	const bool debug = (17233802 == the_ram_boundary_start && 17234376 == the_ram_boundary_end);
//	const bool debug = false;
	BamTools::BamAlignment local_alignment_entry;
	auto& ref_vec = local_reader.GetReferenceData();
	while (local_reader.GetNextAlignment(local_alignment_entry)) {
//		if(!debug) {
//			debug = (string::npos != local_alignment_entry.Name.find("HSQ700642:208:D1D6WACXX:1:1313:4398:90371") || string::npos != local_alignment_entry.Name.find("HSQ700642:191:D18JJACXX:1:2111:17855:64424"));
//		}
		if(debug) {
			cout << "[TEA.get_clipped_entries] entry: " << BamWriter::GetSAMAlignment(local_alignment_entry, ref_vec) << "\n";
		}
		auto cigars = local_alignment_entry.CigarData;
		auto cigars_clipped_pos = cigars;
		vector<string> the_cigar;
		for (auto c : cigars) {
			the_cigar.push_back(boost::lexical_cast<string>(c.Length) + boost::lexical_cast<string>(c.Type));
		}
		string cigar_original = castle::StringUtils::join(the_cigar, "");
		if(debug) {
			cout << "[TEA.get_clipped_entries] cigar org.: " << cigar_original << "\n";
		}
		cigars.erase(remove_if(cigars.begin(), cigars.end(), [](const BamTools::CigarOp& c) {return 'D' == c.Type;}), cigars.end());

		for (int64_t c_id = 0; c_id < (static_cast<int64_t>(cigars.size()) - 1);) {
			auto& cur_cigar = cigars[c_id];
			auto& next_cigar = cigars[c_id + 1];
			if (('M' == cur_cigar.Type && 'M' == next_cigar.Type) || ('M' == cur_cigar.Type && 'I' == next_cigar.Type)) {
				cur_cigar.Length += next_cigar.Length;
				next_cigar.Length = 0;
				cigars.erase(cigars.begin() + (c_id + 1));
				continue;
			}
			++c_id;
		}
		if(debug) {
			cout << "[TEA.get_clipped_entries] cigar chg.: " << castle::StringUtils::join(the_cigar, "") << "\n";
		}

		for (int64_t c_id = 0; c_id < (static_cast<int64_t>(cigars_clipped_pos.size()) - 1);) {
			auto& cur_cigar = cigars_clipped_pos[c_id];
			auto& next_cigar = cigars_clipped_pos[c_id + 1];
			if (('M' == cur_cigar.Type && 'S' != next_cigar.Type)) {
				cur_cigar.Length += next_cigar.Length;
				next_cigar.Length = 0;
				cigars_clipped_pos.erase(cigars_clipped_pos.begin() + (c_id + 1));
				continue;
			}
			++c_id;
		}
		the_cigar.clear();
		for (auto& c : cigars) {
			the_cigar.push_back(boost::lexical_cast<string>(c.Length) + boost::lexical_cast<string>(c.Type));
		}
		string cigar_corrected = castle::StringUtils::join(the_cigar, "");
		string ref_seq;
		string clipped_seq_qual_trimmed;
		string clipped_qual_qual_trimmed;
		string clipped_seq;
		string clipped_qual;
		int64_t clipped_pos = local_alignment_entry.Position;
		int64_t clipped_pos_qual_trimmed = local_alignment_entry.Position;
		int64_t delta = 0;
		if (local_alignment_entry.IsReverseStrand()) {
			if(debug) {
				cout << "[TEA.get_clipped_entries] Rev-0\n";
			}
			int64_t n_del = 0;
			int64_t n_in = 0;
			for (auto cigar : local_alignment_entry.CigarData) {
				if ('D' == cigar.Type) {
					n_del += cigar.Length;
				} else if ('I' == cigar.Type) {
					n_in += cigar.Length;
				}
			}
			delta = n_del - n_in;
//				if('S' == cigars.front().Type) {
//					cout << "Reverse(S1): " << local_alignment_entry.QueryBases.substr(0, cigars.front().Length) << "\n";
//				} else
			if ('S' == cigars.back().Type || 'H' == cigars.back().Type) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] Rev-1\n";
				}
				int64_t match_idx = cigars.size();
				match_idx -= 2;
				if (match_idx >= 0 && 'M' == cigars[match_idx].Type) {
					if(debug) {
						cout << "[TEA.get_clipped_entries] Rev-2\n";
					}
					int64_t match_start = local_alignment_entry.QueryBases.size() - (cigars[match_idx].Length + cigars.back().Length);
//					cout << "Reverse(M2): " << local_alignment_entry.QueryBases.substr(match_start, cigars[match_idx].Length) << "\n";
					ref_seq = local_alignment_entry.QueryBases.substr(match_start, cigars[match_idx].Length);
//						clipped_pos += cigars[match_idx].Length;
				}
//					int64_t match_clipped_idx = cigars_clipped_pos.size();
//					match_clipped_idx -= 2;
//					if(match_clipped_idx >= 0 && 'M' == cigars_clipped_pos[match_clipped_idx].Type) {
//						clipped_pos += cigars_clipped_pos[match_idx].Length;
//					}
				clipped_pos += local_alignment_entry.QueryBases.size() - cigars.back().Length + delta;
				if ('S' == cigars.front().Type) {
					if(debug) {
						cout << "[TEA.get_clipped_entries] Rev-3\n";
					}
					clipped_pos -= cigars.front().Length;
				}
//				cout << "Reverse(S2): " << local_alignment_entry.QueryBases.substr(local_alignment_entry.QueryBases.size() - cigars.back().Length) << "\n";
				if('H' != cigars.back().Type) {
					clipped_seq = local_alignment_entry.QueryBases.substr(local_alignment_entry.QueryBases.size() - cigars.back().Length);
					clipped_qual = local_alignment_entry.Qualities.substr(local_alignment_entry.Qualities.size() - cigars.back().Length);

					int64_t n_low = get_number_of_low_qualities_at_begin(clipped_qual, options.qcutoff);
					++clipped_pos;
					clipped_pos_qual_trimmed = clipped_pos;
					if (n_low > 0) {
						if(debug) {
							cout << "[TEA.get_clipped_entries] Rev-4\n";
						}
	//					cout << "# low(0): " << n_low << "/" << clipped_qual << "\n";
						int64_t n_quality_trimmed = local_alignment_entry.QueryBases.size() - cigars.back().Length + n_low;
						clipped_seq_qual_trimmed = local_alignment_entry.QueryBases.substr(n_quality_trimmed);
						clipped_qual_qual_trimmed = local_alignment_entry.Qualities.substr(n_quality_trimmed);
						clipped_pos_qual_trimmed += n_low;
					}
				}
			}
		} else {
			if(debug) {
				cout << "[TEA.get_clipped_entries] For-0\n";
			}
			if ('S' == cigars.front().Type || 'H' == cigars.front().Type) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] For-1\n";
				}
				uint64_t match_idx = 1;
				if (match_idx < cigars.size() && 'M' == cigars[match_idx].Type) {
					if(debug) {
						cout << "[TEA.get_clipped_entries] For-2\n";
					}
					int64_t match_start = cigars.front().Length;
//					cout << "Forward(M3): " << local_alignment_entry.QueryBases.substr(match_start, cigars[match_idx].Length) << "\n";
					ref_seq = local_alignment_entry.QueryBases.substr(match_start, cigars[match_idx].Length);
				}
//					cout << "Forward(S3): " << local_alignment_entry.QueryBases.substr(0, cigars.front().Length) << "\n";

				if('H' != cigars.front().Type) {
					clipped_seq = local_alignment_entry.QueryBases.substr(0, cigars.front().Length);
					clipped_qual = local_alignment_entry.Qualities.substr(0, cigars.front().Length);
					int64_t n_low = get_number_of_low_qualities_at_end(clipped_qual, options.qcutoff);
					clipped_pos_qual_trimmed = clipped_pos;
					if (n_low > 0) {
						if(debug) {
							cout << "[TEA.get_clipped_entries] For-3\n";
						}
	//					cout << "# low(3): " << n_low << "\n";
						int64_t n_quality_trimmed = cigars.front().Length - n_low;
						clipped_seq_qual_trimmed = local_alignment_entry.QueryBases.substr(0, n_quality_trimmed);
						clipped_qual_qual_trimmed = local_alignment_entry.Qualities.substr(0, n_quality_trimmed);
						clipped_pos_qual_trimmed -= n_low;
					}
				}

			}
//				else if('S' == cigars.back().Type) {
//					cout << "Forward(S4): " << local_alignment_entry.QueryBases.substr(local_alignment_entry.QueryBases.size() - cigars.back().Length) << "\n";
//				}
		}
		ClippedEntry an_entry;
		an_entry.chr = tmp_chr_name;
		an_entry.ram_start = the_ram_boundary_start;
		an_entry.ram_end = the_ram_boundary_end;
		an_entry.strand = local_alignment_entry.IsReverseStrand() ? -1 : 1;
		set<string> alt_rep;
		alt_rep.insert(positive_entry.value.rep_repeat.begin(), positive_entry.value.rep_repeat.end());
		alt_rep.insert(negative_entry.value.rep_repeat.begin(), negative_entry.value.rep_repeat.end());
		an_entry.rep_repeat = castle::StringUtils::join(alt_rep, ",");
		an_entry.negative_pos = local_alignment_entry.Position;
		an_entry.clipped_pos = clipped_pos;
		an_entry.clipped_pos_qual_trimmed = clipped_pos_qual_trimmed;
		an_entry.cigar_original = cigar_original;
		an_entry.cigar_corrected = cigar_corrected;
		an_entry.read_name = local_alignment_entry.Name;
		an_entry.ref_seq = ref_seq;
		an_entry.clipped_seq = clipped_seq;
		an_entry.clipped_qual = clipped_qual;
		an_entry.clipped_seq_qual_trimmed = clipped_seq_qual_trimmed;
		an_entry.clipped_qual_qual_trimmed = clipped_qual_qual_trimmed;
		clipped_entries.push_back(an_entry);
	}

	map<int64_t, int64_t> pos_frequency_positive;
	map<int64_t, int64_t> pos_frequency_negative;
	n_positive_clipped_reads = 0;
	n_negative_clipped_reads = 0;
	for (auto& an_entry : clipped_entries) {
		if (1 == an_entry.strand) {
			++pos_frequency_positive[an_entry.clipped_pos];
			if(an_entry.clipped_pos != an_entry.clipped_pos_qual_trimmed) {
				++pos_frequency_positive[an_entry.clipped_pos_qual_trimmed];
			}
			++n_positive_clipped_reads;
		} else if (-1 == an_entry.strand) {
			++pos_frequency_negative[an_entry.clipped_pos];
			if(an_entry.clipped_pos != an_entry.clipped_pos_qual_trimmed) {
				++pos_frequency_negative[an_entry.clipped_pos_qual_trimmed];
			}
			++n_negative_clipped_reads;
		}
	}

	if(debug) {
		cout << (boost::format("[TEA.get_clipped_entries] Mid: %d\n") % mid_point).str();
		for (auto& an_entry : pos_frequency_positive) {
			cout << (boost::format("[TEA.get_clipped_entries] Positive pos: %d, freq: %d\n") % an_entry.first % an_entry.second).str();
		}
		for (auto& an_entry : pos_frequency_negative) {
			cout << (boost::format("[TEA.get_clipped_entries] Negative pos : %d, freq: %d\n") % an_entry.first % an_entry.second).str();
		}
	}
	max_pos_positive = -1;
	int64_t max_pos_freq_positive = 0;
	max_pos_negative = -1;
	int64_t max_pos_freq_negative = 0;
	for (auto& a_freq : pos_frequency_positive) {
		if (max_pos_freq_positive < a_freq.second) {
			max_pos_positive = a_freq.first;
			max_pos_freq_positive = a_freq.second;
		}
	}
	for (auto& a_freq : pos_frequency_negative) {
		if (max_pos_freq_negative < a_freq.second) {
			max_pos_negative = a_freq.first;
			max_pos_freq_negative = a_freq.second;
		}
	}

	int32_t n_positive_ties = 0;
	for (auto& a_freq : pos_frequency_positive) {
		if (max_pos_freq_positive == a_freq.second) {
			++n_positive_ties;
			if(n_positive_ties > 1) {
				break;
			}
		}
	}

	int32_t n_negative_ties = 0;
	for (auto& a_freq : pos_frequency_negative) {
		if (max_pos_freq_negative == a_freq.second) {
			++n_negative_ties;
			if(n_negative_ties > 1) {
				break;
			}
		}
	}
	if(debug) {
		cout << (boost::format("[TEA.get_clipped_entries] # pos ties: %d, # neg ties: %d\n") % n_positive_ties % n_negative_ties).str();
	}

	// the bp in forward strand is strongly supported, but not the bp in reverse strand
	if(max_pos_freq_positive > max_pos_freq_negative || 1 == n_positive_ties) {
		if(debug) {
			cout << "[TEA.get_clipped_entries] strong pos, weak neg - 1\n";
		}
		bool found_negative_candidate = false;
		// check the positive bp position which is within the right bp_margin
		int64_t freq_pos_negative = 0;
		for (auto& a_freq : pos_frequency_negative) {
			int64_t candidate_pos = a_freq.first;
//			if (max_pos_freq_negative == a_freq.second) {
				if(max_pos_positive <= candidate_pos && candidate_pos < (max_pos_positive + options.bp_margin)) {
					if(freq_pos_negative < a_freq.second) {
						max_pos_negative = candidate_pos;
						freq_pos_negative = a_freq.second;
						found_negative_candidate = true;
					}
				}
//			}
		}
		if(!found_negative_candidate) {
			if(debug) {
				cout << "[TEA.get_clipped_entries] strong pos, weak neg - 2\n";
			}
			freq_pos_negative = 0;
			// check the negative bp position which is within the left bp_margin
			for (auto& a_freq : pos_frequency_negative) {
				int64_t candidate_pos = a_freq.first;
				if(debug) {
					cout << "[TEA.get_clipped_entries] strong pos, weak neg - 2: " << candidate_pos << "\n";
				}
//				if (max_pos_freq_negative == a_freq.second) {
					if((max_pos_positive - options.bp_margin) <= candidate_pos && candidate_pos < max_pos_positive) {
						if(freq_pos_negative < a_freq.second) {
							found_negative_candidate = true;
							freq_pos_negative = a_freq.second;
							max_pos_negative = candidate_pos;
						}
					}
//				}
			}
			if(!found_negative_candidate) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] strong pos, weak neg - 3\n";
				}
				int64_t min_delta = numeric_limits<int64_t>::max();
				// the closest entry to the positive breakpoint
				for (auto& a_freq : pos_frequency_negative) {
					int64_t candidate_pos = a_freq.first;
					int64_t delta = abs(max_pos_positive - candidate_pos);
					if(delta < min_delta) {
						min_delta = delta;
						max_pos_negative = candidate_pos;
					}
				}
			}
		}
	}
	// the bp in reverse strand is strongly supported, but not the bp in forward strand
	else if(max_pos_freq_positive < max_pos_freq_negative || 1 == n_negative_ties) {
		if(debug) {
			cout << "[TEA.get_clipped_entries] strong neg, weak pos - 1\n";
		}
		bool found_positive_candidate = false;
		// check the negative bp position which is within the left bp_margin
		int64_t freq_pos_positive = 0;
		for (auto& a_freq : pos_frequency_positive) {
			int64_t candidate_pos = a_freq.first;
//			if (max_pos_freq_positive == a_freq.second) {
			if((max_pos_negative - options.bp_margin) < candidate_pos && candidate_pos < max_pos_negative) {
				if(freq_pos_positive < a_freq.second) {
					found_positive_candidate = true;
					freq_pos_positive = a_freq.second;
					max_pos_positive = candidate_pos;
				}
			}
//			}
		}

		if(!found_positive_candidate) {
			if(debug) {
				cout << "[TEA.get_clipped_entries] strong neg, weak pos - 2\n";
			}
			// check the positive bp position which is within the right bp_margin
			for (auto& a_freq : pos_frequency_positive) {
				int64_t candidate_pos = a_freq.first;
//				if (max_pos_freq_positive == a_freq.second) {
					if(max_pos_negative < candidate_pos && candidate_pos < (max_pos_negative + options.bp_margin)) {
						if(freq_pos_positive < a_freq.second) {
							found_positive_candidate = true;
							freq_pos_positive = a_freq.second;
							max_pos_positive = candidate_pos;
						}
					}
//				}
			}
			if(!found_positive_candidate) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] strong neg, weak pos - 3\n";
				}
				int64_t min_delta = numeric_limits<int64_t>::max();
				// the closest entry to the positive breakpoint
				for (auto& a_freq : pos_frequency_positive) {
					int64_t candidate_pos = a_freq.first;
					int64_t delta = abs(max_pos_negative - candidate_pos);
					if(delta < min_delta) {
						min_delta = delta;
						max_pos_positive = candidate_pos;
					}
				}
			}
		}
	}
	// neither the bp in forward strand nor in reverse strand is strongly supported
	else {
		if(debug) {
			cout << "[TEA.get_clipped_entries] weak pos, weak neg - 1\n";
		}
		vector<int64_t> positive_candidates;
		vector<int64_t> negative_candidates;
		int64_t pos_delta = numeric_limits<int64_t>::max();

		for (auto& a_freq : pos_frequency_positive) {
			int64_t candidate_pos = a_freq.first;
			int64_t cur_delta = abs(candidate_pos - mid_point);
			if (cur_delta < pos_delta) {
				pos_delta = cur_delta;
				max_pos_freq_positive = a_freq.second;
				max_pos_positive = candidate_pos;
			}
		}
		int64_t neg_delta = numeric_limits<int64_t>::max();
		for (auto& a_freq : pos_frequency_negative) {
			int64_t candidate_pos = a_freq.first;
			int64_t cur_delta = abs(candidate_pos - mid_point);
			if (cur_delta < neg_delta) {
				neg_delta = cur_delta;
				max_pos_freq_positive = a_freq.second;
				max_pos_negative = candidate_pos;
			}
		}
	}
	// if it were not caught by the above conditional statements, then both signals would be strongly supported.

	if(debug) {
		cout << (boost::format("[TEA.count_clipped] Positive pos: %d, freq: %d\n") % max_pos_positive % max_pos_freq_positive).str();
		cout << (boost::format("[TEA.count_clipped] Negative pos: %d, freq: %d\n") % max_pos_negative % max_pos_freq_negative).str();
	}

	for (auto& an_entry : clipped_entries) {
		an_entry.clipped_pos_rep = an_entry.clipped_pos;
		an_entry.clipped_seq_rep = an_entry.clipped_seq;
		an_entry.clipped_qual_rep = an_entry.clipped_qual;
		if (an_entry.clipped_pos != an_entry.clipped_pos_qual_trimmed) {
			an_entry.clipped_pos_rep = an_entry.clipped_pos_qual_trimmed;
			an_entry.clipped_seq_rep = an_entry.clipped_seq_qual_trimmed;
			an_entry.clipped_qual_rep = an_entry.clipped_qual_qual_trimmed;
		}

		if (1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos - max_pos_positive);
			int64_t delta_qual = abs(an_entry.clipped_pos_qual_trimmed - max_pos_positive);
			int64_t selected_delta = min(delta, delta_qual);
			if(selected_delta > options.jittering) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] jittering fail - 1\n";
				}
				continue;
			}
			if (max_pos_freq_positive >= 1) {
				int64_t clipped_freq = 0;
				int64_t clipped_freq_trimmed = 0;
				for (int64_t pos = an_entry.clipped_pos - options.jittering; pos < an_entry.clipped_pos + options.jittering; ++pos) {
					auto an_entry = pos_frequency_positive.find(pos);
					if (pos_frequency_positive.end() == an_entry) {
						continue;
					}
					clipped_freq += an_entry->second;
				}

				for (int64_t pos = an_entry.clipped_pos_qual_trimmed - options.jittering; pos < an_entry.clipped_pos_qual_trimmed + options.jittering; ++pos) {
					auto an_entry = pos_frequency_positive.find(pos);
					if (pos_frequency_positive.end() == an_entry) {
						continue;
					}
					clipped_freq_trimmed += an_entry->second;
				}
//				cout << "d(0): " << clipped_freq << "(" << an_entry.clipped_pos << ")" << "/" << clipped_freq_trimmed << "(" << an_entry.clipped_pos_qual_trimmed << ")" << "\n";
//				if (clipped_freq > clipped_freq_trimmed || (1 == clipped_freq && 1 == clipped_freq_trimmed)) {
//				if (clipped_freq >= clipped_freq_trimmed) {
					an_entry.clipped_pos_rep = an_entry.clipped_pos;
					an_entry.clipped_seq_rep = an_entry.clipped_seq;
					an_entry.clipped_qual_rep = an_entry.clipped_qual;
//				} else {
//					if(debug) {
//						cout << (boost::format("[TEA.get_clipped_entries] freq fail - 1: %d <=> %d\n") % clipped_freq % clipped_freq_trimmed).str();
//					}
//				}
			} else {
				if(debug) {
					cout << "[TEA.get_clipped_entries] single entry - 1\n";
				}
			}
		} else if (-1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos - max_pos_negative);
			int64_t delta_qual = abs(an_entry.clipped_pos_qual_trimmed - max_pos_negative);
			int64_t selected_delta = min(delta, delta_qual);
			if(selected_delta > options.jittering) {
				if(debug) {
					cout << "[TEA.get_clipped_entries] jittering fail - 2\n";
				}
				continue;
			}
			if (max_pos_freq_negative >= 1) {
				int64_t clipped_freq = 0;
				int64_t clipped_freq_trimmed = 0;

				for (int64_t pos = an_entry.clipped_pos - options.jittering; pos < an_entry.clipped_pos + options.jittering; ++pos) {
					auto an_entry = pos_frequency_negative.find(pos);
					if (pos_frequency_negative.end() == an_entry) {
						continue;
					}
					clipped_freq += an_entry->second;
				}

				for (int64_t pos = an_entry.clipped_pos_qual_trimmed - options.jittering; pos < an_entry.clipped_pos_qual_trimmed + options.jittering; ++pos) {
					auto an_entry = pos_frequency_negative.find(pos);
					if (pos_frequency_negative.end() == an_entry) {
						continue;
					}
					clipped_freq_trimmed += an_entry->second;
				}
//				cout << "d(1): " << clipped_freq << "(" << an_entry.clipped_pos << ")" << "/" << clipped_freq_trimmed << "(" << an_entry.clipped_pos_qual_trimmed << ")" << "\n";
//				if (clipped_freq > clipped_freq_trimmed || (1 == clipped_freq && 1 == clipped_freq_trimmed)) {
//				if (clipped_freq >= clipped_freq_trimmed) {
					an_entry.clipped_pos_rep = an_entry.clipped_pos;
					an_entry.clipped_seq_rep = an_entry.clipped_seq;
					an_entry.clipped_qual_rep = an_entry.clipped_qual;
//				} else {
//					if(debug) {
//						cout << (boost::format("[TEA.get_clipped_entries] freq fail - 2: %d <=> %d\n") % clipped_freq % clipped_freq_trimmed).str();
//					}
//				}
			} else {
				if(debug) {
					cout << "[TEA.get_clipped_entries] single entry - 2\n";
				}
			}
		}
	}
//	for (auto& an_entry : clipped_entries) {
//		if (1 == an_entry.strand) {
//			if (an_entry.clipped_pos_rep == max_pos_positive) {
//				++n_aligned_clipped_positive;
//			}
//		} else if (-1 == an_entry.strand) {
//			if (an_entry.clipped_pos_rep == max_pos_negative) {
//				++n_aligned_clipped_negative;
//			}
//		}
//	}
	n_aligned_clipped_positive = max_pos_freq_positive;
	n_aligned_clipped_negative = max_pos_freq_negative;
}

void TEA::output_clipped_stat(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_germline, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree,
		RefRepeatIntervalVector& stat_results, GeneIntervalTree& gene_interval_tree, GeneIntervalVector& gene_results, BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry,
		const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin, const int64_t mid_point) {
	vector<ClippedEntry> clipped_entries;
	int64_t max_pos_positive = 0;
	int64_t max_pos_negative = 0;
	int64_t n_positive_clipped_reads = 0;
	int64_t n_negative_clipped_reads = 0;
	int64_t n_aligned_clipped_positive = 0;
	int64_t n_aligned_clipped_negative = 0;
	get_clipped_entries(clipped_entries, max_pos_positive, max_pos_negative, n_positive_clipped_reads, n_negative_clipped_reads, n_aligned_clipped_positive, n_aligned_clipped_negative, local_reader, the_ram_boundary_start, the_ram_boundary_end, positive_entry, negative_entry, chr, read_length, mid_point);
	string tmp_chr_name(chr);
	if (string::npos == tmp_chr_name.find("chr")) {
		tmp_chr_name = "chr" + chr;
	}
	ClippedStatEntry a_stat_entry;
	a_stat_entry.ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
	a_stat_entry.chr = tmp_chr_name;
	a_stat_entry.s = the_ram_boundary_start;
	a_stat_entry.e = the_ram_boundary_end;
	a_stat_entry.size = the_ram_boundary_end - the_ram_boundary_start + 1;
	a_stat_entry.pbp = max_pos_positive;
	a_stat_entry.nbp = max_pos_negative;
	if (-1 != max_pos_positive && -1 != max_pos_negative) {
		a_stat_entry.tsd = max_pos_negative - max_pos_positive - 1;
	}

	set<string> alt_rep;
	alt_rep.insert(positive_entry.value.rep_repeat.begin(), positive_entry.value.rep_repeat.end());
	alt_rep.insert(negative_entry.value.rep_repeat.begin(), negative_entry.value.rep_repeat.end());
	a_stat_entry.rep_repeat = castle::StringUtils::join(alt_rep, ",");
	set<string> rep_family;
	rep_family.insert(positive_entry.value.family.begin(), positive_entry.value.family.end());
	rep_family.insert(negative_entry.value.family.begin(), negative_entry.value.family.end());
	auto family_itr = rep_family.find("PolyA");
	if (rep_family.end() != family_itr) {
		rep_family.erase(family_itr);
	}
	a_stat_entry.repeat_family = castle::StringUtils::join(rep_family, ",");
	a_stat_entry.rep_suffix = *rep_family.begin();
	set<string> rep_class;
	rep_class.insert(positive_entry.value.repeat_class.begin(), positive_entry.value.repeat_class.end());
	rep_class.insert(negative_entry.value.repeat_class.begin(), negative_entry.value.repeat_class.end());
	auto class_itr = rep_class.find("PolyA");
	if (rep_class.end() != class_itr) {
		rep_class.erase(class_itr);
	}
	a_stat_entry.repeat_class = castle::StringUtils::join(rep_class, ",");

	a_stat_entry.pram = positive_entry.value.pos.size();
	a_stat_entry.nram = negative_entry.value.pos.size();
	a_stat_entry.cr = n_positive_clipped_reads + n_negative_clipped_reads;
	a_stat_entry.pcr = n_positive_clipped_reads;
	a_stat_entry.ncr = n_negative_clipped_reads;
	a_stat_entry.acr = n_aligned_clipped_positive + n_aligned_clipped_negative;
	a_stat_entry.pacr = n_aligned_clipped_positive;
	a_stat_entry.nacr = n_aligned_clipped_negative;
	if (0 != a_stat_entry.cr) {
		a_stat_entry.acrr = a_stat_entry.acr / (double) a_stat_entry.cr;
	}
	if (!positive_entry.value.pos.empty()) {
		a_stat_entry.pram_start = positive_entry.value.pos[0];
		a_stat_entry.pram_end = positive_entry.value.pos.back();
	}
	if (!negative_entry.value.pos.empty()) {
		a_stat_entry.nram_start = negative_entry.value.pos[0];
		a_stat_entry.nram_end = negative_entry.value.pos.back();
	}
	{
		int64_t max_ram = 20;
		int64_t max_cr = 10;
		double s_r1 = min(a_stat_entry.pram, max_ram) / (double) max_ram;
		double s_r2 = min(a_stat_entry.nram, max_ram) / (double) max_ram;
		double s_cl1 = min(a_stat_entry.pacr, max_cr) / (double) max_cr;
		double s_cl2 = min(a_stat_entry.nacr, max_cr) / (double) max_cr;
		a_stat_entry.score = s_r1 + s_r2 + s_cl1 + s_cl2;
	}
	{
		int64_t n_pram_down_nbp = 0;
		int64_t n_pram_up_nbp = 0;
		int64_t n_nram_up_pbp = 0;
		int64_t n_nram_down_pbp = 0;
		for (auto a_pram_pos : positive_entry.value.pos) {
			if (a_pram_pos <= a_stat_entry.nbp) {
				++n_pram_down_nbp;
			} else if (a_pram_pos > a_stat_entry.nbp) {
				++n_pram_up_nbp;
			}
		}
		for (auto an_nram_pos : negative_entry.value.pos) {
			if (an_nram_pos >= a_stat_entry.pbp) {
				++n_nram_up_pbp;
			} else if (an_nram_pos < a_stat_entry.pbp) {
				++n_nram_down_pbp;
			}
		}
		a_stat_entry.s2n = 2 * sqrt(n_pram_down_nbp * n_nram_up_pbp) - (n_pram_up_nbp + n_nram_down_pbp);
		bool clipped = false;
		if (clipped) {
			int64_t pacr = a_stat_entry.pacr;
			int64_t nacr = a_stat_entry.nacr;
			int64_t p_non_acr = a_stat_entry.pcr - pacr;
			int64_t n_non_acr = a_stat_entry.ncr - nacr;
			double s2n_clipped = 2 * sqrt(pacr * nacr) - (p_non_acr + n_non_acr);
			if (s2n_clipped < 0) {
				s2n_clipped = 0;
			}

			a_stat_entry.s2n += 3 * s2n_clipped;
		}
	}
	bool debug = false;
	for (auto& an_entry : clipped_entries) {
//		if("chr1" == an_entry.chr && 47262411 == an_entry.ram_start && 47263181 == an_entry.ram_end) {
//			debug = true;
//		}
		if (1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos_rep - a_stat_entry.pbp);
			if (delta <= options.jittering) {
				an_entry.aligned = 1;
			}
//				df1$cpos[i] - ccnt$pbp
		} else if (-1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos_rep - a_stat_entry.nbp);
			if (delta <= options.jittering) {
				an_entry.aligned = 1;
			}
		}
		out_clipped << an_entry.str() << "\n";
	}
	int64_t s = 0;
	int64_t e = 0;
	if (a_stat_entry.pbp > a_stat_entry.nbp) {
		s = a_stat_entry.nbp;
		e = a_stat_entry.pbp;
	} else {
		e = a_stat_entry.nbp;
		s = a_stat_entry.pbp;
	}

	if (-1 == s || -1 == e) {
		s = e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}
	s -= rmasker_filter_margin;
	e += rmasker_filter_margin;

//		cout << "s: " << s << ", e: " << e << "\n";
	stat_results.clear();
	ref_repeat_interval_tree.find_overlap(s, e, stat_results);
	bool found_class = false;
	for (auto& a_stat : stat_results) {
		if (a_stat.value.chromosome != chr && a_stat.value.chromosome != prefixed_chr) {
			continue;
		}

		auto rep_itr = alt_rep.find(a_stat.value.repeat_name);
		if (alt_rep.end() != rep_itr) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] found class: " << a_stat.value.repeat_name << "\n";
			}

			found_class = true;
			break;
		}
	}
//		cout << "# matched repeat chromosome: " << stat_results.size() << "\n";
	a_stat_entry.oi = found_class ? 0 : 1;

	int64_t gene_positive_s = a_stat_entry.pbp;
	int64_t gene_positive_e = a_stat_entry.pbp;
	int64_t gene_negative_s = a_stat_entry.nbp;
	int64_t gene_negative_e = a_stat_entry.nbp;
	if (-1 == gene_positive_s) {
		gene_positive_s = gene_positive_e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}
	if (-1 == gene_negative_s) {
		gene_negative_s = gene_negative_e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}

	gene_positive_s -= gene_margin;
	gene_positive_e += gene_margin;
	gene_negative_s -= gene_margin;
	gene_negative_e += gene_margin;

	gene_results.clear();
	if(debug) {
		cout << "[TEA.output_clipped_stat] gene_positive_s: " << gene_positive_s << ", gene_positive_e: " << gene_positive_e << "\n";
		cout << "[TEA.output_clipped_stat] gene_negative_s: " << gene_negative_s << ", gene_negative_e: " << gene_negative_e << "\n";
	}
	gene_interval_tree.find_overlap(gene_positive_s, gene_positive_e, gene_results);
	for (auto& a_gene_positive_entry : gene_results) {
		a_stat_entry.pgene = (boost::format("%s_%s") % a_gene_positive_entry.value.name % a_gene_positive_entry.value.type).str();
		if(debug) {
			cout << "[TEA.output_clipped_stat] positive gene: " << a_stat_entry.pgene << "\n";
		}
		break;
	}
	gene_results.clear();
	gene_interval_tree.find_overlap(gene_negative_s, gene_negative_e, gene_results);
	for (auto& a_gene_negative_entry : gene_results) {
		a_stat_entry.ngene = (boost::format("%s_%s") % a_gene_negative_entry.value.name % a_gene_negative_entry.value.type).str();
		if(debug) {
			cout << "[TEA.output_clipped_stat] negative gene: " << a_stat_entry.ngene << "\n";
		}
		break;
	}
	const int64_t half_min_acr = options.min_acr / (double) 2;
	if (!positive_entry.value.pos.empty() && !negative_entry.value.pos.empty()) {
		if (a_stat_entry.acr >= options.min_acr && a_stat_entry.pacr >= half_min_acr && a_stat_entry.nacr >= half_min_acr) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] conf 2-1\n";
			}
			a_stat_entry.conf = 2;
		}
	} else {
		if (a_stat_entry.acr >= options.min_acr) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] conf 2-2\n";
			}
			a_stat_entry.conf = 2;
		}
	}
	if (2 == a_stat_entry.conf) {
		if (a_stat_entry.acrr >= options.min_acrr) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] conf 3\n";
			}
			a_stat_entry.conf = 3;
		}
	}
	if (3 == a_stat_entry.conf) {
		if (a_stat_entry.tsd >= options.min_tsd && a_stat_entry.tsd <= options.max_tsd) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] conf 4\n";
			}
			a_stat_entry.conf = 4;
		}
	}

	if (4 == a_stat_entry.conf) {
		if (a_stat_entry.ram >= options.ram_cutoff) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] conf 5\n";
			}
			a_stat_entry.conf = 5;
		}
	}
	if(options.no_oi && 0 == a_stat_entry.oi) {
		if(debug) {
			cout << "[TEA.output_clipped_stat] conf 0 due to oi == 0\n";
		}
		a_stat_entry.conf = 0;
	}
	bool has_written_germline = false;
	if(options.oneside_ram) {
		if (a_stat_entry.conf >= options.min_out_conf) {
			if(debug) {
				cout << "[TEA.output_clipped_stat] write to .germline (oneside-ram)\n";
			}
			out_germline << options.naive_prefix << "\t" << a_stat_entry.str() << "\n";
			has_written_germline = true;
		}
	} else {
		if (!options.oneside_ram && !positive_entry.value.pos.empty() && !negative_entry.value.pos.empty()) {
			if (a_stat_entry.conf >= options.min_out_conf) {
				if(debug) {
					cout << "[TEA.output_clipped_stat] write to .germline (not oneside-ram)\n";
				}
				out_germline << options.naive_prefix << "\t" << a_stat_entry.str() << "\n";
				has_written_germline = true;
			}
		}
	}
	if(has_written_germline) {
		int64_t n_valid_entries = 0;
		for (auto& an_entry : clipped_entries) {
			if (1 != an_entry.aligned || 1 != an_entry.strand) {
				continue;
			}
			++n_valid_entries;
		}
		if (n_valid_entries > 0) {
			string fq_prefix = (boost::format("%s.%s.%s.%s.%s") % options.naive_prefix % tmp_chr_name % a_stat_entry.s % a_stat_entry.e % a_stat_entry.rep_suffix).str();
			string fq_name = (boost::format("%s/%s.pos.fa") % contig_dir % fq_prefix).str();
			out_p_clipped_filename << fq_prefix << "\n";
			ofstream out_fq(fq_name, ios::binary);
			int64_t n_cnt = 0;
			for (auto& an_entry : clipped_entries) {
				if (1 != an_entry.aligned || 1 != an_entry.strand) {
					continue;
				}
				out_fq << (boost::format(">cr%d\n%s\n") % n_cnt % an_entry.clipped_seq).str();
				++n_cnt;
			}
			auto& positions = positive_entry.value.pos;
			auto& read_names = positive_entry.value.rname;
			int64_t max_pos = min(positions.size(), read_names.size());
			for (int64_t p_id = 0; p_id < max_pos; ++p_id) {
				out_p_mate_rname << a_stat_entry.s << "\t" << a_stat_entry.e << "\t" << a_stat_entry.rep_suffix << "\t" << positions[p_id] << "\t" << read_names[p_id] << "\n";
			}
		}
	}
	if(has_written_germline) {
		int64_t n_valid_entries = 0;
		for (auto& an_entry : clipped_entries) {
			if (1 != an_entry.aligned || -1 != an_entry.strand) {
				continue;
			}
			++n_valid_entries;
		}
		if (n_valid_entries > 0) {
			string fq_prefix = (boost::format("%s.%s.%s.%s.%s") % options.naive_prefix % tmp_chr_name % a_stat_entry.s % a_stat_entry.e % a_stat_entry.rep_suffix).str();
			string fq_name = (boost::format("%s/%s.neg.fa") % contig_dir % fq_prefix).str();
			out_n_clipped_filename << fq_prefix << "\n";
			ofstream out_fq(fq_name, ios::binary);
			int64_t n_cnt = 0;
			for (auto& an_entry : clipped_entries) {
				if (1 != an_entry.aligned || -1 != an_entry.strand) {
					continue;
				}
				out_fq << (boost::format(">cr%d\n%s\n") % n_cnt % an_entry.clipped_seq).str();
				++n_cnt;
			}
			auto& positions = negative_entry.value.pos;
			auto& read_names = negative_entry.value.rname;
			int64_t max_pos = min(positions.size(), read_names.size());
			for (int64_t p_id = 0; p_id < max_pos; ++p_id) {
				out_n_mate_rname << a_stat_entry.s << "\t" << a_stat_entry.e << "\t" << a_stat_entry.rep_suffix << "\t" << positions[p_id] << "\t" << read_names[p_id] << "\n";
			}
		}
	}

	out_cl << a_stat_entry.str() << "\n";
}

void TEA::output_clipped_stat_v(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_germline, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree, RefRepeatIntervalVector& stat_results, const map<int64_t, string>& vannot,
				BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin) {
	vector<ClippedEntry> clipped_entries;
	int64_t max_pos_positive = 0;
	int64_t max_pos_negative = 0;
	int64_t n_positive_clipped_reads = 0;
	int64_t n_negative_clipped_reads = 0;
	int64_t n_aligned_clipped_positive = 0;
	int64_t n_aligned_clipped_negative = 0;
	get_clipped_entries(clipped_entries, max_pos_positive, max_pos_negative, n_positive_clipped_reads, n_negative_clipped_reads, n_aligned_clipped_positive, n_aligned_clipped_negative, local_reader, the_ram_boundary_start, the_ram_boundary_end, positive_entry, negative_entry, chr, read_length, 0);
	string tmp_chr_name(chr);
	if (string::npos == tmp_chr_name.find("chr")) {
		tmp_chr_name = "chr" + chr;
	}
	ClippedStatEntry a_stat_entry;
	a_stat_entry.ram = positive_entry.value.pos.size() + negative_entry.value.pos.size();
	a_stat_entry.chr = tmp_chr_name;
	a_stat_entry.s = the_ram_boundary_start;
	a_stat_entry.e = the_ram_boundary_end;
	a_stat_entry.size = the_ram_boundary_end - the_ram_boundary_start + 1;
	a_stat_entry.pbp = max_pos_positive;
	a_stat_entry.nbp = max_pos_negative;
	if (-1 != max_pos_positive && -1 != max_pos_negative) {
		a_stat_entry.tsd = max_pos_negative - max_pos_positive - 1;
	}

	set<string> alt_rep;
	alt_rep.insert(positive_entry.value.rep_repeat.begin(), positive_entry.value.rep_repeat.end());
	alt_rep.insert(negative_entry.value.rep_repeat.begin(), negative_entry.value.rep_repeat.end());
	a_stat_entry.rep_repeat = castle::StringUtils::join(alt_rep, ",");
	set<string> rep_family;
	rep_family.insert(positive_entry.value.family.begin(), positive_entry.value.family.end());
	rep_family.insert(negative_entry.value.family.begin(), negative_entry.value.family.end());
	auto family_itr = rep_family.find("PolyA");
	if (rep_family.end() != family_itr) {
		rep_family.erase(family_itr);
	}
	a_stat_entry.repeat_family = castle::StringUtils::join(rep_family, ",");
	set<string> rep_class;
	rep_class.insert(positive_entry.value.repeat_class.begin(), positive_entry.value.repeat_class.end());
	rep_class.insert(negative_entry.value.repeat_class.begin(), negative_entry.value.repeat_class.end());
	auto class_itr = rep_class.find("PolyA");
	if (rep_class.end() != class_itr) {
		rep_class.erase(class_itr);
	}
	if(a_stat_entry.repeat_family.empty()) {
		a_stat_entry.repeat_family = a_stat_entry.rep_repeat;
	}
	a_stat_entry.repeat_class = castle::StringUtils::join(rep_class, ",");
	if(a_stat_entry.repeat_class.empty()) {
		a_stat_entry.repeat_class = "-";
	}
	a_stat_entry.pram = positive_entry.value.pos.size();
	a_stat_entry.nram = negative_entry.value.pos.size();
	a_stat_entry.cr = n_positive_clipped_reads + n_negative_clipped_reads;
	a_stat_entry.pcr = n_positive_clipped_reads;
	a_stat_entry.ncr = n_negative_clipped_reads;
	a_stat_entry.acr = n_aligned_clipped_positive + n_aligned_clipped_negative;
	a_stat_entry.pacr = n_aligned_clipped_positive;
	a_stat_entry.nacr = n_aligned_clipped_negative;
	if (0 != a_stat_entry.cr) {
		a_stat_entry.acrr = a_stat_entry.acr / (double) a_stat_entry.cr;
	}
	if (!positive_entry.value.pos.empty()) {
		a_stat_entry.pram_start = positive_entry.value.pos[0];
		a_stat_entry.pram_end = positive_entry.value.pos.back();
	}
	if (!negative_entry.value.pos.empty()) {
		a_stat_entry.nram_start = negative_entry.value.pos[0];
		a_stat_entry.nram_end = negative_entry.value.pos.back();
	}
	{
		int64_t max_ram = 20;
		int64_t max_cr = 10;
		double s_r1 = min(a_stat_entry.pram, max_ram) / (double) max_ram;
		double s_r2 = min(a_stat_entry.nram, max_ram) / (double) max_ram;
		double s_cl1 = min(a_stat_entry.pacr, max_cr) / (double) max_cr;
		double s_cl2 = min(a_stat_entry.nacr, max_cr) / (double) max_cr;
		a_stat_entry.score = s_r1 + s_r2 + s_cl1 + s_cl2;
	}
	{
		int64_t n_pram_down_nbp = 0;
		int64_t n_pram_up_nbp = 0;
		int64_t n_nram_up_pbp = 0;
		int64_t n_nram_down_pbp = 0;
		for (auto a_pram_pos : positive_entry.value.pos) {
			if (a_pram_pos <= a_stat_entry.nbp) {
				++n_pram_down_nbp;
			} else if (a_pram_pos > a_stat_entry.nbp) {
				++n_pram_up_nbp;
			}
		}
		for (auto an_nram_pos : negative_entry.value.pos) {
			if (an_nram_pos >= a_stat_entry.pbp) {
				++n_nram_up_pbp;
			} else if (an_nram_pos < a_stat_entry.pbp) {
				++n_nram_down_pbp;
			}
		}
		a_stat_entry.s2n = 2 * sqrt(n_pram_down_nbp * n_nram_up_pbp) - (n_pram_up_nbp + n_nram_down_pbp);
		bool clipped = false;
		if (clipped) {
			int64_t pacr = a_stat_entry.pacr;
			int64_t nacr = a_stat_entry.nacr;
			int64_t p_non_acr = a_stat_entry.pcr - pacr;
			int64_t n_non_acr = a_stat_entry.ncr - nacr;
			double s2n_clipped = 2 * sqrt(pacr * nacr) - (p_non_acr + n_non_acr);
			if (s2n_clipped < 0) {
				s2n_clipped = 0;
			}

			a_stat_entry.s2n += 3 * s2n_clipped;
		}
	}
	for (auto& an_entry : clipped_entries) {
		if (1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos_rep - a_stat_entry.pbp);
			if (delta <= options.jittering) {
				an_entry.aligned = 1;
			}
//				df1$cpos[i] - ccnt$pbp
		} else if (-1 == an_entry.strand) {
			int64_t delta = abs(an_entry.clipped_pos_rep - a_stat_entry.nbp);
			if (delta <= options.jittering) {
				an_entry.aligned = 1;
			}
		}
		out_clipped << an_entry.str() << "\n";
	}
	int64_t s = 0;
	int64_t e = 0;
	if (a_stat_entry.pbp > a_stat_entry.nbp) {
		s = a_stat_entry.nbp;
		e = a_stat_entry.pbp;
	} else {
		e = a_stat_entry.nbp;
		s = a_stat_entry.pbp;
	}

	if (-1 == s || -1 == e) {
		s = e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}
	s -= rmasker_filter_margin;
	e += rmasker_filter_margin;

//		cout << "s: " << s << ", e: " << e << "\n";
	stat_results.clear();
	ref_repeat_interval_tree.find_overlap(s, e, stat_results);
	bool found_class = false;
	for (auto& a_stat : stat_results) {
		if (a_stat.value.chromosome != chr && a_stat.value.chromosome != prefixed_chr) {
			continue;
		}

		auto rep_itr = alt_rep.find(a_stat.value.repeat_name);
		if (alt_rep.end() != rep_itr) {
			found_class = true;
			break;
		}
	}
//		cout << "# matched repeat chromosome: " << stat_results.size() << "\n";
	a_stat_entry.oi = found_class ? 0 : 1;

	int64_t gene_positive_s = a_stat_entry.pbp;
	int64_t gene_positive_e = a_stat_entry.pbp;
	int64_t gene_negative_s = a_stat_entry.nbp;
	int64_t gene_negative_e = a_stat_entry.nbp;
	if (-1 == gene_positive_s) {
		gene_positive_s = gene_positive_e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}
	if (-1 == gene_negative_s) {
		gene_negative_s = gene_negative_e = (a_stat_entry.s + a_stat_entry.e) / (double) 2;
	}

	gene_positive_s -= gene_margin;
	gene_positive_e += gene_margin;
	gene_negative_s -= gene_margin;
	gene_negative_e += gene_margin;

	const string delim("_&_");
	vector<string> the_viruses;
	castle::StringUtils::tokenize(a_stat_entry.rep_repeat, delim.c_str(), the_viruses);

	for(auto& a_virus: the_viruses) {
		auto an_ncbi_itr = vannot.find(boost::lexical_cast<int64_t>(a_virus));
		if(vannot.end() == an_ncbi_itr) {
			continue;
		}
		a_stat_entry.desc = an_ncbi_itr->second;
		break;
	}


	const int64_t half_min_acr = options.min_acr / (double) 2;
	if (!positive_entry.value.pos.empty() && !negative_entry.value.pos.empty()) {
		if (a_stat_entry.acr >= options.min_acr && a_stat_entry.pacr >= half_min_acr && a_stat_entry.nacr >= half_min_acr) {
			a_stat_entry.conf = 2;
		}
	} else {
		if (a_stat_entry.acr >= options.min_acr) {
			a_stat_entry.conf = 2;
		}
	}
	if (2 == a_stat_entry.conf) {
		if (a_stat_entry.acrr >= options.min_acrr) {
			a_stat_entry.conf = 3;
		}
	}
	if (3 == a_stat_entry.conf) {
		if (a_stat_entry.tsd >= options.min_tsd && a_stat_entry.tsd <= options.max_tsd) {
			a_stat_entry.conf = 4;
		}
	}

	if (4 == a_stat_entry.conf) {
		if (a_stat_entry.ram >= options.ram_cutoff) {
			a_stat_entry.conf = 5;
		}
	}
	if(!a_stat_entry.desc.empty()) {
		out_germline << options.naive_prefix << "\t" << a_stat_entry.str() << "\n";
	}
	{
		int64_t n_valid_entries = 0;
		for (auto& an_entry : clipped_entries) {
			if (1 != an_entry.aligned || 1 != an_entry.strand) {
				continue;
			}
			++n_valid_entries;
		}
		if (n_valid_entries > 0) {
			string fq_prefix = (boost::format("%s.%s.%s.%s.%s") % options.naive_prefix % tmp_chr_name % a_stat_entry.s % a_stat_entry.e % a_stat_entry.rep_repeat).str();
			string fq_name = (boost::format("%s/%s.pos.fa") % contig_dir % fq_prefix).str();
			out_p_clipped_filename << fq_prefix << "\n";
			ofstream out_fq(fq_name, ios::binary);
			int64_t n_cnt = 0;
			for (auto& an_entry : clipped_entries) {
				if (1 != an_entry.aligned || 1 != an_entry.strand) {
					continue;
				}
				out_fq << (boost::format(">cr%d\n%s\n") % n_cnt % an_entry.clipped_seq).str();
				++n_cnt;
			}
			auto& positions = positive_entry.value.pos;
			auto& read_names = positive_entry.value.rname;
			int64_t max_pos = min(positions.size(), read_names.size());
			for (int64_t p_id = 0; p_id < max_pos; ++p_id) {
				out_p_mate_rname << a_stat_entry.s << "\t" << a_stat_entry.e << "\t" << a_stat_entry.rep_repeat << "\t" << positions[p_id] << "\t" << read_names[p_id] << "\n";
			}
		}
	}
	{
		int64_t n_valid_entries = 0;
		for (auto& an_entry : clipped_entries) {
			if (1 != an_entry.aligned || -1 != an_entry.strand) {
				continue;
			}
			++n_valid_entries;
		}
		if (n_valid_entries > 0) {
			string fq_prefix = (boost::format("%s.%s.%s.%s.%s") % options.naive_prefix % tmp_chr_name % a_stat_entry.s % a_stat_entry.e % a_stat_entry.rep_repeat).str();
			string fq_name = (boost::format("%s/%s.neg.fa") % contig_dir % fq_prefix).str();
			out_n_clipped_filename << fq_prefix << "\n";
			ofstream out_fq(fq_name, ios::binary);
			int64_t n_cnt = 0;
			for (auto& an_entry : clipped_entries) {
				if (1 != an_entry.aligned || -1 != an_entry.strand) {
					continue;
				}
				out_fq << (boost::format(">cr%d\n%s\n") % n_cnt % an_entry.clipped_seq).str();
				++n_cnt;
			}
			auto& positions = negative_entry.value.pos;
			auto& read_names = negative_entry.value.rname;
			int64_t max_pos = min(positions.size(), read_names.size());
			for (int64_t p_id = 0; p_id < max_pos; ++p_id) {
				out_n_mate_rname << a_stat_entry.s << "\t" << a_stat_entry.e << "\t" << a_stat_entry.rep_repeat << "\t" << positions[p_id] << "\t" << read_names[p_id] << "\n";
			}
		}
	}

	out_cl << a_stat_entry.str() << "\n";

}

void TEA::output_mate_fa(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.output_mate_fa");
	checker.start();
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string contig_dir = options.prefix + "/assembly_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		contig_dir = options.working_prefix + "/assembly_" + options.rasym + "m";
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;

	boost::unordered_map<string, AlnPairEntry> a_positive_repeat_map;
	boost::unordered_map<string, AlnPairEntry> a_negative_repeat_map;
	vector<string> a_positive_clipped_prefixes;
	vector<string> a_negative_clipped_prefixes;

	vector<string> chr_names;
	for(auto& c_entry : ram) {
		auto c = c_entry.first;
		chr_names.push_back(c);
	}
	sort(chr_names.begin(), chr_names.end(), [&](const string& lhs, const string& rhs)->bool{
		string local_lhs = lhs;
		string local_rhs = rhs;
		boost::replace_all(local_lhs, "chr", "");
		boost::replace_all(local_rhs, "chr", "");
		int64_t lhs_v = numeric_limits<int64_t>::max();
		int64_t rhs_v = lhs_v;
		try {
			lhs_v = boost::lexical_cast<int64_t>(local_lhs);
		} catch(exception& ex) {}
		try {
			rhs_v = boost::lexical_cast<int64_t>(local_rhs);
		} catch(exception& ex) {}
		if(lhs_v < rhs_v) {
			return true;
		} else if(lhs_v > rhs_v) {
			return false;
		}
		return local_lhs < local_rhs;
	});

	vector<function<void()> > tasks;

	tasks.push_back([&] {
//		string c = "22";
		for(auto& c : chr_names) {
//			auto& c = c_entry.first;


			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			vector<string> data;
			string line;
			string clipped_file = cl_prefix + "." + tmp_chr_name + ".p.mate.rname";
			ifstream in_clipped(clipped_file, ios::binary);
			// to ignore the first line
//			if(!headless) {
//				getline(in_clipped, line);
//			}
			while(getline(in_clipped, line, '\n')) {
				castle::StringUtils::tokenize(line, delim_tab, data);
				string a_key = data[4];
				int64_t a_pos = boost::lexical_cast<int64_t>(data[3]);
				string file_name_prefix = naive_prefix + "." + tmp_chr_name + "." + data[0] + "." + data[1] + "." + data[2];
//				const bool debug = ("chr14" == tmp_chr_name) && ("34139887" == data[0]) && ("34140417" == data[1]);
//				if(debug) {
//					cout << "[TEA.output_mate_fa] positive mate: " << a_key << ":" << data[2] << "\n";
//				}
				a_positive_repeat_map[a_key].file_name_prefix = file_name_prefix;
				a_positive_repeat_map[a_key].pos = a_pos;
				a_positive_repeat_map[a_key].chr = c;
			}
		}
	});

	tasks.push_back([&] {
//		string c = "22";
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			vector<string> data;
			string line;
			string clipped_file = cl_prefix + "." + tmp_chr_name + ".n.mate.rname";
			ifstream in_clipped(clipped_file, ios::binary);
			// to ignore the first line
//			getline(in_clipped, line);
			while(getline(in_clipped, line, '\n')) {
				castle::StringUtils::tokenize(line, delim_tab, data);
				string a_key = data[4];
				int64_t a_pos = boost::lexical_cast<int64_t>(data[3]);
				string file_name_prefix = naive_prefix + "." + tmp_chr_name + "." + data[0] + "." + data[1] + "." + data[2];
//				const bool debug = ("chr14" == tmp_chr_name) && ("34139887" == data[0]) && ("34140417" == data[1]);
//				if(debug) {
//					cout << "[TEA.output_mate_fa] negative mate: " << a_key << ":" << data[2] << "\n";
//				}
				a_negative_repeat_map[a_key].file_name_prefix = file_name_prefix;
				a_negative_repeat_map[a_key].pos = a_pos;
				a_negative_repeat_map[a_key].chr = c;
			}
		}
	});
	tasks.push_back([&] {
		//		string c = "22";
		string line;
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			string p_mate_clipped_file = cl_prefix + "." + tmp_chr_name + ".p.clipped.fname";
			ifstream in_clipped(p_mate_clipped_file, ios::binary);
			while(getline(in_clipped, line, '\n')) {
				a_positive_clipped_prefixes.push_back(line);
			}
		}
	});
	tasks.push_back([&] {
		//		string c = "22";
		string line;
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			string n_mate_clipped_file = cl_prefix + "." + tmp_chr_name + ".n.clipped.fname";
			ifstream in_clipped(n_mate_clipped_file, ios::binary);
			while(getline(in_clipped, line, '\n')) {
				a_negative_clipped_prefixes.push_back(line);
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[TEA.output_mate_fa] # positive mates: " << a_positive_repeat_map.size() << ", # negative mates: " << a_negative_repeat_map.size() << "\n";

	string disc_file_name = options.prefix + ".disc.num.bam";

	string cl_disc_file_name = options.prefix + ".cl.sorted.disc.num.bam";

	if("tea_v" == options.program_name) {
		disc_file_name = options.prefix + ".cl.sorted.bam";
	}

	if (!options.working_dir.empty()) {
		disc_file_name = options.working_prefix + ".disc.num.bam";
		if("tea_v" == options.program_name) {
			disc_file_name = options.working_prefix + ".cl.sorted.bam";
		}
		cl_disc_file_name = options.working_prefix + ".cl.sorted.disc.num.bam";
	}
	string disc_bai_name = disc_file_name + ".bai";
	string disc_bni_name = disc_file_name + ".bni";
	string cl_disc_bai_name = cl_disc_file_name + ".bai";
	string cl_disc_bni_name = cl_disc_file_name + ".bni";

	boost::unordered_map<string, vector<string>> positive_mate_reads;
	boost::unordered_map<string, vector<string>> negative_mate_reads;

	int64_t size_block = 8192000;
//	vector<meerkat::BlockBoundary> fixed_size_blocks;

	cout << "[TEA.output_mate_fa] collect mate reads from disc file\n";
//	_output_mate_fa_serial(positive_mate_reads, negative_mate_reads, disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, disc_file_name, disc_bai_name, disc_bni_name, size_block);
//	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	collect_boundaries(fixed_size_blocks, disc_file_name, size_block);
	_output_mate_fa(positive_mate_reads, negative_mate_reads, local_unmapped_included_blocks, disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
	bool found_clipped_reads = false;
	for(auto& an_entry : a_positive_repeat_map) {
		auto& read_name = an_entry.first;
		if(string::npos != read_name.rfind("mu1") || string::npos != read_name.rfind("mu2") || string::npos != read_name.rfind("sc")) {
			found_clipped_reads = true;
			break;
		}
	}
	if(!found_clipped_reads) {
		for(auto& an_entry : a_negative_repeat_map) {
			auto& read_name = an_entry.first;
			if(string::npos != read_name.rfind("mu1") || string::npos != read_name.rfind("mu2") || string::npos != read_name.rfind("sc")) {
				found_clipped_reads = true;
				break;
			}
		}
	}

	if(found_clipped_reads) {
		cout << "[TEA.output_mate_fa] collect mate reads from cl disc file\n";
//		_output_mate_fa_serial(positive_mate_reads, negative_mate_reads, cl_disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
		vector<meerkat::BlockBoundary> local_fixed_size_blocks;
		vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
		vector<meerkat::BlockBoundary> local_independent_blocks;
		collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, cl_disc_file_name, cl_disc_bai_name, cl_disc_bni_name, size_block);
//		collect_boundaries(fixed_size_blocks, cl_disc_file_name, size_block);
		_output_mate_fa(positive_mate_reads, negative_mate_reads, local_unmapped_included_blocks, cl_disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
	}

	for(auto& an_entry : positive_mate_reads) {
		tasks.push_back([&] {
			auto& file_name_prefix = an_entry.first;
			auto& the_seq_vec = an_entry.second;
			string fa_name = (boost::format("%s/%s.rpos.fa") % contig_dir % file_name_prefix).str();
			ofstream out_fa(fa_name, ios::binary);
			for(int64_t n_id = 0; n_id < static_cast<int64_t>(the_seq_vec.size()); ++n_id) {
				auto& a_seq = the_seq_vec[n_id];
				out_fa << (boost::format(">cr%d\n%s\n") % n_id % a_seq).str();
			}
		});
	}

	for(auto& an_entry : negative_mate_reads) {
		tasks.push_back([&] {
			auto& file_name_prefix = an_entry.first;
			auto& the_seq_vec = an_entry.second;
			string fa_name = (boost::format("%s/%s.rneg.fa") % contig_dir % file_name_prefix).str();
			ofstream out_fa(fa_name, ios::binary);
			for(int64_t n_id = 0; n_id < static_cast<int64_t>(the_seq_vec.size()); ++n_id) {
				auto& a_seq = the_seq_vec[n_id];
				out_fa << (boost::format(">cr%d\n%s\n") % n_id % a_seq).str();
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[TEA.output_mate_fa] start assembling\n";
	string all_assembly_files = options.prefix + ".assemblies.target";
	if (!options.working_dir.empty()) {
		all_assembly_files = options.working_prefix + ".assemblies.target";
	}
	{
		ofstream all_files(all_assembly_files, ios::binary);
		cout << "[TEA.output_mate_fa] adds positive FASTA files\n";
		for(auto& an_entry : positive_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string pos_fa_name = (boost::format("%s/%s.pos.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(pos_fa_name)) {
					continue;
				}
				all_files << pos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % pos_fa_name % options.assembler_param % pos_fa_name).str();
	//			all_cmds << assembler_cmd << "\n";
	//			system(assembler_cmd.c_str());
	//		});
		}
		cout << "[TEA.output_mate_fa] adds reverse positive FASTA files\n";
		for(auto& an_entry : positive_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string rpos_fa_name = (boost::format("%s/%s.rpos.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(rpos_fa_name)) {
					continue;
				}
				all_files << rpos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % rpos_fa_name % options.assembler_param % rpos_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		cout << "[TEA.output_mate_fa] adds negative FASTA files\n";
		for(auto& an_entry : negative_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string neg_fa_name = (boost::format("%s/%s.neg.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(neg_fa_name)) {
					continue;
				}
				all_files << neg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % neg_fa_name % options.assembler_param % neg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		cout << "[TEA.output_mate_fa] adds reverse negative FASTA files\n";
		for(auto& an_entry : negative_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string rneg_fa_name = (boost::format("%s/%s.rneg.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(rneg_fa_name)) {
					continue;
				}
				all_files << rneg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % rneg_fa_name % options.assembler_param % rneg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		cout << "[TEA.output_mate_fa] adds positive clipped FASTA files\n";
		for(auto& a_prefix: a_positive_clipped_prefixes) {
	//		tasks.push_back([&]{
				string pos_fa_name = (boost::format("%s/%s.pos.fa") % contig_dir % a_prefix).str();
				if(!boost::filesystem::exists(pos_fa_name)) {
					continue;
				}
				all_files << pos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % pos_fa_name % options.assembler_param % pos_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}

		cout << "[TEA.output_mate_fa] adds negative clipped FASTA files\n";
		for(auto& a_prefix: a_negative_clipped_prefixes) {
	//		tasks.push_back([&]{
				string neg_fa_name = (boost::format("%s/%s.neg.fa") % contig_dir % a_prefix).str();
				if(!boost::filesystem::exists(neg_fa_name)) {
					continue;
				}
				all_files << neg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % neg_fa_name % options.assembler_param % neg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		cout << "[TEA.output_mate_fa] adds germline FASTA files\n";
		bool headless = false;
		for (auto& c : chr_names) {
//			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			vector<string> data;
			vector<string> cols;
			string line;
			string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";

			ifstream in_germline(germline_file, ios::binary);
			// ignore the first line;
			if(!headless) {
				getline(in_germline, line, '\n');
			}
			while(getline(in_germline, line, '\n')) {
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[8], delim_comma, cols);
				string key_name = data[0] + "." + data[1] + "." + data[2] + "." + data[3] + "." + cols[0];
				string pclipped = contig_dir + "/" + key_name + ".pos.fa";
				string nclipped = contig_dir + "/" + key_name + ".neg.fa";
				string prammate = contig_dir + "/" + key_name + ".rpos.fa";
				string nrammate = contig_dir + "/" + key_name + ".rneg.fa";
				if(boost::filesystem::exists(pclipped)) {
					all_files << pclipped << "\n";
				}
				if(boost::filesystem::exists(nclipped)) {
					all_files << nclipped << "\n";
				}
				if(boost::filesystem::exists(prammate)) {
					all_files << prammate << "\n";
				}
				if(boost::filesystem::exists(nrammate)) {
					all_files << nrammate << "\n";
				}
			}
			headless = true;
		}
	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[TEA.output_mate_fa] runs assembling\n";
	string out_log_file = options.prefix + ".assemblies.log";
	if(!options.working_dir.empty()) {
		out_log_file = options.working_prefix + ".assemblies.log";
	}
	string assembler_cmd = (boost::format("cat %s | xargs -n 1 -P %d -I {} %s {} %s > %s") % all_assembly_files % n_cores % options.assembler % options.assembler_param % out_log_file).str();
	system(assembler_cmd.c_str());

	cout << "[TEA.output_mate_fa] write germline calls\n";
	bool headless = false;
	for (auto& c : chr_names) {
		tasks.push_back([&, headless] {
//			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			vector<string> data;
			vector<string> cols;
			string line;
			string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";
			string germline_contig_file = cl_prefix + "." + tmp_chr_name + ".germline.contig";

			ifstream in_germline(germline_file, ios::binary);
			ofstream out_germline(germline_contig_file, ios::binary);
			// ignore the first line;
			if(!headless) {
				getline(in_germline, line, '\n');
			}
			while(getline(in_germline, line, '\n')) {
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[9], delim_comma, cols);
				string key_name = data[0] + "." + data[1] + "." + data[2] + "." + data[3] + "." + cols[0];
				string pclipped_prefix = contig_dir + "/" + key_name + ".pos.fa.cap";
//				cout << pclipped_prefix << "\n";
				string pclipped_contig = pclipped_prefix + ".contigs";
				string pclipped_singlets = pclipped_prefix + ".singlets";

				string pclipped("-");
				get_longest_fa(pclipped, pclipped_contig);
				get_longest_fa(pclipped, pclipped_singlets);

				string nclipped_prefix = contig_dir + "/" + key_name + ".neg.fa.cap";
				string nclipped_contig = nclipped_prefix + ".contigs";
				string nclipped_singlets = nclipped_prefix + ".singlets";

				string nclipped("-");
				get_longest_fa(nclipped, nclipped_contig);
				get_longest_fa(nclipped, nclipped_singlets);

				string prammate_prefix = contig_dir + "/" + key_name + ".rpos.fa.cap";
				string prammate_contig = prammate_prefix + ".contigs";
				string prammate_singlets = prammate_prefix + ".singlets";

				string prammate("-");
				get_longest_fa(prammate, prammate_contig);
				get_longest_fa(prammate, prammate_singlets);

				string nrammate_prefix = contig_dir + "/" + key_name + ".rneg.fa.cap";
				string nrammate_contig = nrammate_prefix + ".contigs";
				string nrammate_singlets = nrammate_prefix + ".singlets";

				string nrammate("-");
				get_longest_fa(nrammate, nrammate_contig);
				get_longest_fa(nrammate, nrammate_singlets);

				string orientation("NA");

				string polyA("-");
				string polyT("-");

				if ("-" != pclipped) {
					int64_t the_pos = pclipped.size();
					the_pos -= 6;
					if (the_pos < 0) {
						the_pos = 0;
					}
					string the_suffix = pclipped.substr(the_pos);
					size_t acnt = count(the_suffix.begin(), the_suffix.end(), 'A');
					if(acnt >= 5) {
						polyA = "polyA";
					}
				}
				if ("-" != nclipped) {
					int64_t the_pos = nclipped.size();
					the_pos -= 6;
					if (the_pos < 0) {
						the_pos = 0;
					}
					string the_suffix = nclipped.substr(the_pos);
					size_t tcnt = count(the_suffix.begin(), the_suffix.end(), 'T');
					if(tcnt >= 5) {
						polyT = "polyT";
					}
				}
				if ("polyA" == polyA && "-" == polyT) {
					orientation = "+";
				}
				if ("polyT" == polyT && "-" == polyA) {
					orientation = "-";
				}

				out_germline << line << "\t" << orientation << "\t" << polyA << "\t" << polyT << "\t" << pclipped << "\t" << nclipped << "\t" << prammate << "\t" << nrammate << "\n";
			}
		});
		headless = true;
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void TEA::_output_mate_fa(boost::unordered_map<string, vector<string>>& positive_mate_reads, boost::unordered_map<string, vector<string>>& negative_mate_reads, vector<meerkat::BlockBoundary>& actual_blocks, const string& a_path,
		const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map, const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map) {
	string a_bai_path;
	get_bai_index_path(a_path, a_bai_path);

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	if (verbose) {
//		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
//			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
//		}
//	}
	// map<read name, pair<first read, second read>>
	vector<boost::unordered_map<string, vector<string>>> positive_mate_reads_list(calculated_n_blocks - 1);
	vector<boost::unordered_map<string, vector<string>>> negative_mate_reads_list(calculated_n_blocks - 1);
//	vector<string> read_groups;
	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
				return;
			}
			int64_t num_total = 0;

			BamTools::BamAlignment local_alignment_entry;
//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

//			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA._output_mate_fa] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA._output_mate_fa] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
//			map<string, int64_t> ref_reverse_index;
//			const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
//			for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
//				auto& a_ref = a_ref_vector[ref_id];
//				ref_reverse_index[a_ref.RefName] = ref_id;
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;

//			const bool is_soft_clipped = string::npos != input_BAM_name.find(".cl.disc.sorted.bam");

				auto& local_positive_mate_reads = positive_mate_reads_list[block_id];
				auto& local_negative_mate_reads = negative_mate_reads_list[block_id];

//				const bool debug = (0 == block_id);
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA._output_mate_fa] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA._output_mate_fa] Block-%d (first) %s\n")
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

					auto the_pos_itr = a_positive_repeat_map.find(local_alignment_entry.Name);
					auto the_neg_itr = a_negative_repeat_map.find(local_alignment_entry.Name);
					if (a_positive_repeat_map.end() != the_pos_itr) {
						auto& aln_pair = the_pos_itr->second;
						int64_t the_pos = aln_pair.pos;
						int64_t aln_pos = local_alignment_entry.Position + 1;

//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
//						const bool debug = ("14" == aln_pair.chr) && (34139886 <= aln_pair.pos || 34140417 >= aln_pair.pos);
//						if(debug) {
//							cout << "[TEA._output_mate_fa] prefix: " << aln_pair.file_name_prefix << ", pos: " << aln_pair.pos << "\n";
//						}
						if(aln_pos != the_pos) {
//							cout << "Pos: " << the_pos << ": " << BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData()) << "\n";
							if(local_alignment_entry.IsReverseStrand()) {
								local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
							}
							auto& the_seq_vec = local_positive_mate_reads[aln_pair.file_name_prefix];
							the_seq_vec.push_back(local_alignment_entry.QueryBases);
						}
//					}
					}
					if (a_negative_repeat_map.end() != the_neg_itr) {
						auto& aln_pair = the_neg_itr->second;
						int64_t the_pos = aln_pair.pos;
						int64_t aln_pos = local_alignment_entry.Position + 1;

//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
						if(aln_pos != the_pos) {
//							cout << "Neg: " << the_pos << ": " << BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData()) << "\n";
							if(local_alignment_entry.IsReverseStrand()) {
								local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
							}
							auto& the_seq_vec = local_negative_mate_reads[aln_pair.file_name_prefix];
							the_seq_vec.push_back(local_alignment_entry.QueryBases);
						}
//					}
					}
				}

				local_reader.Close();
//				done_vector[block_id] = 'D';

//				if(verbose) {
//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//							% local_alignment_entry.AlignmentFlag).str();
//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
//								% block_id % a_block_boundary_str).str();
//					} else {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//				} else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	cout << "[TEA._output_mate_fa] gather scattered information\n";
	tasks.push_back([&, calculated_n_blocks] {
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_positive_mate_reads = positive_mate_reads_list[block_id];
			for (auto& an_entry: local_positive_mate_reads) {
				auto& file_name_prefix = an_entry.first;
				auto& the_seq_vec = an_entry.second;
				auto& the_target_vec = positive_mate_reads[file_name_prefix];
				the_target_vec.insert(the_target_vec.end(), the_seq_vec.begin(), the_seq_vec.end());
			}
		}
	});
	tasks.push_back([&, calculated_n_blocks] {
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_negative_mate_reads = negative_mate_reads_list[block_id];
			for (auto& an_entry: local_negative_mate_reads) {
				auto& file_name_prefix = an_entry.first;
				auto& the_seq_vec = an_entry.second;
				auto& the_target_vec = negative_mate_reads[file_name_prefix];
				the_target_vec.insert(the_target_vec.end(), the_seq_vec.begin(), the_seq_vec.end());
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void TEA::output_mate_fa_v(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram) {
	castle::TimeChecker checker;
	checker.setTarget("TEA.output_mate_fa_v");
	checker.start();
	string cl_dir = options.prefix + "/cluster_" + options.rasym + "m";
	string naive_prefix = options.naive_prefix;
	string contig_dir = options.prefix + "/assembly_" + options.rasym + "m";

	if (!options.working_dir.empty()) {
		cl_dir = options.working_prefix + "/cluster_" + options.rasym + "m";
		contig_dir = options.working_prefix + "/assembly_" + options.rasym + "m";
	}
	string cl_prefix = cl_dir + "/" + naive_prefix;

	boost::unordered_map<string, AlnPairEntry> a_positive_repeat_map;
	boost::unordered_map<string, AlnPairEntry> a_negative_repeat_map;
	vector<string> a_positive_clipped_prefixes;
	vector<string> a_negative_clipped_prefixes;

	vector<function<void()> > tasks;

	tasks.push_back([&] {
//		string c = "22";
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;

			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			vector<string> data;
			vector<string> keys;
			string line;

//			string clipped_file = cl_prefix + "." + tmp_chr_name + ".p.mate.rname";
			string clipped_file = cl_prefix + "." + tmp_chr_name + ".cluster.raw";
			ifstream in_clipped(clipped_file, ios::binary);
			// to ignore the first line
			getline(in_clipped, line);
			while(getline(in_clipped, line, '\n')) {
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[data.size() - 2], delim_comma, keys);
				int64_t a_pos = boost::lexical_cast<int64_t>(data[3]);
				string the_first_taxon_id = data[4];
				auto the_delim_pos = the_first_taxon_id.find("_&_");
				if(string::npos != the_delim_pos) {
					the_first_taxon_id = the_first_taxon_id.substr(0, the_delim_pos);
				}
				string file_name_prefix = naive_prefix + "." + tmp_chr_name + "." + data[1] + "." + data[2] + "." + the_first_taxon_id;
				for(auto& a_key : keys) {
					a_positive_repeat_map[a_key].file_name_prefix = file_name_prefix;
					a_positive_repeat_map[a_key].pos = a_pos;
					a_positive_repeat_map[a_key].chr = c;
				}
			}
		}
	});

	tasks.push_back([&] {
//		string c = "22";
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			const char* delim_comma = ",";
			vector<string> data;
			vector<string> keys;
			string line;
//			string clipped_file = cl_prefix + "." + tmp_chr_name + ".n.mate.rname";
			string clipped_file = cl_prefix + "." + tmp_chr_name + ".cluster.raw";
			ifstream in_clipped(clipped_file, ios::binary);
			// to ignore the first line
			getline(in_clipped, line);
			while(getline(in_clipped, line, '\n')) {
				castle::StringUtils::tokenize(line, delim_tab, data);
				castle::StringUtils::c_string_multi_split(data[data.size() - 1], delim_comma, keys);
				int64_t a_pos = boost::lexical_cast<int64_t>(data[3]);
				string the_first_taxon_id = data[4];
				auto the_delim_pos = the_first_taxon_id.find("_&_");
				if(string::npos != the_delim_pos) {
					the_first_taxon_id = the_first_taxon_id.substr(0, the_delim_pos);
				}
				string file_name_prefix = naive_prefix + "." + tmp_chr_name + "." + data[1] + "." + data[2] + "." + the_first_taxon_id;
				for(auto& a_key : keys) {
					a_negative_repeat_map[a_key].file_name_prefix = file_name_prefix;
					a_negative_repeat_map[a_key].pos = a_pos;
					a_negative_repeat_map[a_key].chr = c;
				}
//				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
//				string a_key = data[4];
//				int64_t a_pos = boost::lexical_cast<int64_t>(data[3]);
//				string file_name_prefix = naive_prefix + "." + tmp_chr_name + "." + data[0] + "." + data[1] + "." + data[2];
//				a_negative_repeat_map[a_key].file_name_prefix = file_name_prefix;
//				a_negative_repeat_map[a_key].pos = a_pos;
//				a_negative_repeat_map[a_key].chr = c;
			}
		}
	});
	tasks.push_back([&] {
		//		string c = "22";
		string line;
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			string p_mate_clipped_file = cl_prefix + "." + tmp_chr_name + ".p.clipped.fname";
			ifstream in_clipped(p_mate_clipped_file, ios::binary);
			while(getline(in_clipped, line, '\n')) {
				a_positive_clipped_prefixes.push_back(line);
			}
		}
	});
	tasks.push_back([&] {
		//		string c = "22";
		string line;
		for(auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			string n_mate_clipped_file = cl_prefix + "." + tmp_chr_name + ".n.clipped.fname";
			ifstream in_clipped(n_mate_clipped_file, ios::binary);
			while(getline(in_clipped, line, '\n')) {
				a_negative_clipped_prefixes.push_back(line);
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[TEA.output_mate_fa_v] # positive mates: " << a_positive_repeat_map.size() << ", # negative mates: " << a_negative_repeat_map.size() << "\n";


//	string disc_file_name = options.prefix + ".disc.num.bam";
//	string cl_disc_file_name = options.prefix + ".cl.sorted.disc.num.bam";

//	if("tea_v" == options.program_name) {
//		disc_file_name = options.prefix + ".cl.sorted.bam";
//	}
//
//	if (!options.working_dir.empty()) {
//		disc_file_name = options.working_prefix + ".disc.num.bam";
//		if("tea_v" == options.program_name) {
//			disc_file_name = options.working_prefix + ".cl.sorted.bam";
//		}
//		cl_disc_file_name = options.working_prefix + ".cl.sorted.disc.num.bam";
//	}

	string disc_file_name = options.prefix + ".um.bam";
	if (!options.working_dir.empty()) {
		disc_file_name = options.working_prefix + ".um.bam";
	}
	string disc_bai_name;
	get_bai_index_path(disc_file_name, disc_bai_name);
	string disc_bni_name = disc_file_name + ".bni";

	boost::unordered_map<string, vector<string>> positive_mate_reads;
	boost::unordered_map<string, vector<string>> negative_mate_reads;

	int64_t size_block = 8192000;
//	vector<meerkat::BlockBoundary> fixed_size_blocks;

	cout << "[TEA.output_mate_fa_v] collect mate reads from disc file\n";
//	_output_mate_fa_serial(positive_mate_reads, negative_mate_reads, disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
	vector<meerkat::BlockBoundary> local_fixed_size_blocks;
	vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
	vector<meerkat::BlockBoundary> local_independent_blocks;
	collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, disc_file_name, disc_bai_name, disc_bni_name, size_block);
//	vector<meerkat::BlockBoundary> actual_blocks = local_unmapped_included_blocks;
//	collect_boundaries(fixed_size_blocks, disc_file_name, size_block);
	_output_mate_fa_v(positive_mate_reads, negative_mate_reads, local_unmapped_included_blocks, disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
//	bool found_clipped_reads = false;
//	for(auto& an_entry : a_positive_repeat_map) {
//		auto& read_name = an_entry.first;
//		if(string::npos != read_name.rfind("mu1") || string::npos != read_name.rfind("mu2") || string::npos != read_name.rfind("sc")) {
//			found_clipped_reads = true;
//			break;
//		}
//	}
//	if(!found_clipped_reads) {
//		for(auto& an_entry : a_negative_repeat_map) {
//			auto& read_name = an_entry.first;
//			if(string::npos != read_name.rfind("mu1") || string::npos != read_name.rfind("mu2") || string::npos != read_name.rfind("sc")) {
//				found_clipped_reads = true;
//				break;
//			}
//		}
//	}

//	if(found_clipped_reads) {
//		cout << "[TEA.output_mate_fa_v] collect mate reads from cl disc file\n";
////		_output_mate_fa_serial(positive_mate_reads, negative_mate_reads, cl_disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
//		vector<meerkat::BlockBoundary> local_fixed_size_blocks;
//		vector<meerkat::BlockBoundary> local_unmapped_included_blocks;
//		vector<meerkat::BlockBoundary> local_independent_blocks;
//		collect_boundaries_pos(local_fixed_size_blocks, local_unmapped_included_blocks, local_independent_blocks, cl_disc_file_name, size_block);
////		collect_boundaries(fixed_size_blocks, cl_disc_file_name, size_block);
//		_output_mate_fa_v(positive_mate_reads, negative_mate_reads, local_unmapped_included_blocks, cl_disc_file_name, a_positive_repeat_map, a_negative_repeat_map);
//	}

	for(auto& an_entry : positive_mate_reads) {
		tasks.push_back([&] {
			auto& file_name_prefix = an_entry.first;
			auto& the_seq_vec = an_entry.second;
			string fa_name = (boost::format("%s/%s.rpos.fa") % contig_dir % file_name_prefix).str();
			ofstream out_fa(fa_name, ios::binary);
			for(int64_t n_id = 0; n_id < static_cast<int64_t>(the_seq_vec.size()); ++n_id) {
				auto& a_seq = the_seq_vec[n_id];
				out_fa << (boost::format(">cr%d\n%s\n") % n_id % a_seq).str();
			}
		});
	}

	for(auto& an_entry : negative_mate_reads) {
		tasks.push_back([&] {
			auto& file_name_prefix = an_entry.first;
			auto& the_seq_vec = an_entry.second;
			string fa_name = (boost::format("%s/%s.rneg.fa") % contig_dir % file_name_prefix).str();
			ofstream out_fa(fa_name, ios::binary);
			for(int64_t n_id = 0; n_id < static_cast<int64_t>(the_seq_vec.size()); ++n_id) {
				auto& a_seq = the_seq_vec[n_id];
				out_fa << (boost::format(">cr%d\n%s\n") % n_id % a_seq).str();
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << "[TEA.output_mate_fa_v] start assembling\n";
	string all_assembly_files = options.prefix + ".asssemblies.target";
	if (!options.working_dir.empty()) {
		all_assembly_files = options.working_prefix + ".asssemblies.target";
	}
	{
		ofstream all_files(all_assembly_files, ios::binary);
		for(auto& an_entry : positive_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string pos_fa_name = (boost::format("%s/%s.pos.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(pos_fa_name)) {
					continue;
				}
				all_files << pos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % pos_fa_name % options.assembler_param % pos_fa_name).str();
	//			all_cmds << assembler_cmd << "\n";
	//			system(assembler_cmd.c_str());
	//		});
		}
		for(auto& an_entry : positive_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string rpos_fa_name = (boost::format("%s/%s.rpos.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(rpos_fa_name)) {
					continue;
				}
				all_files << rpos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % rpos_fa_name % options.assembler_param % rpos_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		for(auto& an_entry : negative_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string neg_fa_name = (boost::format("%s/%s.neg.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(neg_fa_name)) {
					continue;
				}
				all_files << neg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % neg_fa_name % options.assembler_param % neg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		for(auto& an_entry : negative_mate_reads) {
	//		tasks.push_back([&]{
				auto& file_name_prefix = an_entry.first;
				string rneg_fa_name = (boost::format("%s/%s.rneg.fa") % contig_dir % file_name_prefix).str();
				if(!boost::filesystem::exists(rneg_fa_name)) {
					continue;
				}
				all_files << rneg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % rneg_fa_name % options.assembler_param % rneg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}

		for(auto& a_prefix: a_positive_clipped_prefixes) {
	//		tasks.push_back([&]{
				string pos_fa_name = (boost::format("%s/%s.pos.fa") % contig_dir % a_prefix).str();
				if(!boost::filesystem::exists(pos_fa_name)) {
					continue;
				}
				all_files << pos_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % pos_fa_name % options.assembler_param % pos_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}

		for(auto& a_prefix: a_negative_clipped_prefixes) {
	//		tasks.push_back([&]{
				string neg_fa_name = (boost::format("%s/%s.neg.fa") % contig_dir % a_prefix).str();
				if(!boost::filesystem::exists(neg_fa_name)) {
					continue;
				}
				all_files << neg_fa_name << "\n";
	//			string assembler_cmd = (boost::format("%s %s %s > %s.log") % options.assembler % neg_fa_name % options.assembler_param % neg_fa_name).str();
	//			system(assembler_cmd.c_str());
	//		});
		}
		for (auto& c_entry : ram) {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			vector<string> data;
			string line;
			string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";

			ifstream in_germline(germline_file, ios::binary);
			// ignore the first line;
			getline(in_germline, line, '\n');
			while(getline(in_germline, line, '\n')) {
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
				string key_name = data[0] + "." + data[1] + "." + data[2] + "." + data[3] + "." + data[8];
				string pclipped = contig_dir + "/" + key_name + ".pos.fa";
				string nclipped = contig_dir + "/" + key_name + ".neg.fa";
				string prammate = contig_dir + "/" + key_name + ".rpos.fa";
				string nrammate = contig_dir + "/" + key_name + ".rneg.fa";
				if(boost::filesystem::exists(pclipped)) {
					all_files << pclipped << "\n";
				}
				if(boost::filesystem::exists(nclipped)) {
					all_files << nclipped << "\n";
				}
				if(boost::filesystem::exists(prammate)) {
					all_files << prammate << "\n";
				}
				if(boost::filesystem::exists(nrammate)) {
					all_files << nrammate << "\n";
				}
			}
		}
	}
//	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	string assembler_cmd = (boost::format("cat %s | xargs -n 1 -P %d -I {} %s {} %s > {}.log") % all_assembly_files % n_cores % options.assembler % options.assembler_param).str();
//	cout << assembler_cmd << "\n";
	system(assembler_cmd.c_str());

	cout << "[TEA.output_mate_fa_v] write contigs\n";
	for (auto& c_entry : ram) {
		tasks.push_back([&] {
			auto& c = c_entry.first;
			string tmp_chr_name(c);
			if(string::npos == tmp_chr_name.find("chr")) {
				tmp_chr_name = "chr" + c;
			}
			const char* delim_tab = "\t";
			vector<string> data;
			string line;
			string germline_file = cl_prefix + "." + tmp_chr_name + ".germline";
			string germline_contig_file = cl_prefix + "." + tmp_chr_name + ".contig";

			ifstream in_germline(germline_file, ios::binary);
			ofstream out_germline(germline_contig_file, ios::binary);
			// ignore the first line;
			getline(in_germline, line, '\n');
			while(getline(in_germline, line, '\n')) {
				castle::StringUtils::c_string_multi_split(line, delim_tab, data);
//				string key_name = options.naive_prefix + "." + data[0] + "." + data[1] + "." + data[2] + "." + data[7];
				string the_first_taxon_id = data[8];
				auto the_delim_pos = the_first_taxon_id.find("_&_");
				if(string::npos != the_delim_pos) {
					the_first_taxon_id = the_first_taxon_id.substr(0, the_delim_pos);
				}

				string key_name = data[0] + "." + data[1] + "." + data[2] + "." + data[3] + "." + the_first_taxon_id;
				string pclipped_prefix = contig_dir + "/" + key_name + ".pos.fa.cap";
//				cout << pclipped_prefix << "\n";
				string pclipped_contig = pclipped_prefix + ".contigs";
				string pclipped_singlets = pclipped_prefix + ".singlets";

				string pclipped("-");
				get_longest_fa(pclipped, pclipped_contig);
				get_longest_fa(pclipped, pclipped_singlets);

				string nclipped_prefix = contig_dir + "/" + key_name + ".neg.fa.cap";
				string nclipped_contig = nclipped_prefix + ".contigs";
				string nclipped_singlets = nclipped_prefix + ".singlets";

				string nclipped("-");
				get_longest_fa(nclipped, nclipped_contig);
				get_longest_fa(nclipped, nclipped_singlets);

				string prammate_prefix = contig_dir + "/" + key_name + ".rpos.fa.cap";
				string prammate_contig = prammate_prefix + ".contigs";
				string prammate_singlets = prammate_prefix + ".singlets";

				string prammate("-");
				get_longest_fa(prammate, prammate_contig);
				get_longest_fa(prammate, prammate_singlets);

				string nrammate_prefix = contig_dir + "/" + key_name + ".rneg.fa.cap";
				string nrammate_contig = nrammate_prefix + ".contigs";
				string nrammate_singlets = nrammate_prefix + ".singlets";

				string nrammate("-");
				get_longest_fa(nrammate, nrammate_contig);
				get_longest_fa(nrammate, nrammate_singlets);

				string orientation("NA");

				string polyA("-");
				string polyT("-");

				if ("-" != pclipped) {
					int64_t the_pos = pclipped.size();
					the_pos -= 6;
					if (the_pos < 0) {
						the_pos = 0;
					}
					string the_suffix = pclipped.substr(the_pos);
					size_t acnt = count(the_suffix.begin(), the_suffix.end(), 'A');
					if(acnt >= 5) {
						polyA = "polyA";
					}
				}
				if ("-" != nclipped) {
					int64_t the_pos = nclipped.size();
					the_pos -= 6;
					if (the_pos < 0) {
						the_pos = 0;
					}
					string the_suffix = nclipped.substr(the_pos);
					size_t tcnt = count(the_suffix.begin(), the_suffix.end(), 'T');
					if(tcnt >= 5) {
						polyT = "polyT";
					}
				}
				if ("polyA" == polyA && "-" == polyT) {
					orientation = "+";
				}
				if ("polyT" == polyT && "-" == polyA) {
					orientation = "-";
				}

				out_germline << line << "\t" << orientation << "\t" << polyA << "\t" << polyT << "\t" << pclipped << "\t" << nclipped << "\t" << prammate << "\t" << nrammate << "\n";
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;

}

void TEA::_output_mate_fa_v(boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			vector<meerkat::BlockBoundary>& actual_blocks, const string& a_path,
			const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
			const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map) {
	string a_bai_path;
	get_bai_index_path(a_path, a_bai_path);

	int64_t calculated_n_blocks = actual_blocks.size();
	vector<function<void()> > tasks;

// only for debugging
//	const bool verbose = false;
	string done_vector(calculated_n_blocks - 1, 'U');
	set<string> block_boundary_strs;
//	if (verbose) {
//		for (auto itr = actual_blocks.begin(); actual_blocks.end() != itr; ++itr) {
//			block_boundary_strs.insert((boost::format("%s %d-%d %d") % itr->read_name % itr->ref_id % itr->pos % itr->aln_flag).str());
//		}
//	}
	// map<read name, pair<first read, second read>>
	vector<boost::unordered_map<string, vector<string>>> positive_mate_reads_list(calculated_n_blocks - 1);
	vector<boost::unordered_map<string, vector<string>>> negative_mate_reads_list(calculated_n_blocks - 1);
//	vector<string> read_groups;
	string target = "FCC1E2JACXX:5:2314:17665:69999#AATAAGATsc,FCC1E2JACXX:6:1113:12930:67887#AATCAGATsc,FCC1E2JACXX:6:2207:10389:96543#AATAAGATsc,FCC1E2JACXX:5:2112:6578:23827#AATAAGATmu2,FCC1E2JACXX:6:2106:9722:19446#AATAAGATmu2,FCD1JLLACXX:7:1209:20291:27741#AATAAGATmu2,FCC1E2JACXX:6:2102:12105:15616#AATAAGATmu2,FCC1E2JACXX:5:1112:3968:2738#AATAAGATmu2,FCC1E2JACXX:6:2101:5254:54827#AATAAGATmu2,FCC1E2JACXX:6:2204:12157:46324#AATAAGATmu2,FCD1JLLACXX:7:2115:13002:26719#AATAAGATmu2,FCD1JLLACXX:7:2302:14919:52085#AATAAGATmu2,FCC1E2JACXX:6:1310:9638:5746#CATAAGATmu2,FCC1E2JACXX:5:2301:18305:35316#AATAAGATmu2,FCC1E2JACXX:5:2206:10456:95452#AATAAGATmu2,FCC1E2JACXX:6:1104:10876:81194#AATAAGATmu2,FCC1E2JACXX:6:1112:13437:99251#AATAAGATmu2,FCC1E2JACXX:6:2107:15598:98676#AATAAGATmu2,FCC1E2JACXX:5:1307:10179:4241#AATAAGATmu2,FCD1JLLACXX:7:2311:17154:7342#AATAAGATmu2,FCC1E2JACXX:5:1214:12082:15240#AGTAAAAGmu2,FCC1E2JACXX:6:1316:13821:89836#AATAAGAGmu2,FCC1E2JACXX:6:2206:18414:41434#AATAAGATmu2";
	vector<string> debug_str_str;
	castle::StringUtils::c_string_multi_split(target, ",", debug_str_str);
	set<string> debug_str(debug_str_str.begin(), debug_str_str.end());

	for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
		tasks.push_back([&, block_id] {
			BamTools::BamReader local_reader;
			if (!local_reader.Open(a_path, a_bai_path)) {
				return;
			}
			int64_t num_total = 0;

			BamTools::BamAlignment local_alignment_entry;
//			int32_t the_current_ref_id = actual_blocks[block_id].ref_id;
//			int32_t the_current_ref_pos = actual_blocks[block_id].pos;
//
//			int32_t the_next_ref_id = actual_blocks[block_id + 1].ref_id;
//			int32_t the_next_ref_pos = actual_blocks[block_id + 1].pos;
//			uint32_t the_next_aln_flag = actual_blocks[block_id + 1].aln_flag;
			string str_block_id = boost::lexical_cast<string>(block_id);

			string the_next_block_read_name = actual_blocks[block_id + 1].read_name;
//			bool jump_success = local_reader.Jump(actual_blocks[block_id].ref_id, actual_blocks[block_id].jump_pos);
//			if(!jump_success) {
//				cout << (boost::format("[TEA._output_mate_fa] block-%d (Jump fail): %d:%d (%d/%d)-(%d/%d)\n")
//						% block_id % actual_blocks[block_id].ref_id % actual_blocks[block_id].jump_pos % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//				local_reader.Close();
//				return;
//			}
//			if(verbose) {
//				cout << (boost::format("[TEA._output_mate_fa] block-%d (start) (%d/%d)-(%d/%d)\n")
//						% block_id % the_current_ref_id % the_current_ref_pos
//						% the_next_ref_id % the_next_ref_pos).str();
//			}
//			map<string, int64_t> ref_reverse_index;
//			const BamTools::RefVector& a_ref_vector = local_reader.GetReferenceData();
//			for (uint64_t ref_id = 0; ref_id < a_ref_vector.size(); ++ref_id) {
//				auto& a_ref = a_ref_vector[ref_id];
//				ref_reverse_index[a_ref.RefName] = ref_id;
//			}
			int64_t the_current_ref_offset = actual_blocks[block_id].offset;
			int64_t the_next_ref_offset = actual_blocks[block_id + 1].offset;
			auto& m_bgzf = local_reader.GetBGZF();
			if(0 != block_id) {
				if(!m_bgzf.Seek(the_current_ref_offset)) {
					local_reader.Close();
					return;
				}
			}
			int64_t cur_offset = m_bgzf.Tell();
			int64_t prev_offset = cur_offset;

//			const bool is_soft_clipped = string::npos != input_BAM_name.find(".cl.disc.sorted.bam");

				auto& local_positive_mate_reads = positive_mate_reads_list[block_id];
				auto& local_negative_mate_reads = negative_mate_reads_list[block_id];

//				const bool debug = (0 == block_id);
				while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
//					if(verbose && 0 == num_total) {
//						string a_block_boundary_str = (boost::format("%s %d-%d %d")
//								% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//								% local_alignment_entry.AlignmentFlag).str();
//						if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//							cout << (boost::format("[TEA._output_mate_fa] Block-%d (first-wrong) %s\n")
//									% block_id % a_block_boundary_str).str();
//						} else {
//							cout << (boost::format("[TEA._output_mate_fa] Block-%d (first) %s\n")
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
					bool debug = debug_str.end() != debug_str.find(local_alignment_entry.Name);
					auto the_pos_itr = a_positive_repeat_map.find(local_alignment_entry.Name);
					auto the_neg_itr = a_negative_repeat_map.find(local_alignment_entry.Name);
					if (a_positive_repeat_map.end() != the_pos_itr) {
						auto& aln_pair = the_pos_itr->second;
//						int64_t the_pos = aln_pair.pos;
//						int64_t aln_pos = local_alignment_entry.Position + 1;

//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
						if(!local_alignment_entry.IsMapped()) {

//						if(aln_pos != the_pos) {
//							cout << "Pos: " << the_pos << ": " << BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData()) << "\n";
//							if(local_alignment_entry.IsReverseStrand()) {
							local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
//							}
							if(debug) {
								cout << (boost::format("[TEA._output_mate_fa_v] %s\n") % BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData())).str();
								cout << "[TEA._output_mate_fa_v] positive\n";
							}
							auto& the_seq_vec = local_positive_mate_reads[aln_pair.file_name_prefix];
							the_seq_vec.push_back(local_alignment_entry.QueryBases);
						}
//					}
					}
					if (a_negative_repeat_map.end() != the_neg_itr) {

						auto& aln_pair = the_neg_itr->second;
//						int64_t the_pos = aln_pair.pos;
//						int64_t aln_pos = local_alignment_entry.Position + 1;

//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
						if(!local_alignment_entry.IsMapped()) {
							if(debug) {
								cout << (boost::format("[TEA._output_mate_fa_v] %s\n") % BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData())).str();
								cout << "[TEA._output_mate_fa_v] negative\n";
							}
//						if(aln_pos != the_pos) {
//							cout << "Neg: " << the_pos << ": " << BamWriter::GetSAMAlignment(local_alignment_entry, local_reader.GetReferenceData()) << "\n";
							if(local_alignment_entry.IsReverseStrand()) {
								local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
							}
							auto& the_seq_vec = local_negative_mate_reads[aln_pair.file_name_prefix];
							the_seq_vec.push_back(local_alignment_entry.QueryBases);
						}
//					}
					}
				}

				local_reader.Close();
//				done_vector[block_id] = 'D';

//				if(verbose) {
//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//							% local_alignment_entry.AlignmentFlag).str();
//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
//								% block_id % a_block_boundary_str).str();
//					} else {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//				} else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	cout << "[TEA._output_mate_fa] gattter scattered information\n";
	tasks.push_back([&, calculated_n_blocks] {
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_positive_mate_reads = positive_mate_reads_list[block_id];
			for (auto& an_entry: local_positive_mate_reads) {
				auto& file_name_prefix = an_entry.first;
				auto& the_seq_vec = an_entry.second;
				auto& the_target_vec = positive_mate_reads[file_name_prefix];
				the_target_vec.insert(the_target_vec.end(), the_seq_vec.begin(), the_seq_vec.end());
			}
		}
	});
	tasks.push_back([&, calculated_n_blocks] {
		for (int64_t block_id = 0; block_id < calculated_n_blocks - 1; ++block_id) {
			auto& local_negative_mate_reads = negative_mate_reads_list[block_id];
			for (auto& an_entry: local_negative_mate_reads) {
				auto& file_name_prefix = an_entry.first;
				auto& the_seq_vec = an_entry.second;
				auto& the_target_vec = negative_mate_reads[file_name_prefix];
				the_target_vec.insert(the_target_vec.end(), the_seq_vec.begin(), the_seq_vec.end());
			}
		}
	});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);

}

void TEA::_output_mate_fa_serial(boost::unordered_map<string, vector<string>>& positive_mate_reads, boost::unordered_map<string, vector<string>>& negative_mate_reads, const string& input_BAM_name, const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
		const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map) {
	string a_path(input_BAM_name);
	string an_index_path(a_path);
	an_index_path += ".bai";

// only for debugging
	// map<read name, pair<first read, second read>>
//	vector<string> read_groups;
	BamTools::BamReader local_reader;
	if (!local_reader.Open(a_path, an_index_path)) {
		return;
	}
	int64_t num_total = 0;

	BamTools::BamAlignment local_alignment_entry;

//		const bool is_soft_clipped = string::npos != input_BAM_name.find(".cl.disc.sorted.bam");

	string target_str =
			"DHFC08P1:370:C1784ACXX:2:1114:18204:57892,DHFC08P1:370:C1784ACXX:1:1212:15037:61907mu1,DHFC08P1:370:C1784ACXX:1:1113:17042:78731mu1,DHFC08P1:370:C1784ACXX:1:2210:3273:91802mu1,DHFC08P1:370:C1784ACXX:2:1214:6860:23842,DHFC08P1:370:C1784ACXX:1:1204:3054:21186mu1,DHFC08P1:370:C1784ACXX:1:2203:20922:42380,DHFC08P1:370:C1784ACXX:1:2306:20447:97121mu1,DHFC08P1:370:C1784ACXX:2:2215:4090:25856mu1,DHFC08P1:370:C1784ACXX:1:2115:10819:58552,HWI-ST115:415:C170UACXX:6:2202:19413:98859,HWI-ST115:415:C170UACXX:6:2204:19696:50780mu1,DHFC08P1:370:C1784ACXX:1:2314:18589:6896,HWI-ST115:415:C170UACXX:6:1312:5194:66873,DHFC08P1:370:C1784ACXX:1:2102:3399:50694,DHFC08P1:370:C1784ACXX:1:2210:3534:20988,HWI-ST115:415:C170UACXX:6:1111:12602:96456,HWI-ST115:415:C170UACXX:6:2311:15847:27338,DHFC08P1:370:C1784ACXX:2:2308:18049:12744,DHFC08P1:370:C1784ACXX:1:1306:14988:50848,DHFC08P1:370:C1784ACXX:1:1206:21023:31009,DHFC08P1:370:C1784ACXX:1:2302:6153:22406,DHFC08P1:370:C1784ACXX:1:1313:15436:14374,DHFC08P1:370:C1784ACXX:2:2314:12830:49664,DHFC08P1:370:C1784ACXX:1:2307:12504:27551,DHFC08P1:370:C1784ACXX:1:1115:10579:10756sc,DHFC08P1:370:C1784ACXX:1:1214:9958:99422,HWI-ST115:415:C170UACXX:6:2104:9922:86649sc,DHFC08P1:370:C1784ACXX:2:1109:9418:21136mu1";
	vector<string> target_vec;
	castle::StringUtils::c_string_multi_split(target_str, ",", target_vec);
	set<string> target_set(target_vec.begin(), target_vec.end());
	cout << "[TEA._output_mate_fa_serial] # target_set: " << target_set.size() << "\n";

	while (local_reader.LoadNextAlignmentCore(local_alignment_entry)) {
		++num_total;
//			if(is_soft_clipped) {
//				if(string::npos == local_alignment_entry.Name.rfind("mu1") &&
//						string::npos == local_alignment_entry.Name.rfind("mu2") &&
//						string::npos == local_alignment_entry.Name.rfind("sc")) {
//					continue;
//				}
//			}
		auto the_pos_itr = a_positive_repeat_map.find(local_alignment_entry.Name);
		auto the_neg_itr = a_negative_repeat_map.find(local_alignment_entry.Name);
		const bool debug = target_set.end() != target_set.find(local_alignment_entry.Name);
		if (a_positive_repeat_map.end() != the_pos_itr) {
			auto& aln_pair = the_pos_itr->second;
			int64_t the_pos = aln_pair.pos;
			int64_t aln_pos = local_alignment_entry.Position + 1;
			if (debug) {
				cout << (boost::format("%s\t%d\t%d<=>%d\t%s\n") % local_alignment_entry.Name % local_alignment_entry.AlignmentFlag % the_pos % aln_pos % local_alignment_entry.QueryBases).str();
			}

//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
			if (aln_pos != the_pos) {
				if (local_alignment_entry.IsReverseStrand()) {
					local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
				}
				auto& the_seq_vec = positive_mate_reads[aln_pair.file_name_prefix];
				the_seq_vec.push_back(local_alignment_entry.QueryBases);
			}
//					}
		}
		if (a_negative_repeat_map.end() != the_neg_itr) {
			auto& aln_pair = the_neg_itr->second;
			int64_t the_pos = aln_pair.pos;
			int64_t aln_pos = local_alignment_entry.Position + 1;
			if (debug) {
				cout << (boost::format("%s\t%d\t%d<=>%d\t%s\n") % local_alignment_entry.Name % local_alignment_entry.AlignmentFlag % the_pos % aln_pos % local_alignment_entry.QueryBases).str();
			}
//					auto the_ref_id_itr = ref_reverse_index.find(aln_pair.chr);
//					if(ref_reverse_index.end() != the_ref_id_itr) {
//						int64_t the_ref_id = the_ref_id_itr->second;
//						if(local_alignment_entry.Position != the_pos && local_alignment_entry.RefID == the_ref_id) {
			if (aln_pos != the_pos) {
				if (local_alignment_entry.IsReverseStrand()) {
					local_alignment_entry.QueryBases = castle::StringUtils::get_reverse_complement(local_alignment_entry.QueryBases);
				}
				auto& the_seq_vec = negative_mate_reads[aln_pair.file_name_prefix];
				the_seq_vec.push_back(local_alignment_entry.QueryBases);
			}
//					}
		}
	}

	local_reader.Close();
//				done_vector[block_id] = 'D';

//				if(verbose) {
//					string a_block_boundary_str = (boost::format("%s %d-%d %d")
//							% local_alignment_entry.Name % local_alignment_entry.RefID % local_alignment_entry.Position
//							% local_alignment_entry.AlignmentFlag).str();
//					if(block_boundary_strs.end() == block_boundary_strs.find(a_block_boundary_str)) {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last-wrong) %s\n")
//								% block_id % a_block_boundary_str).str();
//					} else {
//						cout << (boost::format("[TEA.BAM_to_FASTQ] Block-%d (last) %s\n")
//								% block_id % a_block_boundary_str).str();
//					}
//				} else {
//					size_t n = count(done_vector.begin(), done_vector.end(), 'D');
//					double processed = n/(double)done_vector.size() * 100.0;
//					cout << (boost::format("%.2f %%\n") % processed).str();
//				}
}

}

/* namespace tea */
