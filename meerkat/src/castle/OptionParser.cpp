/*
 * OptionParser.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "OptionParser.hpp"

namespace castle {

void OptionParser::show_help() {
	if ("dre" == program_name) {
		show_help_dre();
	}
	cout << "Meerkat Preprocess: Parallelized bamreader (Ver. 0.0.1) \n"
			"==============================================================================\n";
	cout << "Syntax:    pmeerkat preprocess [options] <BAM input file>\n\n" << "Options:\n"
			"  -b cutoff  Generate a blacklist of genomic positions with coverage\n"
			"             greater than 'cutoff'.  The input BAM must be sorted if\n"
			"             -b is specified with cutoff > 0.  Set the cutoff to <= 0\n"
			"             to disable blacklist construction.  [0]\n"
			"             WARNING: currently, the code does not take CIGAR strings\n"
			"             into account when computing coverage.  The 'calDepth'\n"
			"             example in samtools does consider the CIGAR string.\n"
			"  -u         Do not process UU pairs.\n"
			"  -x prefix  All output files are prepended with 'prefix'.  [output]\n"
			"  -g isize   Maximum allowed insert size for computing the insert\n"
			"             size distribution.  This is useful to filter out large\n"
			"             inserts that result from SVs.  [1,000]\n"
			"  -c bps     Defines the cutoff for calling soft clips.  A read is\n"
			"             soft clipped if it contains an 'S' CIGAR operation of\n"
			"             at least 'bps' basepairs.  [20]\n"
			"  -s bps     Fragment lengths for the unmapped and softclipped files.\n"
			"             Defaults to the size specified by '-c'.\n"
			"  -n nreads  Output the insert sizes for the first 'nreads' reads\n"
			"             from each read group to determine the insert size\n"
			"             distribution.  If nreads < 0, all reads in each group\n"
			"             will be used.  (Very large values can cause significant\n"
			"             loss of speed!)  [10,000]\n"
			"  -q qv      Read trimming parameter.  Equivalent to BWA's -q option.  [0]\n"
			"  -N integer Directly before writing a read to any output file (except\n"
			"             .unmapped.fq and .softclips.fq), require that the read\n"
			"             have < 'integer' Ns in its sequence.  [5]\n"
			"  -r file    Ignore reads from any read group in this file.  The file\n"
			"             should specify one read group per line.\n"
			"  -f         Force immediate run.  Do not prompt the user to verify\n"
			"             run-time options.\n"
			"  -v         Verbose mode.\n"
			"  -h         Print this help message.\n\n";
	cout << "Contact:   Euncheon Lim <abysslover@gmail.com>\n";
	exit(1);
}
void OptionParser::show_help_dre() {
	cout << "Meerkat dre: Parallelized discordant read extractor (Ver. 0.0.1) \n"
			"==============================================================================\n";
	cout << "Syntax:    pmeerkat dre [options] <BAM input file>\n\n"
			"Options:\n" //
			"  -P         Remove PCR duplicates.  A set of reads are considered\n"//
			"             PCR duplicates if both mapped positions of the two mates\n"//
			"             are the same.  Currently, the first read in each duplicate\n"//
			"             set is kept and the rest are discarded--we could do better\n"//
			"             here by selecting the highest quality read from each set.\n"//
			"  -m         Ignore reads marked as duplicates in the BAM.  Do NOT mix\n"//
			"             this option with -P.\n"//
			"  -s integer Insert size cutoff for discordant read pairs that map to\n"//
			"             the same chromosome and have the correct orientation.\n"//
			"             The value is the number of standard deviations from the\n"//
			"             median insert size.  [3]\n"//
			"  -b file    A file of blacklisted genomic coordinates.  Reads mapping\n"//
			"             to any coordinate in this file will not be considered\n"//
			"             discordant, even if they are.\n"//
			"  -i file    Specify the file containing insert size distribution info.\n"//
			"  -o file    Output file.  Written in BAM format.\n"//
			"  -d file    PCR duplicate file.  Written in BAM format.\n"//
			"  -r file    Ignore reads from any read group in this file.  The file\n"//
			"             should specify one read group per line.\n"//
			"  -f         Force immediate run.  Do not prompt the user to verify\n"//
			"             run-time options.\n"//
			"  -v         Verbose mode.\n"//
			"  -h         Print this help message.\n\n"//
			"Contact:   Euncheon Lim <abysslover@gmail.com>\n";
	exit(1);
}
OptionParser::OptionParser() :
		prefix("output"), no_prompt(true), verbose(false), processUU(true), max_isize(1000), big_s_bps(20), frag_size(-1), cut_sr(-1), nstdevs(3), isize_samples(1000), q(0), min_read_len(0), coverage_cutoff(0), n_cutoff(5), plotrange(1000), sd_cutoff_disc(3), sd_cutoff_cl(-1), support_mps(-1), support_mpf(-1), min_mapq(
				0), alt_map_max(-1), alt_map_max_clip(numeric_limits<int32_t>::max()), sv_size_cutoff(1000000000), support_reads(1), te_size_max(100000), del_ins_size_cutoff_d(0.8), del_ins_size_cutoff_u(1.2), ovl(0.8), filter_dups(false), filter_dups_by_flag(false), remove_dup(false), clip(true), ad_align(true), use_all_align(false), include_other(false), generate_mapped_sc_um_file(false), is_cleaning(false) {
	program_names.insert("preprocess");
	program_names.insert("dre");
	program_names.insert("r_alt");
	program_names.insert("s_alt");
	program_names.insert("sclus");
	program_names.insert("mpd");
	program_names.insert("alg");
	program_names.insert("srd");
	program_names.insert("filter");
	program_names.insert("blast_bp");
	program_names.insert("mechanism");
	program_names.insert("selection");
	program_names.insert("cns");
	program_names.insert("svs");
	program_names.insert("test");

	tmp_directory = "/tmp";
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

OptionParser::OptionParser(int argc, char **argv) :
		prefix("-"), no_prompt(true), verbose(false), processUU(true), max_isize(1000), big_s_bps(20), frag_size(-1), cut_sr(-1), nstdevs(3), isize_samples(1000), q(0), min_read_len(0), coverage_cutoff(0), n_cutoff(5), plotrange(1000), sd_cutoff_disc(3), sd_cutoff_cl(-1), support_mps(-1), support_mpf(-1), min_mapq(
				0), alt_map_max(-1), alt_map_max_clip(numeric_limits<int32_t>::max()), sv_size_cutoff(1000000000), support_reads(1), te_size_max(100000), del_ins_size_cutoff_d(0.8), del_ins_size_cutoff_u(1.2), ovl(0.8), filter_dups(false), filter_dups_by_flag(false), remove_dup(false), clip(true), ad_align(true), use_all_align(false), include_other(false), generate_mapped_sc_um_file(false), is_cleaning(false) {
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();
	tmp_directory = "/tmp";
	string line;
	vector<string> cols;
//	const char* delim = " ";
	// the first line

	program_names.insert("preprocess");
	program_names.insert("dre");
	program_names.insert("r_alt");
	program_names.insert("s_alt");
	program_names.insert("sclus");
	program_names.insert("mpd");
	program_names.insert("alg");
	program_names.insert("srd");
	program_names.insert("filter");
	program_names.insert("blast_bp");
	program_names.insert("mechanism");
	program_names.insert("selection");

	program_names.insert("cns");
	program_names.insert("svs");
	program_names.insert("test");

	if (argc > 1) {
		string local_program_name = argv[1];
		program_name = local_program_name;
		if (program_names.end() == program_names.find(local_program_name)) {
			program_name = "default";
		}
		cout << (boost::format("[OptionParser.main] program: %s\n") % program_name).str();
	}

	for (int i = 1; i < argc; ++i) {
		try {
			string argument = argv[i];
			if ("-a" == argument) {
				ad_align = false;
				if("dre" == program_name) {
					string value(argv[i + 1]);
					infile_name = value;
				}
			} else if ("-b" == argument) {
				string value(argv[i + 1]);
				if ("s_alt" == program_name || "r_alt" == program_name) {
					fname = value;
				} else if ("dre" == program_name) {
					bfile_name = value;
					read_blacklist();
				} else {
					coverage_cutoff = boost::lexical_cast<int>(value);
				}
			} else if ("-c" == argument) {
				string value(argv[i + 1]);
				if("cns" == program_name) {
					input_cnv_filename = value;
					continue;
				}
				big_s_bps = boost::lexical_cast<int>(value);
				sd_cutoff_cl = big_s_bps;
			} else if ("--clean" == argument) {
				is_cleaning = true;
			} else if ("-d" == argument) {
				string value(argv[i + 1]);
				if ("dre" == program_name) {
					dupfile_name = value;
				} else {
					sd_cutoff_disc = boost::lexical_cast<int32_t>(value);
				}
			} else if ("-D" == argument) {
				string value(argv[i + 1]);
				if("null" != value) {
					working_dir = value;
					expand_path(working_dir);
					if('/' != working_dir[working_dir.size() - 1]) {
						working_dir += "/";
					}
				}
			} else if ("-f" == argument) {
				if("dre" == program_name) {
					no_prompt = true;
				} else {
					string value(argv[i + 1]);
					alt_map_max = boost::lexical_cast<int64_t>(value);
				}
			} else if ("-F" == argument) {
				string value(argv[i + 1]);
				reference_path = value;
			} else if ("-g" == argument) {
				string value(argv[i + 1]);
				alt_map_max_clip = boost::lexical_cast<int64_t>(value);
				max_isize = boost::lexical_cast<int>(value);
			} else if ("-i" == argument) {
				input_filename = string(argv[i + 1]);
				if ("dre" == program_name) {
					insfile_name = input_filename;
					read_isinfo_dre();
				}
			} else if ("-l" == argument) {
				string value(argv[i + 1]);
				if ("1" == value) {
					clip = true;
				}
			} else if ("-L" == argument) {
				string value(argv[i + 1]);
				plotrange = boost::lexical_cast<int32_t>(value);
			} else if ("-m" == argument) {
				string value(argv[i + 1]);
				if ("dre" == program_name) {
					filter_dups_by_flag = true;
				}
				if ("1" == value) {
					remove_dup = true;
				} else if("0" == value) {
					filter_dups = false;
				}
			} else if ("-M" == argument) {
				generate_mapped_sc_um_file = true;
			} else if ("-n" == argument) {
				string value(argv[i + 1]);
				isize_samples = boost::lexical_cast<int>(value);
			} else if ("-N" == argument) {
				string value(argv[i + 1]);
				n_cutoff = boost::lexical_cast<int>(value);
			} else if ("-o" == argument) {
				string value(argv[i + 1]);
				if ("dre" == program_name || "cns" == program_name || "svs" == program_name) {
					outfile_name = value;
				} else if("mechanism" == program_name) {
					include_other = true;
				} else {
					support_mpf = boost::lexical_cast<int32_t>(value);
				}
			} else if ("-p" == argument) {
				string value(argv[i + 1]);
				support_mps = boost::lexical_cast<int32_t>(value);
			} else if ("-P" == argument) {
				filter_dups = true;
			} else if ("-q" == argument) {
				string value(argv[i + 1]);
				q = boost::lexical_cast<int>(value);
				support_reads = boost::lexical_cast<int64_t>(value);
			} else if ("-Q" == argument) {
				string value(argv[i + 1]);
				min_mapq = boost::lexical_cast<int32_t>(value);
			} else if ("-r" == argument) {
				if("tea" == program_name || "tea_v" == program_name) {
					continue;
				}
				string value(argv[i + 1]);
				rg_blacklist_fname = value;
				read_rg_blacklist();
			} else if ("-R" == argument) {
				rmskfile = string(argv[i + 1]);
			} else if ("-s" == argument) {
				string value(argv[i + 1]);
				if("cns" == program_name) {
					input_sv_filename = value;
					continue;
				}
				if("tea" == program_name || "tea_v" == program_name) {
					continue;
				}
				frag_size = boost::lexical_cast<int>(value);
				cut_sr = frag_size;
				nstdevs = frag_size;
			} else if("-s1" == argument) {
				string value(argv[i + 1]);
				input_sv1_filename = value;
			} else if("-s2" == argument) {
				string value(argv[i + 1]);
				input_sv2_filename = value;
			} else if("-t" == argument) {
				string value(argv[i + 1]);
				n_cores = min(n_cores, boost::lexical_cast<int32_t>(value));
			} else if ("-T" == argument) {
				string value(argv[i + 1]);
				te_size_max = boost::lexical_cast<int32_t>(value);
			} else if ("-tmp" == argument) {
				boost::filesystem::path value(argv[i + 1]);
				tmp_directory = boost::filesystem::absolute(value).string();
			} else if ("-u" == argument) {
				processUU = false;
				use_all_align = true;
			} else if ("-v" == argument) {
				verbose = true;
			} else if ("-x" == argument) {
				prefix = string(argv[i + 1]);
			} else if ("-z" == argument) {
				output_filename = string(argv[i + 1]);
			} else if ("-sv" == argument) {
				string value(argv[i + 1]);
				sv_size_cutoff = boost::lexical_cast<int64_t>(value);
			} else if ("-h" == argument) {
				show_help();
			}
		} catch (exception& ex) {
			cout << ex.what() << "\n";
			cout << argv[i] << "/" << argv[i + 1] << "\n";
			exit(1);
		}
	}

	/* No longer an option */
	min_read_len = big_s_bps << 1;
	if (frag_size == -1) {
		frag_size = big_s_bps; /* -s defaults to -c */
	}

	if (1 != argc) {
		fname = argv[argc - 1];
		if(!boost::filesystem::exists(fname)) {
			string tmp_name = working_dir + fname;
			if(boost::filesystem::exists(tmp_name)) {
				fname = tmp_name;
			}
		}
		infile_name = fname;
	}
	if (-1 == sd_cutoff_cl) {
		sd_cutoff_cl = sd_cutoff_disc;
	}
	if (use_all_align) {
		ad_align = false;
	}

	/* Forcefully prevent mixing of -P and -m flags */
	if (filter_dups_by_flag) {
		if (filter_dups) {
			cerr << "WARNING: you specified both -P and -m.  Ignoring -P.\n";
		}
		filter_dups = false;
	}
	if (!filter_dups) {
		dupfile_name = ""; /* discard whatever they might've input */
	}
	if (!is_complete()) {
		show_help();
	}

	if("dre" == program_name) {
		auto the_bam_pos = infile_name.rfind(".bam");
		if(string::npos == the_bam_pos) {
			prefix = infile_name;
		} else {
			prefix = infile_name.substr(0, the_bam_pos);
		}
	}

	if(working_dir.empty()) {
		if("tea" != program_name && "tea_v" != program_name) {
			cout << (boost::format("[OptionParser.OptionParser] working directory: %s\n") % prefix).str();
		}
	string sep = prefix.empty() ? "" : ".";
	umfname = prefix + sep + "unmapped.fq.gz";
	clipname = prefix + sep + "softclips.fq.gz";
	isinfoname = prefix + sep + "isinfo";
	umrdistname = prefix + sep + "unmapped.rdist";
	scrdistname = prefix + sep + "softclips.rdist";
	blistname = prefix + sep + "blacklist.gz";
	split1name = prefix + sep + "sr.1.fq.gz";
	split2name = prefix + sep + "sr.2.fq.gz";

	} else {
		auto delim_pos = prefix.rfind("/");

		if((prefix.size() - 1) == delim_pos) {
			delim_pos = prefix.rfind("/", prefix.size() - 2);
			working_prefix = prefix.substr(delim_pos + 1);
		} else {
			if(string::npos == delim_pos) {
				working_prefix = working_dir + prefix.substr(0);
			} else {
				working_prefix = working_dir + prefix.substr(delim_pos + 1);
			}
		}
		if("-" != prefix) {
			if(!boost::filesystem::exists(working_prefix)) {
				boost::filesystem::create_directories(working_prefix);
			}
		}
		cout << (boost::format("[OptionParser.OptionParser] working directory: %s\n") % working_dir).str();
		cout << (boost::format("[OptionParser.OptionParser] working prefix: %s\n") % working_prefix).str();
		string sep = prefix.empty() ? "" : ".";
		umfname = working_prefix + sep + "unmapped.fq.gz";
		clipname = working_prefix + sep + "softclips.fq.gz";
		isinfoname = working_prefix + sep + "isinfo";
		umrdistname = working_prefix + sep + "unmapped.rdist";
		scrdistname = working_prefix + sep + "softclips.rdist";
		blistname = working_prefix + sep + "blacklist.gz";
		split1name = working_prefix + sep + "sr.1.fq.gz";
		split2name = working_prefix + sep + "sr.2.fq.gz";
	}
	if ("preprocess" == program_name) {
		print_options();
	} else if ("dre" == program_name) {
		print_options_dre();
	} else {
		read_isinfo();
		read_point_black_lists();
	}

	if("-" != prefix && "dre" != program_name) {
		if (!boost::filesystem::exists(prefix) && !boost::filesystem::exists(working_prefix)) {
			cout << "ERROR: could not create directory '" << prefix << "'\n";
			exit(1);
		}
	}
	if("/tmp" == tmp_directory && !working_dir.empty()) {
		tmp_directory = working_prefix + "/tmp";
		if (!boost::filesystem::exists(tmp_directory)) {
			if (!boost::filesystem::create_directories(tmp_directory)) {
				cout << "ERROR: could not create temporary directory '" << prefix << "'\n";
			}
		}
	}
}

set<string> OptionParser::read_rg_blacklist() {
	if (!boost::filesystem::exists(rg_blacklist_fname)) {
		cout << "warning: cannot open readgroup blacklist '" << rg_blacklist_fname << "'\n";
		return rg_blacklist;
	}
	ifstream f(rg_blacklist_fname, ios::binary);
	if (f.is_open()) {
		string line;
		while (getline(f, line)) {
			rg_blacklist.insert(line);
		}
	} else {
		cout << "warning: cannot open readgroup blacklist '" << rg_blacklist_fname << "'\n";
	}
	return rg_blacklist;
}

void OptionParser::read_blacklist() {
	if (bfile_name.empty()) {
		return;
	}
	istream *bfile;
	igzstream gzf;
	ifstream f;
	bool is_gz = false;
	if (".gz" == bfile_name.substr(bfile_name.size() - 3, 3)) {
		is_gz = true;
	}

	if (is_gz) {
		gzf.open(bfile_name.c_str());
		bfile = &gzf;
	} else {
		f.open(bfile_name.c_str(), ios::binary);
		bfile = &f;
	}

	char buf[100];
	long pos, cov;

	while (bfile->good()) {
		*bfile >> buf;
		*bfile >> pos;
		*bfile >> cov;

		string chr(buf);
		memset(&buf[0], 0, 100);
		blacklist[chr][pos] = true;
	}

	if (is_gz)
		((igzstream *) bfile)->close();
	else
		((ifstream *) bfile)->close();

}

void OptionParser::read_point_black_lists() {
	string line;
	igzstream in(blistname.c_str());
	const char* delims = "\t";
	vector<string> a_cols;
	while (getline(in, line, '\n')) {
		StringUtils::c_string_multi_split(line, delims, a_cols);
		point_black_lists[make_pair(a_cols[0], boost::lexical_cast<int64_t>(a_cols[1]))] = boost::lexical_cast<int64_t>(a_cols[2]);
	}
}
void OptionParser::read_isinfo() {
	string line;
	ifstream in(isinfoname, ios::binary);
	const char* delims = "\t";
	vector<string> a_cols;
	while (getline(in, line, '\n')) {
		if (string::npos != line.find("Read length")) {
			StringUtils::c_string_multi_split(line, delims, a_cols);
			string rg = a_cols[a_cols.size() - 1];
			getline(in, line, '\n');
			double tmp = boost::lexical_cast<double>(line);
			if (is["rlu"].end() == is["rlu"].find("selected")) {
				is["rlu"]["selected"] = tmp;
			} else {
				if (tmp > is["rlu"]["selected"]) {
					is["rlu"]["selected"] = tmp;
				}
			}
			if (is["rld"].end() == is["rld"].find("selected")) {
				is["rld"]["selected"] = tmp;
			} else {
				if (is[rg]["rl"] < is["rld"]["selected"]) {
					is["rld"]["selected"] = tmp;
				}
			}
		} else if (string::npos != line.find("Median")) {
			StringUtils::c_string_multi_split(line, delims, a_cols);
			string rg = a_cols[a_cols.size() - 1];
			getline(in, line, '\n');
			double tmp = boost::lexical_cast<double>(line);
			is[rg]["median"] = tmp;
		} else if (string::npos != line.find("Standard deviation")) {
			StringUtils::c_string_multi_split(line, delims, a_cols);
			string rg = a_cols[a_cols.size() - 1];
			getline(in, line, '\n');
			double tmp_sd = boost::lexical_cast<double>(line);
			is[rg]["sd"] = tmp_sd;
			double tmp_isu = is[rg]["median"] + tmp_sd * sd_cutoff_cl + 1;
			is[rg]["isu"] = tmp_isu;
			double tmp_isd = is[rg]["median"] - tmp_sd * sd_cutoff_cl - 1;
			is[rg]["isd"] = tmp_isd;
			if (is["isu"].end() == is["isu"].find("selected")) {
				is["isu"]["selected"] = tmp_isu;
			} else {
				if (is[rg]["isu"] > is["isu"]["selected"]) {
					is["isu"]["selected"] = tmp_isu;
				}
			}
		}
	}
}

void OptionParser::read_isinfo_dre() {
	ifstream isf(insfile_name, ios::binary);
	if (isf.fail()) {
		cout << "read_isinfo: failed to open insert info file '" << insfile_name << "'\n";
		exit(1);
	}

	string line;
	string rg;
	while (true) {
		getline(isf, line);
		if (isf.eof()) {
			break;
		}

		if (string::npos != line.find("Median insert size")) {
			rg = get_rg(line);
			getline(isf, line);
			double val = strtod(line.c_str(), NULL);
			auto iter = isinfo.find(rg);
			if (iter != isinfo.end()) {
				isinfo[rg].first = val;
			} else {
				isinfo[rg].first = val;
				isinfo[rg].second = 0;
			}
		} else if (string::npos != line.find("Standard deviation of insert size")) {
			rg = get_rg(line);
			getline(isf, line);
			double val = strtod(line.c_str(), NULL);
			auto iter = isinfo.find(rg);
			if (isinfo.end() != iter) {
				isinfo[rg].second = val;
			} else {
				isinfo[rg].first = 0;
				isinfo[rg].second = val;
			}
		}
	}

}
void OptionParser::print_options() {
	cout << "Options:\n   input BAM: " << fname << "\n   prefix: " << prefix << "\n";
	if (coverage_cutoff > 0) {
		cout << "   generating blacklist: coverage >= " << coverage_cutoff << " reads\n";
	} else {
		cout << "   no blacklist will be generated\n";
	}
	cout << "   max insert size: " << max_isize << " bases\n" << "   soft clip threshold: " << big_s_bps << " bases\n" << "   fragment size: " << frag_size << " bases\n" << "   insert size samples: ";
	if (isize_samples < 0) {
		cout << "all";
	} else {
		cout << isize_samples;
	}
	cout << " reads\n";
	cout << "   blacklisted read groups: ";

	for (auto it = rg_blacklist.begin(); it != rg_blacklist.end(); ++it) {
		cout << *it << " ";
	}
	cout << "\n";
	cout << "   trim reads: " << q << "\n   max N bps per read: " << n_cutoff << "\n";
	cout << "   output files: \n      " << umfname << "\n      " << clipname << "\n      " << isinfoname << "\n      " << umrdistname << "\n      " << scrdistname << "\n      " << split1name << "\n      " << split2name << "\n";
	if (coverage_cutoff > 0) {
		cout << "      " << blistname << "\n";
	}
	cout << "      " << prefix << "/ (directory)" << "\n";

	if (!no_prompt) {
		cout << "Continue? (y/N) ";
		char c;
		cin >> noskipws >> c;
		if (tolower(c) != 'y') {
			cout << "Terminating...\n";
			exit(1);
		}
	} else {
		cout << "Forcing immediate execution.\n";
	}
}

void OptionParser::print_options_dre() {
	/* Options blurb */
	cout << "Options:\n" << "   input BAM: " << f(infile_name) << "\n" << "   insert file: " << f(insfile_name) << "\n";
	map<string, pair<double, double> >::iterator iter = isinfo.begin();
	while (iter != isinfo.end()) {
		cout << "      " << iter->first << ": median " << iter->second.first << ", stdev " << iter->second.second << "\n";
		++iter;
	}
	cout << "   blacklist file: " << f(bfile_name) << "\n";
	if (blacklist.size() > 0) {
		int npos = 0;
		blacklist_t::iterator iter = blacklist.begin();
		while (iter != blacklist.end()) {
			npos += iter->second.size();
			++iter;
		}

		cout << "      " << blacklist.size() << " chromosomes present in blacklist\n" << "      " << npos << " blacklisted positions\n";
	}
	cout << "   ignoring read groups: ";
	for (auto it = rg_blacklist.begin(); it != rg_blacklist.end(); ++it) {
		cout << *it << " ";
	}
	cout << "\n";

	cout << "   output file: " << f(outfile_name) << "\n   PCR dup filter: " << boolalpha << filter_dups << "\n   PCR duplicate file: " << f(dupfile_name) << "\n   using PCR dup flag: " << boolalpha << filter_dups_by_flag << "\n   n stdevs: " << nstdevs << "\n";

	/* User confirmation before run unless -f was specified */
	if (!no_prompt) {
		cout << "Continue? (y/N) ";
		char c;
		cin >> noskipws >> c;
		if (tolower(c) != 'y') {
			cout << "Terminating...\n";
			exit(1);
		}
	} else {
		cout << "Forcing immediate execution.\n";
	}

}
OptionParser::~OptionParser() {
}

bool OptionParser::is_complete() {
	if ("preprocess" == program_name) {
		return !fname.empty();
	} else if ("dre" == program_name) {
		return !infile_name.empty() && !outfile_name.empty();
	} else if("mechanism" == program_name) {
		return !prefix.empty() && !rmskfile.empty();
	} else {
		return !prefix.empty();
	}
	return true;
}
void OptionParser::expand_home_path(string& a_path) {
	if (a_path.empty() || '~' != a_path[0]) {
		return;
	}
	char const* home = getenv("HOME");
	if (home || ((home = getenv("USERPROFILE")))) {
		a_path.replace(0, 1, home);
	} else {
		char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
		a_path.replace(0, 1, string(hdrive) + hpath);
	}
}

void OptionParser::expand_path(string& a_path) {
	boost::filesystem::path current_path(boost::filesystem::current_path());
	if ('/' != a_path[0] && '~' != a_path[0]) {
		a_path = current_path.string() + "/" + a_path;
	} else if ('~' == a_path[0]) {
		expand_home_path(a_path);
	}
}
string OptionParser::get_working_path(const string& a_path) {
	string copied_path(a_path);
	boost::filesystem::path current_path(boost::filesystem::current_path());
	string file_prefix;
	if ('/' != copied_path[0] && '~' != copied_path[0]) {
		file_prefix = current_path.string();
	} else if ('~' == copied_path[0]) {
		expand_home_path(copied_path);
	}
	file_prefix = copied_path.substr(0, copied_path.rfind('/'));
	return file_prefix;
}

} /* namespace castle */
