/*
 * TEAOptionParser.cpp
 *
 *  Created on: Aug 17, 2016
 *      Author: el174
 */

#include "TEAOptionParser.hpp"

namespace tea {

void TEAOptionParser::show_help() {
	if(prefix.empty() || "-" == prefix) {
		prefix = "NONE";
	}
	if(working_dir.empty()) {
		working_dir = "NONE";
	}

	cout << "==============================================================================\n"
			"TEA: Transposable Element Analysis (Ver. 2.0.2)\n"
            "==============================================================================\n"
			"Usage:                       tea <Command> [Options]\n"
			"Commands:                    tea, tea_u, tea_v, transduction, orphan, post\n"
			"       tea                   basic module\n"
			"       tea_u                 module for unmapped [experimental]\n"
			"       tea_v                 module for viral integration [experimental]\n"
			"       transduction          module for detecting transduction [experimental]\n"
			"       orphan                module for detecting orphan [experimental]\n"
			"       post                  post-processing for the transduction, and orphan modules [experimental]\n"
			"       comp                  compare results\n\n"
			"Options:\n\n"
			"       -x [STR]              the absolute path prefix of the INPUT BAM file excluding \".bam\" suffix. e.g.) /home/test.bam -> /home/test [" << prefix << "]\n";
	cout << boolalpha << "       -D [STR]              the path to output folder [" << working_dir << "]\n"
			"       -t [INT]              # cores for parallelization [" << n_cores << "]\n\n"
			"       -M [STR]              the path to a repeat masker text file [NONE]\n"
			"       -N [STR]              the path to a repeat annotation file [NONE]\n"
			"       -V [STR]              the path to a virus annotation file [NONE]\n"
			"       -r [STR]              the path to a repeat reference file [NONE]\n"
			"       -hr [STR]             the path to a human reference file [NONE]\n"

			"       --step [STR]          the name of the start step [fasta, rabam, ram, cbam, rid]\n"
			"       --sub [STR]           the name of a single step to run [fasta, rabam, ram, cbam, rid]\n"
			"       --force               override the previous results [" << is_force << "]\n"
			"       --clean               remove some intermediate files after calling all the steps [" << is_cleaning<< "]\n"
			"       --sampe               the BAM is mapped by bwa sampe [" << is_sampe << "]\n\n"

			"       -a [INT]              the minimum number of bases defining the polyAT [" << min_polyAT << "]\n"
			"       -c [INT]              the quality cutoff [" << qcutoff << "]\n"
			"       -l [STR]              bwa aln parameter [" << aln_param << "]\n"
			"       -m [INT]              the minimum number of matches while collecting a clipped read if the MD tag exists [" << min_matches << "]\n"

			"       -n [INT]              the maximum number of mismatches [" << max_mismatches << "]\n"
			"       -R [STR]              the reference name [" << ref << "]\n"
			"       -s [STR]              samse parameter [" << samse_param << "]\n\n"

			"       --min_ram [INT]       the minimum number of supporting ram [" << min_ram << "]\n"
			"       --oneside_ram         include one-sided ram clusters [" << oneside_ram << "]\n"
			"       --exo                 turn on exogeneous mode [" << exo << "]\n"
			"       --oi                  yields out of instances entries [" << no_oi << "]\n"
			"       --min_out_conf [INT]  the minimum value of confidence score [" << min_out_conf << "]\n"
			"       --ram_cutoff [INT]    the number of minimum ram cutoff [" << ram_cutoff << "]\n"
			"       --jittering [INT]     the number of bases to adjust clipped position [" << jittering << "]\n"
			"       --bp_margin [INT]     the number of bases to define a window within which a partner breakpoint (clipped position) would be sought in either forward or reverse strand [" << bp_margin << "]\n"
			"       --min_acr [INT]       the minimum value of acr [" << min_acr << "]\n"
			"       --min_acrr [INT]      the minimum value of acrr [" << min_acrr << "]\n"
			"       --max_tsd [INT]       the maximum value of tsd [" << max_tsd << "]\n"
			"       --rasym               the name of module [" << rasym << "]\n"
			"       --transduction        the post analysis is transduction [" << is_transduction << "]\n"
			"       --orphan              the post analysis is orphan [" << is_orphan << "]\n"
			"       --l1                  the left .germline file name [" << l1 << "]\n"
			"       --l2                  the right .germlime file name [" << l2 << "]\n"
			"       --comp                the comparison results of the .germlime files [" << comp_out << "]\n\n"
			"Note:                        Please set the tea_base environment to which the TEA is installed. [" << string(getenv("tea_base")) << "]\n"
			"Contact:                     Euncheon Lim <abysslover@gmail.com>\n";
	exit(0);
}
TEAOptionParser::TEAOptionParser() :
		prefix("-"), aln_param("-l 40 -k 2 -n 3"), samse_param("-n 1000"), ref("hg19"), rasym("ra"), assembler("cap3"),
		assembler_param("-i 21 -j 31 -o 16 -s 251 -p 70"), start_step("FASTA"),
		n_cores(1), qcutoff(2), max_mismatches(9), min_matches(25), min_polyAT(10),
		min_ram(3), ram_cutoff(6), jittering(2), bp_margin(50), min_acr(2), min_acrr(0.4), min_tsd(-20), max_tsd(50), min_out_conf(5), min_clipped_len(25),
		is_force(false), no_clipped(false), oneside_ram(false), exo(false), bam_chr(true), merge_family(true), annot_oi(true), stringent_pair(false), annot_gene(true), no_oi(true), is_transduction(false), is_orphan(false), is_cleaning(false), is_sampe(false) {
	ram_cutoff = min_ram;
	sub_name_map["FASTA"] = 0;
	sub_name_map["fasta"] = 0;
	sub_name_map["RABAM"] = 1;
	sub_name_map["rabam"] = 1;
	sub_name_map["RAM"] = 2;
	sub_name_map["ram"] = 2;
	sub_name_map["CBAM"] = 3;
	sub_name_map["cbam"] = 3;
	sub_name_map["RID"] = 4;
	sub_name_map["rid"] = 4;
}

TEAOptionParser::TEAOptionParser(int argc, char **argv) : prefix("-"), aln_param("-l 40 -k 2 -n 3"), samse_param("-n 1000"), ref("hg19"), rasym("ra"), assembler("cap3"),
		assembler_param("-i 21 -j 31 -o 16 -s 251 -p 70"), start_step("FASTA"),
		qcutoff(2), max_mismatches(9), min_matches(25), min_polyAT(10),
		min_ram(3), ram_cutoff(6), jittering(2), bp_margin(50), min_acr(2), min_acrr(0.4), min_tsd(-20), max_tsd(50), min_out_conf(5), min_clipped_len(25),
		is_force(false), no_clipped(false), oneside_ram(false), exo(false), bam_chr(true), merge_family(true), annot_oi(true), stringent_pair(false), annot_gene(true), no_oi(true), is_transduction(false), is_orphan(false), is_cleaning(false), is_sampe(false) {
	string cmd = "[TEAOptionParser.TEAOptionParser] cmd: ";
	for(int i = 0; i < argc; ++i) {
		cmd += string(argv[i]) + " ";
	}
	cout << cmd << "\n";
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
	sub_name_map["FASTA"] = 0;
	sub_name_map["fasta"] = 0;
	sub_name_map["RABAM"] = 1;
	sub_name_map["rabam"] = 1;
	sub_name_map["RAM"] = 2;
	sub_name_map["ram"] = 2;
	sub_name_map["CBAM"] = 3;
	sub_name_map["cbam"] = 3;
	sub_name_map["RID"] = 4;
	sub_name_map["rid"] = 4;
	string line;
	vector<string> cols;
	if(argc > 1) {
		program_name = argv[1];
	}
	sub_module = "all";
	if("tea_v" == program_name) {
		rasym = "va";
		no_clipped = true;
		exo = true;
	}
	for (int i = 1; i < argc; ++i) {
		try {
			string argument = argv[i];
			if ("-a" == argument) {
				string value(argv[i + 1]);
				min_polyAT = boost::lexical_cast<int32_t>(value);
			} else if ("-c" == argument) {
				string value(argv[i + 1]);
				qcutoff = boost::lexical_cast<int32_t>(value);
			} else if ("-D" == argument) {
				string value(argv[i + 1]);
				if("null" != value) {
					working_dir = value;
					expand_path(working_dir);
					if('/' != working_dir[working_dir.size() - 1]) {
						working_dir += "/";
					}
				}
			} else if("-l" == argument) {
				string value(argv[i + 1]);
				aln_param = value;
			} else if ("-m" == argument) {
				string value(argv[i + 1]);
				min_matches = boost::lexical_cast<int32_t>(value);
			} else if("-M" == argument) {
				string value(argv[i + 1]);
				rmasker_rfile = value;
			} else if ("-n" == argument) {
				string value(argv[i + 1]);
				max_mismatches = boost::lexical_cast<int32_t>(value);
			} else if("-N" == argument) {
				string value(argv[i + 1]);
				rannot_file = value;
			} else if("-r" == argument) {
				string value(argv[i + 1]);
				repeat_reference = value;
			} else if("-hr" == argument) {
				string value(argv[i + 1]);
				human_reference = value;
			}
			// hg18, hg19 etc
			else if("-R" == argument) {
				string value(argv[i + 1]);
				ref = value;
			} else if("-s" == argument) {
				string value(argv[i + 1]);
				samse_param = value;
			} else if ("-t" == argument) {
				string value(argv[i + 1]);
				n_cores = min(n_cores, boost::lexical_cast<int32_t>(value));
			} else if("-V" == argument) {
				string value(argv[i + 1]);
				vannot_file = value;
			}  else if ("-x" == argument) {
				string value(argv[i + 1]);
				prefix = value;
			} else if ("-h" == argument) {
				show_help();
			} else if("--force" == argument) {
				is_force = true;
			} else if("--min_ram" == argument) {
				string value(argv[i + 1]);
				min_ram = boost::lexical_cast<int32_t>(value);
			} else if("--ram_cutoff" == argument) {
				string value(argv[i + 1]);
				ram_cutoff = boost::lexical_cast<int32_t>(value);
			} else if("--jittering" == argument) {
				string value(argv[i + 1]);
				jittering = boost::lexical_cast<int32_t>(value);
			} else if("--bp_margin" == argument) {
				string value(argv[i + 1]);
				bp_margin = boost::lexical_cast<int32_t>(value);
			} else if("--min_acr" == argument) {
				string value(argv[i + 1]);
				min_acr = boost::lexical_cast<int32_t>(value);
			} else if("--min_acrr" == argument) {
				string value(argv[i + 1]);
				min_acrr = boost::lexical_cast<double>(value);
			} else if("--max_tsd" == argument) {
				string value(argv[i + 1]);
				max_tsd = boost::lexical_cast<int32_t>(value);
			} else if("--min_out_conf" == argument) {
				string value(argv[i + 1]);
				min_out_conf = boost::lexical_cast<int32_t>(value);
			} else if("--oneside_ram" == argument) {
				oneside_ram = true;
			} else if("--exo" == argument) {
				exo = true;
			} else if("--oi" == argument) {
				no_oi = false;
			} else if("--rasym" == argument) {
				string value(argv[i + 1]);
				rasym = value;
			} else if("--sub" == argument) {
				string value(argv[i + 1]);
				sub_module = value;
			}  else if("--step" == argument) {
				string value(argv[i + 1]);
				start_step = value;
			} else if("--clean" == argument) {
				is_cleaning = true;
			} else if("--sampe" == argument) {
				is_sampe = true;
			} else if("--transduction" == argument) {
				is_transduction = true;
			} else if("--orphan" == argument) {
				is_orphan = true;
			} else if("--l1" == argument) {
				string value(argv[i + 1]);
				l1 = value;
			} else if("--l2" == argument) {
				string value(argv[i + 1]);
				l2 = value;
			} else if("--comp" == argument) {
				string value(argv[i + 1]);
				comp_out = value;
			}
		} catch (exception& ex) {
			cout << ex.what() << "\n";
			cout << argv[i] << ": " << argv[i + 1] << "\n";
			exit(1);
		}
	}
	if (working_dir.empty()) {
		string sep = prefix.empty() ? "" : ".";
		isinfoname = prefix + sep + "isinfo";
	} else {
		auto delim_pos = prefix.rfind("/");
		if ((prefix.size() - 1) == delim_pos) {
			delim_pos = prefix.rfind("/", prefix.size() - 2);
			working_prefix = prefix.substr(delim_pos + 1);
		} else {
			if (string::npos == delim_pos) {
				working_prefix = working_dir + prefix.substr(0);
			} else {
				working_prefix = working_dir + prefix.substr(delim_pos + 1);
			}
		}
		if ("-" != prefix) {
			if (!boost::filesystem::exists(working_prefix)) {
				boost::filesystem::create_directories(working_prefix);
			}
		}
		string sep = prefix.empty() ? "" : ".";
		isinfoname = working_prefix + sep + "isinfo";
	}

	naive_prefix = prefix;
	if (!working_dir.empty()) {
		naive_prefix = working_prefix;
	}
	size_t the_slash_pos = naive_prefix.rfind("/");
	if (string::npos != the_slash_pos) {
		naive_prefix = naive_prefix.substr(the_slash_pos + 1);
	}

	if (exo) {
		oneside_ram = true;
		merge_family = false;
		annot_oi = false;
	}

	if ("um" == rasym) {
		stringent_pair = true;
	} else {
		stringent_pair = false;
	}

	if (!is_complete()) {
		show_help();
	}
}

TEAOptionParser::~TEAOptionParser() {
}

bool TEAOptionParser::is_complete() const {
	if("comp" == program_name) {
		return !l1.empty() && !l2.empty() && !comp_out.empty();
	} else {
		return !prefix.empty() && "-" != prefix;
	}
}
void TEAOptionParser::expand_home_path(string& a_path) {
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

void TEAOptionParser::expand_path(string& a_path) {
	boost::filesystem::path current_path(boost::filesystem::current_path());
	if ('/' != a_path[0] && '~' != a_path[0]) {
		a_path = current_path.string() + "/" + a_path;
	} else if ('~' == a_path[0]) {
		expand_home_path(a_path);
	}
}
string TEAOptionParser::get_working_path(const string& a_path) {
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
} /* namespace tea */
