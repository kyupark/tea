/*
 * TEAOptionParser.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: el174
 */

#ifndef TEA_TEAOPTIONPARSER_HPP_
#define TEA_TEAOPTIONPARSER_HPP_
#include <string>
#include <vector>
#include <iostream>
#include <map>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "../castle/TimeChecker.hpp"

namespace tea {
using namespace std;

class TEAOptionParser {
public:
	void show_help();
	TEAOptionParser();
	TEAOptionParser(int argc, char **argv);
	~TEAOptionParser();
	bool is_complete() const;
	void expand_path(string& a_path);
	void expand_home_path(string& a_path);
	string get_working_path(const string& a_path);
public:
	map<string, int64_t> sub_name_map;
	string program_name;
	string prefix;
	string naive_prefix;
	string working_dir;
	string working_prefix;

	string isinfoname;
	string repeat_reference;
	string human_reference;
	string aln_param;
	string samse_param;

	// hg18, hg19 etc
	string ref;
	string rannot_file;
	string vannot_file;
	string rmasker_rfile;
	string rasym;

	string assembler;
	string assembler_param;
	string start_step;
	string sub_module;

	// comparison
	string l1;
	string l2;
	string comp_out;

	int32_t n_cores;
	int32_t qcutoff;

	int32_t max_mismatches;
	int32_t min_matches;
	int32_t min_polyAT;
	int32_t min_ram;
	int32_t ram_cutoff;
	int32_t jittering;
	int32_t bp_margin;

	int32_t min_acr;
	double min_acrr;
	int32_t min_tsd;
	int32_t max_tsd;
	int32_t min_out_conf;
	int64_t min_clipped_len;

	bool is_force;
	bool no_clipped;
	bool oneside_ram;
	bool exo;
	bool bam_chr;
	bool merge_family;
	bool annot_oi;
	bool stringent_pair;
	bool annot_gene;
	bool no_oi;
	bool is_transduction;
	bool is_orphan;
	bool is_cleaning;
	bool is_sampe;

};

} /* namespace tea */

#endif /* TEA_TEAOPTIONPARSER_HPP_ */
