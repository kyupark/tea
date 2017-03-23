/*
 * OptionParser.hpp
 *
 *  Created on: Jun 4, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef OPTIONPARSER_HPP_
#define OPTIONPARSER_HPP_
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include "TimeChecker.hpp"
#include "StringUtils.hpp"
#include "../third/gzstream.h"

namespace castle {
	using namespace std;
	typedef map< string, map<long, bool> > blacklist_t;
	class OptionParser {
		public:
			string prefix;
			string working_prefix;
			string rg_blacklist_fname;
			string rmskfile;

			string bfile_name;
			string fname;

			string umfname;
			string clipname;
			string isinfoname;
			string umrdistname;
			string scrdistname;
			string blistname;
			string split1name;
			string split2name;
			string program_name;
			string reference_path;

			string input_filename;
			string output_filename;
			string insfile_name;
			string infile_name;
			string outfile_name;
			string input_cnv_filename;
			string input_sv_filename;
			string input_sv1_filename;
			string input_sv2_filename;

			string dupfile_name;
			string working_dir;
			string tmp_directory;

			blacklist_t blacklist;
			set<string> rg_blacklist;
			map<pair<string, int64_t>, int64_t> point_black_lists;
			unordered_map<string, unordered_map<string, double>> is;
			map<string, pair<double, double> > isinfo;
			set<string> program_names;
			bool no_prompt;
			bool verbose;
			bool processUU;

			int max_isize;
			int big_s_bps;
			int frag_size;
			int cut_sr;
			int nstdevs;
			int isize_samples;
			int q;
			int min_read_len;
			int coverage_cutoff;
			int n_cutoff;
			int32_t plotrange;
			int32_t sd_cutoff_disc;
			int32_t sd_cutoff_cl;
			int32_t support_mps;
			int32_t support_mpf;
			int32_t min_mapq;
			int64_t alt_map_max;
			int64_t alt_map_max_clip;
			int64_t sv_size_cutoff;
			int64_t support_reads;
			int32_t te_size_max;
			double del_ins_size_cutoff_d;
			double del_ins_size_cutoff_u;
			double ovl;

			int32_t n_cores;
			bool filter_dups;
			bool filter_dups_by_flag;
			bool remove_dup;
			bool clip;
			bool ad_align;
			bool use_all_align;
			bool include_other;
			bool generate_mapped_sc_um_file;
			bool is_cleaning;
		public:
			void show_help();
			void show_help_dre();
			OptionParser();
			OptionParser(int argc, char **argv);
			~OptionParser();
		public:
			bool is_complete();
			void expand_path(string& a_path);
			void expand_home_path(string& a_path);
			set<string> read_rg_blacklist();
			void read_blacklist();
			void read_point_black_lists();
			void read_isinfo();
			void read_isinfo_dre();
			void print_options();
			void print_options_dre();
			string get_working_path(const string& a_path);
			string get_working_folder_path() const;
		private:
					/* Utility function for read_isinfo. */
					string get_rg(string& s);
					string f(string& s);
	};
	inline string OptionParser::get_rg(string& s) {
		int i = s.find_last_of(" \r\n\t", s.find(":") + 1);
		return s.substr(i+1, s.size());
	}
	/*
	 * This "none" string actually has to match the "none" readgroup
	 * from bamreader or else the correct insert size statistics
	 * won't be used.
	 */
	inline string OptionParser::f(string& s) { return s.size() == 0 ? "none" : s; }
} /* namespace castle */
#endif /* OPTIONPARSER_HPP_ */
