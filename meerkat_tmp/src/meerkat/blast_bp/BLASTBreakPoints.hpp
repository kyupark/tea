/*
 * BLASTBreakPoints.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: el174
 */

#ifndef MEERKAT_BLASTBREAKPOINTS_HPP_
#define MEERKAT_BLASTBREAKPOINTS_HPP_
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/filter/gzip.hpp>

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"
#include "../../third/xfaidx.h"
#include "../ClusterEntry.hpp"
#include "../BWACaller.hpp"
#include "../../third/gzstream.h"


namespace meerkat {
using namespace std;
class BLASTBreakPoints {
public:
	BLASTBreakPoints();
	~BLASTBreakPoints();
	void set_option_parser(const castle::OptionParser& the_options);
	void find_break_points();
	void collect_cluster_name_read_name_map(boost::unordered_map<string, string>& a_read_map);
	void collect_all_local_reads_detail(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist);
	void collect_read_name_sequence_map(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist, const string& a_gz_file_name);
	void collect_read_name_sequence_map_alt(boost::unordered_map<string, boost::unordered_map<string, string>>& a_reads_detail, const boost::unordered_map<string, string>& readslist, const string& a_gz_file_name);
	// intra
	void process_intra_chromosomal_events(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void process_intra_chromosomal_events_alt(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_del(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_inssu(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_inssd(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_insou(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_insod(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_invers(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_del_invers(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_tandem_dup(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_invers_f(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_invers_r(const string& blast_dir, ofstream& INTRAREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);

	// inter
	void process_inter_chromosomal_events(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void process_inter_chromosomal_events_alt(const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_inss(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_inso(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	void write_transl_inter(const string& blast_dir, ofstream& INTERREFINE, const vector<string>& data, const vector<string>& mpd_id, IndexedFasta& fai, const boost::unordered_map<string, string>& readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& readsdetail, const string& n);
	string run_blast(const string blastdir, const vector<string>& ref_data, const string& mpd_id, const string& name, const string& seq, const boost::unordered_map<string, boost::unordered_map<string, string>>& ref_readsdetail, bool forward);

	void parse_blast(vector<int64_t>& ref_bps, vector<int64_t>& ref_bps2, const string& bltoutfile, const boost::unordered_map<string, string>& ref_readslist, const boost::unordered_map<string, boost::unordered_map<string, string>>& ref_readsdetail);

	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};
inline bool BLASTBreakPoints::covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
//		return 0 if ( $a1 =~ /\D/ or $a2 =~ /\D/ or $b1 =~ /\D/ or $b2 =~ /\D/ );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( $a1, $a2 ) = sort { $a <=> $b } ( $a1, $a2 );
//		( $b1, $b2 ) = sort { $a <=> $b } ( $b1, $b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		return true;
	}
	return false;
}
} /* namespace meerkat */

#endif /* MEERKAT_BLASTBREAKPOINTS_HPP_ */
