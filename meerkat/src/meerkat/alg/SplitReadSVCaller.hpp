/*
 * SplitReadSVCaller.hpp
 *
 *  Created on: Jul 8, 2016
 *      Author: el174
 */

#ifndef MEERKAT_SPLITREADSVCALLER_HPP_
#define MEERKAT_SPLITREADSVCALLER_HPP_

#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <string>
#include <limits>

#include <boost/range/adaptor/reversed.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"
#include "../ClusterEntry.hpp"

namespace meerkat {
using namespace std;
class SplitReadSVCaller {
public:
	SplitReadSVCaller();
	~SplitReadSVCaller();
	void set_option_parser(const castle::OptionParser& the_options);
	void call_structural_variants();

	void collect_bp_data(map<string, CoordinateEntry>& cluster_region);
	void detect_intra_chromosomal_events(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region);
	void detect_intra_chromosomal_events_alt(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region);

	// intra
	void detect_del(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
					vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds, vector<string>& temp_cols);
	void detect_inssd(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
				vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_inssu(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
			vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_insod(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
				vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_insou(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
					vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_invers(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
						vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_tandem_dup(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
							vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_invers_f(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
								vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_invers_r(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
									vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	// inter
	void detect_inter_chromosomal_events(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region);
	void detect_inter_chromosomal_events_alt(ofstream& BPREAD, vector<EventEntry>& result_sr, vector<EventEntry>& unsupport_del, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region);
	void detect_inss(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
					vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_inso(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
					vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void detect_transl_inter(ofstream& BPREAD, vector<EventEntry>& result_sr, BamTools::BamReader& reader, const map<string, pair<int64_t, int64_t>>& reverse_index_ref_id, const map<string, CoordinateEntry>& cluster_region, map<int64_t, int32_t>& bp_window,
						vector<string>& data, vector<string>& cl, vector<string>& cluster_ids, vector<string>& mpds);
	void sr_cluster(vector<int64_t>& boundary1, vector<int64_t>& boundary2, vector<int64_t>& support_sr, map<int64_t, vector<string>>& bpread,
			vector<BamTools::BamAlignment>& ref_discord_sr, int64_t cl2, int32_t internal_break);
	void sr_cluster_no_cl(vector<int64_t>& src_return, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& discord_sr);
	void sr_cluster_1(vector<int64_t>& boundary1, vector<int64_t>& boundary2, vector<int64_t>& support_sr, map<int64_t, vector<string>>& bpread,
			vector<BamTools::BamAlignment>& discord_sr, int64_t cl2);
	void sr_cluster_1_no_cl(vector<int64_t>& src_return, map<int64_t, vector<string>>& bpread, vector<BamTools::BamAlignment>& discord_sr);

	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};

inline bool SplitReadSVCaller::covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
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
#endif /* MEERKAT_SPLITREADSVCALLER_HPP_ */
