/*
 * SplitReadReAlign.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#ifndef MEERKAT_SPLITREADREALIGN_HPP_
#define MEERKAT_SPLITREADREALIGN_HPP_

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

#include "../../third/xfaidx.h"

#include "../ClusterEntry.hpp"
#include "../BWACaller.hpp"

namespace meerkat {
using namespace std;

class SplitReadReAlign {
public:
	SplitReadReAlign();
	~SplitReadReAlign();
	void set_option_parser(const castle::OptionParser& the_options);

	void align_split_reads();
	void align_split_reads_alt();
	void align_split_reads_par();

	void collect_intra_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions, boost::unordered_set<string>& cluster_exist);
	void collect_intra_regions(map<string, map<int8_t, int32_t>>& bp_weight, map<string, map<string, string>>& regions, set<string>& cluster_exist);

	void collect_inter_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& regions, boost::unordered_set<string>& cluster_exist);
	void collect_inter_regions(map<string, map<int8_t, int32_t>>& bp_weight, map<string, map<string, map<string, string>>>& regions, set<string>& cluster_exist);

	void write_intra_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions, boost::unordered_set<string>& cluster_exist);
	void write_intra_regions(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, string>>& intra_regions, boost::unordered_set<string>& cluster_exist);

	void write_inter_regions_serial(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& inter_regions, boost::unordered_set<string>& cluster_exist);
	void write_inter_regions(boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight, boost::unordered_map<string, map<string, map<string, string>>>& inter_regions, boost::unordered_set<string>& cluster_exist);
	void prepare_split_alignment();
	void prepare_split_alignment_single();

	void remove_and_merge_temporary_alignment_files();
	void collect_misaligned_reads(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const map<string, map<int8_t, int32_t>>& bp_weight);
	void collect_misaligned_reads_serial(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const map<string, map<int8_t, int32_t>>& bp_weight);
	void collect_misaligned_reads_serial_alt(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight);
	void collect_misaligned_reads_par(boost::unordered_map<string, map<int8_t, string>>& alg_mis, const boost::unordered_map<string, map<int8_t, int32_t>>& bp_weight);

	void adjust_misaligned_reads(boost::unordered_map<string, map<int8_t, string>>& alg_mis);
	void adjust_misaligned_reads_serial(boost::unordered_map<string, map<int8_t, string>>& alg_mis);
	void adjust_misaligned_reads_serial_alt(boost::unordered_map<string, map<int8_t, string>>& alg_mis);
	void adjust_misaligned_reads_par(const boost::unordered_map<string, map<int8_t, string>>& alg_mis);
	void create_misalignment_adjusted_bam();

	pair<int64_t, int64_t> overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
private:
	static const int64_t BLOCK_SIZE = 4 * 1024 * 1024;
	int64_t n_blocks;
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;

};

//TODO: if I change change the statement, if (-1 == b1 && -1 == b2) to if (0 == b1 && 0 == b2), the algorithm finds more SV events.
inline pair<int64_t, int64_t> SplitReadReAlign::overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {

//	if (-1 == b1 && -1 == b2) {
	if (0 == b1 && 0 == b2) {
		return make_pair(a1, a2);
	}
//		( $a1, $a2 ) = sort { $a <=> $b } ( $a1, $a2 );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( $b1, $b2 ) = sort { $a <=> $b } ( $b1, $b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		int64_t a = (a1 > b1) ? a1 : b1;
		int64_t b = (a2 > b2) ? b2 : a2;
		return make_pair(a, b);
	} else {
		// the original version returns (0, 0)
//		return make_pair(b1, b2);
		return make_pair(0, 0);
	}
}

inline bool SplitReadReAlign::covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
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

#endif /* MEERKAT_SPLITREADREALIGN_HPP_ */
