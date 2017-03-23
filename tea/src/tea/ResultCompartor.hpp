/*
 * ResultCompartor.hpp
 *
 *  Created on: Feb 3, 2017
 *      Author: el174
 */

#ifndef TEA_RESULTCOMPARTOR_HPP_
#define TEA_RESULTCOMPARTOR_HPP_
#include <boost/unordered_map.hpp>
#include "TEAOptionParser.hpp"
#include "../castle/TimeChecker.hpp"
#include "../castle/StringUtils.hpp"
#include "../castle/IOUtils.hpp"
#include "../emma/IntervalTree.hpp"

namespace tea {

struct GermlineString {
	int64_t pbp;
	int64_t nbp;
	string value;

	bool operator<(const GermlineString& other) const {
		if(pbp < other.pbp) {
			return true;
		} else if(pbp > other.pbp) {
			return false;
		}
		if(nbp < other.nbp) {
			return true;
		} else if(nbp > other.nbp) {
			return false;
		}
		return value < other.value;
	}
	bool operator==(const GermlineString& other) const {
		return pbp == other.pbp && nbp == other.nbp && value == other.value;
	}
};

typedef emma::IntervalType<GermlineString> GermlineStringIntervalEntry;
typedef std::vector<GermlineStringIntervalEntry> GermlineStringIntervalEntryVector;
typedef emma::IntervalTreeType<GermlineString> GermlineStringIntervalClusterTree;

class ResultCompartor {
public:
	ResultCompartor();
	~ResultCompartor();
	void set_option_parser(const TEAOptionParser& the_options);
	void find_overlaps(const string& comp_suffix);
private:
	int32_t n_cores;
	TEAOptionParser options;
};

} /* namespace tea */

#endif /* TEA_RESULTCOMPARTOR_HPP_ */
