/*
 * CNVOverlapper.hpp
 *
 *  Created on: Oct 14, 2016
 *      Author: el174
 */

#ifndef MEERKAT_CNS_CNVOVERLAPPER_HPP_
#define MEERKAT_CNS_CNVOVERLAPPER_HPP_

#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "../../castle/OptionParser.hpp"
#include "../../castle/TimeChecker.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

#include "../../emma/IntervalTree.hpp"


namespace meerkat {

using namespace std;

struct ColoredString {
	int64_t group_id;
	string type;
	string value;
	int64_t the_size;
	bool operator<(const ColoredString& other) const {
		if (group_id < other.group_id) {
			return true;
		} else if (group_id > other.group_id) {
			return false;
		}
		if (type < other.type) {
			return true;
		} else if (type > other.type) {
			return false;
		}
		return false;
	}
	bool operator==(const ColoredString& other) const {
		if (group_id != other.group_id) {
			return false;
		}
		if (type != other.type) {
			return false;
		}
		return value == other.value && the_size == other.the_size;
	}
};

typedef emma::IntervalType<ColoredString> ColoredStringIntervalEntry;
typedef std::vector<ColoredStringIntervalEntry> ColoredStringIntervalEntryVector;
typedef emma::IntervalTreeType<ColoredString> ColoredStringIntervalClusterTree;

typedef emma::IntervalType<string> StringIntervalEntry;
typedef std::vector<StringIntervalEntry> StringIntervalEntryVector;
typedef emma::IntervalTreeType<string> StringIntervalClusterTree;

class CNVOverlapper {
public:
	CNVOverlapper();
	~CNVOverlapper();
	void set_option_parser(const castle::OptionParser& the_options);
	void find_overlaps();
	void find_sv_overlaps();
private:
	int64_t n_cores;
	castle::OptionParser options;
};

} /* namespace meerkat */

#endif /* MEERKAT_CNS_CNVOVERLAPPER_HPP_ */
