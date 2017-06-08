/*
 * SelectionFilter.hpp
 *
 *  Created on: Jan 12, 2017
 *      Author: el174
 */

#ifndef MEERKAT_FILTER_SELECTIONFILTER_HPP_
#define MEERKAT_FILTER_SELECTIONFILTER_HPP_

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

namespace meerkat {
using namespace std;

class SelectionFilter {
public:
	SelectionFilter();
	~SelectionFilter();
	void set_option_parser(const castle::OptionParser& the_options);
	void filter_variants();
private:
	int64_t n_cores;
	castle::OptionParser options;
};

} /* namespace meerkat */

#endif /* MEERKAT_FILTER_SELECTIONFILTER_HPP_ */
