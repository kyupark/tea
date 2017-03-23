/*
 * OutputFilter.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#ifndef MEERKAT_OUTPUTFILTER_HPP_
#define MEERKAT_OUTPUTFILTER_HPP_

#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <string>

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

#include "../ClusterEntry.hpp"

namespace meerkat {
using namespace std;
class OutputFilter {
public:
	OutputFilter();
	~OutputFilter();
	void set_option_parser(const castle::OptionParser& the_options);
	void filter();
	void filter_intra_chromosomal_events();
	void filter_inter_chromosomal_events();
private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};

} /* namespace meerkat */

#endif /* MEERKAT_OUTPUTFILTER_HPP_ */
