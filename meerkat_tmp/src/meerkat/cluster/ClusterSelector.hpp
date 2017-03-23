/*
 * ClusterSelector.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#ifndef MEERKAT_CLUSTER_CLUSTERSELECTOR_HPP_
#define MEERKAT_CLUSTER_CLUSTERSELECTOR_HPP_

#include <unordered_map>
#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

#include "ClusterContainer.h"
#include "Readpair.h"
#include "Cluster.h"
#include "ReadInfo.h"

namespace meerkat {
using namespace std;
class ClusterSelector {
public:
	ClusterSelector();
	~ClusterSelector();
	void set_option_parser(const castle::OptionParser& the_options);
	void select_cluster();
	list<ReadInfo*> select2(ClusterContainer &con, unordered_map<string, ReadInfo> &readinfo, Cluster &cl, ostream &ostr, int64_t id, int64_t sid);
private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};

} /* namespace meerkat */

#endif /* MEERKAT_CLUSTER_CLUSTERSELECTOR_HPP_ */
