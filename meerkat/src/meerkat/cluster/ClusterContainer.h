#ifndef CLUSTERCONTAINER_H
#define CLUSTERCONTAINER_H

#include <vector>
#include <list>
#include <map>

#include "Readpair.h"

class Cluster;
class ReadInfo;

class cluscomp {
public:
	bool operator() (const Cluster *a, const Cluster *b);
};

using namespace std;

class ClusterContainer {
public:
	ClusterContainer();
	void add(Cluster *cl);
	void remove(Cluster *cl);
	void set_dirty(Cluster *cl);
	void unset_dirty(Cluster *cl);
	void reset_dirty();
	size_t get_size() const;

	list<Cluster*>& get_dirty_clusters();

	/* Remove and reinsert all dirty clusters.  Must be run after
	 * every select(), because select() may change the weight of
	 * some clusters
	 */
	void fix(list<ReadInfo*> R);

	list<Cluster*> best();    /* Find the "best" cluster in the container */

	typedef multimap<Cluster*,Cluster*,cluscomp>::iterator iterator;
	typedef multimap<Cluster*,Cluster*,cluscomp>::value_type value_type;

private:
	/* fix() will not reinsert clusters with weight <= weight_threshold */
	static constexpr const double weight_threshold = 0;
	static constexpr const double accuracy_threshold = 0.0000000001;

	long max_id;
	multimap<Cluster*,Cluster*,cluscomp> clusters;

	/* The following data structures are to speed up the fix() step */
	vector<bool> dcl_filter;  /* bit filter for O(1) insert into dirty_clusters */
	list<Cluster*> dirty_clusters;
};

#endif
