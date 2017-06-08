#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "Cluster.h"
#include "ClusterContainer.h"

/* Relies on a globally declared ReadInfo container.  The ReadInfo
 * container is not really a good idea to begin with, but since it
 * is already implemented and working AND doesn't cause any severe
 * performance problems, we'll go along with it.
 *
 * XXX: ASSUMES THAT READS IN THE INPUT APPEAR IN THE SAME ORDER
 * IN BOTH CLUSTERS.  This assumption may need to be enforced by
 * sorting each cluster's read list by name in the container build
 * phase.
 */

/*
 * A simple helper function for handling cluster ties.  'expr' is
 * an evaluated boolean expression comparing simple O(1) properties
 * of cluster 'a' and 'b'.  If 'expr' is true, then 'a' and 'b'
 * _might_ be equal, but only if they contain exactly the same set
 * of reads.  This function does the slow, read by read comparison
 * of the clusters to determine if they are actually equal.  If
 * they are not equal, return a lexicographic ordering on the first
 * read that is not shared between the clusters.
 *
 * Remember: a == b is defined with < as !(a < b) && !(b < a).
 * So if a == b, then op(a,b) = op(b,a) = false.
 */
/*
 * XXX: this doesn't work if the read lists are not sorted; and
 * they are not currently sorted by code, we only assume the input
 * is correctly sorted.  this needs to be fixed.
 */
bool eq(const Cluster *a, const Cluster *b) {
	if (a->get_size() != b->get_size())
		//return a->get_size() < b->get_size();
		return false;

	/* If our quick test failed, we have to do the very */
	/* long slow test */
	const list<Readpair> &areads = a->get_reads();
	const list<Readpair> &breads = b->get_reads();

	list<Readpair>::const_iterator ait = areads.begin(), bit = breads.begin();
	for (; ait != areads.end() && bit != breads.end(); ++ait, ++bit)
		/* If the reads are used, they don't matter */
		if (ait->get_name() != bit->get_name())
			return false;

	/* All the reads are shared, the clusters are equal */
	return true;
}

bool cluscomp::operator()(const Cluster *a, const Cluster *b) {
	if (abs(a->get_weight() - b->get_weight()) > 0.0000000001)
		return a->get_weight() < b->get_weight();

	bool a_same = a->is_same_chr();
	bool b_same = b->is_same_chr();
	if (a_same != b_same) {
		/* If exactly one of the two clusters is confined to a single
		 * chromosome, then that cluster is greater */
		return b_same;
	}

	/* Now we have two clusters with the same weight */
	/* If the clusters do not share the same read set, then the
	 * ordering is arbitrary.  We choose to order by read name
	 * of the first read in the cluster. */
	if (!eq(a, b)) {
		const list<Readpair> &areads = a->get_reads();
		const list<Readpair> &breads = b->get_reads();
		return areads.front().get_name() < breads.front().get_name();
	}

	/* Now we know both a and b are clusters that contain the same
	 * reads, just mapped to different locations. */
	return false; /* so we want a<b = F and b<a = F */

	/*
	 bool a_same = a->is_same_chr();
	 bool b_same = b->is_same_chr();
	 if (a_same && b_same) {
	 // both clusters are on the same chr
	 // favor the smallest insert size
	 if (b->get_insert_size() == a->get_insert_size())
	 return op(a, b);
	 else
	 return b->get_insert_size() < a->get_insert_size();
	 } else if (a_same || b_same) {
	 // exactly one of the clusters is the same chr
	 return b_same;
	 } else {
	 if (b->get_mismatches() == a->get_mismatches())
	 return op(a, b);
	 else
	 return b->get_mismatches() < a->get_mismatches();
	 }
	 */
}

ClusterContainer::ClusterContainer() :
		max_id(-1) {
}

list<Cluster*>&
ClusterContainer::get_dirty_clusters() {
	return dirty_clusters;
}

size_t ClusterContainer::get_size() const {
	return clusters.size();
}

/* O(log(n)) insert */
void ClusterContainer::add(Cluster *cl) {
	/* Silently reject low weight or empty clusters */
	if (cl->get_weight() <= weight_threshold || cl->get_size() == 0) {
#ifdef DEBUG
		cerr << "add: rejected cluster [" << *cl
		<< "] due to weight or nreads\n";
#endif
		return;
	}

	if (cl->get_id() > max_id) {
		dcl_filter.resize(cl->get_id());
		max_id = cl->get_id();
	}

	ClusterContainer::value_type rec(cl, cl);
	/* track the multimap iterator in the cluster */
	cl->set_it(clusters.insert(rec));
}

/* O(log(n)) removal */
void ClusterContainer::remove(Cluster *cl) {
	clusters.erase(cl->get_it());
}

/* O(1) tracking list */
void ClusterContainer::set_dirty(Cluster *cl) {
	if (!dcl_filter[cl->get_id()]) {
#ifdef DEBUG
		cout << "set cluster id=" << cl->get_id() << " dirty\n";
#endif
		dcl_filter[cl->get_id()] = true;
		dirty_clusters.push_back(cl);
	}
}

void ClusterContainer::unset_dirty(Cluster *cl) {
	dcl_filter[cl->get_id()] = false;
}

void ClusterContainer::reset_dirty() {
	dirty_clusters.clear();
}

#include "ReadInfo.h"

void ClusterContainer::fix(list<ReadInfo*> R) {
	/* Remove all affected clusters from the container */
	for (list<Cluster*>::iterator it = dirty_clusters.begin();
			it != dirty_clusters.end(); ++it) {
		remove(*it);
		unset_dirty(*it);
	}

	/* Remove the selected reads from all clusters containing them */
	for (list<ReadInfo*>::iterator it = R.begin(); it != R.end(); ++it) {
		ReadInfo *ri = *it;
		for (ReadInfo::iterator it2 = ri->begin(); it2 != ri->end(); ++it2) {
			Cluster *c = it2->first;
			list<Readpair>::iterator it3 = it2->second;
			Readpair &r = *it3;
			c->remove(r, it3);
		}
	}

	/* Reinsert the affected clusters */
	for (list<Cluster*>::iterator it = dirty_clusters.begin();
			it != dirty_clusters.end(); ++it)
		add(*it);

	reset_dirty();
	/* XXX: remove get_dirty and reset_dirty since no one needs these */
}

bool sortfn(const Cluster *a, const Cluster *b) {
	bool a_same = a->is_same_chr();
	bool b_same = b->is_same_chr();
	if (a_same && b_same)
		/* both clusters are on the same chr */
		/* favor the smallest insert size */
		return a->get_insert_size() < b->get_insert_size();
	else if (a_same || b_same)
		/* exactly one of the clusters is the same chr */
		/* XXX: this should never happen */
		return b_same;
	else
		return a->get_mismatches() < b->get_mismatches();
}

/*
 * All of the tie-breaking logic is handled in the comparison
 * class.  Here we just return the "greatest" object in the
 * container.
 */
list<Cluster *> ClusterContainer::best() {
	Cluster *best = clusters.rbegin()->second;
	std::pair<ClusterContainer::iterator, ClusterContainer::iterator> range =
			clusters.equal_range(best);

	list<Cluster*> ret;
	for (ClusterContainer::iterator it = range.first; it != range.second; ++it)
		ret.push_back(it->second);

	ret.sort(sortfn);

	return ret;

}
