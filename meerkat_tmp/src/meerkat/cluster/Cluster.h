#ifndef CLUSTER_H
#define CLUSTER_H

#include <ostream>
#include <list>

#include "ClusterContainer.h"

class ClusterContainer;
class Readpair;
using namespace std;
class Cluster {
public:
	Cluster(long cid);

	long get_id() const;
	size_t get_size() const;
	double get_weight() const;
	void recalc_weight();
	double get_insert_size() const;
	size_t get_mismatches() const;
	const list<Readpair>& get_reads() const;
	list<Readpair>& get_reads();

	bool is_same_chr() const;

	list<Readpair>::iterator add(Readpair &ref);
	void remove(Readpair &ref, list<Readpair>::iterator pos);
	ostream& write(ostream &ostr);

	void set_it(ClusterContainer::iterator it);
	ClusterContainer::iterator get_it();


private:
	long id;
	size_t nreads;
	double weight;
	size_t issum;
	size_t mmsum;
	list<Readpair> reads;
	ClusterContainer::iterator it;
};

ostream& operator<< (ostream &ostr, Cluster &cl);

#endif
