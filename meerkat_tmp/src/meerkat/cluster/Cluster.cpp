#include <iostream>
#include <cstdlib>
#include <ostream>
#include <list>

using namespace std;

#include "Cluster.h"
#include "Readpair.h"

Cluster::Cluster(long cid) :
	id(cid),
	nreads(0),
	weight(0),
	issum(0),
	mmsum(0)
{ }

void Cluster::set_it(ClusterContainer::iterator it) { this->it = it; }
ClusterContainer::iterator Cluster::get_it() { return it; }

long Cluster::get_id() const { return id; }

/* NOT reads.size(), because reads are never removed from that container */
/* XXX: not true anymore */
size_t Cluster::get_size() const { return nreads; }

double Cluster::get_weight() const { return weight; }

list<Readpair>& Cluster::get_reads() { return reads; }

const list<Readpair>& Cluster::get_reads() const { return reads; }

/* The average insert size for reads in this cluster */
double Cluster::get_insert_size() const {
	return (double)issum / (double)nreads;
}

size_t Cluster::get_mismatches() const { return mmsum; }

/* Return true if all of the reads in this cluster (and
 * all of their mates) belong to the same chromosome.
 */
bool
Cluster::is_same_chr() const
{
	list<Readpair>::const_iterator it = reads.begin();
	const string &chr = it->get_read().get_chr();
	for ( ; it != reads.end(); ++it) {
		const string &r = it->get_read().get_chr();
		const string &m = it->get_mate().get_chr();
		if (r != chr || m != chr)
			return false;
	}

	return true;
}

list<Readpair>::iterator
Cluster::add(Readpair &ref)
{
	weight += ref.get_weight();
	++nreads;
	int64_t tmp_insert_size = ref.get_insert_size();
	if(-1 != tmp_insert_size) {
		issum += ref.get_insert_size();
	} else {
		// no insert size penalty
		issum += 500000;
	}
	mmsum += ref.get_mismatches();
	list<Readpair>::iterator it;
	for (it = reads.begin(); it != reads.end(); ++it)
		if (it->get_name() > ref.get_name())
			break;
	reads.insert(it, ref);
	--it;
	return it;
}

void Cluster::remove(Readpair &ref, list<Readpair>::iterator pos)
{
	weight += ref.get_weight();
	issum -= ref.get_insert_size();
	mmsum -= ref.get_mismatches();
	--nreads;  /* XXX: nreads has become pointless since weight adjusting is instant */
	if (nreads < 0) {
		cerr << "FATAL: nreads in cluster ["
		     << *this << "] became negative" << endl;
		exit(1);
	}

	/* Clusters own Readpair memory, so this deletes the Readpair */
	reads.erase(pos);
}

/* Does not write the cluster's reads to ostr; just summary info
 * about the cluster */
ostream&
Cluster::write(ostream &ostr)
{
	return ostr << "cluster id=" << id << " weight=" << weight
	     << " #reads=" << reads.size() << " (nreads=" << nreads << ")";
}

ostream&
operator<< (ostream &ostr, Cluster &cl)
{
	return cl.write(ostr);
}
