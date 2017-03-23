#ifndef READINFO_H
#define READINFO_H 1

#include <list>
#include "Cluster.h"

/*
 * Each read is contained in >= 1 cluster.  Each such cluster
 * maintains its own list of reads.  So each read<->cluster
 * relation needs to be tracked by two things:
 *   1. the cluster containing the read
 *   2. the read's place in the cluster's readlist
 */
/* XXX: could probably just typedef rather than have a new class */
class ReadInfo : public list< std::pair<Cluster*,list<Readpair>::iterator> > {
public:
	ReadInfo();
	bool is_used() const;
	void set_used();

	/* XXX: deleteme */
	//typedef std::pair< list<Cluster*>,list<Readpair>::iterator >::iterator
		//iterator;

private:
	bool used; /* XXX: deleteme and related functions */
};

#endif
