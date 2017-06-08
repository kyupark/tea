#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "Readpair.h"

double Readpair::get_weight() const { return weight; }

const string& Readpair::get_name() const { return name; }

const Read& Readpair::get_read() const { return read; }

const Read& Readpair::get_mate() const { return mate; }

size_t Readpair::get_mismatches() const { return mismatch_read + mismatch_mate; }

int64_t Readpair::get_insert_size() const {
	return insert_size;
}


ostream&
Readpair::operator<< (ostream &ostr)
{
	return ostr << name << "\t" << read_group << "\t"
	            << read << "\t" << mate << "\t" << insert_size;
}

/* Global insertion operator */
ostream &
operator<<(ostream &ostr, Readpair &rp)
{
	return rp.operator<<(ostr);
}

/*
 * Format:
 *   cid       - cluster ID
 *   cweight   - cluster weight
 *   cmismatch - mismatches for the whole cluster
 *   name      - name of this read (the mate shares the name)
 *   rg        - read group
 *   r         - read (see below)
 *   m         - mate (see below)
 *   insert    - insert size between r and m
 *   mm1       - mismatches for read 'r'
 *   mm2       - mismatches for read 'm'
 *   weight    - weight of this read-cluster assignment
 *
 *   | cid | cweight | cmismatch | name | rg | r | m | insert | mm1 | mm2 | weight |
 *
 * where 'r' and 'm' denote the read (r) and mate (m).  These fields
 * have format:
 * 
 *   | chrom | strand | start | end | length |
 *
 * NOTE: the cid, cweight and cmismatch fields must be stripped off
 * before calling this constructor.
 */
Readpair::Readpair(istream &data_stream) {
	data_stream >> name >> read_group;
	read = Read(data_stream);
	mate = Read(data_stream);
	data_stream >> insert_size >> mismatch_read >> mismatch_mate >> weight;
}
