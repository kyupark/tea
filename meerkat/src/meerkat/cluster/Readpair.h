#ifndef READPAIR_H
#define READPAIR_H

#include <ostream>
#include <list>

#include "Read.h"

class Cluster;

using namespace std;
/* Individual reads can never be separated from their mates */
class Readpair {
public:
	Readpair();
	Readpair(istream &data_stream);  /* Construct from a line of input */

	ostream& operator<< (ostream &ostr);

	double get_weight() const;
	const string &get_name() const;
	size_t get_mismatches() const;
	int64_t get_insert_size() const;
	const Read& get_read() const;
	const Read& get_mate() const;

private:
	double weight;
	string name;
	string read_group;
	Read read;
	Read mate;
	int64_t insert_size;
	int64_t mismatch_read;
	int64_t mismatch_mate;
};

ostream& operator<< (ostream &ostr, Readpair &rp);

#endif
