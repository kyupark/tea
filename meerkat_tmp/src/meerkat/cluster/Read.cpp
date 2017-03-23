#include <cmath>
#include <iostream>

using namespace std;

#include "Read.h"

Read::Read() :
	chr("NA"),
	strand(0),
	start(-1),
	//end(-1),
	length(-1)
{
}

const string& Read::get_chr() const { return chr; }

Read::Read(istream &data)
{
	data >> chr >> strand >> start >> /* end >> */ length;
}

ostream&
Read::operator<< (ostream &ostr)
{
	return ostr << chr << "\t" << strand << "\t"
	     << start << "\t" // << end << "\t"
	     << length;
}

/* global alias */
ostream &
operator<<(ostream &ostr, Read &r)
{
	return r.operator<<(ostr);
}
