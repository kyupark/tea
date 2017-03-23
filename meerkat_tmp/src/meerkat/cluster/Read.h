#ifndef READ_H
#define READ_H

#include <istream>
using namespace std;

class Read {
public:
	Read();
	Read(istream &data);

	const string& get_chr() const;

	ostream& operator<< (ostream &ostr);

private:
	string chr;
	int strand;
	long start;
	//long end;
	size_t length;
};

ostream& operator<< (ostream &ostr, Read &r);

#endif
