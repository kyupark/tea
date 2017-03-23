#ifndef READGROUP_HPP_
#define READGROUP_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <boost/lexical_cast.hpp>
#include "api/BamAlignment.h"
#include "utils/bamtools_utilities.h"
//#include "../gzstream.h"

namespace meerkat {
	using namespace std;
	struct QENCODING {
			const char *name;
			int min;
			int max;
			int offset;
	};
	class ReadGroup {
		public:
			static const QENCODING known_qencodings[4];
		public:
			ReadGroup();
			ReadGroup(BamTools::BamAlignment &al, int max_isize,
					int isize_samples, string prefix, string str_block_id,
					list<string> blacklist);
			~ReadGroup();
			ReadGroup& operator+=(const ReadGroup& b) {
				// check for self-assignment
				if (this == &b)
					return *this;
				if ((name != b.name) || (prefix != b.prefix)) {
					return *this;
				}
				for (auto b_itr = b.readlens.begin(); b.readlens.end() != b_itr;
						++b_itr) {
					if (find(readlens.begin(), readlens.end(), *b_itr)
							== readlens.end()) {
						readlens.push_back(*b_itr);
					}
				}
				copy(b.readlens.begin(), b.readlens.end(),
						back_inserter(readlens));

				if (static_cast<int>(inserts.size()) < isize_samples) {
					inserts.reserve(inserts.size() + b.inserts.size());
					copy(b.inserts.begin(), b.inserts.end(),
							back_inserter(inserts));
				}

				nreads += b.nreads;
				return *this;
			}
			static BamTools::BamAlignment snip(BamTools::BamAlignment &al,
					int start, int len);
			void witness(BamTools::BamAlignment &al);

			/* Also trims the read based on q and min_read_len.  If the read
			 * becomes too short, it is rejected and record() returns false.
			 */
			bool recordSC(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int clippedMate, ofstream &clipfile, ofstream& split1,
					ofstream& split2, int big_s_bps, int frag_size,
					int n_cutoff);

			bool recordSCAlt(ostream& f1, ostream& f2,
					BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int clippedMate, ofstream &clipfile, ofstream& split1,
					ofstream& split2, int big_s_bps, int frag_size,
					int n_cutoff);
			bool recordSCAltNoRG(BamTools::BamAlignment &a,
					BamTools::BamAlignment &b, int clippedMate,
					ostream &clipfile, ostream& split1, ostream& split2,
					int big_s_bps, int frag_size, int n_cutoff);
			bool recordSCAltRG(ostream& f1, ostream& f2, BamTools::BamAlignment &a,
					BamTools::BamAlignment &b, int clippedMate,
					int big_s_bps, int frag_size, int n_cutoff);
			bool recordSCAltOnlyClip_Split(BamTools::BamAlignment &a,
					BamTools::BamAlignment &b, int clippedMate,
					ofstream &clipfile, ofstream& split1, ofstream& split2,
					int big_s_bps, int frag_size, int n_cutoff);

			bool recordUM(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int big_s_bps, int n_cutoff);
			bool recordUMAlt(ostream& f1, ostream& f2,
					BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int big_s_bps, int n_cutoff);
			bool recordUMAltNoRG(BamTools::BamAlignment &a,
					BamTools::BamAlignment &b, int big_s_bps, int n_cutoff);
			bool recordUMAltRG(ostream& f1, ostream& f2, BamTools::BamAlignment &a,
								BamTools::BamAlignment &b, int big_s_bps, int n_cutoff);

			bool recordUU(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int big_s_bps, int n_cutoff);
			bool recordUUAlt(ofstream& f1, ofstream& f2,
					BamTools::BamAlignment &a, BamTools::BamAlignment &b,
					int big_s_bps, int n_cutoff);

			string getName();
			int getNReads();
			vector<int> getReadlens();
			bool is_blacklisted();

			friend ostream & operator<<(ostream &os, const ReadGroup &rg);
		public:
			static void writeFQpair(ostream &f1, BamTools::BamAlignment &a,
					ostream &f2, BamTools::BamAlignment &b, int n_cutoff);
			static void writeFQ(ostream &f, BamTools::BamAlignment &al);
			static void trim_read(BamTools::BamAlignment &al, int q, int qenc);
			static void split_read(BamTools::BamAlignment &al,
					ostream & split1, ostream& split2, int frag_size,
					int n_cutoff);

			static bool isBigS(const BamTools::BamAlignment& al,
					const BamTools::CigarOp& cop, int n);
			static bool isBigH(const BamTools::BamAlignment& al);
			static bool isMatePair(const BamTools::BamAlignment& a,
					const BamTools::BamAlignment& b);
			static int getMateNumber(const BamTools::BamAlignment &al);
		public:
			string name;
			vector<int> readlens;
		private:
			ofstream f1;
			ofstream f2;
			string prefix;
			string block_id;
			vector<int> inserts;

			int nreads;
			int max_isize;
			int isize_samples;
			bool extract_scs_and_mates;
			bool blacklisted;
	};
	inline bool ReadGroup::isBigS(const BamTools::BamAlignment& al,
			const BamTools::CigarOp& cop, int n) {
		/* cop *can* be a NULL ref if the alignment has a "*" CIGAR string */
		return (al.CigarData.size() > 0)
				&& (cop.Type == 'S' && min(al.Length, (signed) cop.Length) >= n);
	}

	/* Not sure if H operations must be first or last, so just check everything. */
	inline bool ReadGroup::isBigH(const BamTools::BamAlignment& al) {
		for (auto iter = al.CigarData.begin(); iter != al.CigarData.end();
				++iter)
			if (iter->Type == 'H')
				return true;

		return false;
	}

	inline bool ReadGroup::isMatePair(const BamTools::BamAlignment& a,
			const BamTools::BamAlignment& b) {
		return (a.IsFirstMate() && b.IsSecondMate())
				|| (a.IsSecondMate() && b.IsFirstMate());
	}

	inline int ReadGroup::getMateNumber(const BamTools::BamAlignment &al) {
		return al.IsFirstMate() ? 1 : 2;
	}
	inline void ReadGroup::writeFQpair(ostream &f1, BamTools::BamAlignment &a,
			ostream &f2, BamTools::BamAlignment &b, int n_cutoff) {
//		bool debug = string::npos != a.Name.find("ST-E00104:502:HFJN5CCXX:1:1108:32826:10890");
		int n = 0;
		for (size_t i = 0; i < a.QueryBases.size(); ++i) {
			if (a.QueryBases[i] == 'N') {
				++n;
			}
		}
//		if(debug) {
//			cout << "[ReadGroup.writeFQpair] n: " << n << "\n";
//			cout << "[ReadGroup.writeFQpair] QueryBases: " << a.QueryBases << "\n";
//		}
		if (n > n_cutoff) {
			return;
		}
		for (size_t i = 0; i < b.QueryBases.size(); ++i) {
			if (b.QueryBases[i] == 'N') {
				++n;
			}
		}
//		if(debug) {
//			cout << "[ReadGroup.writeFQpair] n: " << n << "\n";
//			cout << "[ReadGroup.writeFQpair] QueryBases: " << b.QueryBases << "\n";
//		}
		if (n > n_cutoff) {
			return;
		}

		writeFQ(f1, a);
		writeFQ(f2, b);
	}
	/*
	 * Write out the substring [start, start+length] (inclusive) of the
	 * aligned read 'al' to the file 'f' in FASTQ format.
	 *
	 * writeFQ() quietly aborts if al contains more than n_cutoff N bases.
	 *
	 * NOTE: Quite similar to the PrintFastq function in the bamtools
	 * utils.
	 */

	inline void ReadGroup::writeFQ(ostream &f, BamTools::BamAlignment &al) {
		string quals = al.Qualities;
		string bases = al.QueryBases;
		if (al.IsReverseStrand()) {
			BamTools::Utilities::Reverse(quals);
			BamTools::Utilities::ReverseComplement(bases);
		}

		/* Using \n instead of endl to avoid unnecessary stream flushes */
		f << "@" << al.Name << "\n" << bases << "\n" << "+" << "\n" << quals
				<< "\n";
	}
} // namespace std

#endif /* READGROUP_HPP_ */
