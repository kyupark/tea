#include <iostream>
#include <algorithm>
#include <vector>

#include "ReadGroup.h"

using namespace std;
using namespace BamTools;

namespace meerkat {
	const QENCODING ReadGroup::known_qencodings[] = {
	/* ASCII-33 */
	{ "Sanger", 33, 126 },

	/* All of the Illumina formats have a 64 offset */
	{ "Solexa/Illumina 1.0", 59, 126 },

	/* Illumina 1.3 no longer allows negative values */
	{ "Illumina 1.3", 64, 126 },

	/* The only difference between 1.3 and 1.5 is that 1.5 does
	 * not utilize values 0 or 1, and value 2 is used as a read
	 * quality indicator. */
	{ "Illumina 1.5", 66, 126 }, };

	ReadGroup::ReadGroup() :
			nreads(0), max_isize(0), isize_samples(0), extract_scs_and_mates(
					false), blacklisted(false) {
	}

	ReadGroup::~ReadGroup() {
		if (!blacklisted) {
			f1.close();
			f2.close();
		}

//		/* Write the insert size distribution to file */
//		ofstream f((prefix + "/" + name + ".is").c_str());
//		for (unsigned int i = 0; i < inserts.size(); ++i)
//			f << inserts[i] << "\n";
//		f.close();
	}

	ReadGroup::ReadGroup(BamAlignment &al, int max_isize, int isize_samples,
			string prefix, string str_block_id, list<string> blacklist) :
			prefix(prefix), block_id(str_block_id), max_isize(max_isize), isize_samples(
					isize_samples), blacklisted(false) {
		if (!al.GetReadGroup(name)) {
			name = "none";
		}

		nreads = 0;

		/* Determine if this read group is in the blacklist */
		for (list<string>::iterator it = blacklist.begin();
				it != blacklist.end(); ++it) {
			if (*it == name) {
				blacklisted = true;
				break;
			}
		}

		if (!blacklisted) {
			f1.open((prefix + "/" + name + "_1." + str_block_id + ".fq.gz"),
					ios::binary);
			f2.open((prefix + "/" + name + "_2." + str_block_id + ".fq.gz"),
					ios::binary);
		}

		witness(al);
	}

	void ReadGroup::witness(BamAlignment &al) {
		++nreads;

		/* Selects first 'isize_samples' samples instead of random sampling */
		/* Don't use negative insert sizes.  If it's negative, then the mate  *
		 * has a positive insert of the same size.  If both are recorded, the *
		 * average will always be 0.                                          */
		if (al.InsertSize > 0 && al.InsertSize <= max_isize
				&& static_cast<int>(inserts.size()) < isize_samples)
			inserts.push_back(al.InsertSize);

		if (find(readlens.begin(), readlens.end(), al.Length)
				== readlens.end()) {
			// Removed by request, Mar 17 2015
			//if (readlens.size() > 0) {
			//cerr << "ERROR: differing read lengths in read group '"
			//<< name << "'" << endl;
			//}
			readlens.push_back(al.Length);
		}
	}

	/* Assumes at least one of 'a' and 'b' are 'S' ops */
	bool copcomp(CigarOp &a, CigarOp &b) {
		if (b.Type != 'S')
			return false;
		if (a.Type != 'S')
			return true;

		return a.Length < b.Length;
	}

	void clipAlignment(BamAlignment &al) {
		int offset, length;
		CigarOp cop1 = al.CigarData.front();
		CigarOp cop2 = al.CigarData.back();

		if (copcomp(cop2, cop1)) {
			offset = 0;
			length = min(al.Length, (signed) cop1.Length);
		} else {
			offset = al.Length - min(al.Length, (signed) cop2.Length);
			length = min(al.Length, (signed) cop2.Length);
		}

		try {
			al.Qualities = al.Qualities.substr(offset, length);
			al.QueryBases = al.QueryBases.substr(offset, length);
		} catch (exception &e) {
			cout << "ERROR: substr failed in clipAlignment()" << endl;
			cout << al.Name << " " << (al.IsReverseStrand() ? "(-)" : "(+)");
			cout << " offset: " << offset << " length: " << length
					<< " taglen: " << al.Length << endl;
			cout << "cop1: " << cop1.Length << cop1.Type << endl;
			cout << "cop2: " << cop2.Length << cop2.Type << endl;
			exit(1);
		}
	}

	/*
	 * snip() doesn't leave a valid BamAlignment; it contains
	 * correct FASTQ data.  Handles negative strand alignments:
	 * 'start=0' will always correspond to the 5'-most basepair
	 * in the read.
	 */
	BamAlignment ReadGroup::snip(BamAlignment &a, int start, int len) {
		BamAlignment copy(a);

		/* Handle reverse strand mappings */
		int converted_start =
				copy.IsReverseStrand() ? copy.Length - start - len : start;
		if(converted_start < 0) {
			return copy;
		}
		copy.Length = len;
		try {
			copy.QueryBases = copy.QueryBases.substr(converted_start, len);
			copy.Qualities = copy.Qualities.substr(converted_start, len);
		} catch (exception &e) {
			cout << "ERROR: substr failed in snip(" << a.Name << ", " << start
					<< ", " << len << ")" << endl;
			cout << (a.IsReverseStrand() ? "(-)" : "(+)")
					<< ", converted_start: " << converted_start << endl;
			cout << a.QueryBases << endl;
			cout << a.Qualities << endl;
			exit(1);
		}
		return copy;
	}

	/*
	 * The alignments a and b are already trimmed by quality score.
	 * Create read pairs out of the front and back 'big_s_bps' bps
	 * of each read, provided the reads are big enough to supply
	 * 2*big_s_bps nonoverlapping bases.  This can produce up to 4
	 * pairs: FF, FB, BF, BB (F - front fragment, B - back fragment).
	 * Lixing wants these named uu1-uu4.
	 */
	bool ReadGroup::recordUU(BamAlignment &a, BamAlignment &b, int big_s_bps,
			int n_cutoff) {
		/* Not even an FF pair can be produced */
		if (a.Length < big_s_bps || b.Length < big_s_bps)
			return false;

		/* Produce the FF pair */
		BamAlignment af = snip(a, 0, big_s_bps);
		af.Name += "uu1";
		BamAlignment bf = snip(b, 0, big_s_bps);
		bf.Name += "uu1";
		writeFQpair(f1, af, f2, bf, n_cutoff);

		/* If the read is long enough, try back fragments */
		BamAlignment ab;
		BamAlignment bb;
		bool aback = false, bback = false;
		if (a.Length >= 2 * big_s_bps) {
			/* BF pair */
			aback = true;
			ab = snip(a, a.Length - big_s_bps, big_s_bps);
			ab.Name += "uu2";
			bf.Name = b.Name + "uu2";
			writeFQpair(f1, ab, f2, bf, n_cutoff);
		}
		if (b.Length >= 2 * big_s_bps) {
			/* FB pair */
			bback = true;
			bb = snip(b, b.Length - big_s_bps, big_s_bps);
			bb.Name += "uu3";
			af.Name = a.Name + "uu3";
			writeFQpair(f1, af, f2, bb, n_cutoff);
		}
		if (aback && bback) {
			/* BB pair */
			ab.Name = a.Name + "uu4";
			bb.Name = b.Name + "uu4";
			writeFQpair(f1, ab, f2, bb, n_cutoff);
		}

		return true;
	}
	bool ReadGroup::recordUUAlt(ofstream& f1, ofstream& f2,
			BamTools::BamAlignment &a, BamTools::BamAlignment &b, int big_s_bps,
			int n_cutoff) {
		/* Not even an FF pair can be produced */
		if (a.Length < big_s_bps || b.Length < big_s_bps)
			return false;

		/* Produce the FF pair */
		BamAlignment af = snip(a, 0, big_s_bps);
		af.Name += "uu1";
		BamAlignment bf = snip(b, 0, big_s_bps);
		bf.Name += "uu1";
		writeFQpair(f1, af, f2, bf, n_cutoff);

		/* If the read is long enough, try back fragments */
		BamAlignment ab;
		BamAlignment bb;
		bool aback = false, bback = false;
		if (a.Length >= 2 * big_s_bps) {
			/* BF pair */
			aback = true;
			ab = snip(a, a.Length - big_s_bps, big_s_bps);
			ab.Name += "uu2";
			bf.Name = b.Name + "uu2";
			writeFQpair(f1, ab, f2, bf, n_cutoff);
		}
		if (b.Length >= 2 * big_s_bps) {
			/* FB pair */
			bback = true;
			bb = snip(b, b.Length - big_s_bps, big_s_bps);
			bb.Name += "uu3";
			af.Name = a.Name + "uu3";
			writeFQpair(f1, af, f2, bb, n_cutoff);
		}
		if (aback && bback) {
			/* BB pair */
			ab.Name = a.Name + "uu4";
			bb.Name = b.Name + "uu4";
			writeFQpair(f1, ab, f2, bb, n_cutoff);
		}

		return true;
	}

	bool ReadGroup::recordUM(BamAlignment &a, BamAlignment &b, int big_s_bps,
			int n_cutoff) {
		BamAlignment mapped = a.IsMapped() ? a : b;
		BamAlignment unmapped = a.IsMapped() ? b : a;

		if (unmapped.Length < big_s_bps || mapped.Length < 2 * big_s_bps)
			return false;

		mapped.Name += "mu";
		unmapped.Name += "mu";

		/* unmapped >= big_s_bps.  Make a fragment from the 5' end */
		BamAlignment cp1 = snip(unmapped, 0, big_s_bps);
		cp1.Name += "1";
		mapped.Name += "1";
		writeFQpair(f1, mapped, f2, cp1, n_cutoff);

		if (unmapped.Length < 2 * big_s_bps)
			return true;

		/* unmapped >= 2*big_s_bps.  Make another fragment from the 3' end */
		BamAlignment cp2 = snip(unmapped, unmapped.Length - big_s_bps,
				big_s_bps);
		cp2.Name += "2";
		mapped.Name[mapped.Name.size() - 1] = '2';
		writeFQpair(f1, mapped, f2, cp2, n_cutoff);

		return true;
	}
	bool ReadGroup::recordUMAlt(ostream& f1, ostream& f2,
			BamTools::BamAlignment &a, BamTools::BamAlignment &b, int big_s_bps,
			int n_cutoff) {
		BamAlignment mapped = a.IsMapped() ? a : b;
		BamAlignment unmapped = a.IsMapped() ? b : a;

		if (unmapped.Length < big_s_bps || mapped.Length < 2 * big_s_bps)
			return false;

		mapped.Name += "mu";
		unmapped.Name += "mu";

		/* unmapped >= big_s_bps.  Make a fragment from the 5' end */
		BamAlignment cp1 = snip(unmapped, 0, big_s_bps);
		cp1.Name += "1";
		mapped.Name += "1";
		writeFQpair(f1, mapped, f2, cp1, n_cutoff);

		if (unmapped.Length < 2 * big_s_bps)
			return true;

		/* unmapped >= 2*big_s_bps.  Make another fragment from the 3' end */
		BamAlignment cp2 = snip(unmapped, unmapped.Length - big_s_bps,
				big_s_bps);
		cp2.Name += "2";
		mapped.Name[mapped.Name.size() - 1] = '2';
		writeFQpair(f1, mapped, f2, cp2, n_cutoff);

		return true;
	}
	bool ReadGroup::recordUMAltNoRG(BamTools::BamAlignment &a,
			BamTools::BamAlignment &b, int big_s_bps, int n_cutoff) {
		BamAlignment mapped = a.IsMapped() ? a : b;
		BamAlignment unmapped = a.IsMapped() ? b : a;

		if (unmapped.Length < big_s_bps || mapped.Length < 2 * big_s_bps) {
			return false;
		}

		/* unmapped >= big_s_bps.  Make a fragment from the 5' end */
		snip(unmapped, 0, big_s_bps);

		if (unmapped.Length < 2 * big_s_bps) {
			return true;
		}

		/* unmapped >= 2*big_s_bps.  Make another fragment from the 3' end */
		snip(unmapped, unmapped.Length - big_s_bps, big_s_bps);

		return true;
//		if (a.IsMapped()) {
//			// mapped=> a
//			// unmapped => b
//			if (b.Length < big_s_bps || a.Length < 2 * big_s_bps)
//				return false;
//			snip(b, 0, big_s_bps);
//			if (b.Length < 2 * big_s_bps)
//				return true;
//		} else {
//			// mapped => b
//			// unmapped => a
//			if (a.Length < big_s_bps || b.Length < 2 * big_s_bps)
//				return false;
//			snip(a, 0, big_s_bps);
//			if (a.Length < 2 * big_s_bps)
//				return true;
//		}
//		return true;
	}

	bool ReadGroup::recordUMAltRG(ostream& f1, ostream& f2,
			BamTools::BamAlignment &a, BamTools::BamAlignment &b, int big_s_bps,
			int n_cutoff) {
		BamAlignment mapped = a.IsMapped() ? a : b;
		BamAlignment unmapped = a.IsMapped() ? b : a;

		if (unmapped.Length < big_s_bps || mapped.Length < 2 * big_s_bps)
			return false;

		mapped.Name += "mu";
		unmapped.Name += "mu";

		/* unmapped >= big_s_bps.  Make a fragment from the 5' end */
		BamAlignment cp1 = snip(unmapped, 0, big_s_bps);
		cp1.Name += "1";
		mapped.Name += "1";
		writeFQpair(f1, mapped, f2, cp1, n_cutoff);

		if (unmapped.Length < 2 * big_s_bps)
			return true;

		/* unmapped >= 2*big_s_bps.  Make another fragment from the 3' end */
		BamAlignment cp2 = snip(unmapped, unmapped.Length - big_s_bps,
				big_s_bps);
		cp2.Name += "2";
		mapped.Name[mapped.Name.size() - 1] = '2';
		writeFQpair(f1, mapped, f2, cp2, n_cutoff);

		return true;
	}

	bool ReadGroup::recordSC(BamAlignment &a, BamAlignment &b, int clippedMate,
			ofstream &clipfile, ofstream& split1, ofstream& split2,
			int big_s_bps, int frag_size, int n_cutoff) {
		/* Make a copy of each alignment so we don't modify the base data */
		BamAlignment clipped, mate;
		if (clippedMate == getMateNumber(a)) {
			clipped = a;
			mate = b;
		} else {
			clipped = b;
			mate = a;
		}

		/* The trim could've ruined the big S op */
		if (clipped.CigarData.empty()) {
			return false;
		}

		/* Handle two completely different filesets and parameters. */
		/* The softclips.fq file uses 'frag_size' to determine if a
		 * read is a valid softclip. */
		if ((isBigS(clipped, clipped.CigarData.front(), frag_size)
				|| isBigS(clipped, clipped.CigarData.back(), frag_size))
				&& clipped.Length >= 2 * frag_size
				&& mate.Length >= 2 * frag_size) {
			writeFQ(clipfile, clipped);
			split_read(clipped, split1, split2, frag_size, n_cutoff);
		}

		/* The _1 and _2 files use 'big_s_bps' to determine if a read
		 * is a valid softclip. */
		if ((isBigS(clipped, clipped.CigarData.front(), big_s_bps)
				|| isBigS(clipped, clipped.CigarData.back(), big_s_bps))
				&& mate.Length >= 2 * big_s_bps) {
			BamAlignment copy(clipped);
			clipAlignment(copy);
			copy.Name += "sc";
			mate.Name += "sc";
			writeFQpair(f1, copy, f2, mate, n_cutoff);
		}

		return true;
	}

	bool ReadGroup::recordSCAlt(ostream& f1, ostream& f2,
			BamTools::BamAlignment &a, BamTools::BamAlignment &b,
			int clippedMate, ofstream &clipfile, ofstream& split1,
			ofstream& split2, int big_s_bps, int frag_size, int n_cutoff) {
		/* Make a copy of each alignment so we don't modify the base data */

		if (clippedMate == getMateNumber(a)) {
			/* The trim could've ruined the big S op */
			if (a.CigarData.empty()) {
				return false;
			}
		} else {
			/* The trim could've ruined the big S op */
			if (b.CigarData.empty()) {
				return false;
			}
		}
		BamAlignment clipped, mate;
		if (clippedMate == getMateNumber(a)) {
			clipped = a;
			mate = b;
		} else {
			clipped = b;
			mate = a;
		}

		/* Handle two completely different filesets and parameters. */
		/* The softclips.fq file uses 'frag_size' to determine if a
		 * read is a valid softclip. */
		if ((isBigS(clipped, clipped.CigarData.front(), frag_size)
				|| isBigS(clipped, clipped.CigarData.back(), frag_size))
				&& clipped.Length >= 2 * frag_size
				&& mate.Length >= 2 * frag_size) {
			writeFQ(clipfile, clipped);
			split_read(clipped, split1, split2, frag_size, n_cutoff);
		}

		/* The _1 and _2 files use 'big_s_bps' to determine if a read
		 * is a valid softclip. */
		if ((isBigS(clipped, clipped.CigarData.front(), big_s_bps)
				|| isBigS(clipped, clipped.CigarData.back(), big_s_bps))
				&& mate.Length >= 2 * big_s_bps) {
			BamAlignment copy(clipped);
			clipAlignment(copy);
			copy.Name += "sc";
			mate.Name += "sc";
			writeFQpair(f1, copy, f2, mate, n_cutoff);
		}

		return true;
	}

	bool ReadGroup::recordSCAltNoRG(BamTools::BamAlignment &clipped,
			BamTools::BamAlignment &mate, int clippedMate, ostream &clipfile,
			ostream& split1, ostream& split2, int big_s_bps, int frag_size,
			int n_cutoff) {
		/* Make a copy of each alignment so we don't modify the base data */
//		BamAlignment clipped, mate;
//		if (clippedMate == getMateNumber(a)) {
//			clipped = a;
//			mate = b;
//		} else {
//			clipped = b;
//			mate = a;
//		}

		/* The trim could've ruined the big S op */
		if (clipped.CigarData.empty())
			return false;

		/* Handle two completely different filesets and parameters. */
		/* The softclips.fq file uses 'frag_size' to determine if a
		 * read is a valid softclip. */
//		bool debug = clipped.Name == "ST-E00104:502:HFJN5CCXX:7:2122:32400:13527";
//		if(debug) {
//			cout << "[ReadGroup.recordSCAltNoRG] frag_size: " << frag_size << "\n";
//			cout << "[ReadGroup.recordSCAltNoRG] clipped.Length: " << clipped.Length << "\n";
//			cout << "[ReadGroup.recordSCAltNoRG] mate.Length: " << mate.Length << "\n";
//			cout << "[ReadGroup.recordSCAltNoRG] clipped_front: " << clipped.CigarData.front().Length << clipped.CigarData.front().Type << "\n";
//			cout << "[ReadGroup.recordSCAltNoRG] clipped_back: " << clipped.CigarData.back().Length << clipped.CigarData.back().Type << "\n";
//		}
		if ((isBigS(clipped, clipped.CigarData.front(), frag_size)
				|| isBigS(clipped, clipped.CigarData.back(), frag_size))
				&& clipped.Length >= 2 * frag_size
				&& mate.Length >= 2 * frag_size) {
//			if(debug) {
//				cout << "[ReadGroup.recordSCAltNoRG] has written\n";
//			}
			writeFQ(clipfile, clipped);
			split_read(clipped, split1, split2, frag_size, n_cutoff);
		}
//		if ((isBigS(mate, mate.CigarData.front(), frag_size)
//				|| isBigS(mate, mate.CigarData.back(), frag_size))
//				&& mate.Length >= 2 * frag_size
//				&& clipped.Length >= 2 * frag_size) {
//			writeFQ(clipfile, mate);
//			split_read(mate, split1, split2, frag_size, n_cutoff);
//		}

		return true;
	}
	bool ReadGroup::recordSCAltRG(ostream& f1, ostream& f2,
			BamTools::BamAlignment &clipped, BamTools::BamAlignment &mate,
			int clippedMate, int big_s_bps, int frag_size, int n_cutoff) {
		/* Make a copy of each alignment so we don't modify the base data */

//		if (clippedMate == getMateNumber(a)) {
//			/* The trim could've ruined the big S op */
//			if (a.CigarData.empty()) {
//				return false;
//			}
//		} else {
//			/* The trim could've ruined the big S op */
//			if (b.CigarData.empty()) {
//				return false;
//			}
//		}
		if (clipped.CigarData.empty()) {
			return false;
		}

//		BamAlignment clipped, mate;
//		if (clippedMate == getMateNumber(a)) {
//			clipped = a;
//			mate = b;
//		} else {
//			clipped = b;
//			mate = a;
//		}

		// we do not write clip and split files in this stage
		/* Handle two completely different filesets and parameters. */
		/* The softclips.fq file uses 'frag_size' to determine if a
		 * read is a valid softclip. */
//		if ((isBigS(clipped, clipped.CigarData.front(), frag_size)
//				|| isBigS(clipped, clipped.CigarData.back(), frag_size))
//				&& clipped.Length >= 2 * frag_size
//				&& mate.Length >= 2 * frag_size) {
//			writeFQ(clipfile, clipped);
//			split_read(clipped, split1, split2, frag_size, n_cutoff);
//		}
		/* The _1 and _2 files use 'big_s_bps' to determine if a read
		 * is a valid softclip. */
//		bool debug = clipped.Name == "ST-E00104:502:HFJN5CCXX:1:1108:32826:10890";
//		if(debug) {
//			cout << "[ReadGroup.recordSCAltRG] frag_size: " << frag_size << "\n";
//			cout << "[ReadGroup.recordSCAltRG] clipped.Length: " << clipped.Length << "\n";
//			cout << "[ReadGroup.recordSCAltRG] mate.Length: " << mate.Length << "\n";
//			cout << "[ReadGroup.recordSCAltRG] clipped_front: " << clipped.CigarData.front().Length << clipped.CigarData.front().Type << "\n";
//			cout << "[ReadGroup.recordSCAltRG] clipped_back: " << clipped.CigarData.back().Length << clipped.CigarData.back().Type << "\n";
//		}
		if ((isBigS(clipped, clipped.CigarData.front(), big_s_bps)
				|| isBigS(clipped, clipped.CigarData.back(), big_s_bps))
				&& mate.Length >= 2 * big_s_bps) {
			BamAlignment copy(clipped);
			clipAlignment(copy);
			size_t n_amb = count(copy.QueryBases.begin(), copy.QueryBases.end(), 'N');
			n_amb += count(mate.QueryBases.begin(), mate.QueryBases.end(), 'N');
			if(n_amb > static_cast<size_t>(n_cutoff)) {
				BamAlignment copy(mate);
				clipAlignment(copy);
				copy.Name += "sc";
				clipped.Name += "sc";
//				if(debug) {
//					cout << "[ReadGroup.recordSCAltRG] has written(mate)\n";
//				}
				writeFQpair(f1, copy, f2, clipped, n_cutoff);
			} else {
				copy.Name += "sc";
				mate.Name += "sc";
//				if(debug) {
//					cout << "[ReadGroup.recordSCAltRG] has written(clipped)\n";
//				}
				writeFQpair(f1, copy, f2, mate, n_cutoff);
			}
		} else {
			if ((isBigS(mate, mate.CigarData.front(), big_s_bps)
							|| isBigS(mate, mate.CigarData.back(), big_s_bps))
							&& clipped.Length >= 2 * big_s_bps) {
				BamAlignment copy(mate);
				clipAlignment(copy);
				copy.Name += "sc";
				clipped.Name += "sc";
//				if(debug) {
//					cout << "[ReadGroup.recordSCAltRG] has written(mate)\n";
//				}
				writeFQpair(f1, copy, f2, clipped, n_cutoff);
			}
		}

		return true;
	}

	bool ReadGroup::recordSCAltOnlyClip_Split(BamTools::BamAlignment &a,
			BamTools::BamAlignment &b, int clippedMate, ofstream &clipfile,
			ofstream& split1, ofstream& split2, int big_s_bps, int frag_size,
			int n_cutoff) {
		/* Make a copy of each alignment so we don't modify the base data */
//		int the_mate_number = getMateNumber(a);
//		BamAlignment clipped, mate;
		if (clippedMate == getMateNumber(a)) {
//			clipped = a;
//			mate = b;
			/* The trim could've ruined the big S op */
			if (a.CigarData.empty()) {
				return false;
			}
			/* Handle two completely different filesets and parameters. */
			/* The softclips.fq file uses 'frag_size' to determine if a
			 * read is a valid softclip. */
			if ((isBigS(a, a.CigarData.front(), frag_size)
					|| isBigS(a, a.CigarData.back(), frag_size))
					&& a.Length >= 2 * frag_size && b.Length >= 2 * frag_size) {
				writeFQ(clipfile, a);
				split_read(a, split1, split2, frag_size, n_cutoff);
			}
		} else {
//			clipped = b;
//			mate = a;
			/* The trim could've ruined the big S op */
			if (b.CigarData.empty()) {
				return false;
			}

			/* Handle two completely different filesets and parameters. */
			/* The softclips.fq file uses 'frag_size' to determine if a
			 * read is a valid softclip. */
			if ((isBigS(b, b.CigarData.front(), frag_size)
					|| isBigS(b, b.CigarData.back(), frag_size))
					&& b.Length >= 2 * frag_size && a.Length >= 2 * frag_size) {
				writeFQ(clipfile, b);
				split_read(b, split1, split2, frag_size, n_cutoff);
			}
		}

		return true;
	}
	string ReadGroup::getName() {
		return name;
	}

	int ReadGroup::getNReads() {
		return nreads;
	}

	vector<int> ReadGroup::getReadlens() {
		return readlens;
	}

	bool ReadGroup::is_blacklisted() {
		return blacklisted;
	}
	std::ostream & operator<<(ostream &out, const ReadGroup &rg) {
		out << "Readgroup [" + rg.name + "]: " << rg.nreads << " reads";
		return out;
	}

	void ReadGroup::trim_read(BamAlignment &al, int q, int qenc) {
		int i, best_i, s, best_s;

		int start_idx = al.IsReverseStrand() ? 0 : al.Length - 1;
		int increment = al.IsReverseStrand() ? 1 : -1;

		/* substr always throws errors on 0 length strings */
		if ((al.Qualities[start_idx] - known_qencodings[qenc].min) >= q
				|| al.Length <= 0)
			return;

		s = 0;
		best_s = INT_MIN;
		best_i = start_idx;
		for (i = start_idx; i >= 0 && i < al.Length; i += increment) {
			s += q - (al.Qualities[i] - known_qencodings[qenc].min);
			if (s > best_s) {
				best_s = s;
				best_i = i;
			}
		}

		int x = al.IsReverseStrand() ? best_i + 1 : al.Length - best_i;

		/* if x >= al.Length, then the new length is 0 */
		int new_start = al.IsReverseStrand() ? min(x, al.Length - 1) : 0;
		int new_len = al.Length - x;

		al.Length = new_len;
		try {
			al.QueryBases = al.QueryBases.substr(new_start, new_len);
			al.Qualities = al.Qualities.substr(new_start, new_len);
		} catch (exception &e) {
			cout << "ERROR: substr failed in trim_read" << endl;
			cout << "trim_read(" << al.Name << ", " << q << ")" << endl;
			cout << new_start << " " << new_len << " x: " << x << " best_s: "
					<< best_s << " best_i: " << best_i << endl;
			cout << al.IsReverseStrand() << endl;
			cout << al.QueryBases << endl;
			cout << al.Qualities << endl;
			exit(1);
		}

		/* The trim may extend past one CIGAR op, but we only *
		 * care about S ops in the first or last position.    */
		if (al.CigarData.size() == 0)
			return;

		int idx = al.IsReverseStrand() ? 0 : al.CigarData.size() - 1;
		al.CigarData[idx].Length = max(0,
				(signed) al.CigarData[idx].Length - x);
	}
	void ReadGroup::split_read(BamTools::BamAlignment &al, ostream & split1,
			ostream& split2, int frag_size, int n_cutoff) {
		if (al.Length < 2 * frag_size) {
			cout << "ERROR: split_read tried to save alignment " << al.Name
					<< " with length " << al.Length << " when -s=" << frag_size
					<< endl;
			exit(1);
		}
//		bool debug = al.Name == "ST-E00104:502:HFJN5CCXX:7:2122:32400:13527";
		BamAlignment left = ReadGroup::snip(al, 0, frag_size);
		string aln_length = boost::lexical_cast<string>(al.Length);
		left.Name += "_" + aln_length;
		BamAlignment right = ReadGroup::snip(al, al.Length - frag_size,
				frag_size);
		right.Name += "_" + aln_length;
//		if(debug) {
//			cout << "[ReadGroup.split_read] left: " << left.Name << ": " << left.QueryBases << "\n";
//			cout << "[ReadGroup.split_read] right: " << right.Name << ": " << right.QueryBases << "\n";
//		}
		ReadGroup::writeFQpair(split1, left, split2, right, n_cutoff);
	}

}

