/*
 * ReadGroupAlt.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: el174
 */

#include "ReadGroupAlt.hpp"

namespace meerkat {
	ReadGroupAlt::ReadGroupAlt() : nreads(0) {
	}

	ReadGroupAlt::~ReadGroupAlt() {
	}

	void ReadGroupAlt::witness(BamTools::BamAlignment &al, int32_t max_isize,
			int isize_samples) {
		++nreads;
		/* Selects first 'isize_samples' samples instead of random sampling */
		/* Don't use negative insert sizes.  If it's negative, then the mate  *
		 * has a positive insert of the same size.  If both are recorded, the *
		 * average will always be 0.                                          */
		if (al.InsertSize > 0 && al.InsertSize <= max_isize
				&& static_cast<int>(inserts.size()) < isize_samples) {
			inserts.push_back(al.InsertSize);
		}

		if (find(readlens.begin(), readlens.end(), al.Length)
				== readlens.end()) {
			readlens.push_back(al.Length);
		}
	}
} /* namespace meerkat */
