/*
 * ReadGroupAlt.hpp
 *
 *  Created on: Jun 6, 2016
 *      Author: el174
 */

#ifndef MEERKAT_READGROUPALT_HPP_
#define MEERKAT_READGROUPALT_HPP_

#include <vector>
#include <algorithm>
#include "api/BamReader.h"

namespace meerkat {
	using namespace std;
	class ReadGroupAlt {
		public:
			int64_t nreads;
			vector<int32_t> readlens;
			vector<int> inserts;
		public:
			ReadGroupAlt();
			~ReadGroupAlt();

			ReadGroupAlt(const ReadGroupAlt& b) {
				nreads = b.nreads;
				readlens = b.readlens;
				inserts = b.inserts;
			}
			ReadGroupAlt& operator=(const ReadGroupAlt& b) {
				// check for self-assignment
				if (this == &b) {
					return *this;
				}
				nreads = b.nreads;
				readlens = b.readlens;
				inserts = b.inserts;
				return *this;
			}
			ReadGroupAlt& operator+=(const ReadGroupAlt& b) {
				// check for self-assignment
				if (this == &b) {
					return *this;
				}
				nreads += b.nreads;
				readlens.insert(readlens.end(), b.readlens.begin(), b.readlens.end());
				inserts.insert(inserts.end(), b.inserts.begin(), b.inserts.end());
				return *this;
			}

			bool operator==(const ReadGroupAlt& b) const {
				if (this == &b) {
					return true;
				}
				return false;
			}
			bool operator!=(const ReadGroupAlt& b) const {
				return !(*this == b);
			}
			void witness(BamTools::BamAlignment &al, int32_t max_isize,
					int isize_samples);
	};

} /* namespace meerkat */

#endif /* MEERKAT_READGROUPALT_HPP_ */
