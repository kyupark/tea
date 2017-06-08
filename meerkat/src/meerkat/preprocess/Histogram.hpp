/*
 * Histogram.hpp
 *
 *  Created on: Jun 2, 2016
 *      Author: el174
 */

#ifndef MEERKAT_HISTOGRAM_HPP_
#define MEERKAT_HISTOGRAM_HPP_

#include <iostream>
#include <map>
#include <set>
#include "api/BamAlignment.h"

namespace meerkat {

	class Histogram {
		public:
			Histogram();
			~Histogram();
			Histogram(const Histogram& b) {
				len_freqs = b.len_freqs;
				total = b.total;
			}
			Histogram& operator=(const Histogram& b) {
				// check for self-assignment
				if (this == &b) {
					return *this;
				}
				len_freqs = b.len_freqs;
				total = b.total;
				return *this;
			}
			Histogram& operator+=(const Histogram& b) {
				// check for self-assignment
				if (this == &b) {
					return *this;
				}
				for(auto itr = b.len_freqs.begin(); b.len_freqs.end() != itr; ++itr) {
					len_freqs[itr->first] += itr->second;
				}
				total += b.total;
				return *this;
			}
			bool operator==(const Histogram& b) const {
				if (this == &b)
					return true;
				if(total != b.total){
					return false;
				}
				if(len_freqs.size() != b.len_freqs.size()) {
					return false;
				}
				std::set<int> keys;
				for(auto itr = len_freqs.begin(); len_freqs.end() != itr; ++itr) {
					keys.insert(itr->first);
				}
				for(auto itr = b.len_freqs.begin(); b.len_freqs.end() != itr; ++itr) {
					keys.insert(itr->first);
				}
				for(auto itr = keys.begin(); keys.end() != itr; ++itr) {
					int a_key = *itr;
					auto lhs_value = len_freqs.find(a_key);
					auto rhs_value = b.len_freqs.find(a_key);
					if(len_freqs.end() == lhs_value ||
							b.len_freqs.end() == rhs_value) {
						return false;
					}
					if(lhs_value->second != rhs_value->second) {
						return false;
					}
				}
				return true;
			}
			bool operator!=(const Histogram& v) const {
				return !(*this == v);
			}
			void add(const BamTools::BamAlignment &al);
			void print(std::ofstream &f);

		private:
			std::map<int, int> len_freqs;
			int total;
	};

} /* namespace meerkat */

#endif /* MEERKAT_HISTOGRAM_HPP_ */
