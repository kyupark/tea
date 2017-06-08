/*
 * Histogram.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: el174
 */

#include "Histogram.hpp"

namespace meerkat {

	Histogram::Histogram() : total(0) {
	}

	Histogram::~Histogram() {
	}
	void Histogram::add(const BamTools::BamAlignment &al)
	{
		++len_freqs[al.Length];
		++total;
	}

	void Histogram::print(std::ofstream &f)
	{
		auto it = len_freqs.begin();

		while (it != len_freqs.end()) {
			f << it->first << "\t" << it->second << "\t"
			  << (it->second / (double)total) << "\n";
			++it;
		}
	}

} /* namespace meerkat */
