/*
 * StringUtils.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "StringUtils.hpp"

namespace castle {

StringUtils::StringUtils() {
}

StringUtils::~StringUtils() {
}
char StringUtils::get_complement_base(char c) {
	switch (c) {
	case 'A':
	case 'a':
		return 'T';
	case 'T':
	case 't':
		return 'A';
	case 'C':
	case 'c':
		return 'G';
	case 'G':
	case 'g':
		return 'C';
	default:
		return 'N';
	}
}

int64_t StringUtils::count_differences(const std::string& s1, const std::string& s2, size_t n) {
	int64_t num_diff = 0;
	for (size_t i = 0; i < n; ++i) {
		if (s1[i] != s2[i])
			++num_diff;
	}
	return num_diff;
}

} /* namespace castle */
