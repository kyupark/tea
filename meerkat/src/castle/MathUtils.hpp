/*
 * MathUtil.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: el174
 */

#ifndef CASTLE_MATHUTILS_HPP_
#define CASTLE_MATHUTILS_HPP_

#include <vector>
#include <numeric>
#include <algorithm>

namespace castle {

using namespace std;
class MathUtils {
public:
	MathUtils();
	~MathUtils();
	static double mean(const vector<double>& v);
	static double std_dev(const vector<double>& v);
};

inline double MathUtils::mean(const vector<double>& v) {
	double sum = accumulate(v.begin(), v.end(), 0.0);
	return sum / v.size();
}

inline double MathUtils::std_dev(const vector<double>& v) {
	double a_mean = mean(v);
	vector<double> diff(v.size());
	transform(v.begin(), v.end(), diff.begin(), [a_mean](double x) {
		return x - a_mean;
	});
	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = sqrt(sq_sum / v.size());
	return stdev;
}

} /* namespace tea */

#endif /* CASTLE_MATHUTILS_HPP_ */
