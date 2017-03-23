#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

/*
 * Original Copyright
Copyright (c) 2011 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/*
 * The license notification regarding the modified part of source code
 * Modified: since May 9, 2015
 * Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *       : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 * Email : abysslover@gmail.com
 * License : GPLv3
 */

#include <vector>
#include <algorithm>
#include <iostream>

namespace emma {

template <class T, typename K = int64_t>
class IntervalType {
public:
    K start;
    K stop;
    T value;
    IntervalType(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

template <class T, typename K>
K intervalStart(const IntervalType<T,K>& i) {
    return i.start;
}

template <class T, typename K>
K intervalStop(const IntervalType<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
  std::ostream& operator<<(std::ostream& out, IntervalType<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = int64_t>
class IntervalStartSorter {
public:
    bool operator() (const IntervalType<T,K>& a, const IntervalType<T,K>& b) {
        return a.start < b.start;
    }
};

template <class T, typename K = int64_t>
class IntervalTreeType {

public:
    typedef IntervalType<T,K> interval;
    typedef std::vector<interval> intervalVector;
    typedef IntervalTreeType<T,K> intervalTree;

    intervalVector intervals;
    intervalTree* left;
    intervalTree* right;
    K center;

    IntervalTreeType<T,K>(void)
        : left(NULL)
        , right(NULL)
        , center(0)
    { }

    IntervalTreeType<T,K>(const intervalTree& other)
        : left(NULL)
        , right(NULL)
    {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = new intervalTree(*other.left);
        }
        if (other.right) {
            right = new intervalTree(*other.right);
        }
    }

    IntervalTreeType<T,K>& operator=(const intervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = new intervalTree(*other.left);
        } else {
            if (left) delete left;
            left = NULL;
        }
        if (other.right) {
            right = new intervalTree(*other.right);
        } else {
            if (right) delete right;
            right = NULL;
        }
        return *this;
    }

    IntervalTreeType<T,K>(
            intervalVector& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            K leftextent = 0,
            K rightextent = 0,
            std::size_t maxbucket = 512
            )
        : left(NULL)
        , right(NULL)
    {

        --depth;
        IntervalStartSorter<T,K> intervalStartSorter;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            intervals = ivals;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
              std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            K leftp = 0;
            K rightp = 0;
            K centerp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                std::vector<K> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            //centerp = ( leftp + rightp ) / 2;
            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            intervalVector lefts;
            intervalVector rights;

            for (typename intervalVector::iterator i = ivals.begin(); i != ivals.end(); ++i) {
                interval& interval = *i;
                if (interval.stop < center) {
                    lefts.push_back(interval);
                } else if (interval.start > center) {
                    rights.push_back(interval);
                } else {
                    intervals.push_back(interval);
                }
            }

            if (!lefts.empty()) {
                left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
            }
            if (!rights.empty()) {
                right = new intervalTree(rights, depth, minbucket, centerp, rightp);
            }
        }
    }

    // this algorithm expect a right open range such as [A, B).
    void find_point(K start, intervalVector& overlapping) const {
    	if (!intervals.empty() && ! (start < intervals.front().start)) {
			for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
				const interval& interval = *i;
				if (interval.stop > start && interval.start <= start) {
					overlapping.push_back(interval);
				}
			}
		}

		if (left && start <= center) {
			left->find_point(start, overlapping);
		}

		if (right && start > center) {
			right->find_point(start, overlapping);
		}
    }

    void find_overlap(K start, K stop, intervalVector& overlap_vec) const {
		if (!intervals.empty() && ! (stop < intervals.front().start)) {
			for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
				const interval& interval = *i;
				if (interval.stop >= start && interval.start <= stop) {
					overlap_vec.push_back(interval);
				}
			}
		}

		if (left && start <= center) {
			left->find_overlap(start, stop, overlap_vec);
		}

		if (right && stop >= center) {
			right->find_overlap(start, stop, overlap_vec);
		}
	}
//		not working
//    void find_overlap_mt(K start, K stop, K center, intervalVector& overlap_vec, const intervalVector& intervals, intervalTree* left, intervalTree* right) const {
//    		if (!intervals.empty() && ! (stop < intervals.front().start)) {
//    			for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
//    				const interval& interval = *i;
//    				if (interval.stop >= start && interval.start <= stop) {
//    					overlap_vec.push_back(interval);
//    				}
//    			}
//    		}
//
//    		if (left && start <= center) {
//    			left->find_overlap_mt(start, stop, center, overlap_vec, intervals, left, right);
//    		}
//
//    		if (right && stop >= center) {
//    			right->find_overlap_mt(start, stop, center, overlap_vec, intervals, left, right);
//    		}
//    	}
//    void findContained(K start, K stop, intervalVector& contained) const {
//        if (!intervals.empty() && ! (stop < intervals.front().start)) {
//            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
//                const interval& interval = *i;
//                if (interval.start >= start && interval.stop <= stop) {
//                    contained.push_back(interval);
//                }
//            }
//        }
//
//        if (left && start <= center) {
//            left->findContained(start, stop, contained);
//        }
//
//        if (right && stop >= center) {
//            right->findContained(start, stop, contained);
//        }
//
//    }

    ~IntervalTreeType(void) {
        // traverse the left and right
        // delete them all the way down
        if (left) {
            delete left;
        }
        if (right) {
            delete right;
        }
    }

};
//typedef IntervalType<std::string> PopInterval;
//typedef std::vector<PopInterval> PopIntervalVector;
//typedef IntervalTreeType<std::string> PopIntervalTree;
//
//typedef IntervalType<uint64_t> PopNumericInterval;
//typedef std::vector<PopNumericInterval> PopNumericIntervalVector;
//typedef IntervalTreeType<uint64_t> PopNumericIntervalTree;
//
//struct Description {
//	string text;
//	uint64_t size;
//};
} /* namespace emma */
#endif
