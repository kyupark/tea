/*
 * Mechanism.hpp
 *
 *  Created on: Jul 25, 2016
 *      Author: el174
 */

#ifndef MEERKAT_MECHANISM_MECHANISM_HPP_
#define MEERKAT_MECHANISM_MECHANISM_HPP_

#include <map>

#include "../ClusterEntry.hpp"

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

namespace meerkat {
using namespace std;
class Mechanism {
public:
	Mechanism();
	~Mechanism();
	void set_option_parser(const castle::OptionParser& the_options);
	void call_variants();
	void call_variants_alt();
	void collect_repeat_data(map<string, vector<RepeatEntry>>& te, map<string, vector<RepeatEntry>>& sr);
	void collect_variants(vector<deque<string>>& variants);
	void call_mechanism(const map<string, vector<RepeatEntry>>& te, const map<string, vector<RepeatEntry>>& sr, vector<deque<string>>& variants);
	void call_mechanism_alt(const map<string, vector<RepeatEntry>>& te, const map<string, vector<RepeatEntry>>& sr, vector<deque<string>>& variants);
	void clean_temporary_files();

	void call_tei(string& te_class, string& te_name, const map<string, vector<RepeatEntry>>& te_map, const string& chr, int64_t start, int64_t end); //
	bool call_vntr(const map<string, vector<RepeatEntry>>& sr_map, const string& chr, int64_t start, int64_t end);

	bool sr_te_overlap(const map<string, vector<RepeatEntry>>& a_target_map, const string& chr, int64_t bp);
	bool te_overlap(const map<string, vector<RepeatEntry>>& te_map, const string& chr, int64_t bp);

	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
	pair<int64_t, int64_t> overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};
inline void Mechanism::call_tei(string& te_class, string& te_name, const map<string, vector<RepeatEntry>>& te_map, const string& chr, int64_t start, int64_t end) {
	if (abs(end - start) > options.te_size_max) {
		return;
	}

//	my @overlap;
	vector<int64_t> overlaps;
	const auto& te_itr = te_map.find(chr);
	if (te_map.end() != te_itr) {
		for (auto& te : te_itr->second) {
			if (te.start < end && te.end < end && te.start < start && te.end < start) {
				continue;
			}
			if (te.start > end && te.end > end && te.start > start && te.end > start) {
				break;
			}
			if (covered(start, end, te.start, te.end)) {
				pair<int64_t, int64_t> a_pair = overlap(start, end, te.start, te.end);
				int64_t overlap1 = a_pair.first;
				int64_t overlap2 = a_pair.second;
				if (abs(overlap2 - overlap1) / double(abs(te.end - te.start) + 0.5) >= options.ovl) {
					if (abs(overlap2 - overlap1) / double(abs(end - start) + 0.5) >= options.ovl) {
						te_class = te.repeat_class;
						te_name = te.name;
						return;
					} else {
						overlaps.push_back(overlap2 - overlap1);
					}
				}
			}
		}
	}
	if (!overlaps.empty()) {
		int64_t sum_overlap = 0;
		for (auto& an_ovl : overlaps) {
			sum_overlap += an_ovl;
		}
		if (sum_overlap / double(abs(end - start) + 0.5) >= options.ovl) {
			te_class = "complex";
			te_name = "1";
		}
	}
}

inline bool Mechanism::call_vntr(const map<string, vector<RepeatEntry>>& sr_map, const string& chr, int64_t start, int64_t end) {
	if (abs(end - start) > options.te_size_max) {
		return false;
	}
//	my @ overlap;
	const auto sr_itr = sr_map.find(chr);
	if (sr_map.end() != sr_itr) {
		for (auto sr : sr_itr->second) {
			if (sr.start < end && sr.end < end && sr.start < start && sr.end < start) {
				continue;
			}
			if (sr.start > end && sr.end > end && sr.start > start && sr.end > start) {
				break;
			}
			if (covered(start, end, sr.start, sr.end)) {
				pair<int64_t, int64_t> a_pair = overlap(start, end, sr.start, sr.end);
				int64_t overlap1 = a_pair.first;
				int64_t overlap2 = a_pair.second;
				if (abs(overlap2 - overlap1) / (abs(end - start) + 0.5) >= options.ovl) {
					return true;
				}
			}
		}
	}
	return false;
}

//# overlap satellite repeat, simple repeat, low complexity repeat
inline bool Mechanism::sr_te_overlap(const map<string, vector<RepeatEntry>>& a_target_map, const string& chr, int64_t bp) {
	auto a_target_itr = a_target_map.find(chr);
	if (a_target_map.end() == a_target_itr) {
		return false;
	}
	for (auto a_target : a_target_itr->second) {
		if (a_target.start < bp && a_target.end < bp) {
			continue;
		}
		if (a_target.start > bp && a_target.end > bp) {
			break;
		}
		if (bp >= a_target.start && bp <= a_target.end) {
			return true;
		}
	}
	return false;
}

inline pair<int64_t, int64_t> Mechanism::overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {

//	if (-1 == b1 && -1 == b2) {
	if (0 == b1 && 0 == b2) {
		return make_pair(a1, a2);
	}
//		( a1, a2 ) = sort { a <=> b } ( a1, a2 );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( b1, b2 ) = sort { a <=> b } ( b1, b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		int64_t a = (a1 > b1) ? a1 : b1;
		int64_t b = (a2 > b2) ? b2 : a2;
		return make_pair(a, b);
	} else {
		// the original version returns (0, 0)
//		return make_pair(b1, b2);
		return make_pair(0, 0);
	}
}

inline bool Mechanism::covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
//		return 0 if ( a1 =~ /\D/ or a2 =~ /\D/ or b1 =~ /\D/ or b2 =~ /\D/ );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( a1, a2 ) = sort { a <=> b } ( a1, a2 );
//		( b1, b2 ) = sort { a <=> b } ( b1, b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		return true;
	}
	return false;
}

} /* namespace meerkat */

#endif /* MEERKAT_MECHANISM_MECHANISM_HPP_ */
