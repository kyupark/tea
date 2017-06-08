/*
 * StringUtils.hpp
 *
 *  Created on: Feb 18, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef STRINGUTILS_HPP_
#define STRINGUTILS_HPP_

#include <string>
#include <sstream>
#include <iostream>
#include <bitset>
#include <vector>
#include <set>

#include <boost/range.hpp>
#include <boost/lexical_cast.hpp>

namespace castle {
typedef std::string::const_iterator str_c_iter;
typedef std::string::iterator str_iter;
typedef boost::iterator_range<str_c_iter> c_string_view;
typedef boost::iterator_range<str_iter> string_view;
class StringUtils {
public:
	StringUtils();
	~StringUtils();
	template<typename C>
	static void string_multi_split(std::string& str, char const* delimiters,
			C& ret_array);
	template<typename C>
	static void c_string_multi_split(const std::string& str,
			char const* delimiters, C& ret_array);
	// including empty tokens
	template<typename C>
	static void tokenize(const std::string& str, const std::string& delims,
			C& ret_array);
	template<typename C>
	static bool is_not_containing_ambiguous_bases(const C& str);
	template<typename C>
	static bool is_not_containing_ambiguous_bases_in_range(const C& str,
			int64_t start_offset, int64_t end_offset);
	static char get_complement_base(char c);
	template<typename C>
	static std::string get_reverse_complement(const C& str);
	template<typename C>
	static std::string get_reverse_string(const C& str);
	template<typename C>
	static std::string get_reverse_complement(C& str);
	template<typename C>
	static std::string get_reverse_string(C& str);
	template<typename C>
	static std::string get_longest_common_sub_string(const C& str1,
			const C& str2);
	template<typename C>
	static void replace_all_non_ACGT_with_N(C& str);
	static size_t rfind_one_of_multiple_patterns(const std::string& s,
			const std::vector<std::string>& patterns);
	static int64_t count_differences(const std::string& s1,
			const std::string& s2, size_t n);
	// trim from start
	static std::string& ltrim(std::string &s);

	// trim from end
	static std::string& rtrim(std::string &s);

	// trim from both ends
	static std::string& trim(std::string &s);
	template<typename T>
	static std::string join(const T& items, const std::string& comma);

	// code from http://stackoverflow.com/questions/4351371/c-performance-challenge-integer-to-stdstring-conversion
	template<typename T>
	static T reduce2(T v) {
		/* pre: ((short*)&v)[i] < 100 for all i
		 *  post:
		 *     ((char*)&v)[2i] = ((short*)&v)[i] / 10
		 *     ((char*)&v)[2i + 1] = ((short*)&v)[i] % 10
		 *
		 *     That is, we split each short in v into its ones and tens digits
		 */

		/* t < 100 --> (t * 410) >> 12 == t / 10
		 *			&& (t * 410) < 0x10000
		 *
		 * For the least order short that's all we need, for the others the
		 * shift doesn't drop the remainder so we mask those out
		 */
		T k = ((v * 410) >> 12) & 0x000F000F000F000Full;

		/*
		 * Then just subtract out the tens digit to get the ones digit and
		 * shift them into the right place
		 */
		return (((v - k * 10) << 8) + k);
	}

	template<typename T>
	static T reduce4(T v) {
		/* pre: ((unsigned*)&v)[i] < 10000 for all i
		 *
		 *  preReduce2:
		 *     ((short*)&v)[2i] = ((unsigned*)&v)[i] / 100
		 *     ((short*)&v)[2i + 1] = ((unsigned*)&v)[i] % 100
		 *
		 *     That is, we split each int in v into its one/ten and hundred/thousand
		 *     digit pairs. Put them into the corresponding short positions and then
		 *     call reduce2 to finish the splitting
		 */

		/* This is basically the same as reduce2 with different constants
		 */
		T k = ((v * 10486) >> 20) & 0xFF000000FFull;
		return reduce2(((v - k * 100) << 16) + (k));
	}

	typedef unsigned long long ull;
	static inline ull reduce8(ull v) {
		/* pre: v  < 100000000
		 *
		 *  preReduce4:
		 *     ((unsigned*)&v)[0] = v / 10000
		 *     ((unsigned*)&v)[1] = v % 10000
		 *
		 *     This should be familiar now, split v into the two 4-digit segments,
		 *     put them in the right place, and let reduce4 continue the splitting
		 */

		/* Again, use the same method as reduce4 and reduce2 with correctly tailored constants
		 */
		ull k = ((v * 3518437209u) >> 45);
		return reduce4(((v - k * 10000) << 32) + (k));
	}

	template<typename T>
	static std::string itostr(T o) {
		/*
		 * Use of this union is not strictly compliant, but, really,
		 * who cares? This is just for fun :)
		 *
		 * Our ones digit will be in str[15]
		 *
		 * We don't actually need the first 6 bytes, but w/e
		 */
		union {
			char str[16];
			unsigned short u2[8];
			unsigned u4[4];
			unsigned long long u8[2];
		};

		/* Technically should be "... ? unsigned(~0) + 1 ..." to ensure correct behavior I think */
		/* Tends to compile to: v = (o ^ (o >> 31)) - (o >> 31); */
		unsigned v = o < 0 ? ~o + 1 : o;

		/* We want u2[3] = v / 100000000 ... that is, the first 2 bytes of the decimal rep */

		/* This is the same as in reduce8, that is divide by 10000. So u8[0] = v / 10000 */
		u8[0] = (ull(v) * 3518437209u) >> 45;

		/* Now we want u2[3] = u8[0] / 10000.
		 * If we added " >> 48 " to the end of the calculation below we would get u8[0] = u8[0] / 10000
		 * Note though that in little endian byte ordering u2[3] is the top 2 bytes of u8[0]
		 * and 64 - 16 == 48... Aha, We've got what we want, the rest of u8[0] is junk, but we don't care
		 */
		u8[0] = (u8[0] * 28147497672ull);

		/* Then just subtract out those digits from v so that u8[1] now holds
		 * the low 8 decimal digits of v
		 */
		u8[1] = v - u2[3] * 100000000;

		/* Split u8[1] into its 8 decimal digits */
		u8[1] = reduce8(u8[1]);

		/* f will point to the high order non-zero char */
		char* f;

		/* branch post: f is at the correct short (but not necessarily the correct byte) */
		if (u2[3]) {
			/* split the top two digits into their respective chars */
			u2[3] = reduce2(u2[3]);
			f = str + 6;
		} else {
			/* a sort of binary search on first non-zero digit */
			unsigned short* k = u4[2] ? u2 + 4 : u2 + 6;
			f = *k ? (char*) k : (char*) (k + 1);
		}
		/* update f to its final position */
		if (!*f)
			f++;

		/* '0' == 0x30 and i < 10 --> i <= 0xF ... that is, i | 0x30 = 'i' *
		 * Note that we could do u8[0] |= ... u8[1] |= ... but the corresponding
		 * x86-64 operation cannot use a 64 bit immediate value whereas the
		 * 32 bit 'or' can use a 32 bit immediate.
		 */
		u4[1] |= 0x30303030;
		u4[2] |= 0x30303030;
		u4[3] |= 0x30303030;

		/* Add the negative sign... note that o is just the original parameter passed */
		if (o < 0)
			*--f = '-';

		/* gcc basically forwards this to std::string(f, str + 16)
		 * but msvc handles it way more efficiently
		 */
		return std::string(f, (str + 16) - f);
	}
};

template<typename C>
void StringUtils::string_multi_split(std::string& str, char const* delimiters,
		C& ret_array) {
	ret_array.clear();
	C output;
	std::bitset<255> delims;
	while (*delimiters) {
		unsigned char code = *delimiters++;
		delims[code] = true;
	}
	str_iter beg;
	bool in_token = false;
	for (str_iter it = str.begin(), end = str.end(); it != end; ++it) {
		if (delims[*it]) {
			if (in_token) {
				output.push_back(typename C::value_type(beg, it));
				in_token = false;
			}
		} else if (!in_token) {
			beg = it;
			in_token = true;
		}
	}
	if (in_token) {
		output.push_back(typename C::value_type(beg, str.end()));
	}
	output.swap(ret_array);
}

template<typename C>
void StringUtils::c_string_multi_split(const std::string& str,
		char const* delimiters, C& ret_array) {
	ret_array.clear();
	C output;
	std::bitset<255> delims;
	while (*delimiters) {
		unsigned char code = *delimiters++;
		delims[code] = true;
	}
	str_c_iter beg;
	bool in_token = false;
	for (str_c_iter it = str.begin(), end = str.end(); it != end; ++it) {
		if (delims[*it]) {
			if (in_token) {
				output.push_back(typename C::value_type(beg, it));
				in_token = false;
			}
		} else if (!in_token) {
			beg = it;
			in_token = true;
		}
	}
	if (in_token) {
		output.push_back(typename C::value_type(beg, str.end()));
	}
	output.swap(ret_array);
}


template<typename C>
void StringUtils::tokenize(const std::string& str, const std::string& delims,
		C& ret_array) {
	ret_array.clear();
	std::size_t start = 0, end = 0;
	while ((end = str.find(delims, start)) != std::string::npos) {
		ret_array.push_back(str.substr(start, end - start));
		start = end + delims.size();
	}
	ret_array.push_back(str.substr(start));
}
template<typename C>
bool StringUtils::is_not_containing_ambiguous_bases(const C& str) {
	for (auto itr = str.begin(); str.end() != itr; ++itr) {
		if ('N' == *itr) {
			return false;
		}
	}
	return true;
}
template<typename C>
bool StringUtils::is_not_containing_ambiguous_bases_in_range(const C& str,
		int64_t start_offset, int64_t end_offset) {
	const int64_t max_id = str.size();
	for (int64_t base_id = start_offset;
			base_id < end_offset && base_id < max_id; ++base_id) {
		if ('N' == str[base_id]) {
			return false;
		}
	}
	return true;
}
template<typename C>
std::string StringUtils::get_reverse_complement(const C& str) {
	std::string out(str.size(), ' ');
	int64_t last_id = str.size() - 1;
	for (int64_t base_id = last_id; base_id >= 0; --base_id) {
		out[last_id - base_id] = get_complement_base(str[base_id]);
	}
	return out;
}

template<typename C>
std::string StringUtils::get_reverse_string(const C& str) {
	std::string out(str.size(), ' ');
	int64_t last_id = str.size() - 1;
	for (int64_t base_id = last_id; base_id >= 0; --base_id) {
		out[last_id - base_id] = str[base_id];
	}
	return out;
}
template<typename C>
std::string StringUtils::get_reverse_complement(C& str) {
	std::string out(str.size(), ' ');
	int64_t last_id = str.size() - 1;
	for (int64_t base_id = last_id; base_id >= 0; --base_id) {
		out[last_id - base_id] = get_complement_base(str[base_id]);
	}
	return out;
}

template<typename C>
std::string StringUtils::get_reverse_string(C& str) {
	std::string out(str.size(), ' ');
	int64_t last_id = str.size() - 1;
	for (int64_t base_id = last_id; base_id >= 0; --base_id) {
		out[last_id - base_id] = str[base_id];
	}
	return out;
}
template<typename C>
std::string StringUtils::get_longest_common_sub_string(const C& str1,
		const C& str2) {
//	std::cout << str1.size() << "/" << str2.size() << "\n";
	std::set<char *> res;
	std::string res_str;
	int longest = 0;

	int **n = (int **) calloc(str1.length() + 1, sizeof(int *));
	for (uint64_t i = 0; i <= str1.length(); ++i) {
		n[i] = (int *) calloc(str2.length() + 1, sizeof(int));
	}

	for (uint64_t i = 0; i < str1.length(); ++i) {
		for (uint64_t j = 0; j < str2.length(); ++j) {
			if (str1[i] == str2[j]) {
				n[i + 1][j + 1] = n[i][j] + 1;
				if (n[i + 1][j + 1] > longest) {
					longest = n[i + 1][j + 1];
					res.clear();
				}
				if (n[i + 1][j + 1] == longest)
					for (uint64_t it = i - longest + 1; it <= i; it++) {
						res.insert((char *) &str1[it]);
					}
			}
		}

	}
	for (std::set<char *>::const_iterator it = res.begin(); it != res.end();
			it++) {
		res_str.append(1, **it);
	}
	for (uint64_t i = 0; i <= str1.length(); i++) {
		free(n[i]);
	}
	free(n);

	return res_str;
}

template<typename C>
void StringUtils::replace_all_non_ACGT_with_N(C& str) {
	for (uint64_t base_id = 0; base_id < str.size(); ++base_id) {
		char c = str[base_id];
		switch (c) {
		case 'A':
		case 'a':
			c = 'A';
			break;
		case 'C':
		case 'c':
			c = 'C';
			break;
		case 'G':
		case 'g':
			c = 'G';
			break;
		case 'T':
		case 't':
			c = 'T';
			break;
		default:
			c = 'N';
			break;
		}
		str[base_id] = c;
	}
}
// trim from start
inline std::string& StringUtils::ltrim(std::string &s) {
	s.erase(s.begin(),
			std::find_if(s.begin(), s.end(),
					std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
inline std::string& StringUtils::rtrim(std::string &s) {
	s.erase(
			std::find_if(s.rbegin(), s.rend(),
					std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
			s.end());
	return s;
}

// trim from both ends
inline std::string& StringUtils::trim(std::string &s) {
	return ltrim(rtrim(s));
}

inline size_t StringUtils::rfind_one_of_multiple_patterns(const std::string& s,
		const std::vector<std::string>& patterns) {
	for (auto needle : patterns) {
		auto pos = s.rfind(needle);
		if (pos != std::string::npos) {
			return pos;
		}
	}
	return std::string::npos;
}
template<typename T>
std::string StringUtils::join(const T& items, const std::string& comma) {
	std::ostringstream os;
	for(auto& a_str : items) {
		os << a_str << comma;
	}
	auto s = os.str();
	if (!s.empty()) {
		s.resize(s.size()-comma.size());
	}
	return s;
}
} /* namespace castle */
#endif /* STRINGUTILS_HPP_ */
