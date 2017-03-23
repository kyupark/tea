/*
 * BinaryEncoder.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: euncheonlim
 */

#ifndef BINARYENCODER_HPP_
#define BINARYENCODER_HPP_

#include <string>
#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/detail/endian.hpp>
//#include <byteswap.h>
#include <boost/multiprecision/cpp_int.hpp>
//#include "CommonTypes.hpp"
#include "StringUtils.hpp"

namespace castle {
class BinaryEncoder {
public:
	uint32_t readLength;
	uint32_t a;
	uint64_t lowerBitClearMask;
	uint32_t sizeOfEncodeVector;
	uint64_t additional_shifts_for_half_msb;
	int64_t n_shift_premer;
	static const uint64_t INITIAL_SHIFT_BITS_LONGLONG = 123;
	static const uint64_t INITIAL_SHIFT_BITS_LONG = 60;
	static const uint64_t TWO_BITS_INITIAL_SHIFT_BITS_LONG = 62;
	static const uint32_t INITIAL_SHIFT_BITS_INT = 28;
	static const uint32_t MAXIMUM_LENGTH_OF_ENCODE_VECTOR_INT = 10;
	static const uint32_t MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONG = 21;
	static const uint32_t TWO_BITS_MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONG = 32;
	static const uint32_t MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONGLONG = 42;
	static const uint64_t TWO_BITS_UNIT_LONG = 2;
	static const uint64_t TWO_BITS_MASK_LONG = 0b11;
	static const uint64_t ONE_LONG = 1;
	static const uint64_t MAX_LONG = -1;
	static const uint64_t THREE_LONG = 3;
	static const uint64_t INITIAL_CLEAR_MASKS_LONG = (ONE_LONG << INITIAL_SHIFT_BITS_LONG) - 1;
	static const uint64_t HALF_MSB_BITS = 32;
	static const uint64_t TWO_BITS_INITIAL_CLEAR_MASKS_LONG = (ONE_LONG << TWO_BITS_INITIAL_SHIFT_BITS_LONG) - 1;
	static uint64_t nucleotideTwoBitsEncodingLongMap[256];
	static uint64_t nucleotideTwoBitsReverseEncodingLongMap[256];
	static uint64_t nucleotideThreeBitsEncodingLongMap[256];

	static char nucleotideTwoBitsDecodingMap[4];
	static char nucleotideTwoBitsReverseDecodingMap[4];
	static char nucleotideThreeBitsDecodingMap[5];
	static uint64_t twoBitsreverseComplementFirstLookUp[64];
	static uint64_t twoBitsreverseComplementSecondLookUp[64];

	BinaryEncoder();
	~BinaryEncoder();
	BinaryEncoder(const BinaryEncoder& b) {
		readLength = b.readLength;
		a = b.a;
		sizeOfEncodeVector = b.sizeOfEncodeVector;
		lowerBitClearMask = b.lowerBitClearMask;
		additional_shifts_for_half_msb = b.additional_shifts_for_half_msb;
		n_shift_premer = b.n_shift_premer;
	}
	BinaryEncoder& operator=(const BinaryEncoder& b) {
		// check for self-assignment
		if (this == &b)
			return *this;
		readLength = b.readLength;
		a = b.a;
		sizeOfEncodeVector = b.sizeOfEncodeVector;
		lowerBitClearMask = b.lowerBitClearMask;
		additional_shifts_for_half_msb = b.additional_shifts_for_half_msb;
		n_shift_premer = b.n_shift_premer;
		return *this;
	}
	bool operator==(const BinaryEncoder& b) const {
		if (this == &b)
			return true;
		return readLength == b.readLength && a == b.a && sizeOfEncodeVector == b.sizeOfEncodeVector && lowerBitClearMask == b.lowerBitClearMask;
	}
	bool operator!=(const BinaryEncoder & v) const {
		return !(*this == v);
	}
	static uint32_t countSetBitsInBitVector(const uint8_t v) {
		static const int S[] = { 1, 2, 4, 8, 16 }; // Magic Binary Numbers
		static const int B[] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF };
		unsigned int c = v - ((v >> 1) & B[0]);
		c = ((c >> S[1]) & B[1]) + (c & B[1]);
		c = ((c >> S[2]) + c) & B[2];
		return c;
	}

	static uint8_t getDeltaInBitVector(uint8_t vector_old, uint8_t vector_new) {
		return vector_old ^ vector_new;
	}

	static std::string getEncodedBasesFromBitVector(const uint8_t c) {
		std::string result = "";
		if (0 != (1 & c)) {
			result += 'N';
		}
		if (0 != (2 & c)) {
			result += 'A';
		}
		if (0 != (4 & c)) {
			result += 'C';
		}
		if (0 != (8 & c)) {
			result += 'G';
		}
		if (0 != (16 & c)) {
			result += 'T';
		}
		return result;
	}

	static char getReverseBase(const uint8_t c) {
		if ('A' == c)
			return 'T';
		else if ('T' == c)
			return 'A';
		else if ('C' == c)
			return 'G';
		else if ('G' == c)
			return 'C';
		else if ('a' == c)
			return 'T';
		else if ('t' == c)
			return 'A';
		else if ('c' == c)
			return 'G';
		else if ('g' == c)
			return 'C';
		return 'N';
	}

	void setReadLength(uint32_t length);
	void set_a(uint32_t a_a);
	uint32_t getReadLength();
	void setSizeOfEncodeVector(uint32_t size);
	uint32_t getSizeOfEncodeVector();
	template<typename C> uint64_t encodeToTwoBitsLong(const C& cs, int64_t offset);
	template<typename C> int64_t is_contain_N_in_anchor_path(const C& cs, int64_t offset);

	template<typename C> uint64_t encode_to_anchor(uint8_t& anchor_base_id, const C& cs, int64_t offset);
	template<typename C> int64_t is_contain_N_in_anchor(const C& cs, int64_t offset);

	template<typename C> uint64_t encodeToTwoBitsLongReverse(const C& cs, int64_t offset);
	template<typename C> uint64_t encode_to_anchor_reverse(uint8_t& anchor_base_id, const C& cs, int64_t offset);
	template<typename C> int64_t is_contain_N_in_anchor_reverse(const C& cs, int64_t offset);

	template<typename C> uint64_t get_minimizer(const C& cs, const int32_t base_start_id, const int32_t base_end_id, const char target_c);
	template<typename C> uint64_t get_minimizer_with_base_id(int64_t& base_id, const C& cs, const int32_t base_start_id, const int32_t base_end_id,
	        const char target_c);

//	template<typename C> void get_premers(const C& cs, const int64_t base_start_id, const int64_t base_end_id, PremerBlockTuple& premers,
//	        const char target_c);
	template<typename C> uint64_t get_premer(const C& cs, const int64_t base_start_id, const int64_t base_end_id, const char target_c);
	template<typename C> uint64_t get_premer_and_direction(const C& cs, const int64_t base_start_id, const int64_t base_end_id, const char target_c,
	        char& direction);

	template<typename C> uint64_t get_premer_forward(const C& cs, const int64_t base_start_id, const int64_t base_end_id, const char target_c);
	template<typename C> uint64_t get_premer_reverse(const C& cs, const int64_t base_start_id, const int64_t base_end_id, const char target_c);
//	template<typename C> void get_all_premers(const C& cs, const int64_t base_start_id, const int64_t base_end_id, PremerBlockTuple& premers);

	template<typename C> int64_t is_containing_low_complexity(const C& cs, int64_t base_start_id, int64_t max_id);
	char get_a_base(uint64_t encoded, uint64_t id);
	std::string reverse_complement(const std::string& cs);
	std::string reverse(const std::string& cs);

	uint64_t reverseComplementTwoBits(uint64_t e);
	void decodeFromTwoBitsLong(uint64_t d, char* decoded);
	void decode_from_anchor(uint64_t d, char* decoded);
	static uint32_t getLongLongsEncodedLength(uint32_t readLength);
	static uint32_t getLongsEncodedLength(uint32_t readLength);
	static uint32_t getIntsEncodedLength(uint32_t readLength);
	static uint32_t roundUpToNearestMultiple(uint32_t numberToBeRounded, uint32_t multiple);
	uint64_t insertAtFirstTwoBitsLong(uint64_t encoded, uint8_t c);
	uint64_t insertAtLastTwoBitsLong(uint64_t encoded, uint8_t c);
	uint32_t get_half_msb(uint64_t encoded);
	uint64_t get_full_msb(uint32_t encoded);
};

inline void printbits(uint64_t v) {
	int i; // for C89 compatibility
	int size = sizeof(v) * 8 - 1;
//	int size = 2;
	for (i = size; i >= 0; i--)
		std::cerr << (char) ('0' + ((v >> i) & 1));
	std::cerr << std::endl;
}

template<typename C>
inline uint64_t BinaryEncoder::encodeToTwoBitsLong(const C& cs, int64_t offset) {
	uint64_t currentEncoding = 0;
	uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG;
	int64_t limit = offset + readLength;
	uint8_t c;
	for (int64_t sequenceIndex = offset; sequenceIndex < limit; ++sequenceIndex) {
		c = cs[sequenceIndex];
		currentEncoding |= nucleotideTwoBitsEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	return currentEncoding;
}

template<typename C> int64_t BinaryEncoder::is_contain_N_in_anchor_path(const C& cs, int64_t offset) {
	for (int64_t sequenceIndex = offset + a + readLength; sequenceIndex >= offset && sequenceIndex > -1;) {
		if ('N' == cs[sequenceIndex--]) {
			return sequenceIndex;
		}
	}
	return -1;
}

template<typename C>
inline uint64_t BinaryEncoder::encode_to_anchor(uint8_t& anchor_base_id, const C& cs, int64_t offset) {
	uint64_t currentEncoding = 0;
	uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG;
	int64_t limit = offset + readLength;
	uint8_t c;
	for (int64_t sequenceIndex = offset; sequenceIndex < limit; ++sequenceIndex) {
		c = cs[sequenceIndex];
		currentEncoding |= nucleotideTwoBitsEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	anchor_base_id = cs[limit];
	const int64_t max_size = limit + a + 1;
	for (int64_t sequenceIndex = limit + 1; sequenceIndex < max_size; ++sequenceIndex) {
		c = cs[sequenceIndex];
		currentEncoding |= nucleotideTwoBitsEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	return currentEncoding;
}

template<typename C>
inline int64_t BinaryEncoder::is_contain_N_in_anchor(const C& cs, int64_t offset) {
	int64_t limit = offset + readLength;
	const int64_t max_size = limit + a;
	for (int64_t base_id = max_size; base_id > limit; --base_id) {
		if ('N' == cs[base_id]) {
			return base_id;
		}
	}
	for (int64_t base_id = limit - 1; base_id >= offset; --base_id) {
		if ('N' == cs[base_id]) {
			return base_id;
		}
	}
	return -1;
}

template<typename C>
inline uint64_t BinaryEncoder::encodeToTwoBitsLongReverse(const C& cs, int64_t offset) {
	uint64_t currentEncoding = 0;
	uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG;
	uint8_t c;
	for (int64_t sequenceIndex = offset + readLength - 1; sequenceIndex >= offset && sequenceIndex > -1;) {
		c = cs[sequenceIndex--];
		currentEncoding |= nucleotideTwoBitsReverseEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	return currentEncoding;
}

template<typename C>
inline uint64_t BinaryEncoder::encode_to_anchor_reverse(uint8_t& anchor_base_id, const C& cs, int64_t offset) {
	uint64_t currentEncoding = 0;
	uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG;
	uint8_t c;
	const int64_t kmer_min_base_id = offset + a;
	int64_t base_id = kmer_min_base_id + readLength;
	for (; base_id > kmer_min_base_id && base_id > -1;) {
		c = cs[base_id--];
		currentEncoding |= nucleotideTwoBitsReverseEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	anchor_base_id = getReverseBase(cs[base_id--]);
	for (; base_id >= offset && base_id > -1;) {
		c = cs[base_id--];
		currentEncoding |= nucleotideTwoBitsReverseEncodingLongMap[c] << shiftBits;
		shiftBits -= TWO_BITS_UNIT_LONG;
	}
	return currentEncoding;
}

template<typename C>
inline int64_t BinaryEncoder::is_contain_N_in_anchor_reverse(const C& cs, int64_t offset) {
	const int64_t kmer_min_base_id = offset + a;
	int64_t base_id = kmer_min_base_id + readLength;
	for (; base_id > kmer_min_base_id && base_id > -1;) {
		if ('N' == cs[base_id--]) {
			return base_id;
		}
	}
	--base_id;
	for (; base_id >= offset && base_id > -1;) {
		if ('N' == cs[base_id--]) {
			return base_id;
		}
	}
	return -1;
}

template<typename C>
inline uint64_t BinaryEncoder::get_minimizer_with_base_id(int64_t& base_id, const C& cs, const int32_t base_start_id, const int32_t base_end_id,
        const char target_c) {
	uint64_t minimizer = INT64_MAX;
	//	uint64_t temp_minimizer = 0;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	for (base_id = base_start_id; base_id <= max_id; ++base_id) {
		if (target_c != cs[base_id]) {
			continue;
		}
		temp_base_id = is_containing_low_complexity(cs, base_id, max_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		minimizer = std::min(minimizer, encodeToTwoBitsLong(cs, base_id));
		minimizer = std::min(minimizer, encodeToTwoBitsLongReverse(cs, base_id));
	}
	if (INT64_MAX == minimizer) {
		minimizer = 0;
	}
	return minimizer;
}
template<typename C>
inline uint64_t BinaryEncoder::get_minimizer(const C& cs, const int32_t base_start_id, const int32_t base_end_id, const char target_c) {
	uint64_t minimizer = INT64_MAX;
//	uint64_t temp_minimizer = 0;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		if (target_c != cs[base_id]) {
			continue;
		}
		temp_base_id = is_containing_low_complexity(cs, base_id, max_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		minimizer = std::min(minimizer, encodeToTwoBitsLong(cs, base_id));
		minimizer = std::min(minimizer, encodeToTwoBitsLongReverse(cs, base_id));
//		temp_minimizer = encodeToTwoBitsLong(cs, base_id);
//		if(0 == temp_minimizer) {
//			base_id += readLength;
//			continue;
//		}
//		minimizer = std::min(temp_minimizer, minimizer);
//		temp_minimizer = encodeToTwoBitsLongReverse(cs, base_id);
//		if(0 == temp_minimizer) {
//			base_id += readLength;
//			continue;
//		}
//		minimizer = std::min(temp_minimizer, minimizer);
	}
	if (INT64_MAX == minimizer) {
		minimizer = 0;
	}
	return minimizer;
}
/*

template<typename C>
inline void BinaryEncoder::get_premers(const C& cs, const int64_t base_start_id, const int64_t base_end_id, PremerBlockTuple& premers,
        const char target_c) {
	uint64_t pre_temp;
	uint64_t pre_temp_reverse;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
//	const int64_t n_shift_premer = readLength >> 2;
//	const int32_t last_base_id = readLength - 1;
	char c;
	char c_reverse;
	premers.set_to_max_values();
//	if("TCACTGGCTCATCCATCATAATCAC" == cs) {
//		std::cout << (boost::format("[BinaryEncoder.get_premers] All: %s(%d)\n") % cs % max_id).str();
//	}
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c = cs[base_id];
		c_reverse = getReverseBase(cs[base_id + readLength - 1]);
		pre_temp = encodeToTwoBitsLong(cs, base_id);
		pre_temp_reverse = encodeToTwoBitsLongReverse(cs, base_id);
		c = std::min(c, c_reverse);
		pre_temp = std::min(pre_temp, pre_temp_reverse);
//		if("TCACTGGCTCATCCATCATAATCAC" == cs) {
//			std::cout << (boost::format("[BinaryEncoder.get_premers] Inner Loop: %s\n") % cs.substr(base_id, readLength)).str();
//		}
		if (target_c == c) {
			if (premers.values[2] > pre_temp) {
				premers.values[3] = premers.values[2];
				premers.values[2] = pre_temp;
				//				premers.base_ids[1] = base_id;
			}
		}
//		switch (c) {
//		case 'A':
//			if (premers.values[0] > pre_temp) {
//				premers.values[1] = premers.values[0];
//				premers.values[0] = pre_temp;
////				premers.base_ids[0] = base_id;
//			}
//			break;
//		case 'C':
//			if (premers.values[2] > pre_temp) {
//				premers.values[3] = premers.values[2];
//				premers.values[2] = pre_temp;
////				premers.base_ids[1] = base_id;
//			}
//			break;
//		case 'G':
//			if (premers.values[4] > pre_temp) {
//				premers.values[5] = premers.values[4];
//				premers.values[4] = pre_temp;
////				premers.base_ids[4] = base_id;
//			}
//			break;
//		case 'T':
//			if (premers.values[6] > pre_temp) {
//				premers.values[7] = premers.values[6];
//				premers.values[6] = pre_temp;
////				premers.base_ids[6] = base_id;
//			}
//			break;
//		default:
//			break;
//		}
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[1] = premers.values[0];
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[3] = premers.values[2];
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		} else if ('G' == c) {
//			if (premers.values[4] > pre_temp) {
//				premers.values[5] = premers.values[4];
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[6] > pre_temp) {
//				premers.values[7] = premers.values[6];
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
//		else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}

//		if (c != get_a_base(pre_temp, 0)) {
//			char decoded[readLength + 1];
//			decodeFromTwoBitsLong(pre_temp, decoded);
//			std::cout
//			        << (boost::format("[BinaryEncoder.get_premers] %s %s %s(%c %c)\n") % cs
//			                % std::string(cs.begin() + base_id, cs.begin() + base_id + readLength) % std::string(decoded) % c
//			                % get_a_base(pre_temp, 0)).str();
//			exit(-1);
//		}
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else
//			if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
//		else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}
	}
//	int64_t base_id = max_id;
//	temp_base_id = is_containing_low_complexity(cs, base_id, max_id);
////	std::cout << cs.substr(base_id) << "\n";
//	if (-1 == temp_base_id) {
//		c = cs[base_id];
//		pre_temp = encodeToTwoBitsLong(cs, base_id);
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
////		else if ('G' == c) {
////			if (premers.values[2] > pre_temp) {
////				premers.values[2] = pre_temp;
////				premers.base_ids[2] = base_id;
////			}
////		} else if ('T' == c) {
////			if (premers.values[3] > pre_temp) {
////				premers.values[3] = pre_temp;
////				premers.base_ids[3] = base_id;
////			}
////		}
//		c = getReverseBase(cs[base_id + readLength - 1]);
//		pre_temp = encodeToTwoBitsLongReverse(cs, base_id);
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}
//	}
	premers.invalidate_max_values();
}

*/
template<typename C> uint64_t BinaryEncoder::get_premer(const C& cs, const int64_t base_start_id, const int64_t base_end_id, const char target_c) {
	uint64_t pre_temp;
	uint64_t a_premer = UINT64_MAX;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	char c;
	char c_reverse;
	const char target_c_rerverse = getReverseBase(target_c);
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c = cs[base_id];
		c_reverse = cs[base_id + readLength - 1];
//		if (target_c == c || target_c_rerverse == c_reverse) {
//			std::cout << cs.substr(base_id, readLength) << "/" << StringUtils::get_reverse_complement(cs.substr(base_id, readLength)) << "\n";
//		}
		if (target_c == c) {
			pre_temp = encodeToTwoBitsLong(cs, base_id);
			a_premer = std::min(a_premer, pre_temp);
		}
		if (target_c_rerverse == c_reverse) {
			pre_temp = encodeToTwoBitsLongReverse(cs, base_id);
			a_premer = std::min(a_premer, pre_temp);
		}
	}
	if (UINT64_MAX == a_premer) {
		a_premer = 0;
	}
	return a_premer;
}

template<typename C> uint64_t BinaryEncoder::get_premer_and_direction(const C& cs, const int64_t base_start_id, const int64_t base_end_id,
        const char target_c, char& direction) {
	uint64_t pre_temp;
	uint64_t pre_temp_reverse;
	uint64_t a_premer = UINT64_MAX;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	char c;
	char c_reverse;
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c = cs[base_id];
		c_reverse = getReverseBase(cs[base_id + readLength - 1]);
		if (target_c == c) {
			pre_temp = encodeToTwoBitsLong(cs, base_id);
			if(a_premer > pre_temp) {
				a_premer = pre_temp;
				direction = 'F';
			}
		} else if (target_c == c_reverse) {
			pre_temp_reverse = encodeToTwoBitsLongReverse(cs, base_id);
			if(a_premer > pre_temp_reverse) {
				a_premer = pre_temp_reverse;
				direction = 'R';
			}
		}
	}
	if (UINT64_MAX == a_premer) {
		direction = 'F';
		a_premer = 0;
	}
	return a_premer;
}

template<typename C> uint64_t BinaryEncoder::get_premer_forward(const C& cs, const int64_t base_start_id, const int64_t base_end_id,
        const char target_c) {
	uint64_t pre_temp;
	uint64_t a_premer = UINT64_MAX;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	char c;
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c = cs[base_id];
		if (target_c == c) {
			pre_temp = encodeToTwoBitsLong(cs, base_id);
			a_premer = std::min(a_premer, pre_temp);
		}
	}
	if (UINT64_MAX == a_premer) {
		a_premer = 0;
	}
	return a_premer;
}

template<typename C> uint64_t BinaryEncoder::get_premer_reverse(const C& cs, const int64_t base_start_id, const int64_t base_end_id,
        const char target_c) {
	uint64_t pre_temp_reverse;
	uint64_t a_premer = UINT64_MAX;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
	char c_reverse;
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c_reverse = getReverseBase(cs[base_id + readLength - 1]);
		if (target_c == c_reverse) {
			pre_temp_reverse = encodeToTwoBitsLongReverse(cs, base_id);
			a_premer = std::min(a_premer, pre_temp_reverse);
		}
	}
	if (UINT64_MAX == a_premer) {
		a_premer = 0;
	}
	return a_premer;
}
/*

template<typename C>
inline void BinaryEncoder::get_all_premers(const C& cs, const int64_t base_start_id, const int64_t base_end_id, PremerBlockTuple& premers) {
	uint64_t pre_temp;
//	uint64_t pre_temp_reverse;
	int64_t temp_base_id = 0;
	const int64_t max_id = base_end_id - readLength;
//	const int64_t n_shift_premer = readLength >> 2;
//	const int32_t last_base_id = readLength - 1;
	char c;
//	char c_reverse;
	premers.set_to_max_values();
//	if("TCACTGGCTCATCCATCATAATCAC" == cs) {
//		std::cout << (boost::format("[BinaryEncoder.get_premers] All: %s(%d)\n") % cs % max_id).str();
//	}
	for (int64_t base_id = base_start_id; base_id <= max_id; ++base_id) {
		temp_base_id = is_containing_low_complexity(cs, base_id, base_end_id);
		if (-1 != temp_base_id) {
			base_id = temp_base_id;
			continue;
		}
		c = cs[base_id];

		pre_temp = encodeToTwoBitsLong(cs, base_id);
//		if("TCACTGGCTCATCCATCATAATCAC" == cs) {
//			std::cout << (boost::format("[BinaryEncoder.get_premers] Inner Loop: %s\n") % cs.substr(base_id, readLength)).str();
//		}
		switch (c) {
		case 'A':
			if (premers.values[0] > pre_temp) {
				premers.values[1] = premers.values[0];
				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
			}
			break;
		case 'C':
			if (premers.values[2] > pre_temp) {
				premers.values[3] = premers.values[2];
				premers.values[2] = pre_temp;
//				premers.base_ids[1] = base_id;
			}
			break;
		case 'G':
			if (premers.values[4] > pre_temp) {
				premers.values[5] = premers.values[4];
				premers.values[4] = pre_temp;
//				premers.base_ids[4] = base_id;
			}
			break;
		case 'T':
			if (premers.values[6] > pre_temp) {
				premers.values[7] = premers.values[6];
				premers.values[6] = pre_temp;
//				premers.base_ids[6] = base_id;
			}
			break;
		default:
			break;
		}
		c = getReverseBase(cs[base_id + readLength - 1]);
		pre_temp = encodeToTwoBitsLongReverse(cs, base_id);
		switch (c) {
		case 'A':
			if (premers.values[0] > pre_temp) {
				premers.values[1] = premers.values[0];
				premers.values[0] = pre_temp;
				//				premers.base_ids[0] = base_id;
			}
			break;
		case 'C':
			if (premers.values[2] > pre_temp) {
				premers.values[3] = premers.values[2];
				premers.values[2] = pre_temp;
				//				premers.base_ids[1] = base_id;
			}
			break;
		case 'G':
			if (premers.values[4] > pre_temp) {
				premers.values[5] = premers.values[4];
				premers.values[4] = pre_temp;
				//				premers.base_ids[4] = base_id;
			}
			break;
		case 'T':
			if (premers.values[6] > pre_temp) {
				premers.values[7] = premers.values[6];
				premers.values[6] = pre_temp;
				//				premers.base_ids[6] = base_id;
			}
			break;
		default:
			break;
		}
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[1] = premers.values[0];
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[3] = premers.values[2];
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		} else if ('G' == c) {
//			if (premers.values[4] > pre_temp) {
//				premers.values[5] = premers.values[4];
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[6] > pre_temp) {
//				premers.values[7] = premers.values[6];
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
//		else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}

//		if (c != get_a_base(pre_temp, 0)) {
//			char decoded[readLength + 1];
//			decodeFromTwoBitsLong(pre_temp, decoded);
//			std::cout
//			        << (boost::format("[BinaryEncoder.get_premers] %s %s %s(%c %c)\n") % cs
//			                % std::string(cs.begin() + base_id, cs.begin() + base_id + readLength) % std::string(decoded) % c
//			                % get_a_base(pre_temp, 0)).str();
//			exit(-1);
//		}
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else
//			if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
//		else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}
	}
//	int64_t base_id = max_id;
//	temp_base_id = is_containing_low_complexity(cs, base_id, max_id);
////	std::cout << cs.substr(base_id) << "\n";
//	if (-1 == temp_base_id) {
//		c = cs[base_id];
//		pre_temp = encodeToTwoBitsLong(cs, base_id);
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}
////		else if ('G' == c) {
////			if (premers.values[2] > pre_temp) {
////				premers.values[2] = pre_temp;
////				premers.base_ids[2] = base_id;
////			}
////		} else if ('T' == c) {
////			if (premers.values[3] > pre_temp) {
////				premers.values[3] = pre_temp;
////				premers.base_ids[3] = base_id;
////			}
////		}
//		c = getReverseBase(cs[base_id + readLength - 1]);
//		pre_temp = encodeToTwoBitsLongReverse(cs, base_id);
//		if ('A' == c) {
//			if (premers.values[0] > pre_temp) {
//				premers.values[0] = pre_temp;
//				premers.base_ids[0] = base_id;
//			}
//		} else if ('C' == c) {
//			if (premers.values[1] > pre_temp) {
//				premers.values[1] = pre_temp;
//				premers.base_ids[1] = base_id;
//			}
//		}else if ('G' == c) {
//			if (premers.values[2] > pre_temp) {
//				premers.values[2] = pre_temp;
//				premers.base_ids[2] = base_id;
//			}
//		} else if ('T' == c) {
//			if (premers.values[3] > pre_temp) {
//				premers.values[3] = pre_temp;
//				premers.base_ids[3] = base_id;
//			}
//		}
//	}
	premers.invalidate_max_values();
}

*/
template<typename C>
inline int64_t BinaryEncoder::is_containing_low_complexity(const C& cs, int64_t base_start_id, int64_t max_id) {
	if(base_start_id + readLength > max_id) {
		return max_id;
	}
//	std::cout << "is_containing " << cs.substr(base_start_id) << "\n";
	max_id = base_start_id + readLength - 1;
//	if("TCACTGGCTCATCCATCATAATCAC" == cs) {
//		std::cout << (boost::format("[BinaryEncoder.is_containing_low_complexity] %s, %d-%d\n") % cs % base_start_id % max_id).str();
//	}
//	if (int64_t(readLength) > max_id - base_start_id) {
//		return max_id;
//	}
	for (int64_t base_id = max_id; base_id >= base_start_id; --base_id) {
		if ('N' == cs[base_id]) {
			return base_id;
		}
	}
	// we check the initial bases that contains the same characters in a length (k >> 2)
//	const int64_t max_start_run_id = std::min(base_start_id + (readLength >> 2), max_id);
////	char start_c = cs[max_start_run_id - 1];
//	uint64_t n_same = 0;
//	for (int64_t base_id = base_start_id; base_id < max_start_run_id; ++base_id) {
//		if (cs[base_id] == cs[base_id + 1]) {
//			++n_same;
//		}
//	}
//	if(n_same == max_start_run_id) {
////		std::cout << (boost::format("[BinaryEncoder.is_containing_low_complexity] encountered: %s(%d-%d)\n") % cs.substr(base_start_id, readLength) % base_start_id % max_start_run_id).str();
//		return max_start_run_id;
//	}
	for (int64_t base_id = max_id; base_id > base_start_id + 2; --base_id) {
		if (cs[base_id - 1] == cs[base_id] && cs[base_id - 2] == cs[base_id - 1] && cs[base_id - 3] == cs[base_id - 2]) {
			return base_id;
		}
	}
	return -1;
}

inline char BinaryEncoder::get_a_base(uint64_t encoded, uint64_t id) {
	uint64_t n_shifts = TWO_BITS_INITIAL_SHIFT_BITS_LONG - (id << 1);
	return nucleotideTwoBitsDecodingMap[(encoded >> n_shifts) & TWO_BITS_MASK_LONG];
}

inline void BinaryEncoder::decodeFromTwoBitsLong(uint64_t d, char* decoded) {
	uint32_t processed = 0;
	for (uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG; processed < readLength && shiftBits >= 0; shiftBits -= TWO_BITS_UNIT_LONG) {
		decoded[processed++] = nucleotideTwoBitsDecodingMap[(int) ((d & (TWO_BITS_MASK_LONG << shiftBits)) >> shiftBits)];
	}
	decoded[readLength] = 0;
}

inline void BinaryEncoder::decode_from_anchor(uint64_t d, char* decoded) {
	uint32_t processed = 0;
	for (uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG; processed < readLength + a && shiftBits >= 0; shiftBits -= TWO_BITS_UNIT_LONG) {
		decoded[processed++] = nucleotideTwoBitsDecodingMap[(int) ((d & (TWO_BITS_MASK_LONG << shiftBits)) >> shiftBits)];
	}
	decoded[readLength + a] = 0;
}

inline uint64_t BinaryEncoder::insertAtFirstTwoBitsLong(uint64_t encoded, uint8_t c) {
	encoded >>= TWO_BITS_UNIT_LONG;
	encoded |= nucleotideTwoBitsEncodingLongMap[c] << TWO_BITS_INITIAL_SHIFT_BITS_LONG;
	return encoded;
}

inline uint64_t BinaryEncoder::insertAtLastTwoBitsLong(uint64_t encoded, uint8_t c) {
	encoded &= TWO_BITS_INITIAL_CLEAR_MASKS_LONG;
	encoded <<= TWO_BITS_UNIT_LONG;
	uint64_t shiftBits = TWO_BITS_INITIAL_SHIFT_BITS_LONG - (TWO_BITS_UNIT_LONG * (readLength - 1));
	encoded |= nucleotideTwoBitsEncodingLongMap[c] << shiftBits;
	return encoded;
}

inline uint32_t BinaryEncoder::get_half_msb(uint64_t encoded) {
	return encoded >> (HALF_MSB_BITS + additional_shifts_for_half_msb);
}

inline uint64_t BinaryEncoder::get_full_msb(uint32_t encoded) {
	uint64_t encoded_temp = encoded;
	return encoded_temp << (HALF_MSB_BITS + additional_shifts_for_half_msb);
}

inline uint64_t BinaryEncoder::reverseComplementTwoBits(uint64_t e) {
	// swap the ranged bit(2bit) contrasting the center
	// i: start positions of bit sequences to swap
	// j: target positions of bit sequences to swap
	// x: XOR temporary
	uint64_t i, j, x;

	for (uint32_t pos = 0; pos < (readLength >> 1); ++pos) {
		i = twoBitsreverseComplementFirstLookUp[pos];
		j = twoBitsreverseComplementSecondLookUp[i];
		x = ((e >> i) ^ (e >> j)) & THREE_LONG;
		e = e ^ ((x << i) | (x << j));
	}
	e = ~e & lowerBitClearMask;
	return e;
}
} /* namespace castle */
#endif /* BINARYENCODER_HPP_ */
