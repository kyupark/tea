/*
 * BinaryEncoder.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: euncheonlim
 */

#include "BinaryEncoder.hpp"
//#include <byteswap.h>

using namespace std;

namespace castle {

uint64_t BinaryEncoder::nucleotideTwoBitsEncodingLongMap[] = {
// 0 - 15
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 16 - 31
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 32 - 07
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /*'-'*/, 0, 0,
        // 08 - 63
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 60 - 79
        0, 0 /* 'A' */, 0, 1 /* 'C' */, 0, 0, 0, 2 /* 'G' */, 0, 0, 0, 0, 0, 0, 0, 0,
        // 80 - 95
        0, 0, 0, 0, 3 /* 'T' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 96 - 111
        0, 0 /* 'a' */, 0, 1 /* 'c' */, 0, 0, 0, 2 /* 'g' */, 0, 0, 0, 0, 0, 0, 0, 0,
        // 112 - 127
        0, 0, 0, 0, 3 /* 't' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 128 - 103
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 100 - 159
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 160 - 175
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 176 - 191
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 192 - 207
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 208 - 223
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 220 - 200
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 201 - 255
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

uint64_t BinaryEncoder::nucleotideTwoBitsReverseEncodingLongMap[] = {
// 0 - 15
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 16 - 31
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 32 - 07
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /*'-'*/, 3, 3,
        // 08 - 63
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 60 - 79
        3, 3 /* 'A' */, 3, 2 /* 'C' */, 3, 3, 3, 1 /* 'G' */, 3, 3, 3, 3, 3, 3, 3, 3,
        // 80 - 95
        3, 3, 3, 3, 0 /* 'T' */, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 96 - 111
        3, 3 /* 'a' */, 3, 2 /* 'c' */, 3, 3, 3, 1 /* 'g' */, 3, 3, 3, 3, 3, 3, 3, 3,
        // 112 - 127
        3, 3, 3, 3, 0 /* 't' */, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 128 - 103
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 100 - 159
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 160 - 175
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 176 - 191
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 192 - 207
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 208 - 223
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 220 - 200
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        // 201 - 255
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

uint64_t BinaryEncoder::nucleotideThreeBitsEncodingLongMap[] = {
// 0 - 15
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 16 - 31
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 32 - 07
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /*'-'*/, 0, 0,
        // 08 - 63
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 60 - 79
        0, 1 /* 'A' */, 0, 2 /* 'C' */, 0, 0, 0, 3 /* 'G' */, 0, 0, 0, 0, 0, 0, 0, 0,
        // 80 - 95
        0, 0, 0, 0, 4 /* 'T' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 96 - 111
        0, 1 /* 'a' */, 0, 2 /* 'c' */, 0, 0, 0, 3 /* 'g' */, 0, 0, 0, 0, 0, 0, 0, 0,
        // 112 - 127
        0, 0, 0, 0, 4 /* 't' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 128 - 103
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 100 - 159
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 160 - 175
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 176 - 191
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 192 - 207
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 208 - 223
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 220 - 200
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // 201 - 255
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
char BinaryEncoder::nucleotideTwoBitsDecodingMap[] = { 'A', 'C', 'G', 'T' };
char BinaryEncoder::nucleotideTwoBitsReverseDecodingMap[] = { 'A', 'C', 'G', 'T' };
char BinaryEncoder::nucleotideThreeBitsDecodingMap[] = { 'N', 'A', 'C', 'G', 'T' };
uint64_t BinaryEncoder::twoBitsreverseComplementFirstLookUp[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
uint64_t BinaryEncoder::twoBitsreverseComplementSecondLookUp[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

uint32_t BinaryEncoder::roundUpToNearestMultiple(uint32_t numberToBeRounded, uint32_t multiple) {
	if (0 == multiple)
		return numberToBeRounded;
	int remainder = numberToBeRounded % multiple;
	if (0 == remainder)
		return numberToBeRounded;
	return numberToBeRounded + multiple - remainder;
}

uint32_t BinaryEncoder::getLongLongsEncodedLength(uint32_t readLength) {
	return roundUpToNearestMultiple(readLength, MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONGLONG) / MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONGLONG;
}

uint32_t BinaryEncoder::getLongsEncodedLength(uint32_t readLength) {
	return roundUpToNearestMultiple(readLength, MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONG) / MAXIMUM_LENGTH_OF_ENCODE_VECTOR_LONG;
}

uint32_t BinaryEncoder::getIntsEncodedLength(uint32_t readLength) {
	return roundUpToNearestMultiple(readLength, MAXIMUM_LENGTH_OF_ENCODE_VECTOR_INT) / MAXIMUM_LENGTH_OF_ENCODE_VECTOR_INT;
}

BinaryEncoder::BinaryEncoder() :
		readLength(0), a(0), lowerBitClearMask(MAX_LONG), sizeOfEncodeVector(0), additional_shifts_for_half_msb(0), n_shift_premer(0) {
}

BinaryEncoder::~BinaryEncoder() {
}

void BinaryEncoder::setReadLength(uint32_t length) {
	readLength = length;
	lowerBitClearMask = MAX_LONG << (TWO_BITS_INITIAL_SHIFT_BITS_LONG - ((readLength - ONE_LONG) << ONE_LONG));
	uint64_t i, j;
	for (uint32_t pos = 0; pos < (readLength >> 1); ++pos) {
		i = (TWO_BITS_INITIAL_SHIFT_BITS_LONG - (pos << ONE_LONG));
		j = TWO_BITS_INITIAL_SHIFT_BITS_LONG - ((readLength - pos - 1) << ONE_LONG);
		twoBitsreverseComplementFirstLookUp[pos] = i;
		twoBitsreverseComplementSecondLookUp[i] = j;
	}
	additional_shifts_for_half_msb = (16 - length) << 1;
	n_shift_premer = readLength >> 2;
}

void BinaryEncoder::set_a(uint32_t a_a) {
	a = a_a;
	std::cout << (boost::format("[BinaryEncoder.set_a] Anchor size: %d\n") % a).str();
}

uint32_t BinaryEncoder::getReadLength() {
	return readLength;
}

void BinaryEncoder::setSizeOfEncodeVector(uint32_t size) {
	sizeOfEncodeVector = size;
}

uint32_t BinaryEncoder::getSizeOfEncodeVector() {
	return sizeOfEncodeVector;
}

string BinaryEncoder::reverse_complement(const string& cs) {
	string copied(cs.rbegin(), cs.rend());
	char c = '\0';
	for(uint32_t base_id = 0; base_id < copied.size(); ++base_id) {
		c = copied[base_id];
		copied[base_id] = getReverseBase(c);
	}
	return copied;
}
string BinaryEncoder::reverse(const string& cs) {
	string copied(cs.rbegin(), cs.rend());
	return copied;
}

} /* namespace castle */
