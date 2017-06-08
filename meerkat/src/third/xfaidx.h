/*
 * xfaidx.cpp
 *
 *  Created on: Oct 10, 2011
 *      Author: lindenb
 */
#ifndef XFAIDX_H
#define XFAIDX_H
#include <string>
#include <memory>
#include <iostream>
#include <stdint.h>
#include "../castle/StringUtils.hpp"

class IndexedFasta {
private:
	void* ptr;
public:
	IndexedFasta(const char* fasta);
	~IndexedFasta();
	int build(std::string fa_name);
	std::string fetch(const std::string& chrom, int32_t start0, int32_t end0);
	std::string fetch_1_system(const std::string& chrom, int32_t start0, int32_t end0);
	int32_t size();
private:
	std::string _fetch(const std::string& chrom, int32_t start0, int32_t end0);
};
#endif
