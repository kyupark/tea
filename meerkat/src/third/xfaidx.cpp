/*
 * xfaidx.cpp
 *
 *  Created on: Oct 10, 2011
 *      Author: lindenb
 */
#include <cstdio>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include "faidx.h"
#include "xfaidx.h"

using namespace std;

IndexedFasta::IndexedFasta(const char* fasta) {
	faidx_t *idx = fai_load(fasta);
	if (idx == NULL) {
		std::cerr << "Cannot load indexed fasta " << fasta << "(" << strerror(errno) << ")";
		exit(1);
	}
	ptr = idx;
}

IndexedFasta::~IndexedFasta() {
	if (ptr != NULL) {
		fai_destroy(static_cast<faidx_t*>(ptr));
	}
}

int32_t IndexedFasta::size() {
	return ::faidx_fetch_nseq(static_cast<faidx_t*>(ptr));
}

int IndexedFasta::build(std::string fa_name) {
	return ::fai_build(fa_name.c_str());
}
std::string IndexedFasta::fetch(const std::string& chrom, int32_t start0, int32_t end0) {
	if(start0 > end0) {
		return castle::StringUtils::get_reverse_complement(_fetch(chrom, end0, start0));
	}
	return _fetch(chrom, start0, end0);
}

std::string IndexedFasta::fetch_1_system(const std::string& chrom, int32_t start0, int32_t end0) {
	if(start0 > end0) {
		return castle::StringUtils::get_reverse_complement(_fetch(chrom, end0 - 1, start0 - 1));
	}
	return _fetch(chrom, start0 - 1, end0 - 1);
}

std::string IndexedFasta::_fetch(const std::string& chrom, int32_t start0, int32_t end0) {
	int len;
	char *s = ::faidx_fetch_seq(static_cast<faidx_t*>(ptr), (char *)chrom.c_str(), start0, end0, &len);
	if (s == NULL) {
		return std::string();
	}
	std::string result(s, len);
	free(s);
	return result;
}

