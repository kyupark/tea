/*
 * AlternativeMap.hpp
 *
 *  Created on: Jun 21, 2016
 *      Author: el174
 */

#ifndef MEERKAT_ALTERNATIVEMAP_HPP_
#define MEERKAT_ALTERNATIVEMAP_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "../../castle/StringUtils.hpp"

namespace meerkat {
using namespace std;
struct AlternativeMap {
	string read_name;
	string ref_id;
	int64_t read_start;
	int32_t read_strand;
	string mate_ref_id;
	int64_t mate_start;
	int32_t mate_strand;
	int64_t insert_size;
	int64_t nm;
	int64_t mate_nm;
	string read_group;
	int64_t aln_len;
	int64_t mate_aln_len;
	double weight;
	AlternativeMap() :
			read_start(-1), read_strand(1), mate_start(-1), mate_strand(1), insert_size(
					numeric_limits<int64_t>::max()), nm(-1), mate_nm(-1), aln_len(-1), mate_aln_len(-1), weight(
					1) {
	}

	bool read_from_string(const string& a_str, vector<string>& data) {
		castle::StringUtils::tokenize(a_str, "\t", data);
		if (data.size() < 14) {
			return false;
		}
		read_name = data[0];
		ref_id = data[1];
		read_start = boost::lexical_cast<int64_t>(data[2]);
		read_strand = boost::lexical_cast<int32_t>(data[3]);
		mate_ref_id = data[4];
		mate_start = boost::lexical_cast<int64_t>(data[5]);
		mate_strand = boost::lexical_cast<int32_t>(data[6]);
		if (data[7].empty()) {
			insert_size = -1;
		} else {
			insert_size = boost::lexical_cast<int64_t>(data[7]);
		}
		nm = boost::lexical_cast<int64_t>(data[8]);
		mate_nm = boost::lexical_cast<int64_t>(data[9]);
		read_group = data[10];
		aln_len = boost::lexical_cast<int64_t>(data[11]);
		mate_aln_len = boost::lexical_cast<int64_t>(data[12]);
		weight = boost::lexical_cast<double>(data[13]);
		return true;
	}
	bool operator ==(const AlternativeMap& other) const {
		return mate_ref_id == other.mate_ref_id
				&& read_start == other.read_start
				&& mate_start == other.mate_start
				&& read_strand == other.read_strand
				&& mate_strand == other.mate_strand
				&& ref_id == other.ref_id
				&& insert_size == other.insert_size
				&& nm == other.nm
				&& mate_nm == other.nm
				&& aln_len == other.aln_len
				&& mate_aln_len == other.mate_aln_len
				&& weight == other.weight
				&& read_group == other.read_group
				&& read_name == other.read_name;
	}
	bool operator !=(const AlternativeMap& other) const {
		return !(*this == other);
	}

	bool operator <(const AlternativeMap& other) const {
		if (ref_id > other.ref_id) {
			return false;
		} else if (ref_id < other.ref_id) {
			return true;
		}
		// lhs.ref_id == rhs.ref_id
		if (mate_ref_id > other.mate_ref_id) {
			return false;
		} else if (mate_ref_id < other.mate_ref_id) {
			return true;
		}
		// lhs.mate_ref_id == rhs.mate_ref_id

		if (read_strand > other.read_strand) {
			return false;
		} else if (read_strand < other.read_strand) {
			return true;
		}
		// lhs.read_strand == rhs.read_strand

		if (mate_strand > other.mate_strand) {
			return false;
		} else if (mate_strand < other.mate_strand) {
			return true;
		}
		// lhs.mate_strand == rhs.mate_strand

		if (read_start > other.read_start) {
			return false;
		} else if (read_start < other.read_start) {
			return true;
		}
		// lhs.read_start == rhs.read_start

		if (mate_start > other.mate_start) {
			return false;
		} else if (mate_start < other.mate_start) {
			return true;
		}
		// lhs.mate_start == rhs.mate_start

		if (insert_size > other.insert_size) {
			return false;
		} else if (insert_size < other.insert_size) {
			return true;
		}
		// lhs.mate_start == rhs.mate_start
		if (read_group > other.read_group) {
			return false;
		} else if (read_group < other.read_group) {
			return true;
		}
		return false;
	}
	string str() const {
		return (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s")
				% read_name % ref_id % read_start % read_strand % mate_ref_id
				% mate_start % mate_strand % insert_size % nm % mate_nm
				% read_group % aln_len % mate_aln_len % weight).str();
	}
	string str_cl() const {
		return (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s")
				% read_name % read_group % ref_id % read_strand % read_start
				% aln_len % mate_ref_id % mate_strand % mate_start
				% mate_aln_len % insert_size % nm % mate_nm % weight).str();
	}
};

}
;

#endif /* MEERKAT_ALTERNATIVEMAP_HPP_ */
