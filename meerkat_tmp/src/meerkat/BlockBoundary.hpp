/*
 * BlockBoundary.hpp
 *
 *  Created on: Jun 8, 2016
 *      Author: el174
 */

#ifndef MEERKAT_BLOCKBOUNDARY_HPP_
#define MEERKAT_BLOCKBOUNDARY_HPP_
#include <string>
#include <boost/lexical_cast.hpp>

namespace meerkat {

	using namespace std;
	struct BlockBoundary {
			string read_name;
			int32_t ref_id;
			int64_t offset;
			int32_t pos;
			uint32_t aln_flag;
			int64_t jump_pos;
			BlockBoundary() :
					ref_id(0), offset(0), pos(0), aln_flag(0), jump_pos(0) {
			}
			BlockBoundary(const BlockBoundary& b) {
				read_name = b.read_name;
				ref_id = b.ref_id;
				offset = b.offset;
				pos = b.pos;
				aln_flag = b.aln_flag;
				jump_pos = b.jump_pos;
			}
			BlockBoundary& operator=(const BlockBoundary& b) {
				// check for self-assignment
				if (this == &b) {
					return *this;
				}
				read_name = b.read_name;
				ref_id = b.ref_id;
				offset = b.offset;
				pos = b.pos;
				aln_flag = b.aln_flag;
				jump_pos = b.jump_pos;
				return *this;
			}
			bool operator<(const BlockBoundary& b) const {
				return offset < b.offset;
//				if (ref_id < b.ref_id) {
//					return true;
//				} else if (ref_id > b.ref_id) {
//					return false;
//				}
//				if (pos < b.pos) {
//					return true;
//				} else if (pos > b.pos) {
//					return false;
//				}
////				if (jump_pos < b.jump_pos) {
////					return true;
////				} else if (jump_pos > b.jump_pos) {
////					return false;
////				}
//				if (aln_flag < b.aln_flag) {
//					return true;
//				} else if (aln_flag > b.aln_flag) {
//					return false;
//				}
//				return read_name < b.read_name;
			}
			bool operator==(const BlockBoundary& b) const {
				if (this == &b) {
					return true;
				}
				return offset == b.offset;
//				return ref_id == b.ref_id && pos == b.pos
//						&& aln_flag == b.aln_flag && read_name == b.read_name;
			}
			bool operator!=(const BlockBoundary& b) const {
				return !(*this == b);
			}
			string str() const {
				string result(read_name);
				result += "\t" + boost::lexical_cast<string>(ref_id) + "\t"
						+ boost::lexical_cast<string>(pos) + "\t"
						+ boost::lexical_cast<string>(aln_flag);
				return result;
			}
	};

}

#endif /* MEERKAT_BLOCKBOUNDARY_HPP_ */
