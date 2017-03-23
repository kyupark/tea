/*
 * Hashs.hpp
 *
 *  Created on: Jun 14, 2016
 *      Author: el174
 */

#ifndef CASTLE_HASHS_HPP_
#define CASTLE_HASHS_HPP_

#include "MurmurHash3.hpp"
#include <unordered_map>
#include <unordered_set>
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>
#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace castle {
	struct UInt64Hash {
			inline uint32_t operator()(const uint64_t& v) const {
				uint32_t h;
				uint64_t* value = const_cast<uint64_t*>(&v);
				MurmurHash3_x86_32(value, sizeof(uint64_t), 0, &h);
				return h;
			}
	};
	struct Int64Hash {
			inline uint32_t operator()(const int64_t& v) const {
				uint32_t h;
				int64_t* value = const_cast<int64_t*>(&v);
				MurmurHash3_x86_32(value, sizeof(int64_t), 0, &h);
				return h;
			}
	};
	struct UInt32Hash {
			inline uint32_t operator()(const uint32_t& v) const {
				uint32_t h;
				uint32_t* value = const_cast<uint32_t*>(&v);
				MurmurHash3_x86_32(value, sizeof(uint32_t), 0, &h);
				return h;
			}
	};
	struct Int32Hash {
			inline uint32_t operator()(const int32_t& v) const {
				int32_t h;
				int32_t* value = const_cast<int32_t*>(&v);
				MurmurHash3_x86_32(value, sizeof(int32_t), 0, &h);
				return h;
			}
	};
	struct StringHash {
			inline uint32_t operator()(const std::string& v) const {
				uint32_t h;
				MurmurHash3_x86_32(&v[0], v.size(), 0, &h);
				return h;
			}
	};
}
;

#endif /* CASTLE_HASHS_HPP_ */
