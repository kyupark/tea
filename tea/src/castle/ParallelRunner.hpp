/*
 * ParallelRunner.hpp
 *
 *  Created on: Feb 20, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef PARALLELRUNNER_HPP_
#define PARALLELRUNNER_HPP_

#include <vector>
#include <algorithm>
#include <functional>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "../meerkat/BlockBoundary.hpp"
#include "StringUtils.hpp"
#include "TimeChecker.hpp"
#include "concurrent_queue.h"

namespace castle {
using namespace std;
using namespace moodycamel;
struct ConcurrentQueueTrait: public ConcurrentQueueDefaultTraits {
	static const size_t BLOCKS_SIZE = 256;
};

class ParallelRunner {

public:
	ParallelRunner();
	~ParallelRunner();
	static function<void(size_t, size_t)> empty_ranged_func;
	static void run_tasks(vector<function<void()> >& a_tasks, function<void(size_t, size_t)>& block_callback);
	static void run_step_wise(const size_t n_cores, vector<function<void()> >& a_tasks, function<void(size_t, size_t)>& block_callback);
	static void run_unbalanced_load(const uint32_t n_cores, vector<function<void()> >& a_tasks);
	static int64_t get_next_start_pos(uint32_t id, int64_t chunk_size, int32_t delta, int64_t already_processed_size);
	static int64_t get_next_end_pos(uint32_t id, uint32_t max_id, int64_t chunk_size, int64_t max_size, int64_t already_processed_size);
	// a bfi index contains the offsets for the BAM sorted by read name
	static void create_bfi_index(vector<int64_t>& block_boundary, const string& input_path, const string& output_index_path, const int64_t block_size);
	static void load_bfi_index(vector<int64_t>& block_boundary, const string& input_path);
	static void create_bni_index(const string& input_path, const string& output_index_path);
	static void load_bni_index(vector<BamTools::BlockOffset>& offset_blocks, const string& input_BAM_name, const string& output_bni_index_path);
	static void load_bai_bni_index(vector<int64_t>& block_boundary, const string& input_path, const string& output_bni_index_path, const int64_t block_size, const int64_t n_cores);
	template<typename T> static T cas(volatile T *ptr, T oval, T nval);
	template<typename T> static T aaf(volatile T *ptr, T amount);
};

template<typename T>
inline T ParallelRunner::cas(volatile T *ptr, T oval, T nval) {
	return __sync_val_compare_and_swap(ptr, oval, nval);
}

template<typename T>
inline T ParallelRunner::aaf(volatile T* ptr, T amount) {
	return __sync_add_and_fetch(ptr, amount);
}

} /* namespace castle */
#endif /* PARALLELRUNNER_HPP_ */
