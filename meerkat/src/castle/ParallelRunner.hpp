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

inline void ParallelRunner::run_tasks(vector<function<void()> >& a_tasks, function<void(size_t, size_t)>& block_callback) {
	boost::thread_group group;
	for (uint32_t start = 0; start < a_tasks.size(); ++start) {
		group.create_thread([&a_tasks, start]() {
			a_tasks[start]();
		});
	}
	group.join_all();
	block_callback(0, a_tasks.size());
	a_tasks.clear();
}

inline void ParallelRunner::run_step_wise(const size_t n_cores, vector<function<void()> >& a_tasks,
        function<void(size_t, size_t)>& block_callback) {

//	const size_t steps = (a_tasks.size() + n_cores - 1) / n_cores;

//	uint64_t block_id = 0;
	for (size_t start = 0; start < a_tasks.size(); start += n_cores) {
		boost::thread_group group;
		for (size_t i = start; i < start + n_cores && i < a_tasks.size(); ++i) {
			group.create_thread([&a_tasks, i]() {
//				const size_t end = min(a_tasks.size(), start + steps);
				    a_tasks[i]();
			    });
		}
		group.join_all();
		block_callback(start, std::min(start + n_cores, a_tasks.size()));
//		cout << (boost::format("[ParallelRunner.run_step_wise] Block: %d\n") % block_id++).str();
	}
	a_tasks.clear();
}

inline void ParallelRunner::run_unbalanced_load(const uint32_t n_cores, vector<function<void()> >& a_tasks) {
	ConcurrentQueue<function<void() >, ConcurrentQueueTrait> tasks_queue;
	tasks_queue.enqueue_bulk(a_tasks.begin(), a_tasks.size());
	const uint64_t n_tasks = a_tasks.size();
	if(n_tasks < n_cores) {
		boost::atomic_uint64_t done_consumers(0);
		boost::thread threads[n_tasks];
		for (uint32_t core_id = 0; core_id != n_tasks; ++core_id) {
			threads[core_id] = boost::thread([&]() {
				function<void() > item;
				do {
					while (tasks_queue.try_dequeue(item)) {
						item();
					}
				} while (done_consumers.fetch_add(1, boost::memory_order_acq_rel) + 1 == n_tasks);
			});
		}
		for (uint32_t core_id = 0; core_id != n_tasks; ++core_id) {
			threads[core_id].join();
		}
	} else {
		boost::atomic_uint64_t done_consumers(0);
		boost::thread threads[n_cores];
		for (uint32_t core_id = 0; core_id != n_cores; ++core_id) {
			threads[core_id] = boost::thread([&]() {
				function<void() > item;
				do {
					while (tasks_queue.try_dequeue(item)) {
						item();
					}
				} while (done_consumers.fetch_add(1, boost::memory_order_acq_rel) + 1 == n_tasks);
			});
		}
		for (uint32_t core_id = 0; core_id != n_cores; ++core_id) {
			threads[core_id].join();
		}
	}
	a_tasks.clear();
}
inline int64_t ParallelRunner::get_next_start_pos(uint32_t id, int64_t chunk_size, int32_t delta, int64_t already_processed_size) {
	return already_processed_size + max<int64_t>(0, id * chunk_size - delta);
}
inline int64_t ParallelRunner::get_next_end_pos(uint32_t id, uint32_t max_id, int64_t chunk_size, int64_t max_size, int64_t already_processed_size) {
	if (id == max_id)
		return max_size;
	return already_processed_size + (id + 1) * chunk_size;
}
} /* namespace castle */
#endif /* PARALLELRUNNER_HPP_ */
