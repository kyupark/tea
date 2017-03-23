/*
 * BWACaller.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: el174
 */

#ifndef MEERKAT_BWACALLER_HPP_
#define MEERKAT_BWACALLER_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <thread>
#include "../castle/TimeChecker.hpp"
#include "../castle/IOUtils.hpp"
#include "../castle/StringUtils.hpp"
#include "../castle/ParallelRunner.hpp"
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

namespace meerkat {
using namespace std;
class BWACaller {
public:
	BWACaller();
	~BWACaller();
	void set_n_cores(const int32_t n_cores);
	void align_single_reads(const string& sam_path, const string& reference_path, const string& aln_param, const string& samse_param, const string& file_path);
	int64_t split_FASTQ(const string& file_path_1, const string& file_path_2, int64_t n_lines = 262144);
	int64_t split_FASTQ_alt(const string& file_path_1, const string& file_path_2, int64_t n_lines = 262144);
	int64_t split_FASTQ_single(const string& file_path_1, int64_t n_lines = 262144);
	void align_reads(const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2);
	void collect_align_tasks(vector<function<void()> >& tasks, const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2);
	void collect_align_tasks_alt(vector<function<void()> >& tasks, const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2, const int64_t n_blocks_1);
	void collect_single_align_tasks(vector<function<void()> >& tasks, const int64_t n_blocks, const string& sam_path, const string& reference_path, const string& aln_param, const string& samse_param, const string& file_path);
	void merge_SAM(const string& file_path, const string& bam_path, const string& sam_path, int64_t n_blocks);
private:
	int64_t n_blocks_1;
	int64_t n_blocks_2;
	int32_t n_cores;
};

} /* namespace meerkat */

#endif /* MEERKAT_BWACALLER_HPP_ */
