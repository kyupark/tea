/*
 * ParallelBamReader.hpp
 *
 *  Created on: Jun 2, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#ifndef PARALLELBAMREADER_HPP_
#define PARALLELBAMREADER_HPP_
#include <string>
#include <sstream>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <boost/regex.hpp>
#include <parallel/algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/ParallelRunner.hpp"
#include "../../castle/IOUtils.hpp"
#include "../../castle/Hashs.hpp"
#include "../../emma/IntervalTree.hpp"
#include "api/BamReader.h"
#include "api/BamWriter.h"

#include "Histogram.hpp"
#include "ReadGroup.h"
#include "ReadGroupAlt.hpp"

#include "PairedAlignment.hpp"
#include "../BlockBoundary.hpp"
#include "../../third/gzstream.h"
#include "../BWACaller.hpp"

namespace meerkat {

struct BAMStatistics {
	BAMStatistics() :
			num_total(0), num_unmapped(0), num_clipped(0), max_read_len(0) {
	}
	int64_t num_total;
	int64_t num_unmapped;
	int64_t num_clipped;
	int64_t max_read_len;
};

using namespace std;
using namespace BamTools;

typedef google::dense_hash_map<string, int8_t, castle::StringHash> StringInt8Map;
typedef google::dense_hash_set<string, castle::StringHash> StringDenseSet;
typedef google::dense_hash_map<string, BamAlignment, castle::StringHash> StringAlignmentMap;
typedef google::sparse_hash_map<string, BamAlignment, castle::StringHash> SparseStringAlignmentMap;
typedef google::dense_hash_map<string, PairedAlignmentWithTwoIds, castle::StringHash> StringPairedAlignmentWithTwoIdsMap;
class ParallelBamReader {
	static const bool silent = false;
	static const bool verbose = false;
public:
	ParallelBamReader();
	~ParallelBamReader();
	void set_option_parser(const castle::OptionParser& the_options);
	void collect_boundaries(const int64_t size_block);
	void collect_chromosome_wide_boundaries();
	void collect_boundaries_alt(vector<BlockOffset>& offset_blocks);
	void preprocess();
	void collect_basic_statistics();
	void load_first_stat();
	void store_first_stat();
	void collect_basic_statistics_alt();
	void output_blacklist_file();
	void find_quality_standard_and_max_read_length(vector<uint64_t>& quality_sample_space);

	void output_unmapped();

	void collect_second_statistics();
	void output_softclips();
	void output_read_groups_alt();
	void output_unprocessed_read_groups(boost::unordered_map<string, PairedAlignmentWithTwoIds>& paired_alns, vector<string>& read_groups, map<string, int64_t>& read_groups_reverse_index);
	void create_reports();
	void align_clipped_reads();
	void report_qualities(vector<uint64_t>& quality_sample_space);

	void create_bni_even_index(const string& a_path, string bni_index_path);

//	pair<int64_t, int64_t> overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
//	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);

//			void output_softclips();
//			void output_softclips_and_mates();

private:
	int64_t n_cores;
	int64_t num_total;
	int64_t max_read_len;
	int64_t num_unmapped;
	int64_t num_clipped;
	int64_t unmapped_rejected;
	int32_t n_read_groups;
	int qenc;
	vector<BlockBoundary> actual_blocks;
	vector<BlockBoundary> fixed_size_blocks;
	vector<BlockBoundary> unmapped_included_blocks;
	vector<uint64_t> offset_blocks;
	vector<BlockBoundary> chromosome_wide_blocks;
	castle::OptionParser options;
//	unordered_map<string, int8_t> softclips;
//	StringDenseSet softclips;
	boost::unordered_set<string> softclips;
	set<string> black_listed;
	set<string> uu_names;

private:
	static const int64_t size_block = 4096000;
	void witness_qualities(const BamAlignment &al, vector<uint64_t>& local_sample_space);
	void get_bai_index_path(const string& parent_path, string& bai_path);
	void get_bni_index_path(const string& parent_path, string& bni_path);
};
inline void ParallelBamReader::witness_qualities(const BamAlignment &al, vector<uint64_t>& local_sample_space) {
	for (int i = 0; i < al.Length; ++i)
		++local_sample_space[static_cast<int>(al.Qualities[i])];
}

inline void ParallelBamReader::get_bai_index_path(const string& parent_path, string& bai_path) {
	bai_path = parent_path;
	bai_path += ".bai";
	if(boost::filesystem::exists(bai_path)) {
		return;
	}
	string target_path(parent_path);
	target_path.back() = 'i';
	if(boost::filesystem::exists(target_path)) {
		bai_path = target_path;
		return;
	}
	bai_path = options.prefix + ".tmp.bai";
	if (!options.working_dir.empty()) {
		bai_path = options.working_prefix + ".tmp.bai";
	}
	if(boost::filesystem::exists(bai_path)) {
		return;
	}
	string sambamba_cmd = (boost::format("sambamba index -t %d %s %s") % n_cores % parent_path % bai_path).str();
	system(sambamba_cmd.c_str());
}

inline void ParallelBamReader::get_bni_index_path(const string& parent_path, string& bni_path) {

	string target_path(parent_path);
	target_path[target_path.size() - 2] = 'n';
	target_path[target_path.size() - 1] = 'i';
	if(boost::filesystem::exists(target_path)) {
		bni_path = target_path;
		return;
	}
	bni_path = parent_path;
	bni_path += ".bni";
	if(boost::filesystem::exists(bni_path)) {
		return;
	}
//	bni_path = options.prefix + ".bni";
//	if (!options.working_dir.empty()) {
//		bni_path = options.working_prefix + ".bni";
//	}
}

//TODO: if I change change the statement, if (-1 == b1 && -1 == b2) to if (0 == b1 && 0 == b2), the algorithm finds more SV events.
//inline pair<int64_t, int64_t> ParallelBamReader::overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
//
////	if (-1 == b1 && -1 == b2) {
//	if (0 == b1 && 0 == b2) {
//		return make_pair(a1, a2);
//	}
////		( $a1, $a2 ) = sort { $a <=> $b } ( $a1, $a2 );
//	if (a1 > a2) {
//		swap(a1, a2);
//	}
//	if (b1 > b2) {
//		swap(b1, b2);
//	}
////		( $b1, $b2 ) = sort { $a <=> $b } ( $b1, $b2 );
//	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
//		int64_t a = (a1 > b1) ? a1 : b1;
//		int64_t b = (a2 > b2) ? b2 : a2;
//		return make_pair(a, b);
//	} else {
//		// the original version returns (0, 0)
////		return make_pair(b1, b2);
//		return make_pair(0, 0);
//	}
//}
}
/* namespace meerkat */

#endif /* PARALLELBAMREADER_HPP_ */
