/*
 * ParallelDiscordExtractor.hpp
 *
 *  Created on: Jun 14, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#ifndef MEERKAT_PARALLELDISCORDEXTRACTOR_HPP_
#define MEERKAT_PARALLELDISCORDEXTRACTOR_HPP_

#include <functional>
#include <algorithm>
#include <parallel/algorithm>
#include <set>
#include <boost/format.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamMultiReader.h>
#include "../../castle/OptionParser.hpp"
#include "../../castle/TimeChecker.hpp"
#include "../../castle/IOUtils.hpp"
#include "../../castle/ParallelRunner.hpp"
#include "../../castle/Hashs.hpp"

#include "../BlockBoundary.hpp"
#include "Group.hpp"

namespace meerkat {
	using namespace std;
	using namespace BamTools;
	typedef google::dense_hash_map<string, bool, castle::StringHash> StringBoolMap;
	class ParallelDiscordExtractor {
		public:
			ParallelDiscordExtractor();
			~ParallelDiscordExtractor();
			void set_option_parser(const castle::OptionParser& the_options);
			void extract_discordant_reads();
			void extract_discordant_reads_serial();
			void collect_boundaries(const int64_t size_block);
			void collect_boundaries_alt(vector<BlockOffset>& offset_blocks);
			void collect_independent_boundaries();
			void collect_chromosome_wide_boundaries();
			void collect_clipped_discordants();
			void collect_clipped_discordants_serial();
			void reset_stats();
			void write_discordants_alt();
			void write_discordants_serial();
			void collect_discordants();
			void collect_discordants_alt();
			void collect_discordants_sambamba();
			void create_bni_even_index(const string& a_path, string bni_index_path);
		private:
			static const int64_t size_block = 8192000;
			static const bool verbose = false;
			static const bool silent = false;
			int64_t n_cores;

			int64_t not_paired;
			int64_t not_both_mapped;
			int64_t not_same_chr;
			int64_t same_pos;
			int64_t not_diff_strands;
			int64_t weird_insert_size;
			int64_t big_insert_size;
			int64_t num_memorized_clipped_discs;

			castle::OptionParser options;
			set<string> black_listed;
			vector<BlockBoundary> fixed_size_blocks;
			vector<BlockBoundary> independent_blocks;
			vector<BlockBoundary> unmapped_included_blocks;
			vector<BlockBoundary> chromosome_wide_blocks;
			unordered_map<string, bool> memo_clipped_disc;
		private:
			string f(const string& s);
			bool is_discordant(const BamAlignment &al,
				double median, double stdev, int nstdevs,
				int64_t& not_paired, int64_t& not_both_mapped,
				int64_t& same_pos, int64_t& not_same_chr,
				int64_t& not_diff_strands, int64_t& weird_insert_size,
				int64_t& big_insert_size);
			bool is_clipped(const BamAlignment &al);
			bool check_blacklist(castle::blacklist_t &blacklist,
				const RefVector &refs, long ref_id, long pos);
			bool is_blacklisted(const BamAlignment &al, const RefVector &refs,
				castle::blacklist_t &blist);
			void get_bai_index_path(const string& parent_path, string& bai_path);
			void get_bni_index_path(const string& parent_path, string& bni_path);
	};
	inline string ParallelDiscordExtractor::f(const string& s) {
		return s.size() == 0 ? "none" : s;
	}
	inline bool ParallelDiscordExtractor::is_discordant(const BamAlignment &al,
		double median, double stdev, int nstdevs, int64_t& not_paired,
		int64_t& not_both_mapped, int64_t& same_pos, int64_t& not_same_chr,
		int64_t& not_diff_strands, int64_t& weird_insert_size,
		int64_t& big_insert_size) {
		/* Must be a paired end read */
		if (!al.IsPaired()) {
			++not_paired;
			return false;
		}

		/* Both reads must be mapped */
		if (!al.IsMapped() || !al.IsMateMapped()) {
			++not_both_mapped;
			return false;
		}

		if (al.Position == al.MatePosition && al.RefID == al.MateRefID) {
			++same_pos;
			return false;
		}

		/* If the mates map to different chroms, then discordant */
		if (al.RefID != al.MateRefID) {
			++not_same_chr;
			return true;
		}

		/* If both reads map to the same strand, then discordant */
		if (al.IsReverseStrand() == al.IsMateReverseStrand()) {
			++not_diff_strands;
			return true;
		}

		/* If the insert size is the wrong sign, then discordant */
		if ((al.IsReverseStrand() && al.InsertSize > 0)
				|| (!al.IsReverseStrand() && al.InsertSize < 0)) {
			++weird_insert_size;
			return true;
		}

		/* If the insert size is too large, then discordant */
		if (abs(al.InsertSize) > ceil(median + nstdevs * stdev)) {
			++big_insert_size;
			return true;
		}

		return false;
	}

	/* Check for either S or H CIGAR operations */
	inline bool ParallelDiscordExtractor::is_clipped(const BamAlignment &al) {
		vector<CigarOp>::const_iterator iter;
		for (iter = al.CigarData.begin(); iter != al.CigarData.end(); ++iter) {
			if (iter->Type == 'H' || iter->Type == 'S') {
				return true;
			}
		}
		return false;
	}
	/* vector refs */
	inline bool ParallelDiscordExtractor::check_blacklist(
		castle::blacklist_t &blacklist, const RefVector &refs, long ref_id,
		long pos) {
		auto iter = blacklist.find(refs[ref_id].RefName);
		if (blacklist.end() != iter) {
			auto iter2 = iter->second.find(pos);
			if (iter->second.end() != iter2) {
				return true;
			}
		}
		return false;
	}
	inline bool ParallelDiscordExtractor::is_blacklisted(const BamAlignment &al,
		const RefVector &refs, castle::blacklist_t &blist) {
		/* If either mate is blacklisted, both mates are discarded */
		return check_blacklist(blist, refs, al.RefID, al.Position) ||
				check_blacklist(blist, refs, al.MateRefID, al.MatePosition);
	}

	inline void ParallelDiscordExtractor::get_bai_index_path(const string& parent_path, string& bai_path) {
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
	inline void ParallelDiscordExtractor::get_bni_index_path(const string& parent_path, string& bni_path) {
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

//		bni_path = options.prefix + ".bni";
//		if (!options.working_dir.empty()) {
//			bni_path = options.working_prefix + ".bni";
//		}
	}

} /* namespace meerkat */

#endif /* MEERKAT_PARALLELDISCORDEXTRACTOR_HPP_ */
