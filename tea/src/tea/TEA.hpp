/*
 * TEA.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: el174
 */

#ifndef TEA_TEA_HPP_
#define TEA_TEA_HPP_

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/range/adaptor/reversed.hpp>

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "../castle/ParallelRunner.hpp"
#include "../castle/StringUtils.hpp"
#include "../castle/IOUtils.hpp"
#include "../castle/MathUtils.hpp"
#include "../meerkat/BWACaller.hpp"
#include "../meerkat/BlockBoundary.hpp"
#include "../meerkat/preprocess/ReadGroup.h"
//#include "../meerkat/dre/ParallelDiscordExtractor.hpp"
#include "../meerkat/BlockBoundary.hpp"
#include "TEAOptionParser.hpp"
#include "RAMEntry.hpp"

namespace tea {
using namespace BamTools;
class TEA {
public:
	TEA();
	~TEA();
	void set_option_parser(const TEAOptionParser& options);
	void preprocess();
	void preprocess_v();
	void preprocess_u();
	void run_rid();
	void run_vid();
	void run_uid();
	// transduction
	void run_transduction();
	void create_germline_transduction(const string& out_path, const string& in_path);
	void create_discord_bed(const string& out_path, const string& in_path);
	void read_ram_read_ids(boost::unordered_set<string>& ram_read_ids);
	void read_ram_ids(boost::unordered_set<string>& ram_ids, const string& in_path, boost::unordered_map<string, int32_t>& cluster_entries);
	void create_umm_tmp_from_cluster_raw(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_set<string>& ram_id);
	void create_umm_tmp_from_cluster(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_map<string, int32_t> cluster_entries, const boost::unordered_set<string>& ram_read_ids);
	void write_non_dup_umm(const string& out_path, const string& in_path, const boost::unordered_map<string, int32_t>& dup_cnt_map);
	// transduction contig
	void run_transduction_contig();
	void create_contig_two_ram(const string& out_path, const string& in_path);
	void create_fa_from_germline_contig(const string& out_path_1, const string& out_path_2, const string& in_path);
	int64_t rfind_the_last_of(int64_t& n_len, const string& str, const char chr);
	int64_t find_the_last_of(int64_t& n_len, const string& str, const char chr);

	void collect_two_ram_map(boost::unordered_map<int64_t, string>& two_ram_map, const string& in_path);
	void collect_two_ram_seq_id_ref_set(boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map);
	void collect_aln_sam_repeat(set<string>& repeat_selected_seq_id, boost::unordered_set<string>& repeat_two_ram_id, const string& in_path);
	void collect_aln_ram_sam_repeat(boost::unordered_map<string, string>& repeat_ram_aligned_map, const string& in_path);
	void create_germline_transduction(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& repeat_ram_aligned_map, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);

	// orphan
	void run_orphan();
	void collect_gene_sets(set<string>& oo, set<string>& gene, const string& in_path);
	void collect_rescue_sets(set<string>& rescue, set<string>& gene, const string& in_path);
	void create_bed_filtered_intersected(const string& out_path, const set<string>& ooo, const string& in_path);
	void collect_cluster_id_pos_map(boost::unordered_map<string, string>& cluster_id_pos_map, const string& in_path);
	void create_intersected_filtered_insertion(const string& out_path, const boost::unordered_map<string, string> cluster_id_pos_map, const string& in_path);
	void collect_cluster_id_pair_map(boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path);
	void create_orphan_umm_from_cluster(const string& out_path, const boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path);

	// orphan contig
	void run_orphan_contig();
	void create_orphan_fa_from_germline_contig(const string& out_path, const string& in_path);
	void create_germline_orphan(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);

	void clean();

	// post processing
	void post_process();
	void create_orphan_list(const string& out_path, set<string>& orphan, const string& in_path);
	void create_transduction_filtered(const string& out_path, set<string>& trans, const string& in_path, const set<string>& orphan);
	void create_contig_filtered_fa(const string& out_path, const string& in_path);
	void collect_aln_repeat_selected_seq_id(set<string>& repeat_selected_seq_id, const string& in_path);
	void create_germline_contig_tmp(const string& out_path, const string& in_path, const set<string>& o);
	void create_germline_contig_tmp_tmp(const string& out_path, const string& in_path, const set<string>& allset);
	void refine_germline_contig_tmp_tmp(const string& out_path, const string& in_path);
	void create_germline_contig_tmp_tmp_filtered_fa(const string& out_path, const string& in_path);
	void collect_refined_aln_sam_ref(set<string>& rname, boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map);
	void collect_refined_aln_sam_repeat(set<string>& rrname, set<string>& oo, set<string>& ooo, const string& in_path, const set<string>& o);
	void create_short_transduction_list(const string& out_path, const set<string>& candidate, const set<string>& o, const set<string>& overlap, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);
	void collect_transduction_set(set<string>& o, const string& in_path);
	void create_post_contig_list(const string& out_path, const set<string>& o, const string& in_path);

	void output_raw_file(const string& chr, const string& cl_prefix,
			const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl, const map<int64_t, int64_t>& pm_cl, const boost::unordered_set<int64_t>& positive_only, const boost::unordered_set<int64_t>& negative_only, const int64_t read_length, const int64_t fragment_size, const bool headless);
	void BAM_to_FASTQ_serial(const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void BAM_to_FASTQ(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);

	void MEMBAM_to_FASTQ(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void collect_boundaries_un(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block);
	void collect_boundaries_pos(vector<meerkat::BlockBoundary>& fixed_size_blocks, vector<meerkat::BlockBoundary>& unmapped_included_blocks, vector<meerkat::BlockBoundary>& independent_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block);
	void collect_boundaries_alt(vector<BlockOffset>& offset_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path);
	void create_bni_even_index(const string& a_path, const string& a_bai_path, const string& a_bni_path);
//	void collect_boundaries(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& input_BAM_name, int64_t size_block, bool verbose = false);
private:
	void _BAM_to_FASTQ(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _MEMBAM_to_FASTQ(vector<int64_t>& block_boundary, const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _MEMBAM_to_FASTQ_serial(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void remove_entry_enclosed_with_large_H(vector<BamAlignment>& alns);
	void convert_H_to_S(vector<BamAlignment>& alns);
	void fill_H_to_S(BamAlignment& aln, const AlnSeqQualEntry& aln_seq_entry);
	void split_query_to_segments(vector<string>& seqs, vector<string>& quals, vector<BamAlignment>& alns);
	void format_isize();
	void create_disc_FASTQs();
	void create_um_FASTQs();
	void generate_va_bams();
	void generate_vam_files();
	void generate_um_bams(const string& a_path, const string& a_bai_path, const string& a_bni_path);
	void generate_um_raw_bam_serial();
	void generate_um_raw_bam(vector<meerkat::BlockBoundary>& actual_blocks);
	void remove_duplicates_serial(const string& in_file_name, const string& out_file_name);
	void remove_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_bam_file_name, const string& out_sam_file_name);
	void generate_umm(const string& in_file_name, const string& out_file_name);
	void write_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_file_name);
	void generate_ra_bams();
	void generate_ram_files();
	void generate_ram_file(const string& refbam, const string& rbamf, const string& ramf, const string& disc_1_ra_bam, const string& disc_2_ra_bam, bool exo = false, bool headless = false);
	void load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, vector<string>& rabam_files);
	void _load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, const int32_t end, const string& rabam);
	// The function, write_ram_and_bam, is not fully implemented, incorrect and inefficient
	void write_ram_and_bam(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo);
	void write_ram_and_bam_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless = false);
	void write_ram_and_bam_mem_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless = false);

	void generate_cbam_files();
	void _generate_cbam_files_mem();
	void _generate_cbam_files_mem_alt();
	void _generate_cbam_files_mem_alt2();
	void _generate_cbam_files_mem_org();
	void _generate_cbam_files_sampe();
	void generate_cbam_files_serial();
	void is_containing_S_or_H(bool& has_S, bool& has_H, vector<BamAlignment>& alns);
	bool has_S_or_H(const BamAlignment& aln);
	bool has_tag(const BamAlignment& aln, const char tag);

	void find_quality_standard(BamTools::BamReader& a_reader);
	void load_repeat_annotation(boost::unordered_map<string, pair<string, string>>& rannot);
	void load_virus_annotation(map<int64_t, string>& vannot);
	void load_ref_annotation(set<string>& chrl,
			boost::unordered_map<string, pair<string, string>>& rannot,
//			boost::unordered_map<string, boost::unordered_map<string, vector< pair<int64_t, int64_t> > > >& ril_annot,
			boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
			boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot,
			boost::unordered_map<string, GeneIntervalVector>& gene_annot, bool out_chrl, bool out_gap);
	void load_chr(set<string>& chrl) const;
	void read_gap_rfile(boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot, const string& file_name, const set<string>& chrl);
	void read_gene_rfile(boost::unordered_map<string, GeneIntervalVector>& gene_annot, const string& file_name, const set<string>& chrl);
	void load_read_length(boost::unordered_map<string, int32_t>& rl);
	void load_insert_size(boost::unordered_map<string, boost::unordered_map<string, double>>& is, boost::unordered_map<string, int32_t>& rl);
	void load_ram(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram, boost::unordered_map<string, pair<string, string>>& rannot, const bool rm_dup = false);
	void get_cluster_alt(const string& chr, RAMIntervalVector& cl, vector<RAMRepeatEntry>& sram, boost::unordered_map<string, pair<string, string>>& rannot, const int32_t strand, int64_t gap_cutoff);
	void pair_cluster_alt(map<int64_t, int64_t>& pm_cl, RAMIntervalVector& p_cl, RAMIntervalVector& n_cl, const int64_t gap_cutoff, const int64_t read_length, bool stringent_pair = false);
	void count_clipped(
			boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
			boost::unordered_map<string, GeneIntervalVector>& gene_annot, const string& chr, const string& cl_prefix, const string& contig_dir,
			const map<int64_t, int64_t>& pm_cl, const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl,
			boost::unordered_set<int64_t>& positive_only, boost::unordered_set<int64_t>& negative_only,
			const int64_t read_length, const int64_t fragment_size,
			const int64_t rmasker_filter_margin, const int64_t gene_margin, const bool headless);
	void count_clipped_v(
				boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
				map<int64_t, string>& vannot,
				const string& chr, const string& cl_prefix, const string& contig_dir,
				const map<int64_t, int64_t>& pm_cl, const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl,
				boost::unordered_set<int64_t>& positive_only, boost::unordered_set<int64_t>& negative_only,
				const int64_t read_length, const int64_t fragment_size,
				const int64_t rmasker_filter_margin, const int64_t gene_margin);
	void get_clipped_entries(vector<ClippedEntry>& clipped_entries, int64_t& max_pos_positive, int64_t& max_pos_negative, int64_t& n_positive_clipped_reads, int64_t& n_negative_clipped_reads,
			int64_t& n_aligned_clipped_positive, int64_t& n_aligned_clipped_negative,
			BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const int64_t read_length, const int64_t mid_point);
	void output_clipped_stat(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_germline, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree, RefRepeatIntervalVector& stat_results, GeneIntervalTree& gene_interval_tree, GeneIntervalVector& gene_results,
			BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin, const int64_t mid_point);
	void output_clipped_stat_v(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_germline, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree, RefRepeatIntervalVector& stat_results, const map<int64_t, string>& vannot,
				BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin);
	void output_mate_fa(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram);

	void _output_mate_fa(boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name,
			const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
			const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map);

	void output_mate_fa_v(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram);
	void _output_mate_fa_v(boost::unordered_map<string, vector<string>>& positive_mate_reads,
				boost::unordered_map<string, vector<string>>& negative_mate_reads,
				vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name,
				const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
				const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map);
	void _output_mate_fa_serial(boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			const string& input_BAM_name,
			const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
			const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map);

private:
	int64_t get_number_of_good_qualities(const string& a_qual, const int64_t qcutoff);
	int64_t get_number_of_low_qualities_at_begin(const string& a_qual, const int64_t qcutoff);
	int64_t get_number_of_low_qualities_at_end(const string& a_qual, const int64_t qcutoff);
	int64_t get_cpos(int32_t pos, std::vector<BamTools::CigarOp>& cigar, const string& qual, int8_t strand);
	void get_longest_fa(string& a_contig, const string& fa_path);
	void get_bai_index_path(const string& parent_path, string& bai_path);
	void get_bni_index_path(const string& parent_path, string& bni_path);
private:
	int32_t n_cores;
	int32_t qenc;
	uint8_t min_qual;
	TEAOptionParser options;
	vector<meerkat::BlockBoundary> fixed_size_blocks;
	vector<meerkat::BlockBoundary> independent_blocks;
	vector<meerkat::BlockBoundary> unmapped_included_blocks;
};

// maximal length of a good quality value
inline int64_t TEA::get_number_of_good_qualities(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	int64_t max_quals = 0;
	for(uint64_t qual_id = 0; qual_id < a_qual.size(); ++qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value > limit_val) {
			++n_quals;
			max_quals = max(n_quals, max_quals);
		} else {
			n_quals = 0;
		}
	}
	return max_quals;
}

inline int64_t TEA::get_number_of_low_qualities_at_begin(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	for(uint64_t qual_id = 0; qual_id < a_qual.size(); ++qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value == limit_val) {
			++n_quals;
		} else {
			break;
		}
	}
	return n_quals;
}
inline int64_t TEA::get_number_of_low_qualities_at_end(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	int64_t min_id = a_qual.size();
	--min_id;
	for(int64_t qual_id = min_id; qual_id >=0; --qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value == limit_val) {
			++n_quals;
		} else {
			break;
		}
	}
	return n_quals;
}

inline int64_t TEA::get_cpos(int32_t pos, std::vector<BamTools::CigarOp>& cigar, const string& qual, int8_t strand) {
	int32_t cpos = pos;
//	#print "pos: $pos, cigar: $cigar, qual: $qual, strand: $strand\n";

	if (1 == strand) {
		return cpos;
	}
	if(cigar.empty()) {
		return cpos;
	}

//	# for the negative strand clipped reads, count the number of base pairs (M & I not D)
	int32_t gap = 0;
	uint64_t c_id = 0;
	if('S' == cigar.front().Type) {
		c_id = 1;
	}
	for(; c_id < cigar.size(); ++c_id) {
		auto& a_cigar = cigar[c_id];
		if('I' == a_cigar.Type || 'S' == a_cigar.Type) {
			continue;
		}
		gap += a_cigar.Length;
	}
//	$cigar =~ s/^\d+S(.*)/$1/;
//	my @cnt = split(/[MINHPSD]/, $cigar);
//	my @chr = split(/\d+/, $cigar);
//	shift @chr;
//	for (my $i=0; $i<@chr; $i++) {
//		if ($chr[$i] ne "I" && $chr[$i] ne "S") { $gap += $cnt[$i] }
//	}
	cpos += gap;
//	#print "gap:$gap\n";
//	#print "cpos: $cpos\n";
	return -cpos;
}

inline void TEA::get_longest_fa(string& a_contig, const string& fa_path) {
	stringstream tmp_lines;
	string line;
	ifstream in(fa_path, ios::binary);
	while(getline(in, line, '\n')) {
		if('>' == line[0]) {
			string the_line = tmp_lines.str();
			if(the_line.size() > a_contig.size()) {
				a_contig = the_line;
			}
			tmp_lines.str(string());
			continue;
		}
		tmp_lines << line;
	}
	string the_line = tmp_lines.str();
	if(the_line.size() > a_contig.size()) {
		a_contig = the_line;
	}
}

inline void TEA::get_bai_index_path(const string& parent_path, string& bai_path) {
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
inline void TEA::get_bni_index_path(const string& parent_path, string& bni_path) {
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

} /* namespace tea */

#endif /* TEA_TEA_HPP_ */
