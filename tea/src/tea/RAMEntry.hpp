/*
 * RAMEntry.hpp
 *
 *  Created on: Aug 18, 2016
 *      Author: el174
 */

#ifndef TEA_RAMENTRY_HPP_
#define TEA_RAMENTRY_HPP_
#include <string>
#include <vector>
#include <set>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "../emma/IntervalTree.hpp"
namespace tea {
using namespace std;
	struct RAMEntry {
		BamTools::BamAlignment bam;
		string map;
//		BamTools::BamAlignment map;
	};
	struct RefRepeatEntry {
		string chromosome;
		int64_t start;
		int64_t end;
		string repeat_name;
		string repeat_class;
		string repeat_family;
	};
	typedef emma::IntervalType<RefRepeatEntry> RefRepeatIntervalEntry;
	typedef std::vector<RefRepeatIntervalEntry> RefRepeatIntervalVector;
	typedef emma::IntervalTreeType<RefRepeatEntry> RefRepeatIntervalTree;

	struct GeneEntry {
		char strand;
		int64_t s;
		int64_t e;
		string name;
		string type;
		string str() const {
			return (boost::format("%s\t%d\t%d\t%s\t%s") % strand % s % e % name % type).str();
		}
	};

	typedef emma::IntervalType<GeneEntry> GeneIntervalEntry;
	typedef std::vector<GeneIntervalEntry> GeneIntervalVector;
	typedef emma::IntervalTreeType<GeneEntry> GeneIntervalTree;

	struct RAMRepeatEntry {
		string repeat_name;
		string repeat_class;
		string repeat_family;
		string read_name;
		int64_t pos;
		RAMRepeatEntry() : pos(0){
		}
		RAMRepeatEntry(const RAMRepeatEntry& other) {
			repeat_name = other.repeat_name;
			repeat_class = other.repeat_class;
			repeat_family = other.repeat_family;
			read_name = other.read_name;
			pos = other.pos;
		}
		RAMRepeatEntry& operator=(const RAMRepeatEntry& other) {
			// check for self-assignment
			if (this == &other) {
				return *this;
			}
			repeat_name = other.repeat_name;
			repeat_class = other.repeat_class;
			repeat_family = other.repeat_family;
			read_name = other.read_name;
			pos = other.pos;
			return *this;
		}
		bool operator==(const RAMRepeatEntry& other) const {
			return pos == other.pos && read_name == other.read_name && repeat_name == other.repeat_name && repeat_class == other.repeat_class && repeat_family == other.repeat_family;
		}
		bool operator!=(const RAMRepeatEntry& other) const {
			return !(*this == other);
		}
		bool operator<(const RAMRepeatEntry& other) const {
			if(pos < other.pos) {
				return true;
			} else if(pos > other.pos) {
				return false;
			}
			if(repeat_name < other.repeat_name) {
				return true;
			} else if(repeat_name > other.repeat_name) {
				return false;
			}
			if(repeat_family < other.repeat_family) {
				return true;
			} else if(repeat_family > other.repeat_family) {
				return false;
			}
			if(repeat_class < other.repeat_class) {
				return true;
			} else if(repeat_class > other.repeat_class) {
				return false;
			}
			return (read_name < other.read_name);
		}
		string str() const {
			return (boost::format("%s, %s, %s, %s, %s") % pos % repeat_name % repeat_class % repeat_family % read_name).str();
		}
	};
	struct RepeatClusterEntry {
		RepeatClusterEntry() :
			 s(0), e(0), global_cluster_id(0), ram(0) {}
		int64_t s;
		int64_t e;
		int64_t global_cluster_id;
		int64_t ram;
		set<string> rep_repeat;
		vector<string> family;
		vector<string> repeat_class;
		vector<string> repeat_name;
		vector<int64_t> pos;
		vector<string> rname;
		bool operator<(const RepeatClusterEntry& other) const {
			if(s < other.s) {
				return true;
			} else if(s > other.s) {
				return false;
			}
			if(e < other.e) {
				return true;
			} else if(e > other.e) {
				return false;
			}
			return false;
		}
		string str() const {
			int64_t the_size = e - s + 1;
			return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s") % s % e % the_size %
					castle::StringUtils::join(rep_repeat, ",") %
					castle::StringUtils::join(repeat_class, ",") %
					castle::StringUtils::join(family, ",") %
					ram %
					castle::StringUtils::join(pos, ",") %
					castle::StringUtils::join(rname, ",")
					).str();
		}
	};
	typedef emma::IntervalType<RepeatClusterEntry> RAMIntervalEntry;
	typedef std::vector<RAMIntervalEntry> RAMIntervalVector;
	typedef emma::IntervalTreeType<RepeatClusterEntry> RAMIntervalTree;
	struct ClippedEntry {
		ClippedEntry() : ram_start(0), ram_end(0), strand(0), negative_pos(0), clipped_pos(0), clipped_pos_qual_trimmed(0), clipped_pos_rep(0), aligned(0) {

		}
		string chr;
		int64_t ram_start;
		int64_t ram_end;
		int8_t strand;
		string rep_repeat;
		int64_t negative_pos;
		int64_t clipped_pos;
		int64_t clipped_pos_qual_trimmed;
		int64_t clipped_pos_rep;
		int8_t aligned;
		string cigar_original;
		string cigar_corrected;
		string read_name;
		string ref_seq;
		string clipped_seq;
		string clipped_qual;
		string clipped_seq_qual_trimmed;
		string clipped_qual_qual_trimmed;
		string clipped_seq_rep;
		string clipped_qual_rep;
		string str() const {

//			if(-1 == strand) {
//				return (boost::format("%s\t%d\t%d\t%s\t-%d(-%d)\t%s(%s)\t%s\t%s\t%s\t%s") % chr % ram_start %
//						ram_end % rep_repeat % negative_pos % clipped_pos_rep % cigar_original % cigar_corrected %
//						read_name % ref_seq % clipped_seq_rep % clipped_qual_rep).str();
//			}
//			return (boost::format("%s\t%d\t%d\t%s\t%d(%d)\t%s(%s)\t%s\t%s\t%s\t%s") % chr % ram_start %
//									ram_end % rep_repeat % negative_pos % clipped_pos_rep % cigar_original % cigar_corrected %
//									read_name % ref_seq % clipped_seq_rep % clipped_qual_rep).str();
			if(-1 == strand) {
				return (boost::format("%s\t%d\t%d\t%s\t-%d\t%s\t%s\t%s\t%s\t%s\t%s") % chr % ram_start %
						ram_end % rep_repeat % clipped_pos_rep % static_cast<int32_t>(aligned) % cigar_original %
						read_name % ref_seq % clipped_seq_rep % clipped_qual_rep).str();
			}
			return (boost::format("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s") % chr % ram_start %
									ram_end % rep_repeat % clipped_pos_rep % static_cast<int32_t>(aligned) % cigar_original %
									read_name % ref_seq % clipped_seq_rep % clipped_qual_rep).str();
		}
	};
	struct ClippedStatEntry {
		ClippedStatEntry() : s(0), e(0), size(0), tsd(-9999), pbp(0), nbp(0), ram(0), pram(0), nram(0), cr(0), pcr(0), ncr(0), acr(0), pacr(0), nacr(0), acrr(0.0),
				pram_start(0), pram_end(0), nram_start(0), nram_end(0), pgene("-"), ngene("-"), score(0.0), s2n(0.0), oi(0), desc("-"), conf(1){

		}
		string chr;
		int64_t s;
		int64_t e;
		int64_t size;
		int64_t tsd;
		int64_t pbp;
		int64_t nbp;
		string rep_repeat;
		string repeat_family;
		string repeat_class;
		string rep_suffix;
		int64_t ram;
		int64_t pram;
		int64_t nram;
		int64_t cr;
		int64_t pcr;
		int64_t ncr;
		int64_t acr;
		int64_t pacr;
		int64_t nacr;
		double acrr;
		int64_t pram_start;
		int64_t pram_end;
		int64_t nram_start;
		int64_t nram_end;
		string pgene;
		string ngene;
		double score;
		double s2n;
		int8_t oi;
		string desc;
		int32_t conf;
		string str() const {
			return (boost::format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\t%s\t%.2f\t%.2f\t%d\t%s\t%d") % chr % s %
					e % size % tsd % pbp % nbp % rep_repeat % repeat_family % repeat_class % ram % pram % nram % cr % pcr % ncr % acr % pacr % nacr
					% acrr % pram_start % pram_end % nram_start % nram_end % pgene % ngene % score % s2n % static_cast<int32_t>(oi) % desc % conf).str();
		}
	};
	struct AlnPairEntry {
		string file_name_prefix;
		int64_t pos;
		string chr;
	};

	struct AlnSeqQualEntry {
		string seq1;
		string seq2;
		string qual1;
		string qual2;
	};
}

#endif /* TEA_RAMENTRY_HPP_ */
