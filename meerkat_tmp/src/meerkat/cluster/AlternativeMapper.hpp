/*
 * AlternativeMapper.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

#ifndef MEERKAT_ALTERNATIVEMAPPER_HPP_
#define MEERKAT_ALTERNATIVEMAPPER_HPP_

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "api/BamReader.h"
#include "api/BamWriter.h"

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"

#include "../cluster/AlternativeMap.hpp"

namespace meerkat {
using namespace std;
using namespace BamTools;
class AlternativeMapper {

public:
	AlternativeMapper();
	~AlternativeMapper();
	void set_option_parser(const castle::OptionParser& the_options);
	void select_alternative_mapping();
	void select_alternative_mapping_alt();
	void select_alternative_mapping_serial();

	void create_raw_alternative_mapping_serial();
	void create_raw_alternative_mapping();
	void create_raw_alternative_mapping_alt();
	void create_raw_alternative_mapping_par();
	void _create_raw_alternative_mapping_par(const string& output_path, const string& input_path);
	void create_bfi(vector<int64_t>& block_boundary, const string& a_path);
	bool alt_map(vector<BamAlignment>& disc_alg, ofstream& rawalt_fh, const RefVector& ref_vec);
	bool write_alt_map(ostream& out, const RefVector& ref_vec, const boost::unordered_map<string, int64_t>& ref_reverse_index, vector<BamAlignment>& recent_als);

	void select_most_likely_alignments(vector<BamAlignment>& alignments);
	void select_maximum_likely_alignments(vector<BamAlignment>& alignments);

	bool has_valid_insertion_size(ostream& out, const RefVector& ref_vec, const boost::unordered_map<string, int64_t>& ref_reverse_index, const vector<BamAlignment>& alignments, bool& found_valid_entry);

private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;
};


inline void AlternativeMapper::select_most_likely_alignments(vector<BamAlignment>& alignments) {
	if (3 > alignments.size()) {
		return;
	}

//	int64_t max_match_1 = -1;
//	int64_t max_match_id_1 = -1;
//	int64_t max_match_2 = -1;
//	int64_t max_match_id_2 = -1;
//	for(uint64_t a_id = 0; a_id < alignments.size(); ++ a_id) {
//		auto& al = alignments[a_id];
//		int64_t n_matches = 0;
//		for(auto a_cigar : al.CigarData) {
//			if('M' == a_cigar.Type) {
//				n_matches += a_cigar.Length;
//			}
//		}
//		if(al.IsFirstMate()) {
//			if(n_matches > 0 && max_match_1 < n_matches) {
//				max_match_1 = n_matches;
//				max_match_id_1 = a_id;
//			}
//		} else if(al.IsSecondMate()) {
//			if(n_matches > 0 && max_match_2 < n_matches) {
//				max_match_2 = n_matches;
//				max_match_id_2 = a_id;
//			}
//		}
//	}
//	vector<BamAlignment> max_alignments;
//	if(-1 != max_match_1) {
//		max_alignments.push_back(alignments[max_match_id_1]);
//	}
//	if(-1 != max_match_2) {
//		max_alignments.push_back(alignments[max_match_id_2]);
//	}
//	alignments.swap(max_alignments);
//	const bool debug = "ST-E00104:502:HFJN5CCXX:2:2122:22830:10415" == alignments[0].Name;
	vector<int64_t> pair_1;
	vector<int64_t> pair_2;
	for (uint64_t a_id = 0; a_id < alignments.size(); ++a_id) {
		if (alignments[a_id].IsFirstMate()) {
			pair_1.push_back(a_id);
		} else if (alignments[a_id].IsSecondMate()) {
			pair_2.push_back(a_id);
		}
	}
	//TODO: I think that to increase the sensitivity and precision, the sum of the match should be
	// used instead.
	pair<int64_t, int64_t> the_largest_match_in_pair_1(0, 0);
	pair<int64_t, int64_t> the_largest_match_in_pair_2(0, 0);
	vector<pair<int64_t, int64_t>> match_up_and_down_1(pair_1.size(), make_pair(0, 0));
	vector<pair<int64_t, int64_t>> match_up_and_down_2(pair_2.size(), make_pair(0, 0));
	for (uint64_t p_id = 0; p_id < pair_1.size(); ++p_id) {
		for (auto cigar : alignments[pair_1[p_id]].CigarData) {
			if ('M' != cigar.Type) {
				continue;
			}
			if (cigar.Length > the_largest_match_in_pair_1.second) {
				the_largest_match_in_pair_1.first = p_id;
				the_largest_match_in_pair_1.second = cigar.Length;
			}
		}
	}
	for (uint64_t p_id = 0; p_id < pair_2.size(); ++p_id) {
		for (auto cigar : alignments[pair_2[p_id]].CigarData) {
			if ('M' != cigar.Type) {
				continue;
			}
			if (cigar.Length > the_largest_match_in_pair_2.second) {
				the_largest_match_in_pair_2.first = p_id;
				the_largest_match_in_pair_2.second = cigar.Length;
			}
		}
	}
	{
		int64_t state = 0;
		for (uint64_t p_id = 0; p_id < match_up_and_down_1.size(); ++p_id) {
			for (auto cigar : alignments[pair_1[p_id]].CigarData) {
				if ('D' == cigar.Type) {
					continue;
				}
				if (the_largest_match_in_pair_1.second == cigar.Length) {
					state = 1;
					continue;
				}
				// count the base match ups.
				if (0 == state) {
					match_up_and_down_1[p_id].first += cigar.Length;
				} else {
					// count the base match downs.
					match_up_and_down_1[p_id].second += cigar.Length;
				}
			}
		}
	}
	{
		int64_t state = 0;
		for (uint64_t p_id = 0; p_id < match_up_and_down_2.size(); ++p_id) {
			for (auto cigar : alignments[pair_2[p_id]].CigarData) {
				if ('D' == cigar.Type) {
					continue;
				}
				if (the_largest_match_in_pair_2.second == cigar.Length) {
					state = 1;
					continue;
				}
				// count the base match ups.
				if (0 == state) {
					match_up_and_down_2[p_id].first += cigar.Length;
				} else {
					// count the base match downs.
					match_up_and_down_2[p_id].second += cigar.Length;
				}
			}
		}
	}
	while (alignments.size() > 2) {
		if (1 < pair_1.size()) {
			// select the entry with the smallest match among pair-1 entries, which will be removed.
			int64_t lhs_id = 0;
			for (uint64_t rhs_id = 1; rhs_id < pair_1.size(); ++rhs_id) {
				auto& lhs_aln_1 = alignments[pair_1[lhs_id]];
				auto& rhs_aln_1 = alignments[pair_1[rhs_id]];

				if (lhs_aln_1.IsReverseStrand()) {
					if (rhs_aln_1.IsReverseStrand()) {
						if (match_up_and_down_1[lhs_id].first > match_up_and_down_1[rhs_id].first) {
							lhs_id = rhs_id;
						}
					} else {
						if (match_up_and_down_1[lhs_id].first > match_up_and_down_1[rhs_id].second) {
							lhs_id = rhs_id;
						}
					}
				} else {
					if (rhs_aln_1.IsReverseStrand()) {
						if (match_up_and_down_1[lhs_id].second > match_up_and_down_1[rhs_id].first) {
							lhs_id = rhs_id;
						}
					} else {
						if (match_up_and_down_1[lhs_id].second > match_up_and_down_1[rhs_id].second) {
							lhs_id = rhs_id;
						}
					}
				}
			}
			alignments.erase(alignments.begin() + pair_1[lhs_id]);
		}
		//		bool debug = "ST-E00104:502:HFJN5CCXX:1:1210:6329:16463"
		//				== alignments[0].Name;
		if (2 == alignments.size()) {
			return;
		}
		if (1 < pair_2.size()) {
			// select the entry with the smallest match among pair-1 entries, which will be removed.
			int64_t lhs_id = 0;
			for (uint64_t rhs_id = 1; rhs_id < pair_2.size(); ++rhs_id) {
				auto& lhs_aln_2 = alignments[pair_2[lhs_id]];
				auto& rhs_aln_2 = alignments[pair_2[rhs_id]];

				if (lhs_aln_2.IsReverseStrand()) {
					if (rhs_aln_2.IsReverseStrand()) {
						//						if (debug) {
						//							cout << "here-del-0: " << pair_2[lhs_id] << "/"
						//									<< pair_2[rhs_id] << "\n";
						//						}
						if (match_up_and_down_2[lhs_id].first > match_up_and_down_2[rhs_id].first) {
							//							if (debug) {
							//								cout << "here-del-1: " << pair_2[lhs_id] << "/"
							//										<< pair_2[rhs_id] << "\n";
							//							}
							lhs_id = rhs_id;
						}
					} else {
						//						if (debug) {
						//							cout << "here-del-2: " << pair_2[lhs_id] << "/"
						//									<< pair_2[rhs_id] << "\n";
						//						}
						if (match_up_and_down_2[lhs_id].first > match_up_and_down_2[rhs_id].second) {
							//							if (debug) {
							//								cout << "here-del-3: " << pair_2[lhs_id] << "/"
							//										<< pair_2[rhs_id] << "\n";
							//							}
							lhs_id = rhs_id;
						}
					}
				} else {
					if (rhs_aln_2.IsReverseStrand()) {
						//						if (debug) {
						//							cout << "here-del-4: " << pair_2[lhs_id] << "/"
						//									<< pair_2[rhs_id] << "\n";
						//						}
						if (match_up_and_down_2[lhs_id].second > match_up_and_down_2[rhs_id].first) {
							//							if (debug) {
							//								cout << "here-del-5: " << pair_2[lhs_id] << "/"
							//										<< pair_2[rhs_id] << "\n";
							//							}
							lhs_id = rhs_id;
						}
					} else {
						//						if (debug) {
						//							cout << "here-del-6: " << pair_2[lhs_id] << "/"
						//									<< pair_2[rhs_id] << "\n";
						//						}
						if (match_up_and_down_2[lhs_id].second > match_up_and_down_2[rhs_id].second) {
							//							if (debug) {
							//								cout << "here-del-7: " << pair_2[lhs_id] << "/"
							//										<< pair_2[rhs_id] << "\n";
							//							}
							lhs_id = rhs_id;
						}
					}
				}
			}
			alignments.erase(alignments.begin() + pair_2[lhs_id]);
		}
	}
}
inline void AlternativeMapper::select_maximum_likely_alignments(vector<BamAlignment>& alignments) {
	int64_t max_match_1 = -1;
	int64_t max_match_id_1 = -1;
	int64_t max_match_2 = -1;
	int64_t max_match_id_2 = -1;
	for(uint64_t a_id = 0; a_id < alignments.size(); ++ a_id) {
		auto& al = alignments[a_id];
		int64_t n_matches = 0;
		for(auto a_cigar : al.CigarData) {
			if('M' == a_cigar.Type) {
				n_matches += a_cigar.Length;
			}
		}
		if(al.IsFirstMate()) {
			if(n_matches > 0 && max_match_1 < n_matches) {
				max_match_1 = n_matches;
				max_match_id_1 = a_id;
			}
		} else if(al.IsSecondMate()) {
			if(n_matches > 0 && max_match_2 < n_matches) {
				max_match_2 = n_matches;
				max_match_id_2 = a_id;
			}
		}
	}
	vector<BamAlignment> max_alignments;
	if(-1 != max_match_1) {
		max_alignments.push_back(alignments[max_match_id_1]);
	}
	if(-1 != max_match_2) {
		max_alignments.push_back(alignments[max_match_id_2]);
	}
	alignments.swap(max_alignments);
}

inline bool AlternativeMapper::has_valid_insertion_size(ostream& out, const RefVector& ref_vec, const boost::unordered_map<string, int64_t>& ref_reverse_index, const vector<BamTools::BamAlignment>& alignments, bool& found_valid_entry) {
	if (2 != alignments.size()) {
		return false;
	}
//	const bool debug = string::npos != alignments[0].Name.find("ST-E00104:502:HFJN5CCXX:7:1123:12266:39510");
//		const bool debug = false;
	auto& data1 = alignments[0];
	auto& data2 = alignments[1];
	int32_t strand = data1.IsReverseStrand() ? -1 : 1;
	int32_t mstrand = data2.IsReverseStrand() ? -1 : 1;

	int64_t max_len_1 = 0;
	int64_t max_len_2 = 0;
//	if (debug) {
//		cout << "[AlternativeMapper.has_valid_insertion_size] here-21-1-a: " << BamWriter::GetSAMAlignment(data1, ref_vec) << "\n";
//		cout << "[AlternativeMapper.has_valid_insertion_size] here-21-1-b: " << BamWriter::GetSAMAlignment(data2, ref_vec) << "\n";
//	}
	for (auto cigar : data1.CigarData) {
		if ('M' == cigar.Type && max_len_1 < cigar.Length) {
			max_len_1 = cigar.Length;
		}
	}
//		if (debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-21-2\n";
//		}
	for (auto cigar : data2.CigarData) {
		if ('M' == cigar.Type && max_len_2 < cigar.Length) {
			max_len_2 = cigar.Length;
		}
	}
//		if (debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-21-3\n";
//		}
	//	If both the entry 1 and entry 2 in the same chromosome, the insert size represents the distance between two entries

	int64_t insert_size = -1;
	if (data1.RefID == data2.RefID) {
//			if (debug) {
//				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-4\n";
//			}
		if (strand == mstrand) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-5\n";
//				}
			insert_size = abs(data1.Position - data2.Position) + 1;
		} else if (1 == strand && mstrand == -1) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-6\n";
//				}
			insert_size = abs(data1.Position - data2.Position) + max_len_2;
		} else if (-1 == strand && 1 == mstrand) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-7\n";
//				}
			insert_size = abs(data1.Position - data2.Position) - max_len_1 + 2;
		}
	}
//		if (debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-21-8: " << insert_size << "\n";
//		}
	string rg;
	if (!data1.GetReadGroup(rg)) {
//			if (debug) {
//				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-9\n";
//			}
		rg = "none";
	}

	if (data1.RefID == data2.RefID && 1 == strand && -1 == mstrand && insert_size <= options.is[rg]["isu"]) {
//			if (debug) {
////				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-10: sd_cutoff_cl: " << options.sd_cutoff_cl << "\n";
////				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-10: isd: " << options.is[rg]["isd"] << "\n";
////				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-10: isd_rep: " << options.is["isd"]["selected"] << "\n";
//				cout << (boost::format("[AlternativeMapper.has_valid_insertion_size] here-21-10: insert: %s isu: %s\n") % insert_size % options.is[rg]["isu"]).str();
////				cout << "[AlternativeMapper.has_valid_insertion_size] here-21-10: isu_rep: " << options.is["isu"]["selected"] << "\n";
//
//			}
		return false;
	}
	uint32_t nm = 0;
	uint32_t mnm = 0;

	data1.GetTag("NM", nm);
	data2.GetTag("NM", mnm);
	if (options.ad_align) {
//		if (debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-21-11\n";
//		}
		char xt = 'U';
		char mate_xt = 'U';
		string xt_tag_str;
		string mate_xt_tag_str;
		if (data1.GetTag("XT", xt_tag_str)) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-12\n";
//				}
			xt = xt_tag_str[0];
		}
		if (data2.GetTag("XT", mate_xt_tag_str)) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-13\n";
//				}
			mate_xt = mate_xt_tag_str[0];
		}
		if (!data1.AlignmentFlag) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-14\n";
//				}
			xt = 'R';
		}
		if (!data2.AlignmentFlag) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-15\n";
//				}
			mate_xt = 'R';
		}
		string xa;
		int32_t tr = 0;
		vector<string> alt_mapping;
		int64_t n_alt_mapping = 0;

		const char* delim_semi_colon = ";";
		if ('R' == xt) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-16\n";
//				}
			if (!data1.GetTag("XA", xa)) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-17\n";
//					}
				tr = 1;
			}
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-18\n";
//				}
			castle::StringUtils::c_string_multi_split(xa, delim_semi_colon, alt_mapping);
			n_alt_mapping = alt_mapping.size();
			if (-1 != options.alt_map_max && n_alt_mapping > options.alt_map_max) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-19\n";
//					}
				tr = 1;
			}
		}

		string mate_xa;
		int32_t mate_tr = 0;
		vector<string> mate_alt_mapping;
		int64_t n_mate_alt_mapping = 0;
		if ('R' == mate_xt) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-20\n";
//				}
			if (!data2.GetTag("XA", mate_xa)) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-21\n";
//					}
				mate_tr = 1;
			}
		}
		castle::StringUtils::c_string_multi_split(mate_xa, delim_semi_colon, mate_alt_mapping);
		n_mate_alt_mapping = mate_alt_mapping.size();
		if (-1 != options.alt_map_max && n_mate_alt_mapping > options.alt_map_max) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-22\n";
//				}
			mate_tr = 1;
		}
		if ('R' == xt || 'R' == mate_xt) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-23\n";
//				}
			if (1 == tr || 1 == mate_tr) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-24\n";
//					}
				return false;
			}
		}
		double weight = 1;
		if ('R' == xt || 'R' == mate_xt) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-25\n";
//				}
			weight = (n_alt_mapping > n_mate_alt_mapping) ? n_alt_mapping + 1 : n_mate_alt_mapping + 1;
			weight = 1 / weight;
		}

//		if(debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-out-0\n";
//		}
		found_valid_entry = true;
		if(data1.RefID <= data2.RefID) {
			out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % ref_vec[data1.RefID].RefName
				% (data1.Position + 1) % strand //
				% ref_vec[data2.RefID].RefName % (data2.Position + 1) % mstrand //
				% insert_size % nm % mnm % rg % max_len_1 % max_len_2 % static_cast<int64_t>(weight)).str();
		} else {
			out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data2.Name % ref_vec[data2.RefID].RefName
				% (data2.Position + 1) % mstrand //
				% ref_vec[data1.RefID].RefName % (data1.Position + 1) % strand //
				% insert_size % mnm % nm % rg % max_len_2 % max_len_1 % static_cast<int64_t>(weight)).str();
		}
		vector<string> xa_cols;
		const char* xa_delims = ",";
		auto& ref_blacklist = options.point_black_lists;
		auto& ref_is = options.is;
		for (auto a_map : alt_mapping) {
			castle::StringUtils::c_string_multi_split(a_map, xa_delims, xa_cols);
			string seqid_alt = xa_cols[0];
			int64_t start_alt = boost::lexical_cast<int64_t>(xa_cols[1].substr(1));
			string CIGAR_str = xa_cols[2];
			uint32_t nm_alt = boost::lexical_cast<int64_t>(xa_cols[3]);
			// check coverage based blacklist
			if (ref_blacklist.end() != ref_blacklist.find(make_pair(seqid_alt, start_alt))
					|| ref_blacklist.end() != ref_blacklist.find(make_pair(seqid_alt, start_alt + (options.cut_sr << 1)))) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-28\n";
//					}
				continue;
			}
			string mseqid_alt = ref_vec[data2.RefID].RefName;
			int64_t mstart_alt = data2.Position;
			int32_t mstrand_alt = mstrand;
			uint32_t mnm_alt = mnm;
			int64_t len_alt = max_len_1;
			int64_t mlen_alt = max_len_2;

			if ((ref_vec[data1.RefID].RefName == seqid_alt) && (data1.Position == start_alt)) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-29\n";
//					}
				continue;
			}
			int32_t strand_alt = '+' == xa_cols[1][0] ? 1 : -1;
			int64_t isize_alt = 0;
			if (seqid_alt == mseqid_alt) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-30\n";
//					}
				if (start_alt <= mstart_alt) {
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-31\n";
//						}
					if (strand_alt == mstrand_alt) {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-32\n";
//							}
						isize_alt = abs(mstart_alt - start_alt) + 1;
					} else {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-33\n";
//							}
						if (-1 == strand_alt && 1 == mstrand_alt) {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-34\n";
//								}
							isize_alt = abs(mstart_alt - start_alt) - mlen_alt + 2;
						} else {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-35\n";
//								}
							isize_alt = abs(mstart_alt - start_alt) + mlen_alt;
						}
					}
				} else {
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-36\n";
//						}
					if (strand_alt == mstrand_alt) {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-37\n";
//							}
						isize_alt = abs(start_alt - mstart_alt) + 1;
					} else {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-38\n";
//							}
						if (-1 == strand_alt && 1 == mstrand_alt) {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-39\n";
//								}
							isize_alt = abs(start_alt - mstart_alt) + mlen_alt;
						} else {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-40\n";
//								}
							isize_alt = abs(start_alt - mstart_alt) - mlen_alt + 2;
						}
					}
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-41\n";
//						}
					swap(seqid_alt, mseqid_alt);
					swap(start_alt, mstart_alt);
					swap(strand_alt, mstrand_alt);
					swap(nm_alt, mnm_alt);
					swap(len_alt, mlen_alt);
				}
				if (seqid_alt > mseqid_alt) {
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-42\n";
//						}
					swap(seqid_alt, mseqid_alt);
					swap(start_alt, mstart_alt);
					swap(strand_alt, mstrand_alt);
					swap(nm_alt, mnm_alt);
					swap(len_alt, mlen_alt);
				}
			}
			if (seqid_alt == mseqid_alt && 1 == strand_alt && -1 == mstrand_alt && isize_alt <= ref_is[rg]["isu"]) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-43\n";
//					}
				return false;
			}
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-44\n";
//				}
			auto seqid_itr = ref_reverse_index.find(seqid_alt);
			auto mseqid_itr = ref_reverse_index.find(mseqid_alt);
			int32_t ref_id = -1;
			int32_t mate_ref_id = -1;
			if(ref_reverse_index.end() != seqid_itr) {
				ref_id = seqid_itr->second;
			}
			if(ref_reverse_index.end() != mseqid_itr) {
				mate_ref_id = mseqid_itr->second;
			}
//			if(debug) {
//				cout << "[AlternativeMapper.has_valid_insertion_size] here-out-1\n";
//			}
			found_valid_entry = true;
			if(ref_id <= mate_ref_id) {
				out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % seqid_alt % (start_alt + 1) % strand_alt //
					% mseqid_alt % (mstart_alt + 1) % mstrand_alt //
					% isize_alt % nm_alt % mnm_alt % rg % len_alt % mlen_alt % static_cast<int64_t>(weight)).str();
			} else {
				out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % mseqid_alt % (mstart_alt + 1) % mstrand_alt //
									% seqid_alt % (start_alt + 1) % strand_alt //
									% isize_alt % mnm_alt % nm_alt % rg % mlen_alt % len_alt % static_cast<int64_t>(weight)).str();
			}
			//			out << "$readname\t$seqid_alt\t$start_alt\t$strand_alt\t$mseqid_alt\t$mstart_alt\t$mstrand_alt\t$isize_alt\t$nm_alt\t$mnm_alt\t$rg\t$len_alt\t$mlen_alt\t$weight\n";
		}
		for (auto a_map : mate_alt_mapping) {
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-45\n";
//				}
			castle::StringUtils::c_string_multi_split(a_map, xa_delims, xa_cols);
			string mseqid_alt = xa_cols[0];
			int64_t mstart_alt = boost::lexical_cast<int64_t>(xa_cols[1].substr(1));
			string mCIGAR_str = xa_cols[2];
			uint32_t mnm_alt = boost::lexical_cast<int64_t>(xa_cols[3]);
			if (ref_blacklist.end() != ref_blacklist.find(make_pair(mseqid_alt, mstart_alt))
					|| ref_blacklist.end() != ref_blacklist.find(make_pair(mseqid_alt, mstart_alt + (options.cut_sr << 1)))) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-46\n";
//					}
				continue;
			}
			string seqid_alt = ref_vec[data1.RefID].RefName;
			int64_t start_alt = data1.Position;
			int32_t strand_alt = strand;
			uint32_t nm_alt = nm;
			int64_t len_alt = max_len_1;
			int64_t mlen_alt = max_len_2;
			string mseqid = ref_vec[data2.RefID].RefName;
			int64_t mstart = data2.Position;
			if (mseqid == mseqid_alt && mstart == mstart_alt) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-47\n";
//					}
				continue;
			}
			int32_t mstrand_alt = '+' == xa_cols[1][0] ? 1 : -1;
			int64_t isize_alt = 0;

			if (seqid_alt == mseqid_alt) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-48\n";
//					}
				if (start_alt <= mstart_alt) {
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-49\n";
//						}
					if (strand_alt == mstrand_alt) {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-50\n";
//							}
						isize_alt = abs(mstart_alt - start_alt) + 1;
					} else {
						if (-1 == strand_alt && 1 == mstrand_alt) {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-51\n";
//								}
							isize_alt = abs(mstart_alt - start_alt) - mlen_alt + 2;
						} else {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-52\n";
//								}
							isize_alt = abs(mstart_alt - start_alt) + mlen_alt;
						}
					}
				} else {
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-53\n";
//						}
					if (strand_alt == mstrand_alt) {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-54\n";
//							}
						isize_alt = abs(start_alt - mstart_alt) + 1;
					} else {
//							if (debug) {
//								cout << "[AlternativeMapper.has_valid_insertion_size] here-21-55\n";
//							}
						if (-1 == strand_alt && 1 == mstrand_alt) {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-56\n";
//								}
							isize_alt = abs(start_alt - mstart_alt) + mlen_alt;
						} else {
//								if (debug) {
//									cout << "[AlternativeMapper.has_valid_insertion_size] here-21-57\n";
//								}
							isize_alt = abs(start_alt - mstart_alt) - mlen_alt + 2;
						}
					}
//						if (debug) {
//							cout << "[AlternativeMapper.has_valid_insertion_size] here-21-58\n";
//						}
					swap(seqid_alt, mseqid_alt);
					swap(start_alt, mstart_alt);
					swap(strand_alt, mstrand_alt);
					swap(nm_alt, mnm_alt);
					swap(len_alt, mlen_alt);
				}
			}
			if (seqid_alt > mseqid_alt) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-59\n";
//					}
				swap(seqid_alt, mseqid_alt);
				swap(start_alt, mstart_alt);
				swap(strand_alt, mstrand_alt);
				swap(nm_alt, mnm_alt);
				swap(len_alt, mlen_alt);
			}
			if (seqid_alt == mseqid_alt && 1 == strand_alt && -1 == mstrand_alt && isize_alt <= ref_is[rg]["isu"]) {
//					if (debug) {
//						cout << "[AlternativeMapper.has_valid_insertion_size] here-21-60\n";
//					}
				continue;
			}
//				if (debug) {
//					cout << "[AlternativeMapper.has_valid_insertion_size] here-21-61\n";
//				}
			auto seqid_itr = ref_reverse_index.find(seqid_alt);
			auto mseqid_itr = ref_reverse_index.find(mseqid_alt);
			int32_t ref_id = -1;
			int32_t mate_ref_id = -1;
			if(ref_reverse_index.end() != seqid_itr) {
				ref_id = seqid_itr->second;
			}
			if(ref_reverse_index.end() != mseqid_itr) {
				mate_ref_id = mseqid_itr->second;
			}
//			if(debug) {
//				cout << "[AlternativeMapper.has_valid_insertion_size] here-out-2\n";
//			}
			found_valid_entry = true;
			if(ref_id <= mate_ref_id) {
				out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % seqid_alt % (start_alt + 1) % strand_alt //
					% mseqid_alt % (mstart_alt + 1) % mstrand_alt //
					% isize_alt % nm_alt % mnm_alt % rg % len_alt % mlen_alt % static_cast<int64_t>(weight)).str();
			} else {
				out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % mseqid_alt % (mstart_alt + 1) % mstrand_alt //
					% seqid_alt % (start_alt + 1) % strand_alt //
					% isize_alt % mnm_alt % nm_alt % rg % mlen_alt % len_alt % static_cast<int64_t>(weight)).str();
			}
			//			push @toprint, "$readname\t$seqid_alt\t$start_alt\t$strand_alt\t$mseqid_alt\t$mstart_alt\t$mstrand_alt\t$isize_alt\t$nm_alt\t$mnm_alt\t$rg\t$len_alt\t$mlen_alt\t$weight\n";
		}
	} else {
//		if(debug) {
//			cout << "[AlternativeMapper.has_valid_insertion_size] here-out-3\n";
//		}
		found_valid_entry = true;
		if(data1.RefID <= data2.RefID) {
			out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data1.Name % ref_vec[data1.RefID].RefName % (data1.Position + 1)
				% strand //
				% ref_vec[data2.RefID].RefName % (data2.Position + 1) % mstrand //
				% insert_size % nm % mnm % rg % max_len_1 % max_len_2).str();
		} else {
			out << (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % data2.Name % ref_vec[data2.RefID].RefName % (data2.Position + 1)
				% mstrand //
				% ref_vec[data1.RefID].RefName % (data1.Position + 1) % strand //
				% insert_size % mnm % nm % rg % max_len_2 % max_len_1).str();
		}
		//		print $rawalt_fh "$readname\t$seqid\t$start\t$strand\t$mseqid\t$mstart\t$mstrand\t$isize\t$nm\t$mnm\t$rg\t$len\t$mlen\n";
	}
	return found_valid_entry;
}


} /* namespace meerkat */

#endif /* MEERKAT_ALTERNATIVEMAPPER_HPP_ */
