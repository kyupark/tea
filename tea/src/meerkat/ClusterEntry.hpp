/*
 * ClusterEntry.hpp
 *
 *  Created on: Jun 25, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 */

#ifndef MEERKAT_CLUSTERENTRY_HPP_
#define MEERKAT_CLUSTERENTRY_HPP_

#include <string>
#include <boost/format.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include "../emma/IntervalTree.hpp"
#include "../castle/StringUtils.hpp"

namespace meerkat {
using namespace std;
namespace bmidx = boost::multi_index;

struct IntervalEntryType {
//	IntervalEntryType(const string& a_p_id, const string& a_s_id, int64_t a_start_pos, int64_t an_end_pos, int64_t a_mate_start_pos, int64_t a_mate_end_pos) :
//		p_id(a_p_id), s_id(a_s_id), start_pos(a_start_pos), end_pos(an_end_pos), mate_start_pos(a_mate_start_pos), mate_end_pos(a_mate_end_pos) {
//	}
	string p_id;
	string s_id;
	// 1 denotes the forward strand while -1 represents the reverse strand
	// (first aln.: second aln.)
	// type 0: (1:-1)
	// type 1: (1:1)
	// type 2: (-1:-1)
	// type 3: (-1:1)
//	int8_t the_type;

//	string seq_id;
//	int64_t start_pos;
//	int64_t end_pos;
//
//	string mate_seq_id;
//	int64_t mate_start_pos;
//	int64_t mate_end_pos;
	string str() const {
		return (boost::format("%s_%s")
//				"(%s) %s:(%d-%d), %s:(%d-%d)")
				% p_id % s_id
//				% static_cast<int32_t>(the_type) % seq_id % start_pos % end_pos % mate_seq_id % mate_start_pos % mate_end_pos
				).str();
	}
};

typedef emma::IntervalType<IntervalEntryType> IntervalEntry;
typedef std::vector<IntervalEntry> IntervalEntryVector;
typedef emma::IntervalTreeType<IntervalEntryType> IntervalClusterTree;

struct ClusterFullEntry {
	string readname;
	string rg;
	string ref_id;
	int64_t n_supports;
	int32_t strand;
	int64_t start;
	int64_t end;
	int64_t len;
	string mate_ref_id;
	int64_t n_mate_supports;
	int32_t mate_strand;
	int64_t mate_start;
	int64_t mate_end;
	int64_t mate_len;
	int64_t isize;
	vector<int64_t> cbps;
	ClusterFullEntry() : n_supports(0), strand(0), start(0), end(0), len(0), n_mate_supports(0),
			mate_strand(0), mate_start(0), mate_end(0), mate_len(0), isize(0) {}
};

struct ClusterEntry {
	string readname;
	string rg;
	string seqid;
	int32_t strand;
	int64_t start;
	int64_t end;
	int64_t len;
	string mseqid;
	int32_t mstrand;
	int64_t mstart;
	int64_t mend;
	int64_t mlen;
	int64_t isize;
};


struct OrientationEntry {
	int32_t strand;
	int32_t mate_strand;
	string ref_id;
	string mate_ref_id;
};
struct CoordinateEntry {
	string chr_bp_1;
	int64_t start;
	int64_t end;
	string chr_bp_2;
	int64_t mate_start;
	int64_t mate_end;
	int32_t orientation;
	string str;
	CoordinateEntry() : start(-1), end(-1), mate_start(-1), mate_end(-1), orientation(0) {
	}
	CoordinateEntry(const CoordinateEntry& other) {
		chr_bp_1 = other.chr_bp_1;
		start = other.start;
		end = other.end;
		chr_bp_2 = other.chr_bp_2;
		mate_start = other.mate_start;
		mate_end = other.mate_end;
		orientation = other.orientation;
		str = other.str;
	}
	CoordinateEntry& operator=(const CoordinateEntry& other) {
				// check for self-assignment
			if (this == &other)
				return *this;
			chr_bp_1 = other.chr_bp_1;
					start = other.start;
					end = other.end;
					chr_bp_2 = other.chr_bp_2;
					mate_start = other.mate_start;
					mate_end = other.mate_end;
					orientation = other.orientation;
					str = other.str;
					return *this;
		}
	bool operator==(const CoordinateEntry& other) const {
		return chr_bp_1 == other.chr_bp_1 && start == other.start && end == other.end && chr_bp_2 == other.chr_bp_2
				&& mate_start == other.mate_start && mate_end == other.mate_end && orientation == other.orientation;
	}
};
struct EventEntry {
	string type;
	string cluster_id;
	// for MatePairDiscordantCaller, this value really represent mate_cluster_id
	// for SplitReadSVCaller, this value represents the original n_supports string
	string mate_cluster_id;
	string n_suppports_str;
	string n_mate_supports_str;
	int64_t n_supports;
	int64_t n_mate_support;
	string ref_id;
	int64_t event_start;
	int8_t strand;
	int8_t mate_strand;
	int64_t event_end;
	int64_t event_size_1;
	int64_t event_size_2;
	int64_t distance;
	int64_t distance_1;
	int64_t distance_2;
	string mate_ref_id;
	int64_t mate_event_start;
	int64_t mate_event_end;
	vector<int64_t> cbps;
	EventEntry() : n_supports(0), n_mate_support(0), event_start(-1), strand(0), mate_strand(0), event_end(-1),
			event_size_1(-1), event_size_2(-1), distance(-1), distance_1(-1), distance_2(-1), mate_event_start(-1), mate_event_end(-1){
	}

	EventEntry(const EventEntry& other) {
		type = other.type;
		cluster_id = other.cluster_id;
		mate_cluster_id = other.mate_cluster_id;
		n_suppports_str = other.n_suppports_str;
		n_mate_supports_str = other.n_mate_supports_str;
		n_supports = other.n_supports;
		n_mate_support = other.n_mate_support;
		ref_id = other.ref_id;
		event_start = other.event_start;
		strand = other.strand;
		mate_strand = other.mate_strand;
		event_end = other.event_end ;
		event_size_1 = other.event_size_1;
		event_size_2 = other.event_size_2;
		distance = other.distance;
		distance_1 = other.distance_1;
		distance_2 = other.distance_2;
		mate_ref_id = other.mate_ref_id;
		mate_event_start = other.mate_event_start;
		mate_event_end = other.mate_event_end;
		cbps = other.cbps;
	}
	EventEntry& operator=(const EventEntry& other) {
			// check for self-assignment
		if (this == &other)
			return *this;
		type = other.type;
		cluster_id = other.cluster_id;
		mate_cluster_id = other.mate_cluster_id;
		n_suppports_str = other.n_suppports_str;
		n_mate_supports_str = other.n_mate_supports_str;
		n_supports = other.n_supports;
		n_mate_support = other.n_mate_support;
		ref_id = other.ref_id;
		event_start = other.event_start;
		strand = other.strand;
		mate_strand = other.mate_strand;
		event_end = other.event_end ;
		event_size_1 = other.event_size_1;
		event_size_2 = other.event_size_2;
		distance = other.distance;
		distance_1 = other.distance_1;
		distance_2 = other.distance_2;
		mate_ref_id = other.mate_ref_id;
		mate_event_start = other.mate_event_start;
		mate_event_end = other.mate_event_end;
		cbps = other.cbps;
		return *this;
	}

	bool operator<(const EventEntry& other) const {
		if (cluster_id < other.cluster_id) {
			return true;
		} else if (cluster_id > other.cluster_id) {
			return false;
		}
		if (mate_cluster_id < other.mate_cluster_id) {
			return true;
		} else if (cluster_id > other.cluster_id) {
			return false;
		}
		if (event_start < other.event_start) {
			return true;
		} else if (event_start > other.event_start) {
			return false;
		}
		if (event_end < other.event_end) {
			return true;
		} else if (event_end > other.event_end) {
			return false;
		}

		if (mate_event_start < other.mate_event_start) {
			return true;
		} else if (mate_event_start > other.mate_event_start) {
			return false;
		}
		if (mate_event_end < other.mate_event_end) {
			return true;
		} else if (mate_event_end > other.mate_event_end) {
			return false;
		}
		return false;
	}
	bool operator==(const EventEntry& other) const {
		if (type == other.type && cluster_id == other.cluster_id
				&& mate_cluster_id == other.mate_cluster_id && n_supports == other.n_supports
				&& n_mate_support == other.n_mate_support && ref_id == other.ref_id
				&& event_start == other.event_start
				&& event_end == other.event_end
				&& strand == other.strand
				&& mate_strand == other.mate_strand
				&& event_size_1 == other.event_size_1
				&& event_size_2 == other.event_size_2
				&& distance == other.distance
				&& distance_1 == other.distance_1
				&& distance_2 == other.distance_2
				&& mate_ref_id == other.mate_ref_id
				&& mate_event_start == other.mate_event_start
				&& mate_event_end == other.mate_event_end) {
			if (cbps.size() != other.cbps.size()) {
				return false;
			}
			for (uint64_t c_id = 0; c_id < cbps.size(); ++c_id) {
				if (cbps[c_id] != other.cbps[c_id]) {
					return false;
				}
			}
			return true;
		}
		return false;
	}
	bool operator!=(const EventEntry& other) const {
		return !(*this == other);
	}
	string str() {
		if(type.empty()) {
			return "";
		}
//		if(cbps.size() > 3 && (-1 == cbps[2] || -1 == cbps[3])) {
//			return "";
//		}
		if(cbps.size() < 8) {
			cbps.resize(8);
		}
		vector<int64_t> cbps_tmp_1(cbps.begin(), cbps.begin() + 4);
		vector<int64_t> cbps_tmp_2(cbps.begin() + 4, cbps.begin() + 8);
		string cbps1_str = castle::StringUtils::join(cbps_tmp_1, ":");
		string cbps2_str = castle::StringUtils::join(cbps_tmp_2, ":");

//		if(8 == cbps.size()) {
//			vector<int64_t> cbps_tmp_1(cbps.begin(), cbps.begin() + 4);
//			cbps1_str = castle::StringUtils::join(cbps_tmp_1, ":");
//			vector<int64_t> cbps_tmp_2(cbps.begin() + 4, cbps.begin() + 8);
//			cbps2_str = castle::StringUtils::join(cbps_tmp_2, ":");
//		}
//		if(!ref_id.empty() && !mate_ref_id.empty()) {
////			if(string::npos != type.find("invers")) {
//				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\t%s")
//						% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
//									% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
//									% (-1 == event_size_2 ? "" : boost::lexical_cast<string>(event_size_2) + "\t")
//									% cbps1_str
//									% cbps2_str).str() ;
////			} else {
////				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s")
////					% type % cluster_id % n_supports % n_mate_support % ref_id % event_start
////					% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end % event_size_2 % distance %
////					cbps1_str % cbps2_str).str();
////			}
//		}
//		return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s")
//		% type % cluster_id % n_supports % ref_id % event_start % event_end % event_size_1 % cbps1_str).str();
		// this event may have changed its value, hence discarding mate_cluster_id and n_mate_support
		if(cluster_id == mate_cluster_id && n_supports == n_mate_support) {
			mate_cluster_id = "";
			n_mate_support = 0;
		}
		return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
		% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
											% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
											% event_size_2 % distance % distance_1 % distance_2 % static_cast<int32_t>(strand)
											% static_cast<int32_t>(mate_strand) % cbps1_str % cbps2_str).str();
	}
	string mpd_str() {
			if(type.empty()) {
				return "";
			}
	//		if(cbps.size() > 3 && (-1 == cbps[2] || -1 == cbps[3])) {
	//			return "";
	//		}
			if(cbps.size() < 8) {
				cbps.resize(8);
			}
			if(type == "tandem_dup" || type == "invers_f" || type == "invers_r" || type == "transl_inter" || type == "del") {
				for(int64_t c_id = 0; c_id < 4; ++ c_id) {
					auto a_cbp = cbps[c_id];
					if(0 == a_cbp) {
						return "";
					}
				}
			} else {
				for(auto a_cbp: cbps) {
					if(0 == a_cbp) {
						return "";
					}
				}
			}
			vector<int64_t> cbps_tmp_1(cbps.begin(), cbps.begin() + 4);
			vector<int64_t> cbps_tmp_2(cbps.begin() + 4, cbps.begin() + 8);
			string cbps1_str = castle::StringUtils::join(cbps_tmp_1, ":");
			string cbps2_str = castle::StringUtils::join(cbps_tmp_2, ":");

			if(cluster_id == mate_cluster_id && n_supports == n_mate_support) {
				mate_cluster_id = "";
				n_mate_support = 0;
			}
			if(type == "tandem_dup" || type == "invers_f" || type == "invers_r") {
				return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % type % cluster_id % n_supports % ref_id % event_start
				% event_end % event_size_1 % cbps1_str).str();
			} else if(string::npos != type.find("inssu") || string::npos != type.find("inssd") || string::npos != type.find("insod") || string::npos != type.find("insou")) {
				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2 % distance % cbps1_str % cbps2_str).str();
			} else if(string::npos != type.find("inso") || string::npos != type.find("inss")) {
				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2 % cbps1_str % cbps2_str).str();
			} else if(string::npos != type.find("invers")) {
				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% cbps1_str % cbps2_str).str();
			} else if(type == "transl_inter") {
				if(cbps[2] < mate_event_end && cbps[3] > mate_event_end) {
					return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
						% type % cluster_id % n_supports % ref_id % event_start % static_cast<int32_t>(strand)
						% mate_ref_id % mate_event_end
						% static_cast<int32_t>(mate_strand) % cbps1_str).str();
				} else {
					return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
						% type % cluster_id % n_supports % ref_id % event_start % static_cast<int32_t>(strand)
						% mate_ref_id % mate_event_start
						% static_cast<int32_t>(mate_strand) % cbps1_str).str();
				}
			} else if(type == "del") {
				return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % n_supports % ref_id % event_start
				% event_end % event_size_1 % cbps1_str).str();
			} else {
				return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2 % distance % distance_1 % distance_2 % static_cast<int32_t>(strand)
				% static_cast<int32_t>(mate_strand) % cbps1_str % cbps2_str).str();
			}
		}
	string mpd_pure_str() {
				if(type.empty()) {
					return "";
				}
		//		if(cbps.size() > 3 && (-1 == cbps[2] || -1 == cbps[3])) {
		//			return "";
		//		}
				if(cbps.size() < 8) {
					cbps.resize(8);
				}
				vector<int64_t> cbps_tmp_1(cbps.begin(), cbps.begin() + 4);
				vector<int64_t> cbps_tmp_2(cbps.begin() + 4, cbps.begin() + 8);
				string cbps1_str = castle::StringUtils::join(cbps_tmp_1, ":");
				string cbps2_str = castle::StringUtils::join(cbps_tmp_2, ":");

				if(cluster_id == mate_cluster_id && n_supports == n_mate_support) {
					mate_cluster_id = "";
					n_mate_support = 0;
				}
				if(type == "tandem_dup" || type == "invers_f" || type == "invers_r") {
					return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % type % cluster_id % n_supports % ref_id % event_start
					% event_end % event_size_1 % cbps1_str).str();
				} else if(string::npos != type.find("inssu") || string::npos != type.find("inssd") || string::npos != type.find("insod") || string::npos != type.find("insou")) {
					return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
					% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
					% event_size_2 % distance % cbps1_str % cbps2_str).str();
				} else if(string::npos != type.find("inso") || string::npos != type.find("inss")) {
					return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
					% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
					% event_size_2 % cbps1_str % cbps2_str).str();
				} else if(string::npos != type.find("invers")) {
					return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
					% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
					% cbps1_str % cbps2_str).str();
				} else if(type == "transl_inter") {
					if(cbps[2] < mate_event_end && cbps[3] > mate_event_end) {
						return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
							% type % cluster_id % n_supports % ref_id % event_start % static_cast<int32_t>(strand)
							% mate_ref_id % mate_event_end
							% static_cast<int32_t>(mate_strand) % cbps1_str).str();
					} else {
						return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
							% type % cluster_id % n_supports % ref_id % event_start % static_cast<int32_t>(strand)
							% mate_ref_id % mate_event_start
							% static_cast<int32_t>(mate_strand) % cbps1_str).str();
					}
				} else if(type == "del") {
					return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% type % cluster_id % n_supports % ref_id % event_start
					% event_end % event_size_1 % cbps1_str).str();
				} else {
					return (boost::format("%s\t%s/%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
					% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
					% event_size_2 % distance % distance_1 % distance_2 % static_cast<int32_t>(strand)
					% static_cast<int32_t>(mate_strand) % cbps1_str % cbps2_str).str();
				}
			}
	string sr_str() {
			if(type.empty()) {
				return "";
			}
	//		if(cbps.size() > 3 && (-1 == cbps[2] || -1 == cbps[3])) {
	//			return "";
	//		}

			if("transl_inter" == type) {
				return (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% type % cluster_id % mate_cluster_id % n_supports
				// data[3]:n_supports/data[4]:n_mate_support
				% ref_id % event_start % static_cast<int32_t>(strand) % mate_ref_id % mate_event_start % static_cast<int32_t>(mate_strand)).str();
			} else if("del" == type || "invers_r" == type) {
				return (boost::format("%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\n")
				% type % cluster_id % mate_cluster_id % n_supports % ref_id % event_start
				% event_end % event_size_1).str();
			} else if("tandem_dup" == type|| "invers_f" == type ) {
				return (boost::format("%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\n")
				% type % cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1).str();
			} else if("invers" == type) {
				return (boost::format("%s\t%s\t%s\t%s/%s\t%s\t%d\t%d\t%d\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1).str();
			} else if ("del_inssu" == type || string::npos != type.find("inssd") || string::npos != type.find("inssu") || string::npos != type.find("insod") || string::npos != type.find("insou")) {
				return (boost::format("%s\t%s\t%s\t%d/%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2).str();
			} else if("inss" == type || "inso" == type || "del_inss" == type || "del_inso" == type) {
				return (boost::format("%s\t%s\t%s\t%d/%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2).str();
			} else {
				return (boost::format("%s\t%s\t%s\t%d/%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n")
				% type % cluster_id % mate_cluster_id % n_supports % n_mate_support % ref_id % event_start
				% event_end % event_size_1 % mate_ref_id % mate_event_start % mate_event_end
				% event_size_2 % distance_1 % distance_2).str();
			}

		}
};

struct BLASTEntry {
	// entry[0]: read name
	string query_id;
	// entry[1]:  cluster name
	string subject_name;
	// entry[2]
	double percent_identity;
	// entry[3]
	int64_t aln_length;
	// entry[4]
	int64_t n_mismatches;
	// entry[5]
	int64_t n_gap_opens;
	// entry[6]
	int64_t query_start;
	// entry[7]
	int64_t query_end;
	// entry[8]
	int64_t subject_start;
	// entry[9]
	int64_t subject_end;
	// entry[10]
	double expect_value;
	// entry[11]
	double bit_score;

};

struct RepeatEntry {
	int64_t start;
	int64_t end;
	string repeat_class;
	string name;
};

struct Int64Pair {
	int64_t first;
	int64_t second;
	Int64Pair(int64_t first, int64_t second) :
			first(first), second(second) {
	}
	bool operator<(const Int64Pair& e) const {
		return first < e.first;
	}
};

// define a multiply indexed set with indices by first key and second value
typedef bmidx::multi_index_container<Int64Pair,
		bmidx::indexed_by<
				// sort by greater<int64_t> on second
				bmidx::ordered_non_unique<bmidx::member<Int64Pair, int64_t, &Int64Pair::second> >, //
				bmidx::hashed_unique<bmidx::member<Int64Pair, int64_t, &Int64Pair::first> >//
		>//
> int64_t_value_sortedset;

typedef bmidx::multi_index_container<Int64Pair,
		bmidx::indexed_by<
				// sort by greater<int64_t> on second
				bmidx::ordered_non_unique<bmidx::member<Int64Pair, int64_t, &Int64Pair::second>, std::greater<int64_t>>, //
				bmidx::hashed_unique<bmidx::member<Int64Pair, int64_t, &Int64Pair::first> >//
		>//
> int64_t_value_desc_sortedset;

struct StringInt64Pair {
	string first;
	int64_t second;
	StringInt64Pair(string first, int64_t second) :
			first(first), second(second) {
	}
	bool operator<(const StringInt64Pair& e) const {
		return first < e.first;
	}
};

// define a multiply indexed set with indices by first key and second value
typedef bmidx::multi_index_container<StringInt64Pair,
		bmidx::indexed_by<
				// sort by greater<int64_t> on second
				bmidx::ordered_non_unique<bmidx::member<StringInt64Pair, int64_t, &StringInt64Pair::second> >, //
				bmidx::hashed_unique<bmidx::member<StringInt64Pair, string, &StringInt64Pair::first> >//
		>//
> string_int64_t_value_sortedset;

typedef bmidx::multi_index_container<StringInt64Pair,
		bmidx::indexed_by<
				// sort by greater<int64_t> on second
				bmidx::ordered_non_unique<bmidx::member<StringInt64Pair, int64_t, &StringInt64Pair::second>, std::greater<int64_t> >, //
				bmidx::hashed_unique<bmidx::member<StringInt64Pair, string, &StringInt64Pair::first> >//
		>//
> string_int64_t_value_desc_sortedset;
}

#endif /* MEERKAT_CLUSTERENTRY_HPP_ */
