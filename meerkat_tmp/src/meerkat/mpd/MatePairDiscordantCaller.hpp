/*
 * MatePairDiscordantCaller.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lixing Yang
 */

/*
 * MatePairDiscordantCaller.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: el174
 */

#ifndef MEERKAT_MATEPAIRDISCORDANTCALLER_HPP_
#define MEERKAT_MATEPAIRDISCORDANTCALLER_HPP_

#include "../../castle/TimeChecker.hpp"
#include "../../castle/OptionParser.hpp"
#include "../../castle/StringUtils.hpp"
#include "../../castle/IOUtils.hpp"
#include "../ClusterEntry.hpp"

namespace meerkat {
using namespace std;

class MatePairDiscordantCaller {
	static const bool silent = false;
	static const bool verbose = false;
public:
	MatePairDiscordantCaller();
	~MatePairDiscordantCaller();
	void set_option_parser(const castle::OptionParser& the_options);
	void select_candidate_events_in_clusters();
	void select_candidate_events_in_clusters_alt();

	void select_discordant_clusters();
	void select_discordant_clusters_serial();

	void call_initial_events();
	void call_initial_events_serial();
	void call_intra_chromosomal_events(vector<EventEntry>& result_mp);
	void call_intra_chromosomal_events_alt(vector<EventEntry>& result_mp);

	void call_inter_chromosomal_events(vector<EventEntry>& result_mp);
	void call_inter_chromosomal_events_serial(vector<EventEntry>& result_mp);
	void call_rest_events(vector<EventEntry>& result_mp);
	void call_rest_events_alt(vector<EventEntry>& result_mp);
	void call_smaller_events_among_intra_events(vector<EventEntry>& result_mp);
	void call_smaller_events_among_intra_events_alt(vector<EventEntry>& result_mp);

	void call_smaller_events_among_inter_events(vector<EventEntry>& result_mp);
	void call_smaller_events_among_inter_events_alt(vector<EventEntry>& result_mp);

	void determine_a_cluster(ofstream& out, vector<ClusterEntry>& cluster, map<string, map<string, OrientationEntry>>& orientation,
				map<string, map<string, pair<int64_t, int64_t>>>& starts_map,
		map<string, map < string, pair<int64_t, int64_t>>>& mstarts_map,
		map<string, map < string, vector<int64_t>>>& cbps_map,
		map<string, map<string, map<int64_t, vector<string>>>>& cluster_type,
		map<string, int64_t>& support, map<string, int64_t>& support_f,
		const string& lastpid, const string& lastsid, int64_t l, int64_t j);

	pair<int64_t, int64_t> overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2);
	bool covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2);

private:
	int64_t n_cores;
	castle::OptionParser options;
	set<string> black_listed;

	map<string, map<string, OrientationEntry>> orientation;
	map<string, map<string, pair<int64_t, int64_t>>> starts_map;
	map<string, map<string, pair<int64_t, int64_t>>> mstarts_map;
	map<string, map<string, vector<int64_t>>> cbps_map;
	map<string, map<string, map<int64_t, vector<string>>>> cluster_type;
	map<string, int64_t> support;
	map<string, int64_t> support_f;
	vector<IntervalEntryVector> the_intervals;
	string the_last_pid;
};


inline void MatePairDiscordantCaller::determine_a_cluster(ofstream& out, vector<ClusterEntry>& cluster,
		map<string, map<string, OrientationEntry>>& orientation, map<string, map<string, pair<int64_t, int64_t>>>&starts_map,
		map<string, map < string, pair<int64_t, int64_t>>>& mstarts_map,
		map<string, map < string, vector<int64_t>>>& cbps_map,
		map<string, map<string, map<int64_t, vector<string>>>>& cluster_type,
		map<string, int64_t>& support, map<string, int64_t>& support_f,
		const string& lastpid, const string& lastsid, int64_t l, int64_t j) {
			auto& ref_is = options.is;
			vector<int64_t> starts;
			vector<int64_t> mstarts;
//				vector<int64_t> isizes;
		vector<int64_t> cbps;
		pair<int64_t, int64_t> cbp_0(-1, -1);
		pair<int64_t, int64_t> cbp_1(-1, -1);
//		const bool debug = "8" == lastpid;
				//|| "41657" == lastpid;
		const bool debug = false;
		if(debug) {
			cout << "cbp_test-10: "<< lastpid << "\n";
		}
		for (auto& ce : cluster) {
//					string& readname = ce.readname;
			string& rg = ce.rg;
			string& seqid = ce.seqid;
			int32_t strand = ce.strand;
			int64_t start = ce.start;
			int64_t end = ce.end;
			int64_t len = ce.len;
			string& mseqid = ce.mseqid;
			int32_t mstrand = ce.mstrand;
			int64_t mstart = ce.mstart;
			int64_t mend = ce.mend;
			int64_t mlen = ce.mlen;
//					int64_t isize = ce.isize;

			starts.push_back(start);
			starts.push_back(end);
			mstarts.push_back(mstart);
			mstarts.push_back(mend);
//					if(-1 != isize) {
//						isizes.push_back(isize);
//					}
			if (-1 == end) {
				end = start;
			}
			if (-1 == mend) {
				mend = mstart;
			}

			if (1 == strand) {
				pair<int64_t, int64_t> a_result = overlap(end - 10,
				start + ref_is[rg]["isu"] - mlen, cbp_0.first,
				cbp_0.second);
				cbp_0.first = a_result.first;
				cbp_0.second = a_result.second;
					if(debug) {
						cout << "cbp_test-10-0: "<< cbp_0.first << "," << cbp_0.second << "\n";
					}
//	 -10 or +10: expand the window to avoid one read too close to a break point
			} else {
				pair<int64_t, int64_t> a_result = overlap(
				end - ref_is[rg]["isu"] + mlen, start + 10,
				cbp_0.first, cbp_0.second);
				cbp_0.first = a_result.first;
				cbp_0.second = a_result.second;
					if(debug) {
						cout << "cbp_test-10-1-a: "<< static_cast<int64_t>(end - ref_is[rg]["isu"] + mlen) << "," << (start + 10) << "\n";
						cout << "cbp_test-10-1-b: "<< cbp_0.first << "," << cbp_0.second << "\n";
					}
			}
			if (1 == mstrand) {
				pair<int64_t, int64_t> a_result = overlap(mend - 10,
				mstart + ref_is[rg]["isu"] - len, cbp_1.first,
				cbp_1.second);
				cbp_1.first = a_result.first;
				cbp_1.second = a_result.second;
					if(debug) {
						cout << "cbp_test-10-2: "<< cbp_1.first << "," << cbp_1.second << "\n";
					}
			} else {
				pair<int64_t, int64_t> a_result = overlap(
				mend - ref_is[rg]["isu"] + len, mstart + 10,
				cbp_1.first, cbp_1.second);
				cbp_1.first = a_result.first;
				cbp_1.second = a_result.second;
				if(debug) {
					cout << "cbp_test-10-3: "<< cbp_1.first << "," << cbp_1.second << "\n";
				}
			}
			orientation[lastpid][lastsid].strand = strand;
			orientation[lastpid][lastsid].mate_strand = mstrand;
			orientation[lastpid][lastsid].ref_id = seqid;
			orientation[lastpid][lastsid].mate_ref_id = mseqid;
		}
//				auto starts_min_max = minmax_element(starts.begin(), starts.end());
//				auto mstarts_min_max = minmax_element(mstarts.begin(), mstarts.end());

		sort(starts.begin(), starts.end());
		sort(mstarts.begin(), mstarts.end());
//				sort(isizes.begin(), isizes.end());
		if(debug) {
			cout << (boost::format("2-0: |%s|%s\n2-1: %s|%s|%s|%s\n")
			% lastpid % ref_is["isu"]["selected"] % starts[0] % starts.back() % mstarts[0] % mstarts.back()).str();
		}
		if (!cluster.empty()) {
			if(debug) {
				cout << "the cluster data are inserted: " << lastpid << "/" << lastsid << "\n";
			}
			starts_map[lastpid][lastsid].first = starts[0];
			starts_map[lastpid][lastsid].second = starts.back();
			mstarts_map[lastpid][lastsid].first = mstarts[0];
			mstarts_map[lastpid][lastsid].second = mstarts.back();
			cbps_map[lastpid][lastsid].resize(4);
			cbps_map[lastpid][lastsid][0] = cbp_0.first;
			cbps_map[lastpid][lastsid][1] = cbp_0.second;
			cbps_map[lastpid][lastsid][2] = cbp_1.first;
			cbps_map[lastpid][lastsid][3] = cbp_1.second;
		}
		if(debug) {
			cout << (boost::format("2-2: strand: %d, mstrand: %d\n") % orientation[lastpid][lastsid].strand % orientation[lastpid][lastsid].mate_strand).str();
		}
		// strand type 0
		if (1 == orientation[lastpid][lastsid].strand && -1 == orientation[lastpid][lastsid].mate_strand) {
			if(debug) {
				cout << "strand type : 0\n";
				cout << "(" << starts_map[lastpid][lastsid].first << "," << starts_map[lastpid][lastsid].second
						<< "), (" << mstarts_map[lastpid][lastsid].first << "," << mstarts_map[lastpid][lastsid].second
						<< ")\n";
			}
			int64_t size = -1;
			if (orientation[lastpid][lastsid].ref_id == orientation[lastpid][lastsid].mate_ref_id) {
				size = mstarts_map[lastpid][lastsid].first
				- starts_map[lastpid][lastsid].second;
			}

			string p_s_id = lastpid + "_" + lastsid;
			cluster_type[orientation[lastpid][lastsid].ref_id][orientation[lastpid][lastsid].mate_ref_id][0].push_back(
			p_s_id);
			if(verbose) {
				cout << (boost::format("2-3: size: %d, p_s_id: %s, seq_id: %s, mseq_id: %s\n")
				% size % p_s_id % orientation[lastpid][lastsid].ref_id % orientation[lastpid][lastsid].mate_ref_id).str();
				cout << (boost::format("2-4: starts[0]: %d, starts[-1]: %s, mstarts[0]: %s, mstarts[-1]: %s, cbp_0-0: %d, cbp_0-1: %d, cbp_1-0: %d, cbp_1-1: %d\n")
				% starts[0] % starts.back() % starts[0] % starts.back() % cbp_0.first % cbp_0.second % cbp_1.first % cbp_1.second).str();
			}
			if (size > 0 || orientation[lastpid][lastsid].ref_id != orientation[lastpid][lastsid].mate_ref_id) {
				if(-1 == size) {
					out << (boost::format(
					"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n")
					% lastpid % lastsid % support[lastpid]
					% orientation[lastpid][lastsid].ref_id
					% starts_map[lastpid][lastsid].second
					% orientation[lastpid][lastsid].strand
					% orientation[lastpid][lastsid].mate_ref_id
					% mstarts_map[lastpid][lastsid].first
					% orientation[lastpid][lastsid].mate_strand).str();
				} else {
					out << (boost::format(
					"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
					% lastpid % lastsid % support[lastpid]
					% orientation[lastpid][lastsid].ref_id
					% starts_map[lastpid][lastsid].second
					% orientation[lastpid][lastsid].strand
					% orientation[lastpid][lastsid].mate_ref_id
					% mstarts_map[lastpid][lastsid].first
					% orientation[lastpid][lastsid].mate_strand
					% size).str();
				}
//					DISCCL <<
//"$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
			}
		}
		// strand type 1
		if (1 == orientation[lastpid][lastsid].strand && 1 == orientation[lastpid][lastsid].mate_strand) {
			if(debug) {
				cout << "strand type : 1\n";
			}
			int64_t size = -1;
			if (orientation[lastpid][lastsid].ref_id == orientation[lastpid][lastsid].mate_ref_id) {
				size = mstarts_map[lastpid][lastsid].second - starts_map[lastpid][lastsid].second;
			}
			string p_s_id = lastpid + "_" + lastsid;
			cluster_type[orientation[lastpid][lastsid].ref_id][orientation[lastpid][lastsid].mate_ref_id][1].push_back(p_s_id);
			if(-1 == size) {
				out << (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n")
				% lastpid % lastsid % support[lastpid]
				% orientation[lastpid][lastsid].ref_id
				% starts_map[lastpid][lastsid].second
				% orientation[lastpid][lastsid].strand
				% orientation[lastpid][lastsid].mate_ref_id
				% mstarts_map[lastpid][lastsid].second
				% orientation[lastpid][lastsid].mate_strand).str();
			} else {
				out << (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% lastpid % lastsid % support[lastpid]
				% orientation[lastpid][lastsid].ref_id
				% starts_map[lastpid][lastsid].second
				% orientation[lastpid][lastsid].strand
				% orientation[lastpid][lastsid].mate_ref_id
				% mstarts_map[lastpid][lastsid].second
				% orientation[lastpid][lastsid].mate_strand % size).str();
			}

//					DISCCL << "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		// strand type 2
		if (-1 == orientation[lastpid][lastsid].strand && -1 == orientation[lastpid][lastsid].mate_strand) {
			if(debug) {
				cout << "strand type : 2\n";
			}
			int64_t size = -1;
			if (orientation[lastpid][lastsid].ref_id == orientation[lastpid][lastsid].mate_ref_id) {
				size = mstarts_map[lastpid][lastsid].first - starts_map[lastpid][lastsid].first;
			}
			string p_s_id = lastpid + "_" + lastsid;
			cluster_type[orientation[lastpid][lastsid].ref_id][orientation[lastpid][lastsid].mate_ref_id][2].push_back(p_s_id);
			if(-1 == size) {
				out << (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n")
				% lastpid % lastsid % support[lastpid]
				% orientation[lastpid][lastsid].ref_id
				% starts_map[lastpid][lastsid].first
				% orientation[lastpid][lastsid].strand
				% orientation[lastpid][lastsid].mate_ref_id
				% mstarts_map[lastpid][lastsid].first
				% orientation[lastpid][lastsid].mate_strand).str();
			} else {
				out << (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n")
				% lastpid % lastsid % support[lastpid]
				% orientation[lastpid][lastsid].ref_id
				% starts_map[lastpid][lastsid].first
				% orientation[lastpid][lastsid].strand
				% orientation[lastpid][lastsid].mate_ref_id
				% mstarts_map[lastpid][lastsid].first
				% orientation[lastpid][lastsid].mate_strand % size).str();
			}
//					DISCCL << "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		// strand type 3
		if (-1 == orientation[lastpid][lastsid].strand && 1 == orientation[lastpid][lastsid].mate_strand) {
			int64_t size = -1;
			if (orientation[lastpid][lastsid].ref_id == orientation[lastpid][lastsid].mate_ref_id) {
				size = mstarts_map[lastpid][lastsid].second
				- starts_map[lastpid][lastsid].first;
			}
			string p_s_id = lastpid + "_" + lastsid;
			if(debug) {
				cout << "strand type : 3: " << p_s_id << ":" << orientation[lastpid][lastsid].ref_id << "/" << orientation[lastpid][lastsid].mate_ref_id << "\n";
			}
			cluster_type[orientation[lastpid][lastsid].ref_id][orientation[lastpid][lastsid].mate_ref_id][3].push_back(p_s_id);
			if(-1 == size) {
				out << (boost::format(
					"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n")
					% lastpid % lastsid % support[lastpid]
					% orientation[lastpid][lastsid].ref_id
					% starts_map[lastpid][lastsid].first
					% orientation[lastpid][lastsid].strand
					% orientation[lastpid][lastsid].mate_ref_id
					% mstarts_map[lastpid][lastsid].second
					% orientation[lastpid][lastsid].mate_strand).str();
			} else {
				out << (boost::format(
				"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
				% lastpid % lastsid % support[lastpid]
				% orientation[lastpid][lastsid].ref_id
				% starts_map[lastpid][lastsid].first
				% orientation[lastpid][lastsid].strand
				% orientation[lastpid][lastsid].mate_ref_id
				% mstarts_map[lastpid][lastsid].second
				% orientation[lastpid][lastsid].mate_strand % size).str();
			}
//					DISCCL << "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		if (l < options.support_mps || j < options.support_mpf) {
			if(debug) {
				cout << "erased: " << j << "/" << support[lastpid] << "\n";
			}
			auto s_map_itr = starts_map.find(lastpid);
			starts_map.erase(s_map_itr);
			auto ms_map_itr = mstarts_map.find(lastpid);
			mstarts_map.erase(ms_map_itr);
			auto c_map_itr = cbps_map.find(lastpid);
			cbps_map.erase(c_map_itr);
			auto o_itr = orientation.find(lastpid);
			orientation.erase(o_itr);
			auto s_itr = support.find(lastpid);
			support.erase(s_itr);
			auto s_f_itr = support_f.find(lastpid);
			support_f.erase(s_f_itr);
		}
	}



//TODO: if I change change the statement, if (-1 == b1 && -1 == b2) to if (0 == b1 && 0 == b2), the algorithm finds more SV events.
inline pair<int64_t, int64_t> MatePairDiscordantCaller::overlap(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {

//	if (-1 == b1 && -1 == b2) {
	if ((-1 == b1 && -1 == b2) || (0 == b1 && 0 == b2)) {
		a1 = max(static_cast<int64_t>(0), a1);
		a2 = max(static_cast<int64_t>(0), a2);
		return make_pair(a1, a2);
	}
//		( $a1, $a2 ) = sort { $a <=> $b } ( $a1, $a2 );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( $b1, $b2 ) = sort { $a <=> $b } ( $b1, $b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		int64_t a = (a1 > b1) ? a1 : b1;
		int64_t b = (a2 > b2) ? b2 : a2;
		return make_pair(a, b);
	} else {
		// the original version returns (0, 0)
//		return make_pair(b1, b2);
		return make_pair(0, 0);
	}
}


inline bool MatePairDiscordantCaller::covered(int64_t a1, int64_t a2, int64_t b1, int64_t b2) {
//		return 0 if ( $a1 =~ /\D/ or $a2 =~ /\D/ or $b1 =~ /\D/ or $b2 =~ /\D/ );
	if (a1 > a2) {
		swap(a1, a2);
	}
	if (b1 > b2) {
		swap(b1, b2);
	}
//		( $a1, $a2 ) = sort { $a <=> $b } ( $a1, $a2 );
//		( $b1, $b2 ) = sort { $a <=> $b } ( $b1, $b2 );
	if ((a1 >= b1 && a1 <= b2) || (a2 >= b1 && a2 <= b2) || (b1 >= a1 && b1 <= a2) || (b2 >= a1 && b2 <= a2)) {
		return true;
	}
	return false;
}


} /* namespace meerkat */

#endif /* MEERKAT_MATEPAIRDISCORDANTCALLER_HPP_ */
