/*
 * ClusterSelector.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: Euncheon Lim, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
 *      Email: euncheon_lim@hms.harvard.edu
 *
 *  The original version is written by Lovelace J. Luquette
 */

#include "ClusterSelector.hpp"

namespace meerkat {

ClusterSelector::ClusterSelector() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

ClusterSelector::~ClusterSelector() {
}
void ClusterSelector::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	black_listed = set<string>(options.rg_blacklist.begin(), options.rg_blacklist.end());
}
void ClusterSelector::select_cluster() {
//	if(boost::filesystem::exists(options.output_filename)) {
//		return;
//	}

	int64_t nreads = 0;
	int64_t cid;
	vector<Cluster*> clusters;
	unordered_map<string, ReadInfo> readinfo;

	string input_file(options.input_filename);
	cerr << "Reading input file " << input_file << "\n";

	int64_t n_clusters = 0;
	string line;
	{
		ifstream ifs(input_file, ios::binary);
		const char * delim_tab = "\t";
		vector<string> data;
		while(getline(ifs, line, '\n')) {
			castle::StringUtils::c_string_multi_split(line, delim_tab, data);
			if(data.size() < 1) {
				continue;
			}
			n_clusters = max(n_clusters, boost::lexical_cast<int64_t>(data[0]));
		}
//		string line = castle::IOUtils::get_last_line(input_file);
//		istringstream linestream(line);
//		linestream >> n_clusters;
	}

	++n_clusters;
	clusters.resize(n_clusters);
	for(int64_t cid = 0; cid < n_clusters; ++cid) {
		clusters[cid] = new Cluster(cid);
	}
	ifstream ifs(input_file, ios::binary);
	while (getline(ifs, line, '\n')) {
		++nreads;
		istringstream linestream(line);
		linestream >> cid;
//		if (cid >= clusters.size()) {
////			clusters.push_back(new Cluster(cid));
////			cout << "here-2:" << clusters.size() << "\n";
//		}

		Readpair rp(linestream);
//		cout << "here-3:" << line << "\n";
		/* the memory for each readpair is _owned by the cluster_ */
		ReadInfo::value_type elem(clusters[cid],
				clusters[cid]->add(rp));
		readinfo[rp.get_name()].push_back(elem);
	}

	/* Build the container */
	cout << "Building container [" << clusters.size() << " clusters, " << nreads
				<< " reads (" << readinfo.size() << " unique reads)]\n";
	cerr << "Building container [" << clusters.size() << " clusters, " << nreads
			<< " reads (" << readinfo.size() << " unique reads)]\n";
	ClusterContainer con;
	vector<Cluster*>::iterator it = clusters.begin();
	for (; it != clusters.end(); ++it) {
#ifdef DEBUG
		cout << "DEBUG: " << *it << "\n";
#endif
		con.add(*it);
	}

	/* Select the best cluster, write it, fix the container */
	cerr << "Starting selection process..\n";
	int64_t cluster_id = 0;
	ofstream out(options.output_filename, ios::binary);
	while (con.get_size() > 0) {
		int64_t secondary_id = 0;
		list<ReadInfo*> reads;
		list<Cluster*> cls = con.best();
		for (list<Cluster*>::iterator it = cls.begin(); it != cls.end(); ++it) {
			list<ReadInfo*> new_reads = select2(con, readinfo, **it, out,
					cluster_id, secondary_id++);
			reads.splice(reads.begin(), new_reads);
		}
		con.fix(reads);
		++cluster_id;
	}

	for(int64_t cid = 0; cid < n_clusters; ++cid) {
		Cluster* a_ptr = clusters[cid];
		if(a_ptr) {
			delete a_ptr;
		}
	}
	cerr << "Finished all the selections\n";
}

list<ReadInfo*> ClusterSelector::select2(ClusterContainer &con, unordered_map<string, ReadInfo> &readinfo, Cluster &cl, ostream &ostr, int64_t id, int64_t sid) {
	list<Readpair> &reads = cl.get_reads();

	/* Build a list of readinfo pointers and print out the reads */
	/* XXX: this read list should really be owned by ClusterContainer */
	/* rather than passed by value back to main() */
	list<ReadInfo*> R;
	for (list<Readpair>::iterator rp = reads.begin(); rp != reads.end(); ++rp) {
		ReadInfo &ri = readinfo[rp->get_name()];
		if (!ri.is_used()) {
			R.push_back(&ri);
		}
		ri.set_used();
		ostr << id << "\t" << sid << "\t" << cl.get_id() << "\t" << cl.get_weight() << "\t" << cl.get_mismatches() << "\t" << *rp << "\n";

		/* Mark each cluster containing this read dirty */
		for (ReadInfo::iterator it = ri.begin(); it != ri.end(); ++it) {
			con.set_dirty((*it).first);
		}
	}

	return R;
}
} /* namespace meerkat */
