/*
 * SelectionFilter.cpp
 *
 *  Created on: Jan 12, 2017
 *      Author: el174
 */

#include "SelectionFilter.hpp"

namespace meerkat {

SelectionFilter::SelectionFilter() {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

SelectionFilter::~SelectionFilter() {
}

void SelectionFilter::set_option_parser(const castle::OptionParser& the_options) {
	options = the_options;
	n_cores = options.n_cores;
}
void SelectionFilter::filter_variants() {
	string variant_file = options.prefix + ".variants";
	string variant_filtered_file = options.prefix + ".variants.filtered";
	if(!options.working_dir.empty()) {
		variant_file = options.working_prefix + ".variants";
		variant_filtered_file = options.working_prefix + ".variants.filtered";
	}
	vector<string> data;
	const char* delim_tab = "\t";
	string line;
	ifstream in(variant_file, ios::binary);
	ofstream out(variant_filtered_file, ios::binary);
	while(getline(in, line, '\n')) {
		castle::StringUtils::c_string_multi_split(line, delim_tab, data);
//		auto& type = data[0];
		auto& homology = data[data.size() - 2];
		if("1000" == homology || "1000/1000" == homology) {
			continue;
		}
		out << line << "\n";
	}
//	del_ins event: deletion part/insertion part
//	homology length is too large compared with read length -> mFP
//	homology of 1000 -> FP
//	negative insertion size -> mFP
//	single cluster support -> not reliable
//	both mapped at a repetitive region (based on BP tag) -> mFP
//	in an MDA sample, event with opposite direction -> hFP
//	in an MDA sample, inversion (++, -- signal) -> hFP (to select possible TP, the size should be larger than 1kb)
//	an MDA sample may have too many events of length smaller than 1kb -> mFP
//	del_invers (inversion size, del1_size, del2_size) with negative inversion size -> mFP
//	insod of size less than -100 -> mTP
//
//	invers_r of size larger than 1kb -> mTP
//	bp tag of SR_SR->mFP, since it mapped at unreliable regions -> mFP
//	inso of negative distance  compared with read length -> mFP
//	del_ins: the reported values are for the deletion part, the insertion source is unknown since it would be too short or repeats (ATATAT)
//
//	transl_inter: it has insufficient supports (blast) to be either del_inss or inss type.
}
} /* namespace meerkat */
