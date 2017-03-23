/*
 * bam_reader_test.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: el174
 */

#include "castle/TimeChecker.hpp"
#include "castle/OptionParser.hpp"
#include "sub_modules.hpp"

using namespace std;

int main(int argc, char **argv) {
	setvbuf(stdout, NULL, _IONBF, 0);
	castle::OptionParser option_parser(argc, argv);
	castle::TimeChecker checker;
	checker.setTarget("ParallelBamReader.main");
	checker.start();

	if("preprocess" == option_parser.program_name) {
		meerkat::ParallelBamReader pbr;
		pbr.set_option_parser(option_parser);
		pbr.preprocess();
	} else if("dre" == option_parser.program_name) {
		meerkat::ParallelDiscordExtractor pde;
		pde.set_option_parser(option_parser);
		pde.extract_discordant_reads();
//		pde.extract_discordant_reads_serial();
	} else if("r_alt" == option_parser.program_name) {
		meerkat::AlternativeMapper am;
		am.set_option_parser(option_parser);
		am.create_raw_alternative_mapping_par();
	} else if("s_alt" == option_parser.program_name) {
		meerkat::AlternativeMapper am;
		am.set_option_parser(option_parser);
		am.select_alternative_mapping_alt();
	} else if("sclus" == option_parser.program_name) {
		meerkat::ClusterSelector cs;
		cs.set_option_parser(option_parser);
		cs.select_cluster();
	} else if("mpd" == option_parser.program_name) {
		meerkat::MatePairDiscordantCaller mpdc;
		mpdc.set_option_parser(option_parser);
//		mpdc.select_candidate_events_in_clusters();
		mpdc.select_candidate_events_in_clusters_alt();
	} else if("alg" == option_parser.program_name) {
		meerkat::SplitReadReAlign srra;
		srra.set_option_parser(option_parser);
//		srra.align_split_reads();
//		srra.align_split_reads_alt();
		srra.align_split_reads_par();
	} else if("srd" == option_parser.program_name) {
		meerkat::SplitReadSVCaller srd;
		srd.set_option_parser(option_parser);
		srd.call_structural_variants();
	} else if("filter" == option_parser.program_name) {
		meerkat::OutputFilter of;
		of.set_option_parser(option_parser);
		of.filter();
	} else if("blast_bp" == option_parser.program_name) {
		meerkat::BLASTBreakPoints bbp;
		bbp.set_option_parser(option_parser);
		bbp.find_break_points();
	} else if("mechanism" == option_parser.program_name) {
		meerkat::Mechanism m;
		m.set_option_parser(option_parser);
		m.call_variants();
	} else if("selection" == option_parser.program_name) {
		meerkat::SelectionFilter s;
		s.set_option_parser(option_parser);
		s.filter_variants();
	} else if("cns" == option_parser.program_name) {
		meerkat::CNVOverlapper cnvo;
		cnvo.set_option_parser(option_parser);
		cnvo.find_overlaps();
	} else if("svs" == option_parser.program_name) {
		meerkat::CNVOverlapper cnvo;
		cnvo.set_option_parser(option_parser);
		cnvo.find_sv_overlaps();
	} else if("test" == option_parser.program_name) {
//		meerkat::SplitReadSVCaller srd;
//		srd.set_option_parser(option_parser);
//		srd.call_structural_variants_test();
		meerkat::Mechanism m;
		m.set_option_parser(option_parser);
		m.call_variants_alt();
	}
	cout << checker;
	return 0;
}
