//============================================================================
// Name        : tea.cpp
// Author      : Euncheon Lim
// Version     :
// Copyright   : Park Lab @ Harvard Medical School
// Description : TEA 2.0
//============================================================================

#include "castle/TimeChecker.hpp"
#include "castle/OptionParser.hpp"
#include "sub_modules.hpp"

using namespace std;

int main(int argc, char **argv) {
	setvbuf(stdout, NULL, _IONBF, 0);
	tea::TEAOptionParser option_parser(argc, argv);
	castle::TimeChecker checker;
	checker.setTarget("TEA.main");
	checker.start();

	if("tea" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.preprocess();
		tea.clean();
	} else if("tea_v" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.preprocess_v();
		tea.run_vid();
		tea.clean();
	} else if("tea_u" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.preprocess_u();
		tea.run_uid();
		tea.clean();
	} else if("transduction" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.run_transduction();
		tea.run_transduction_contig();
		tea.clean();
	} else if("orphan" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.run_orphan();
		tea.run_orphan_contig();
		tea.clean();
	} else if("post" == option_parser.program_name) {
		tea::TEA tea;
		tea.set_option_parser(option_parser);
		tea.post_process();
	} else if("comp" == option_parser.program_name) {
		tea::ResultCompartor rc;
		rc.set_option_parser(option_parser);
		string comp_suffix = ".1";
		rc.find_overlaps(comp_suffix);
		swap(option_parser.l1, option_parser.l2);
		rc.set_option_parser(option_parser);
		comp_suffix = ".2";
		rc.find_overlaps(comp_suffix);
	}
	cout << checker;
	return 0;
}
