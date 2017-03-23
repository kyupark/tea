/*
 * BWACaller.cpp
 *
 *  Created on: Aug 6, 2016
 *      Author: el174
 */

#include "BWACaller.hpp"

namespace meerkat {

BWACaller::BWACaller() :
		n_blocks_1(0), n_blocks_2(0) {
	castle::TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

BWACaller::~BWACaller() {
}
void BWACaller::set_n_cores(const int32_t n_core) {
	n_cores = n_core;
}

void BWACaller::align_single_reads(const string& sam_path, const string& reference_path, const string& aln_param, const string& samse_param, const string& file_path) {
	if(0 == castle::IOUtils::get_file_size(file_path)) {
		return;
	}
	string fq_sai = file_path + ".sai";
	string bwa_err = file_path + ".bwa.err";
	string bwa_aln_cmd = (boost::format("bwa aln %s -t %d %s %s >%s 2>%s") % aln_param % n_cores % reference_path % file_path % fq_sai % bwa_err).str();
	string bwa_cmd_samse = (boost::format("bwa samse %s -f %s %s %s %s 2>>%s") % samse_param % sam_path % reference_path % fq_sai % file_path % bwa_err).str();
	system(bwa_aln_cmd.c_str());
	system(bwa_cmd_samse.c_str());
}
int64_t BWACaller::split_FASTQ(const string& file_path_1, const string& file_path_2, int64_t n_lines) {
	int64_t n_file_lines_1 = castle::IOUtils::get_number_of_lines(file_path_1);
	int64_t n_file_lines_2 = castle::IOUtils::get_number_of_lines(file_path_2);
	if (n_file_lines_1 != n_file_lines_2) {
		cerr << "[BWACaller.split_FASTQ] Your input FASTQ files do not have the same number of lines\n";
		exit(1);
	}
	int64_t n_ideal_lines = n_file_lines_1 / (double) n_cores / (double) 4;
	if (n_ideal_lines < n_lines) {
		n_lines = n_ideal_lines;
	}
	cout << (boost::format("[BWACaller.split_FASTQ] # n_lines: %d\n") % n_lines).str();
	vector<function<void()> > tasks;
	vector<int64_t> the_pos_1;
	vector<int64_t> the_pos_2;
	tasks.push_back([&, n_lines] {
		string line;
		ifstream in(file_path_1, ios::binary);
		int64_t cur_lines = 0;
		int64_t cur_pos = 0;
		the_pos_1.push_back(cur_pos);
		// each read entry has 4 lines. n_lines * 4
			const int64_t fq_n_lines = n_lines * 4;
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				if(cur_lines >= fq_n_lines) {
					the_pos_1.push_back(cur_pos);
					cur_lines = 0;
				}
			}
			the_pos_1.push_back(cur_pos + 10);
		});
	tasks.push_back([&, n_lines] {
		string line;
		ifstream in(file_path_2, ios::binary);
		int64_t cur_lines = 0;
		int64_t cur_pos = 0;
		the_pos_2.push_back(cur_pos);
		// each read entry has 4 lines. n_lines * 4
			const int64_t fq_n_lines = n_lines * 4;
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				if(cur_lines >= fq_n_lines) {
					the_pos_2.push_back(cur_pos);
					cur_lines = 0;
				}
			}
			the_pos_2.push_back(cur_pos + 10);
		});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	n_blocks_1 = the_pos_1.size() - 1;

	// split actual files
	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([&, block_id, n_lines] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = the_pos_1[block_id];
			int64_t next_pos = the_pos_1[block_id + 1];
			int64_t cur_lines = 0;
			string line;
			string out_file_path = file_path_1 + "." + str_block_id;
			if(castle::IOUtils::get_file_size(out_file_path) > 0) {
				return;
			}
			ifstream in(file_path_1, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			ofstream out(out_file_path, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				out << line << "\n";
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	n_blocks_2 = the_pos_2.size() - 1;
	for (int64_t block_id = 0; block_id < n_blocks_2; ++block_id) {
		tasks.push_back([&, block_id, n_lines] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = the_pos_2[block_id];
			int64_t next_pos = the_pos_2[block_id + 1];
			int64_t cur_lines = 0;
			string line;
			string out_file_path = file_path_2 + "." + str_block_id;
			if(castle::IOUtils::get_file_size(out_file_path) > 0) {
				return;
			}
			ifstream in(file_path_2, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			ofstream out(out_file_path, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				out << line << "\n";
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	return n_blocks_1;
}

int64_t BWACaller::split_FASTQ_alt(const string& file_path_1, const string& file_path_2, int64_t n_lines) {
	int64_t n_file_lines_1 = castle::IOUtils::get_number_of_lines(file_path_1);
	int64_t n_file_lines_2 = castle::IOUtils::get_number_of_lines(file_path_2);
	if (n_file_lines_1 != n_file_lines_2) {
		cerr << "[BWACaller.split_FASTQ] Your input FASTQ files do not have the same number of lines\n";
		exit(1);
	}
	cout << (boost::format("[BWACaller.split_FASTQ] # n_lines: %d\n") % n_lines).str();
	vector<function<void()> > tasks;
	vector<int64_t> the_pos_1;
	vector<int64_t> the_pos_2;
	tasks.push_back([&, n_lines] {
		string line;
		ifstream in(file_path_1, ios::binary);
		int64_t cur_lines = 0;
		int64_t cur_pos = 0;
		the_pos_1.push_back(cur_pos);
		// each read entry has 4 lines. n_lines * 4
			const int64_t fq_n_lines = n_lines * 4;
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				if(cur_lines >= fq_n_lines) {
					the_pos_1.push_back(cur_pos);
					cur_lines = 0;
				}
			}
			the_pos_1.push_back(cur_pos + 10);
		});
	tasks.push_back([&, n_lines] {
		string line;
		ifstream in(file_path_2, ios::binary);
		int64_t cur_lines = 0;
		int64_t cur_pos = 0;
		the_pos_2.push_back(cur_pos);
		// each read entry has 4 lines. n_lines * 4
			const int64_t fq_n_lines = n_lines * 4;
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				if(cur_lines >= fq_n_lines) {
					the_pos_2.push_back(cur_pos);
					cur_lines = 0;
				}
			}
			the_pos_2.push_back(cur_pos + 10);
		});
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	n_blocks_1 = the_pos_1.size() - 1;

	// split actual files
	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([&, block_id, n_lines] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = the_pos_1[block_id];
			int64_t next_pos = the_pos_1[block_id + 1];
			int64_t cur_lines = 0;
			string line;
			string out_file_path = file_path_1 + "." + str_block_id;
			if(castle::IOUtils::get_file_size(out_file_path) > 0) {
				return;
			}
			ifstream in(file_path_1, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			ofstream out(out_file_path, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				out << line << "\n";
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	n_blocks_2 = the_pos_2.size() - 1;
	for (int64_t block_id = 0; block_id < n_blocks_2; ++block_id) {
		tasks.push_back([&, block_id, n_lines] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = the_pos_2[block_id];
			int64_t next_pos = the_pos_2[block_id + 1];
			int64_t cur_lines = 0;
			string line;
			string out_file_path = file_path_2 + "." + str_block_id;
			if(castle::IOUtils::get_file_size(out_file_path) > 0) {
				return;
			}
			ifstream in(file_path_2, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			ofstream out(out_file_path, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				out << line << "\n";
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	return n_blocks_1;
}


int64_t BWACaller::split_FASTQ_single(const string& file_path_1, int64_t n_lines) {
	cout << (boost::format("[BWACaller.split_FASTQ] # n_lines: %d\n") % n_lines).str();
	vector<function<void()> > tasks;
	vector<int64_t> the_pos_1;
	vector<int64_t> the_pos_2;
	string line;
	ifstream in(file_path_1, ios::binary);
	int64_t cur_lines = 0;
	int64_t cur_pos = 0;
	the_pos_1.push_back(cur_pos);
// each read entry has 4 lines. n_lines * 4
	const int64_t fq_n_lines = n_lines * 4;
	while(getline(in, line, '\n')) {
		cur_pos += line.size() + 1;
		++cur_lines;
		if(cur_lines >= fq_n_lines) {
			the_pos_1.push_back(cur_pos);
			cur_lines = 0;
		}
	}
	the_pos_1.push_back(cur_pos + 10);
	n_blocks_1 = the_pos_1.size() - 1;

	// split actual files
	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([&, block_id, n_lines] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			int64_t cur_pos = the_pos_1[block_id];
			int64_t next_pos = the_pos_1[block_id + 1];
			int64_t cur_lines = 0;
			string line;
			string out_file_path = file_path_1 + "." + str_block_id;
			if(castle::IOUtils::get_file_size(out_file_path) > 0) {
				return;
			}
			ifstream in(file_path_1, ios::binary);
			if(0 != cur_pos) {
				in.seekg(cur_pos, ios::beg);
			}
			ofstream out(out_file_path, ios::binary);
			while(getline(in, line, '\n')) {
				cur_pos += line.size() + 1;
				++cur_lines;
				out << line << "\n";
				if(cur_pos >= next_pos) {
					break;
				}
			}
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	return n_blocks_1;
}

void BWACaller::align_reads(const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2) {
	castle::TimeChecker checker;
	checker.setTarget("BWACaller.align_reads");
	checker.start();
	vector<function<void()> > tasks;
	string reference_fai = reference_path + ".fai";
	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([&, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string fq_1 = file_path_1 + "." + str_block_id;
			string fq_2 = file_path_2 + "." + str_block_id;
			string fq_sai_1 = file_path_1 + ".sai." + str_block_id;
			string fq_sai_2 = file_path_2 + ".sai." + str_block_id;
			string fq_sam = sam_path + "." + str_block_id;
			string bwa_err = file_path_1 + ".bwa.err." + str_block_id;
			string bwa_cmd_1 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>%s") % aln_param % reference_path % fq_1 % fq_sai_1 % bwa_err).str();
			string bwa_cmd_2 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>>%s") % aln_param % reference_path % fq_2 % fq_sai_2 % bwa_err).str();
			string bwa_cmd_sampe = (boost::format("bwa sampe %s -f %s %s %s %s %s %s 2>>%s")
					% sampe_param % fq_sam % reference_path % fq_sai_1 % fq_sai_2 % fq_1 % fq_2 % bwa_err).str();
			system(bwa_cmd_1.c_str());
			system(bwa_cmd_2.c_str());
			system(bwa_cmd_sampe.c_str());
		});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void BWACaller::collect_align_tasks(vector<function<void()> >& tasks, const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2) {
	split_FASTQ(file_path_1, file_path_2);

	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([sam_path, reference_path, aln_param, sampe_param, file_path_1, file_path_2, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string reference_fai = reference_path + ".fai";
			string fq_1 = file_path_1 + "." + str_block_id;
			string fq_2 = file_path_2 + "." + str_block_id;
			string fq_sai_1 = file_path_1 + ".sai." + str_block_id;
			string fq_sai_2 = file_path_2 + ".sai." + str_block_id;
			string fq_sam = sam_path + "." + str_block_id;
			string bwa_err = file_path_1 + ".bwa.err." + str_block_id;
			if(castle::IOUtils::get_file_size(fq_sam) > 0) {
				return;
			}
			string bwa_cmd_1 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>%s") % aln_param % reference_path % fq_1 % fq_sai_1 % bwa_err).str();
			string bwa_cmd_2 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>>%s") % aln_param % reference_path % fq_2 % fq_sai_2 % bwa_err).str();
			string bwa_cmd_sampe = (boost::format("bwa sampe %s -f %s %s %s %s %s %s 2>>%s")
					% sampe_param % fq_sam % reference_path % fq_sai_1 % fq_sai_2 % fq_1 % fq_2 % bwa_err).str();
			system(bwa_cmd_1.c_str());
			system(bwa_cmd_2.c_str());
			system(bwa_cmd_sampe.c_str());
		});
	}
}

void BWACaller::collect_align_tasks_alt(vector<function<void()> >& tasks, const string& sam_path, const string& reference_path, const string& aln_param, const string& sampe_param, const string& file_path_1, const string& file_path_2, const int64_t n_blocks_1) {
	for (int64_t block_id = 0; block_id < n_blocks_1; ++block_id) {
		tasks.push_back([sam_path, reference_path, aln_param, sampe_param, file_path_1, file_path_2, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string reference_fai = reference_path + ".fai";
			string fq_1 = file_path_1 + "." + str_block_id;
			string fq_2 = file_path_2 + "." + str_block_id;
			string fq_sai_1 = file_path_1 + ".sai." + str_block_id;
			string fq_sai_2 = file_path_2 + ".sai." + str_block_id;
			string fq_sam = sam_path + "." + str_block_id;
			string bwa_err = file_path_1 + ".bwa.err." + str_block_id;
			if(castle::IOUtils::get_file_size(fq_sam) > 0) {
				return;
			}
			string bwa_cmd_1 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>%s") % aln_param % reference_path % fq_1 % fq_sai_1 % bwa_err).str();
			string bwa_cmd_2 = (boost::format("bwa aln %s -t 1 %s %s >%s 2>>%s") % aln_param % reference_path % fq_2 % fq_sai_2 % bwa_err).str();
			string bwa_cmd_sampe = (boost::format("bwa sampe %s -f %s %s %s %s %s %s 2>>%s")
					% sampe_param % fq_sam % reference_path % fq_sai_1 % fq_sai_2 % fq_1 % fq_2 % bwa_err).str();
			system(bwa_cmd_1.c_str());
			system(bwa_cmd_2.c_str());
			system(bwa_cmd_sampe.c_str());
		});
	}
}

void BWACaller::collect_single_align_tasks(vector<function<void()> >& tasks, const int64_t n_blocks, const string& sam_path, const string& reference_path, const string& aln_param, const string& samse_param, const string& file_path) {
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([sam_path, reference_path, aln_param, samse_param, file_path, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string reference_fai = reference_path + ".fai";
			string fq = file_path + "." + str_block_id;
			if(0 == castle::IOUtils::get_file_size(fq)) {
				return;
			}
			string fq_sai = file_path + ".sai." + str_block_id;
			string fq_sam = sam_path + "." + str_block_id;
			string bwa_err = file_path + ".bwa.err." + str_block_id;
			if(castle::IOUtils::get_file_size(fq_sam) > 0) {
				return;
			}
			string bwa_aln_cmd = (boost::format("bwa aln %s -t 1 %s %s >%s 2>%s") % aln_param % reference_path % fq % fq_sai % bwa_err).str();
			string bwa_cmd_samse = (boost::format("bwa samse %s -f %s %s %s %s 2>>%s") % samse_param % fq_sam % reference_path % fq_sai % fq % bwa_err).str();
			system(bwa_aln_cmd.c_str());
			system(bwa_cmd_samse.c_str());
		});
	}
}
void BWACaller::merge_SAM(const string& file_path, const string& bam_path, const string& sam_path, int64_t n_blocks) {
	vector<function<void()> > tasks;
	vector<string> merge_file_names;
	vector<string> remove_file_names;
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		string str_block_id = boost::lexical_cast<string>(block_id);
		string fq = file_path + "." + str_block_id;
		string fq_sai = file_path + ".sai." + str_block_id;
		string fq_sam = sam_path + "." + str_block_id;
		string fq_tmp_sam = sam_path + ".tmp.sam." + str_block_id;

		merge_file_names.push_back(fq_tmp_sam);
		remove_file_names.push_back(fq);
		remove_file_names.push_back(fq_sai);
		remove_file_names.push_back(fq_sam);
		remove_file_names.push_back(fq_tmp_sam);
	}
	string first_file_name = sam_path;
	if (n_blocks > 0) {
		first_file_name = merge_file_names[0];
	}
	for (int64_t block_id = 0; block_id < n_blocks; ++block_id) {
		tasks.push_back([&, sam_path, block_id] {
			string str_block_id = boost::lexical_cast<string>(block_id);
			string fq_sam = sam_path + "." + str_block_id;
			string fq_tmp_sam = sam_path + ".tmp.sam." + str_block_id;
			bool is_first_rg = first_file_name == fq_tmp_sam;
			string line;
			ofstream out(fq_tmp_sam, ios::binary);
			ifstream in(fq_sam, ios::binary);
			// write the SAM header.
				if(is_first_rg) {
					while (getline(in, line, '\n')) {
						if(line.empty()) {
							continue;
						}
						out << line << "\n";
					}
				} else {
					while (getline(in, line, '\n')) {
						if(line.empty()) {
							continue;
						}
						if('@' == line[0]) {
							continue;
						}
						out << line << "\n";
					}
				}
			});
	}
	castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
	castle::IOUtils::plain_file_merge(sam_path, merge_file_names, n_cores, false);
	string sam_to_bam_cmd = (boost::format("samtools view -1 -Sb -@ %d -o %s %s") % n_cores % bam_path % sam_path).str();
	system(sam_to_bam_cmd.c_str());
	castle::IOUtils::remove_files(remove_file_names, n_cores);
}
} /* namespace meerkat */
