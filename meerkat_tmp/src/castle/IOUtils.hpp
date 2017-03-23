/*
 * IOUtils.hpp
 *
 *  Created on: Feb 6, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef IOUTILS_HPP_
#define IOUTILS_HPP_
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/range.hpp>
#include "BinaryEncoder.hpp"
#include "StringUtils.hpp"
#include "ParallelRunner.hpp"
//#include "../ropebwt/kseq.h"
#include <zlib.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
//#include "../ropebwt/gzstream.h"
//KSEQ_INIT(gzFile, gzread)

namespace castle {
using namespace std;
namespace bfs = boost::filesystem;
typedef enum {
	FASTQ, FASTA, NON_SEQUENCING
} Type_Of_Format;
typedef enum {
	SINGLE_ENDED, PAIRED_ENDED, DONT_KNOW
} Type_Of_Read;
typedef enum {
	SANGER, ILLUMINA, UNKNOWN
} Quality_Standard;
typedef struct {
	string path_1;
	string path_2;
	int64_t size_1;
	int64_t size_2;
} Paired_End_Info;
typedef struct {
	string read;
	string qual;
	Quality_Standard standard;
	uint8_t min_quality;
	uint8_t max_quality;
} FASTQ_Info;
typedef string::iterator str_iter;
typedef boost::iterator_range<str_iter> string_view;
class IOUtils {
public:
	IOUtils();
	~IOUtils();

	static string expand_home(const string& a_path);
	static Type_Of_Format get_file_format(const string& a_path);
	static void save_FASTQ_with_trimming(const string& a_path);
	static void selective_search(const boost::filesystem::path& root, const string& regex_text, vector<string>& ret);

	// Do not call this function within a thread since this function is already parallelized.
	static void find_skip_points(vector<uint64_t>& skip_points, const string& a_path, const uint64_t chunk_size, const uint64_t max_index,
	        const uint64_t n_blocks, const uint64_t n_cores);
	static void find_fastq_skip_points(vector<uint64_t>& skip_points, const string& a_path, const uint64_t chunk_size, const uint64_t max_index,
	        const uint64_t n_blocks, const uint64_t n_cores);
	static void find_fasta_skip_points(vector<uint64_t>& skip_points, const string& a_path, const uint64_t chunk_size, const uint64_t max_index,
		        const uint64_t n_blocks, const uint64_t n_cores);
	static void replace_unknown_bases(string& a_read);
	static int64_t get_trimmed_indexes(int64_t& base_start_id, int64_t& base_end_id, const string& a_read);
	static void trim_and_replace_unknown_bases(string& a_sequence, string& a_quality);
	static string read_fully(const string& a_path);
	static FASTQ_Info read_fastq(const string& a_path, bool verbose);
	static Quality_Standard get_quality_standard(const string& a_path);
	static FASTQ_Info get_min_max_quality(const string& a_content);
	static string read_fasta(const string& a_path);
	static size_t get_file_size(const string& a_path);
	static void remove_files(const vector<string>& paths, const int32_t n_cores);
	static void remove(const string& a_path);
	static void remove_if_empty(const string& a_path);
	static uint64_t get_number_of_lines(const string& a_path);
	static string get_last_line(const string& a_path);
	static uint64_t read_to_buffer(istream& input_stream, vector<char>& buffer);
	static uint64_t count_line_feeds(const vector<char> & buffer, int64_t size_of_read);
	template<typename C> static void string_multi_split(string& str, char const* delimiters, C& ret_array);
	static int64_t get_next_start_pos(uint32_t id, int64_t chunk_size, int32_t delta, int64_t already_processed_size);
	static int64_t get_next_end_pos(uint32_t id, uint32_t max_id, int64_t chunk_size, int64_t max_size, int64_t already_processed_size);
	static void get_recurisve_paths(const bfs::path& root, const string& ext, vector<bfs::path>& ret);
	static void get_recurisve_multiple_paths(const bfs::path& root, const string& exts, vector<bfs::path>& ret);
	static bool isGzip(const std::string& filename);
	static void plain_file_merge(const string& out_filename, const vector<string>& in_filenames, const int32_t n_cores, bool remove = true);
	static void plain_file_merge_serial(const string& out_filename, const vector<string>& in_filenames, const int32_t n_cores, bool remove = true);
	static void plain_file_compress_and_merge(const string& out_filename, const vector<string>& in_filenames, const int32_t n_cores, bool remove = true);
	// Wrapper function for opening a reader of compressed or uncompressed file
//	static std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
//	static std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);

	static void liftrlimit();
//	static uint64_t kputsn(const char *p, uint64_t l, kstring_t *s);

public:
	static const string illumina_asciis;
	static const string sanger_asciis;
	static const string PBWT_EXT;
	static const string PRBWT_EXT;
	static const string BWT_EXT;
	static const string RBWT_EXT;
	static const string SAI_EXT;
	static const string RSAI_EXT;
	static const string GZIP_EXT;
	static const string POP_PREFIX;
	// SAI bit vector supporting upto 65535 datasets
	static const string BSAI_EXT;
	// SAI bit vector supporting upto 65535 datasets
	static const string RBSAI_EXT;

	static const uint64_t KILO = 1024;
	static const uint64_t MEGA = KILO * KILO;
	static const uint64_t GIGA = MEGA * KILO;
	static const uint64_t READ_BUFFER_SIZE = 8 * KILO;
	static const uint64_t PARTIAL_READ_BUFFER_SIZE = 4 * MEGA;
	static uint64_t BLOCK_SIZE;
	static const uint64_t WRITE_BUFFER_SIZE = 32 * MEGA;
	static const uint64_t SMALL_WRITE_BUFFER_SIZE = 64 * KILO;
	static const uint64_t MEDIUM_WRITE_BUFFER_SIZE = 2 * MEGA;
	static const uint64_t MAX_FILE_SIZE_FOR_CHUNK = 16 * GIGA;
	static const uint64_t MIDDLE_BLOCK_SIZE = 64 * MEGA;
	static const uint64_t BLOCK_READ_SIZE = 16 * MEGA;
	static constexpr uint64_t MAX_OPEN_FILE = 511;
//	static constexpr uint64_t MAX_OPEN_FILE = 255;
};
/*
inline uint64_t IOUtils::kputsn(const char *p, uint64_t l, kstring_t *s) {
	if (s->l + l + 1 >= s->m) {
		char *tmp;
		s->m = s->l + l + 2;
		kroundup32(s->m);
		if ((tmp = (char*) realloc(s->s, s->m)))
			s->s = tmp;
		else {
			std::cout << "[IOUtils.kputsn] An attempt to allocate memory has failed of size " << s->m << "\n";;
			exit(1);
			return EOF;
		}
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}
*/

template<typename C> void IOUtils::string_multi_split(string& str, char const* delimiters, C& ret_array) {
	ret_array.clear();
	C output;
	bitset<255> delims;
	while (*delimiters) {
		unsigned char code = *delimiters++;
		delims[code] = true;
	}
	str_iter beg;
	bool in_token = false;
	for (str_iter it = str.begin(), end = str.end(); it != end; ++it) {
		if (delims[*it]) {
			if (in_token) {
				output.push_back(typename C::value_type(beg, it));
				in_token = false;
			}
		} else if (!in_token) {
			beg = it;
			in_token = true;
		}
	}
	if (in_token) {
		output.push_back(typename C::value_type(beg, str.end()));
	}
	output.swap(ret_array);
}
inline uint64_t IOUtils::read_to_buffer(istream& input_stream, vector<char>& buffer) {
	input_stream.read(&buffer[0], buffer.size());
	return input_stream.gcount();
}
inline uint64_t IOUtils::count_line_feeds(const vector<char> & buffer, int64_t size_of_read) {
	uint64_t n_lines = 0;
	for (int i = 0; i < size_of_read; ++i) {
		if ('\n' == buffer[i]) {
			++n_lines;
		}
	}
	return n_lines;
}
inline int64_t IOUtils::get_next_start_pos(uint32_t id, int64_t chunk_size, int32_t delta, int64_t already_processed_size) {
	return already_processed_size + max<int64_t>(0, id * chunk_size - delta);
}
inline int64_t IOUtils::get_next_end_pos(uint32_t id, uint32_t max_id, int64_t chunk_size, int64_t max_size, int64_t already_processed_size) {
	if (id == max_id)
		return max_size;
	return already_processed_size + (id + 1) * chunk_size;
}
inline void IOUtils::get_recurisve_paths(const bfs::path& root, const string& ext, vector<bfs::path>& ret) {
	if (!bfs::exists(root) || !bfs::is_directory(root))
			return;

	bfs::recursive_directory_iterator it(root);
	bfs::recursive_directory_iterator endit;

	while (it != endit) {
		if (bfs::is_regular_file(*it) && it->path().extension() == ext) {
			ret.push_back(it->path());
		}
		++it;
	}
}
inline void IOUtils::get_recurisve_multiple_paths(const bfs::path& root, const string& exts, vector<bfs::path>& ret) {
	if (!bfs::exists(root) || !bfs::is_directory(root))
			return;

	bfs::recursive_directory_iterator it(root);
	bfs::recursive_directory_iterator endit;
	vector<string> ext_arr;
	StringUtils::c_string_multi_split(exts, "|", ext_arr);
	while (it != endit) {
		if (bfs::is_regular_file(*it)) {
			auto an_extension = it->path().extension();
			for(uint64_t ext_id = 0; ext_id < ext_arr.size(); ++ext_id) {
				if(an_extension == ext_arr[ext_id]) {
					ret.push_back(it->path());
				}
			}
		}
		++it;
	}
}
} /* namespace castle */
#endif /* IOUTILS_HPP_ */
