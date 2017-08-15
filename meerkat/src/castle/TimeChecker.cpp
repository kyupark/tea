/*
 * TimeChecker.cpp
 *
 *  Created on: Dec 15, 2012
 *      Author: Eun-Cheon Lim @ Max-Planck-Institute for Developmental Biology
 *      under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email: abysslover@gmail.com
 */

#include "TimeChecker.hpp"
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "boost/format.hpp"
#include <iostream>
#include <fstream>

namespace castle {
const double TimeChecker::KILO_DOUBLE = 1024;
const double TimeChecker::MEGA_DOUBLE = 1024 * KILO_DOUBLE;
const double TimeChecker::GIGA_DOUBLE = 1024 * MEGA_DOUBLE;

TimeChecker::TimeChecker() :
		start_time(0), start_cpu_time(0), original_bytes(0) {
}

TimeChecker::~TimeChecker() {
}

void TimeChecker::setTarget(std::string target) {
	this->target = target;
}

void TimeChecker::start() {
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	original_bytes = (memInfo.freeram + memInfo.freeswap) * memInfo.mem_unit;
	struct timeval tp;
	struct timezone tzp;
	::gettimeofday(&tp, &tzp);
	start_time = tp.tv_sec + tp.tv_usec * 1e-6;
	struct rusage r;
	::getrusage(RUSAGE_SELF, &r);
	start_cpu_time = r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	std::cout << (boost::format("[%s] starts at %s\n") % target % buf);
}
void TimeChecker::start_without_output() {
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	original_bytes = (memInfo.freeram + memInfo.freeswap) * memInfo.mem_unit;
	struct timeval tp;
	struct timezone tzp;
	::gettimeofday(&tp, &tzp);
	start_time = tp.tv_sec + tp.tv_usec * 1e-6;
	struct rusage r;
	::getrusage(RUSAGE_SELF, &r);
	start_cpu_time = r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
}

std::string TimeChecker::str() {
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	//	unsigned int mem_unit = memInfo.mem_unit;
	//	long total_bytes = (memInfo.totalram + memInfo.totalswap) * mem_unit;
	//	long free_bytes = (memInfo.freeram + memInfo.freeswap) * memInfo.mem_unit;
	struct timeval tp;
	struct timezone tzp;
	::gettimeofday(&tp, &tzp);
	double current_time = tp.tv_sec + tp.tv_usec * 1e-6;
	double real_time = current_time - start_time;

	struct rusage r;
	::getrusage(RUSAGE_SELF, &r);
	double cpu_time = r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	return (boost::format("[%s] Real Time: %.3f sec, CPU: %.3f sec\n") % target % real_time % cpu_time).str();
}

uint32_t TimeChecker::get_number_of_cores() {
	return sysconf(_SC_NPROCESSORS_ONLN);
}

uint32_t TimeChecker::get_delta_in_second() {
	struct timeval tp;
	struct timezone tzp;
	::gettimeofday(&tp, &tzp);
	double current_time = tp.tv_sec + tp.tv_usec * 1e-6;
	uint32_t delta = current_time - start_time;
	return delta;
}

void TimeChecker::process_mem_usage(double& vm_usage, double& resident_set) {
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime
	        >> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}
uint64_t TimeChecker::get_max_available_memory() {
	uint64_t available_memory = 0;
	std::string token;
	std::ifstream file("/proc/meminfo");
	long mem_total = 0;
	long active_memory = 0;
	while (file >> token) {
		if ("MemTotal:" == token) {
			file >> mem_total;
		} else if ("Active:" == token) {
			file >> active_memory;
			available_memory = (mem_total - active_memory) * 1024;
		}
	}
	return available_memory;
}
void TimeChecker::get_available_memory(double& available_memory) {
	std::string token;
	std::ifstream file("/proc/meminfo");
	long mem_total = 0;
	long active_memory = 0;
	while (file >> token) {
		if ("MemTotal:" == token) {
			file >> mem_total;
		} else if ("Active:" == token) {
			file >> active_memory;
			available_memory = (mem_total - active_memory) / MEGA_DOUBLE;
		}
	}
}

std::ostream& operator<<(std::ostream& out, const TimeChecker& x) {
	double available_memory = 0;
	TimeChecker::get_available_memory(available_memory);
//			struct sysinfo memInfo;
//			sysinfo(&memInfo);
//			memInfo.totalram /= x.GIGA_DOUBLE;
//			memInfo.freeram /= x.GIGA_DOUBLE;
//			memInfo.sharedram /= x.GIGA_DOUBLE;
//			memInfo.bufferram /= x.GIGA_DOUBLE;
//			memInfo.totalswap /= x.GIGA_DOUBLE;
//			memInfo.freeswap /= x.GIGA_DOUBLE;
//			memInfo.totalhigh /= x.GIGA_DOUBLE;
//			memInfo.freehigh /= x.GIGA_DOUBLE;
//
//			out << boost::format("[%s] Total: %d(GB), Free: %d(GB), Shared: %d(GB), Buffer: %d(GB), Total Swap: %d(GB), Free Swap: %d(GB)\n")
//						% x.target % memInfo.totalram % memInfo.freeram % memInfo.sharedram % memInfo.bufferram % memInfo.totalswap % memInfo.freeswap;
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	struct timeval tp;
	struct timezone tzp;
	::gettimeofday(&tp, &tzp);
	double current_time = tp.tv_sec + tp.tv_usec * 1e-6;
	double real_time = current_time - x.start_time;

	struct rusage r;
	::getrusage(RUSAGE_SELF, &r);
	double vm, rss;
	TimeChecker::process_mem_usage(vm, rss);
	vm /= x.MEGA_DOUBLE;
	rss /= x.MEGA_DOUBLE;

	double current_cpu_time = r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	double cpu_time = current_cpu_time - x.start_cpu_time;
	int numCPU = sysconf(_SC_NPROCESSORS_ONLN);

	out
	        << boost::format("[%s] Ends at %s, Real: %.3f sec, CPU: %.3f sec, RSS: %.3f(GB)/%.3f(GB), CPU: %d\n") % x.target % buf % real_time
	                % cpu_time % rss % available_memory % numCPU;

//					out << boost::format("Total memory: %.3f(MB), Available memory: %.3f (MB), Used memory: %.3f (MB)\n")
//									% (total_bytes / x.MEGA_DOUBLE) % ((total_bytes - free_bytes) / x.MEGA_DOUBLE)
//									% (free_bytes / x.MEGA_DOUBLE);
	return out;
}
double TimeChecker::cputime() {
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double TimeChecker::realtime() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
} /* namespace castle */
