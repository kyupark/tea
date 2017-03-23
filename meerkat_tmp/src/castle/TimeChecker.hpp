/*
 * TimeChecker.hpp
 *
 *  Created on: Dec 15, 2012
 *      Author: Eun-Cheon Lim @ Max-Planck-Institute for Developmental Biology
 *      under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email: abysslover@gmail.com
 */

#ifndef TIMECHECKER_HPP_
#define TIMECHECKER_HPP_
#include <boost/timer.hpp>
#include <boost/cstdint.hpp>
#include <string>
#include <ctime>

namespace castle {
		class TimeChecker {
				static const double KILO_DOUBLE;
				static const double MEGA_DOUBLE;
				static const double GIGA_DOUBLE;
				double start_time;
				double start_cpu_time;
			private:
				unsigned long long original_bytes;
			public:
				std::string target;
				TimeChecker();
				virtual ~TimeChecker();
				void setTarget(std::string target);
				void start();
				void start_without_output();
				friend std::ostream& operator<<(std::ostream& out, const TimeChecker& x);
				std::string str();
				static void process_mem_usage(double& vm_usage, double& resident_set);
				static void get_available_memory(double& available_memory);
				static uint64_t get_max_available_memory();
				uint32_t get_number_of_cores();
				uint32_t get_delta_in_second();
				static double cputime();
				static double realtime();
		};

} /* namespace castle */
#endif /* TIMECHECKER_HPP_ */
