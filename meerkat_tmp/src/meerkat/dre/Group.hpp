#ifndef MEERKAT_GROUP_HPP_
#define MEERKAT_GROUP_HPP_

#include <vector>
#include <utility>
#include <map>
#include "api/BamAlignment.h"
#include "api/BamWriter.h"

#include "PCRDups.hpp"

namespace meerkat {
	using namespace std;

	class Group {
		public:
			Group();
			void add(BamTools::BamAlignment &al);

			/* Filter out PCR duplicates and write out non-dup reads */
			int64_t write(BamTools::BamWriter &writer, record &prev_dups,
				BamTools::BamWriter& dups, bool write_dup);

			int64_t write_alt(BamTools::BamWriter &writer, record &prev_dups,
							BamTools::BamWriter& dups, bool write_dup);

			/* Empty out this group */
			void clear();

			/* Should we write out the current group and make a new one? */
			bool should_write(BamTools::BamAlignment &al);
			private:
			vector<BamTools::BamAlignment> reads;
			vector<PCRDups> dup_groups;
			int n;
			long pos;
	};

}
#endif
