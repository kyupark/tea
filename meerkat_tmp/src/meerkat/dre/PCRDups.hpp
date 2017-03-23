#ifndef MEERKAT_PCRDUPS_HPP_
#define MEERKAT_PCRDUPS_HPP_

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include "api/BamAlignment.h"

namespace meerkat {
	using namespace std;
	typedef pair<long, long> coords;
	typedef pair<coords, coords> frag_coords;
	typedef map<frag_coords, string> record;

	class PCRDups {
		public:
			PCRDups(BamTools::BamAlignment &a, BamTools::BamAlignment &b) {
				dups.push_back(a);
				dups.push_back(b);
			}
			static frag_coords make_frag_coords(BamTools::BamAlignment &al);

			/* Does 'al' belong to this PCR duplicate group? */
			void add(BamTools::BamAlignment &al) {
				dups.push_back(al);
			}

			/* return a representative of the group */
			/* always exists, because PCRDups must have >= 2 elements */
			BamTools::BamAlignment &rep() {
				return dups[0];
			}

			/* For now, the "best" alignment is just the first one in the list */
			BamTools::BamAlignment best(record &prev_dups) {
				frag_coords x = make_frag_coords(rep());
				record::iterator it = prev_dups.find(x);

				/* No previous matching duplicate block, just return rep */
				if (it == prev_dups.end())
					return rep();

				/* We already saw the matching dup block, pick the mate of
				 * the same read we picked before. */
				string prev_name = it->second;
				for (size_t i = 0; i < dups.size(); ++i) {
					if (dups[i].Name == prev_name)
						return dups[i];
				}

				cout << "ERROR: could not find the best read in dup group at <"
						<< rep().RefID << ", " << rep().Position << ">" << endl;
				cout << "Previous dup block contained mate: " << it->second
						<< endl;
				cout << "This block contains dups:" << endl;
				for (size_t i = 0; i < dups.size(); ++i)
					cout << dups[i].Name << "\t" << dups[i].Position
							<< "\tmate:" << dups[i].MatePosition << endl;
//				exit(1);
				BamTools::BamAlignment empty;
				return empty;
			}

			void write_except(BamTools::BamWriter& out,
				const BamTools::BamAlignment& best) {
				for (size_t i = 0; i < dups.size(); ++i) {
					if (dups[i].Name != best.Name) {
						out.SaveSAMAlignment(dups[i]);
					}
				}
			}

			void write_except_alt(BamTools::BamWriter& out,
				const BamTools::BamAlignment& best) {
				for (size_t i = 0; i < dups.size(); ++i) {
					if (dups[i].Name != best.Name) {
						out.SaveAlignment(dups[i]);
					}
				}
			}
		private:
			vector<BamTools::BamAlignment> dups;
	};
	inline frag_coords PCRDups::make_frag_coords(BamTools::BamAlignment &al) {
		/* The positions need to be tagged with + and - for discordant */
		/* mate pairs that are on the same strand (i.e., ++ pairs and --pairs */
		coords a(al.RefID, (al.IsReverseStrand() ? -1 : 1) * al.Position),
				b(al.MateRefID,
						(al.IsMateReverseStrand() ? -1 : 1) * al.MatePosition);

		/* If the two reads are on different strands */
		if (al.IsReverseStrand() != al.IsMateReverseStrand()) {
			return al.IsReverseStrand() ? frag_coords(b, a) : frag_coords(a, b);
		}
		/* If the two reads are on the same strand... */
		else {
			return a < b ? frag_coords(a, b) : frag_coords(b, a);
		}
	}
}
#endif
