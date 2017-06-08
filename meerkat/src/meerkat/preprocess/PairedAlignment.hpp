/*
 * PairedAlignment.hpp
 *
 *  Created on: Jun 8, 2016
 *      Author: el174
 */

#ifndef MEERKAT_PAIREDALIGNMENT_HPP_
#define MEERKAT_PAIREDALIGNMENT_HPP_

#include "api/BamReader.h"

namespace meerkat {
	using namespace BamTools;
	struct PairedAlignment {
		BamAlignment pair_1;
		BamAlignment pair_2;
		bool should_process;
		PairedAlignment() :
				should_process(false) {
		}
		PairedAlignment(const PairedAlignment& b) {
			pair_1 = b.pair_1;
			pair_2 = b.pair_2;
			should_process = b.should_process;
		}
		PairedAlignment& operator=(const PairedAlignment& b) {
			// check for self-assignment
			if (this == &b) {
				return *this;
			}
			pair_1 = b.pair_1;
			pair_2 = b.pair_2;
			should_process = b.should_process;
			return *this;
		}
	};

	struct PairedAlignmentWithTwoIds {
		BamAlignment pair_1;
		BamAlignment pair_2;
		int64_t id_1;
		int64_t id_2;
		PairedAlignmentWithTwoIds() :
				id_1(-1), id_2(-1) {
		}
		PairedAlignmentWithTwoIds(const PairedAlignmentWithTwoIds& b) {
			pair_1 = b.pair_1;
			pair_2 = b.pair_2;
			id_1 = b.id_1;
			id_2 = b.id_2;
		}
		PairedAlignmentWithTwoIds& operator=(const PairedAlignmentWithTwoIds& b) {
			// check for self-assignment
			if (this == &b) {
				return *this;
			}
			pair_1 = b.pair_1;
			pair_2 = b.pair_2;
			id_1 = b.id_1;
			id_2 = b.id_2;
			return *this;
		}
	};
}

#endif /* MEERKAT_PAIREDALIGNMENT_HPP_ */
