// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/loop_modeling/utilities/LoopFilter.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

using namespace std;

using protocols::filters::FilterOP;

LoopFilter::LoopFilter(FilterOP filter) {
	filter_ = filter;
}

bool LoopFilter::do_apply(Pose & pose) {
	return filter_->apply(pose);
}

} // namespace utilities
} // namespace kinematic_closure
} // namespace protocols

