// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/loophash/FastGapMover.fwd.hh
/// @brief Closes gaps
/// @author


#ifndef INCLUDED_protocols_loophash_FastGapMover_fwd_hh
#define INCLUDED_protocols_loophash_FastGapMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loophash {

// Forward
class FastGapMover;

// Types
typedef utility::pointer::shared_ptr< FastGapMover > FastGapMoverOP;
typedef utility::pointer::shared_ptr< FastGapMover const > FastGapMoverCOP;

}
}

#endif
