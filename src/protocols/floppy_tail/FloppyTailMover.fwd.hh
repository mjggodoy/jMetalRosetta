// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/floppy_tail/FloppyTailMover.fwd.hh
/// @brief FloppyTail forward declaration
/// @author Steven Lewis smlewi@gmail.com
/// @author Jeliazko Jeliazkov jeliazkov@jhu.edu (symmetry)

#ifndef INCLUDED_protocols_floppy_tail_FloppyTailMover_fwd_hh
#define INCLUDED_protocols_floppy_tail_FloppyTailMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace floppy_tail {

class FloppyTailMover;
typedef utility::pointer::shared_ptr< FloppyTailMover > FloppyTailMoverOP;

} //floppy_tail
} //protocols

#endif //INCLUDED_protocols_floppy_tail_FloppyTailMover_fwd_hh
