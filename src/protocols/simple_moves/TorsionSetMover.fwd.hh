// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/TorsionSetMover.fwd.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_simple_moves_TorsionSetMover_FWD_HH
#define INCLUDED_protocols_simple_moves_TorsionSetMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class TorsionSetMover;
typedef utility::pointer::shared_ptr< TorsionSetMover > TorsionSetMoverOP;
typedef utility::pointer::shared_ptr< TorsionSetMover const > TorsionSetMoverCOP;

} //simple_moves
} //protocols

#endif
