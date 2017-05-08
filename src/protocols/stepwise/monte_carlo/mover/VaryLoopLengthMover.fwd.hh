// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.fwd.hh
/// @brief In stepwise design, vary desired loop lengths by updating FullModelParameters
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_stepwise_monte_carlo_mover_VaryLoopLengthMover_fwd_hh
#define INCLUDED_protocols_stepwise_monte_carlo_mover_VaryLoopLengthMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

class VaryLoopLengthMover;

typedef utility::pointer::shared_ptr< VaryLoopLengthMover > VaryLoopLengthMoverOP;
typedef utility::pointer::shared_ptr< VaryLoopLengthMover const > VaryLoopLengthMoverCOP;



} //protocols
} //stepwise
} //monte_carlo
} //mover


#endif //INCLUDED_protocols_stepwise_monte_carlo_mover_VaryLoopLengthMover_fwd_hh





