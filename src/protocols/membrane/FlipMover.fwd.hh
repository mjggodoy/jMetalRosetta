// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/FlipMoverCreator.hh
/// @brief      Flips a span or protein in the membrane (Rosetta Scripts Hook)
/// @details Flips a span, protein or part of a pose in the membrane,
///    depending on the jump number.
///    ONLY FOR FIXED MEMBRANE AND FLEXIBLE PROTEIN
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_FlipMover_fwd_hh
#define INCLUDED_protocols_membrane_FlipMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class FlipMover;
typedef utility::pointer::shared_ptr< FlipMover > FlipMoverOP;
typedef utility::pointer::shared_ptr< FlipMover const > FlipMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_FlipMover_fwd_hh
