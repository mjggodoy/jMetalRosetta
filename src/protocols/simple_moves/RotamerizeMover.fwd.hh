// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/simple_moves/RotamerizeMover.fwd.hh
/// @brief forward declarations for rotamerizing
/// @author Jim Havranek

#ifndef INCLUDED_protocols_simple_moves_RotamerizeMover_FWD_HH
#define INCLUDED_protocols_simple_moves_RotamerizeMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class RotamerizeMover;
typedef utility::pointer::shared_ptr< RotamerizeMover > RotamerizeMoverOP;
typedef utility::pointer::shared_ptr< RotamerizeMover const > RotamerizeMoverCOP;

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_RotamerizeMover_fwd
