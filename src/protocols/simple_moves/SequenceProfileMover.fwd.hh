// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SequenceProfileMover.cc
/// @brief  BS mover to get around a stupid "mover" that was embedded in the parser
/// @author Brian Weitzner brian.weitzner@gmail.com, Steven Lewis smlewi@gmail.com
/// @date   Rebased to next year.


#ifndef INCLUDED_protocols_simple_moves_SequenceProfileMover_fwd_hh
#define INCLUDED_protocols_simple_moves_SequenceProfileMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

// Forward
class SequenceProfileMover;

typedef utility::pointer::shared_ptr< SequenceProfileMover > SequenceProfileMoverOP;
typedef utility::pointer::shared_ptr< SequenceProfileMover const > SequenceProfileMoverCOP;

} // namespace simple_moves
} // namespace protocols

#endif //INCLUDED_protocols_simple_moves_SequenceProfileMover_fwd_hh
