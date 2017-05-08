// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CyclizationMover.hh
/// @brief Forward header file for CyclizationMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_simple_moves_CyclizationMover_fwd_hh
#define INCLUDED_protocols_simple_moves_CyclizationMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class CyclizationMover;
typedef utility::pointer::shared_ptr< CyclizationMover > CyclizationMoverOP;
typedef utility::pointer::shared_ptr< CyclizationMover const > CyclizationMoverCOP;

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_CyclizationMover_FWD_HH
