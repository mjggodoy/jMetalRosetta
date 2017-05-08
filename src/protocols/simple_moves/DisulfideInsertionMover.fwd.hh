// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DisulfideInsertionMover.fwd.hh
/// @brief Forward header file for DisulfideInsertionMover
/// @author Orly Marcu ( orly.marcu@mail.huji.ac.il )
/// @date Jan. 12, 2015

#ifndef INCLUDED_protocols_simple_moves_DisulfideInsertionMover_fwd_hh
#define INCLUDED_protocols_simple_moves_DisulfideInsertionMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class DisulfideInsertionMover;
typedef utility::pointer::shared_ptr< DisulfideInsertionMover > DisulfideInsertionMoverOP;
typedef utility::pointer::shared_ptr< DisulfideInsertionMover const > DisulfideInsertionMoverCOP;

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_DisulfideInsertionMover_FWD_HH
