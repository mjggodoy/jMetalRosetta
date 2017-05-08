// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMover_FWD_HH
#define INCLUDED_protocols_loop_modeling_LoopMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace loop_modeling {

class LoopMover;

typedef utility::pointer::shared_ptr<LoopMover> LoopMoverOP;
typedef utility::pointer::shared_ptr<LoopMover const> LoopMoverCOP;
typedef utility::vector1<LoopMoverOP> LoopMoverOPs;
typedef utility::vector1<LoopMoverCOP> LoopMoverCOPs;

}
}

#endif

