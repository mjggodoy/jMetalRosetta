// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/Loops.fwd.hh
/// @brief  Loops and Loop class forward declarations header
/// @author Chu Wang


#ifndef INCLUDED_protocols_loops_Loops_FWD_HH
#define INCLUDED_protocols_loops_Loops_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {

// Forward
class Loops;

typedef utility::pointer::shared_ptr< Loops > LoopsOP;
typedef utility::pointer::shared_ptr< Loops const > LoopsCOP;

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_Loops_FWD_HH
