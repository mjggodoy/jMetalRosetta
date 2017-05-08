// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleState.fwd.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_multistate_design_SingleState_fwd_hh
#define INCLUDED_protocols_multistate_design_SingleState_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace multistate_design {

class SingleState;
typedef utility::pointer::shared_ptr< SingleState > SingleStateOP;
typedef utility::pointer::shared_ptr< SingleState const > SingleStateCOP;
typedef utility::vector1< SingleStateOP > SingleStateOPs;
typedef utility::vector1< SingleStateCOP > SingleStateCOPs;

}
}

#endif
