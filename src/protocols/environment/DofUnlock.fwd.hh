// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DofUnlock.fwd.hh
/// @brief definition of the DofUnlock class
/// @author

#ifndef INCLUDED_protocols_environment_DofUnlock_fwd_hh
#define INCLUDED_protocols_environment_DofUnlock_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace environment {

class DofUnlock;
typedef utility::pointer::shared_ptr< DofUnlock > DofUnlockOP;
typedef utility::pointer::shared_ptr< DofUnlock const > DofUnlockCOP;

} // environment
} // protocols

#endif //INCLUDED_protocols_moves_DofUnlock_fwd_HH
