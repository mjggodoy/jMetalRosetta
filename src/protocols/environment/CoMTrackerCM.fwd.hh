// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CoMTrackerCM.fwd.hh
/// @brief definition of the CoMTrackerCM class
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_CoMTrackerCM_fwd_hh
#define INCLUDED_protocols_environment_CoMTrackerCM_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace environment {

class CoMTrackerCM;
typedef utility::pointer::shared_ptr< CoMTrackerCM > CoMTrackerCMOP;
typedef utility::pointer::shared_ptr< CoMTrackerCM const > CoMTrackerCMCOP;

} // environment
} // protocols

#endif //INCLUDED_protocols_enviroment_CoMTrackerCM_fwd_HH
