// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/DOF_Node.fwd.hh
/// @brief  Kinematics
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_DOF_Node_fwd_hh
#define INCLUDED_core_optimization_DOF_Node_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace optimization {

class DOF_Node;
typedef utility::pointer::shared_ptr< DOF_Node > DOF_NodeOP;
typedef utility::pointer::shared_ptr< DOF_Node const > DOF_NodeCOP;

} // namespace kinematics
} // namespace core


#endif
