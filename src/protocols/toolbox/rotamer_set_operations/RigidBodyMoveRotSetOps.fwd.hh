// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/RigidBodyMoveRotSetOps.fwd.hh
/// @brief  classes for rigidb body movement during rotamer packing
/// @author Florian Richter, floric@u.washington.edu, sep 2009

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_RigidBodyMoveRotSetOps_fwd_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_RigidBodyMoveRotSetOps_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

class RigidBodyMoveRSO;

typedef utility::pointer::shared_ptr< RigidBodyMoveRSO > RigidBodyMoveRSOOP;
typedef utility::pointer::shared_ptr< RigidBodyMoveRSO const > RigidBodyMoveRSOCOP;

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictOperationsBase_FWD_HH

