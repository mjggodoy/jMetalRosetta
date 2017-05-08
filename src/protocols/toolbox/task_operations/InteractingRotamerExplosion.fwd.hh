// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/InteractingRotamerOperation.fwd.hh
/// @brief
/// @author Florian Richter, florian.richter.1@hu-berlin.de

#ifndef INCLUDED_protocols_toolbox_task_operations_InteractingRotamerOperation_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_InteractingRotamerOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class InteractingRotamerOperation;

typedef utility::pointer::shared_ptr< InteractingRotamerOperation > InteractingRotamerOperationOP;
typedef utility::pointer::shared_ptr< InteractingRotamerOperation const > InteractingRotamerOperationCOP;

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_InteractingRotamerOperation_FWD_HH
