// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterfaceTaskOperation.fwd.hh
/// @brief
/// @author Brian Weitzner (brian.weitzner@jhu.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_InterfaceTaskOperation_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_InterfaceTaskOperation_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {


class InterfaceTaskOperation;
typedef utility::pointer::shared_ptr< InterfaceTaskOperation > InterfaceTaskOperationOP;
typedef utility::pointer::shared_ptr< InterfaceTaskOperation > InterfaceTaskOperationCOP;

}
}
}

#endif
