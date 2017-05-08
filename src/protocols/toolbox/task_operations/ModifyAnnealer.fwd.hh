// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ModifyAnnealer.fwd.hh
///
/// @brief Task operation to set high and low temps for annealer as well as whether or not to do a quench step
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_toolbox_task_operations_ModifyAnnealer_FWD_HH
#define INCLUDED_protocols_toolbox_task_operations_ModifyAnnealer_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class ModifyAnnealer;
typedef utility::pointer::shared_ptr< ModifyAnnealer > ModifyAnnealerOP;
typedef utility::pointer::shared_ptr< ModifyAnnealer const > ModifyAnnealerCOP;

} //task_operations
} //toolbox
} //protocols

#endif


