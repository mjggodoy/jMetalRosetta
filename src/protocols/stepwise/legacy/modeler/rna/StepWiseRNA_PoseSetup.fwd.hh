// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.fwd.hh
/// @author Rhiju Das
/// @author Parin Sripakdeevong

#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_PoseSetup_FWD_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_PoseSetup_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

class StepWiseRNA_PoseSetup;
typedef utility::pointer::shared_ptr< StepWiseRNA_PoseSetup > StepWiseRNA_PoseSetupOP;
typedef utility::pointer::shared_ptr< StepWiseRNA_PoseSetup const > StepWiseRNA_PoseSetupCOP;

} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
