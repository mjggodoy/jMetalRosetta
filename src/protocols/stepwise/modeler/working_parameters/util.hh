// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/working_parameters/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_StepWiseRNA_WorkingParametersUtil_HH
#define INCLUDED_protocols_stepwise_modeler_rna_StepWiseRNA_WorkingParametersUtil_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace working_parameters {

StepWiseWorkingParametersOP
setup_working_parameters_for_swa( utility::vector1< core::Size > const & moving_res_list,
	core::pose::Pose const & pose,
	core::pose::PoseCOP native_pose,
	utility::vector1< core::Size > const & bridge_res,
	utility::vector1< core::Size > const & working_minimize_res );

StepWiseWorkingParametersOP
setup_working_parameters_explicit( core::Size const rebuild_res,
	core::pose::Pose const & pose,
	core::pose::PoseCOP native_pose );

bool
figure_out_rebuild_bulge_mode( core::pose::Pose const & pose, core::Size const rebuild_res );

bool
figure_out_sample_both_sugar_base_rotamer( core::pose::Pose const & pose, bool const floating_base, core::Size const rebuild_suite );

} //working_parameters
} //modeler
} //stepwise
} //protocols

#endif
