// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/legacy/screener/SimplePoseSelection.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_SimplePoseSelection_HH
#define INCLUDED_protocols_stepwise_screener_SimplePoseSelection_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/legacy/screener/SimplePoseSelection.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

class SimplePoseSelection: public stepwise::screener::StepWiseScreener {

public:

	//constructor
	SimplePoseSelection( core::pose::Pose const & pose,
		utility::vector1< Size > const & moving_res_list,
		protocols::stepwise::modeler::options::StepWiseModelerOptionsCOP options,
		bool const full_optimize );

	//destructor
	~SimplePoseSelection();

public:

	std::string
	name() const { return "SimplePoseSelection"; }

	stepwise::screener::StepWiseScreenerType
	type() const { return stepwise::screener::SIMPLE_POSE_SELECTION; }

	bool
	check_screen();

	void
	finalize();

	utility::vector1< core::pose::PoseOP > pose_list() { return pose_list_; }

private:

	core::pose::Pose const & pose_;
	utility::vector1< Size > const moving_res_list_;
	protocols::stepwise::modeler::options::StepWiseModelerOptionsCOP options_;
	bool const full_optimize_;

	utility::vector1< core::pose::PoseOP > pose_list_;

};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
