// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/ResampleMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_ResampleMover_HH
#define INCLUDED_protocols_stepwise_monte_carlo_ResampleMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.fwd.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

class ResampleMover: public protocols::moves::Mover {

public:

	//constructor
	ResampleMover( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler );

	//destructor
	~ResampleMover();

public:

	using moves::Mover::apply;

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	bool
	apply( core::pose::Pose & pose,
		std::string & move_type );

	bool
	apply( core::pose::Pose & pose,
		StepWiseMove const & swa_move );

	bool
	apply( core::pose::Pose & pose,
		StepWiseMove const & swa_move,
		std::string & move_type );

	void set_minimize_single_res( bool const & setting ){ minimize_single_res_ = setting; }
	bool minimize_single_res() const{ return minimize_single_res_; }

	void
	set_options( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options ){ options_ = options; }

	Size
	get_remodel_res( StepWiseMove const & swa_move, core::pose::Pose const & pose ) const;

	void
	slide_jump_randomly( core::pose::Pose & pose, core::Size & remodel_res ) const;

private:

	protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler_;
	StepWiseMoveSelectorOP swa_move_selector_;
	protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options_;

	bool minimize_single_res_;
	bool slide_docking_jumps_;
};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
