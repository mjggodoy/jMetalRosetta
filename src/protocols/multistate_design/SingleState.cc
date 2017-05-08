// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleState.cc
/// @brief
/// @author ashworth

#include <protocols/multistate_design/SingleState.hh>

#include <core/pose/Pose.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

////////////////////////////////////////////////////////////////////////////////////////////////////
SingleState::SingleState()
: utility::pointer::ReferenceCount(),
	pose_p_(/* 0 */), is_positive_state_(false), best_score_(0.), fitness_function_(nullptr)
{
	fitness_function_ = SingleStateFitnessFunctionCOP( SingleStateFitnessFunctionOP( new SingleStateFitnessFunction() ) );
}

SingleState::~SingleState()= default;

SingleState::SingleState( core::pose::Pose const & pose, bool is_positive )
: utility::pointer::ReferenceCount(),
	pose_p_(/* 0 */), is_positive_state_( is_positive ), best_score_(0.), fitness_function_(nullptr)
{
	pose_p_ = core::pose::PoseOP( new core::pose::Pose );
	*pose_p_ = pose;
	fitness_function_ = SingleStateFitnessFunctionCOP( SingleStateFitnessFunctionOP( new SingleStateFitnessFunction() ) );
}

SingleState::SingleState( SingleState const & other )
: utility::pointer::ReferenceCount(), pose_p_(/* 0 */), is_positive_state_(false), best_score_(0.), fitness_function_(nullptr)
{
	pose_p_ = core::pose::PoseOP( new core::pose::Pose );
	*pose_p_ = other.pose();
	is_positive_state_ = other.is_positive_state();
	best_score_ = other.best_score();
	fitness_function_ = other.fitness_function();
}

core::pose::Pose const & SingleState::pose() const { return *pose_p_; }
core::pose::Pose & SingleState::nonconst_pose() { return *pose_p_; }

} // namespace multistate_design
} // namespace protocols
