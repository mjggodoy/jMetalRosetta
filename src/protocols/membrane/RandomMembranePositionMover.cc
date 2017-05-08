// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/RandomMembranePositionMover.cc
///
/// @brief      Random membrane position mover
/// @details Make random perturbations in the position of the membrane
///    as define by some center position and normal coordinate
///    (also of fixed thickness). Rotation mover rotates normal
///    from a random angle, Translation mover picks a randomly
///    translated position from the current position.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

// Unit Headers
#include <protocols/membrane/RandomMembranePositionMover.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/UniformPositionMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.RandomMembranePositionMover" );

namespace protocols {
namespace membrane {

// Random Rotation Mover //////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
RandomPositionRotationMover::RandomPositionRotationMover() :
	protocols::moves::Mover(),
	rot_mag_( 4.0 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify an order of magnitude to rotate
RandomPositionRotationMover::RandomPositionRotationMover(
	core::Real rot_mag,
	core::SSize rb_jump
) :
	protocols::moves::Mover(),
	rot_mag_( rot_mag ),
	rb_jump_( rb_jump )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
RandomPositionRotationMover::RandomPositionRotationMover( RandomPositionRotationMover const & ) = default;

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
RandomPositionRotationMover &
RandomPositionRotationMover::operator=( RandomPositionRotationMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new RandomPositionRotationMover( *this ) );

}

/// @brief Destructor
RandomPositionRotationMover::~RandomPositionRotationMover() = default;

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
RandomPositionRotationMover::get_name() const {
	return "RandomPositionRotationmover";
}

/// @brief Apply Rotation
void
RandomPositionRotationMover::apply( Pose & pose ) {

	using namespace protocols::simple_moves;
	using namespace core::conformation::membrane;

	// Check the pose is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// Check the membrane fold tree is reasonable
	if ( ! pose.conformation().membrane_info()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}

	// Compute random rotation
	core::Vector current_normal( pose.conformation().membrane_info()->membrane_normal(pose.conformation()) );
	current_normal.normalize();
	core::Real theta = 2*numeric::random::rg().uniform() * rot_mag_;

	// Apply Uniform Rotation
	UniformRotationMoverOP rotate( new UniformRotationMover( theta, current_normal, rb_jump_ ) );
	rotate->apply( pose );

}

// Random Translation Mover //////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Random Position Mover
RandomPositionTranslationMover::RandomPositionTranslationMover() :
	protocols::moves::Mover(),
	trans_mag_( 2.0 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify a magnitude of translation
RandomPositionTranslationMover::RandomPositionTranslationMover(
	core::Real trans_mag,
	core::SSize rb_jump
) :
	protocols::moves::Mover(),
	trans_mag_( trans_mag ),
	rb_jump_( rb_jump )
{}


/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
RandomPositionTranslationMover::RandomPositionTranslationMover( RandomPositionTranslationMover const & ) = default;

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
RandomPositionTranslationMover &
RandomPositionTranslationMover::operator=( RandomPositionTranslationMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new RandomPositionTranslationMover( *this ) );

}

/// @brief Destructor
RandomPositionTranslationMover::~RandomPositionTranslationMover() = default;

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
RandomPositionTranslationMover::get_name() const {
	return "RandomPositionTranslationMover";
}

/// @brief Apply Translation to membrane position
/// @brief Translate membrane position to new center
void
RandomPositionTranslationMover::apply( core::pose::Pose & pose ) {

	using namespace core;
	using namespace protocols::simple_moves;

	// Check the pose is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// Check the membrane fold tree is reasonable
	if ( ! pose.conformation().membrane_info()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}

	// Compute new position based on random translation
	core::conformation::Conformation const & conf( pose.conformation() );
	Vector current_center( conf.membrane_info()->membrane_center(conf) );
	Vector current_normal( conf.membrane_info()->membrane_normal(conf) );
	Vector delta_trans = 2*numeric::random::rg().uniform()-1 * trans_mag_ * current_normal;
	Vector new_position = current_center + delta_trans;

	// Apply translation
	UniformTranslationMoverOP translate( new UniformTranslationMover( new_position, rb_jump_ ) );
	translate->apply( pose );


}

} // membrane
} // protocols
