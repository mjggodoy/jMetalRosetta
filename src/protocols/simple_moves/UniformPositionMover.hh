// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/simple_moves/UniformPositionMover.hh
///
/// @brief      Apply a uniform (deterministic move) to a given position
/// @details Generic movers for applying a uniform rotation or translation move
///    along a jump to a given partner in the pose.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_simple_moves_UniformPositionMover_hh
#define INCLUDED_protocols_simple_moves_UniformPositionMover_hh

// Unit Headers
#include <protocols/simple_moves/UniformPositionMover.fwd.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace simple_moves {

/// @brief Uniform Rotation Mover
class UniformRotationMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Custom Constructor
	/// @details Specify a new normal to rotate membranes to
	/// to move this position to
	UniformRotationMover(
		core::Real alpha,
		core::Vector axis,
		core::SSize rb_jump
	);

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	UniformRotationMover( UniformRotationMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	UniformRotationMover &
	operator=( UniformRotationMover const & src );

	/// @brief Destructor
	~UniformRotationMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	std::string get_name() const override;

	/// @brief Apply Rotation
	/// @brief Rotate the membrane to the new normal position
	void apply( Pose & pose ) override;

private: // methods

	/// @brief Construct a Default Membrane Position Mover
	UniformRotationMover();

private:

	// Store new normal axis
	core::Real alpha_;
	core::Vector axis_;

	// Store jump num
	core::SSize rb_jump_;

};

/// @brief Uniform Translation Mover
class UniformTranslationMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Custom Constructor
	/// @details Specify a new center position to translate this stub to
	UniformTranslationMover(
		core::Vector new_position_,
		core::SSize rb_jump
	);

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	UniformTranslationMover( UniformTranslationMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	UniformTranslationMover &
	operator=( UniformTranslationMover const & src );

	/// @brief Destructor
	~UniformTranslationMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	std::string get_name() const override;

	/// @brief Apply Translation to membrane position
	/// @brief Translate membrane position to new center
	void apply( Pose & pose ) override;

private:

	/// @brief Construct a Default Membrane Position Mover
	UniformTranslationMover();

private:

	// Store new center
	core::Vector new_position_;

	// Store jump num
	core::SSize rb_jump_;


};

} // simple_moves
} // protocols

#endif // INCLUDED_protocols_simple_moves_UniformPositionMover_hh
