// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CombinePoseMover.hh
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_simple_moves_CombinePoseMover_hh
#define INCLUDED_protocols_simple_moves_CombinePoseMover_hh

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/simple_moves/CombinePoseMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/io/silent/SilentStruct.hh>
//#include <core/io/silent/SilentStruct.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

class CombinePoseMover : public protocols::moves::Mover

{
public:
	//constructor
	CombinePoseMover( core::scoring::ScoreFunctionCOP sfxn,
		core::pose::Pose const &pose_ref );

	~CombinePoseMover() override;


	void apply( core::pose::Pose &pose2 ) override;


	std::string get_name() const override { return "CombinePoseMover"; }

	void set_default();

	void set_max_struct( core::Size const n ){ max_struct_ = n; }
	void set_max_try( core::Size const n ){ max_struct_try_ = n; }
	void set_minfrac_crossover( core::Real const value ){ minfrac_crossover_ = value; }
	void set_maxfrac_crossover( core::Real const value ){ maxfrac_crossover_ = value; }
	void set_nonideal( bool const setting ){ nonideal_ = setting; }
	void set_minimize( bool const setting ){ do_minimize_ = setting; }
	void set_cartesian( bool const setting ){ cartesian_crossover_ = setting; }
	void set_rmsdcut( core::Real const setting ){ rmsdcut_ = setting; }

	core::Size max_struct() const { return max_struct_; }
	core::Size max_struct_try() const { return max_struct_try_; }

	core::pose::Pose const pose_ref() const { return pose_ref_; }

	std::string pose_tag() const { return pose_tag_; }
	void set_pose_tag( std::string const value ) { pose_tag_ = value; }

	void set_store_silents( bool const value ){ store_silents_ = value; }

	// Structures / history
	std::vector< core::io::silent::SilentStructOP >
	return_silent() const { return sampled_structures_; }

	void clear_combine_history() { combine_history_.resize( 0 ); }
	void append_combine_history( std::vector< core::Size > const v )
	{ combine_history_.push_back( v ); }

	std::vector< std::vector< core::Size > >
	return_combine_history() const { return combine_history_; }

private:
	core::pose::Pose pose_ref_;

	// Parameters/options
	core::Size max_struct_;
	core::Size max_struct_try_;
	core::Real minfrac_crossover_;
	core::Real maxfrac_crossover_;
	bool cartesian_crossover_;
	bool nonideal_;
	bool store_silents_;
	bool do_minimize_;

	// Tag
	std::string pose_tag_;

	// Minimization score
	core::scoring::ScoreFunctionCOP sfxn_;

	// Simple filter for clash
	core::scoring::ScoreFunctionCOP sfxn0_;
	core::Real vdwcut_;
	core::Real rmsdcut_;

	// structure storage
	std::vector< std::vector< core::Size > > combine_history_;
	std::vector< core::io::silent::SilentStructOP > sampled_structures_;

};

} // namespace simple_moves
} // namespace protocols

#endif
