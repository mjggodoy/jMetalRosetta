// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protcols/normalmode/NormalModeMinimizer.hh
/// @brief  High-level atom tree minimizer class
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_normalmode_NormalModeMinimizer_hh
#define INCLUDED_protocols_normalmode_NormalModeMinimizer_hh

// Unit headers
#include <protocols/normalmode/NormalModeMinimizer.fwd.hh>
#include <protocols/normalmode/NormalModeMultiFunc.hh>

// Package headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace normalmode {

/// @brief High-level atom tree minimizer class
class NormalModeMinimizer : public protocols::moves::Mover
{
public:

	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::optimization::MinimizerOptionsOP MinimizerOptionsOP;
	typedef core::optimization::MinimizerOptionsCOP MinimizerOptionsCOP;
	typedef core::Real Real;

public:

	// c-tor
	NormalModeMinimizer();

	~NormalModeMinimizer() override;


	void
	apply( core::pose::Pose & pose ) override;

	/// @brief After minimization has concluded, the user may access the deriv-check result,
	/// assuming that they have run the NormalModeMinimizer with deriv_check = true;
	core::optimization::NumericalDerivCheckResultOP
	deriv_check_result() const;

	void set_modes( utility::vector1< core::Size > const modes_using_in ){
		modes_using_ = modes_using_in;
	}


	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;


	// XRW TEMP  std::string get_name() const override;
	protocols::moves::MoverOP clone() const override { return protocols::moves::MoverOP( new NormalModeMinimizer( *this ) ); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new NormalModeMinimizer ); }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void
	deriv_check_local(
		core::pose::Pose const &pose,
		protocols::normalmode::NormalModeMultifunc f ) const;


private:
	utility::vector1< core::Size > modes_using_;

	core::kinematics::MoveMapOP movemap_;
	ScoreFunctionCOP scorefxn_;
	MinimizerOptionsOP options_;

	Real dampen_;

	core::optimization::NumericalDerivCheckResultOP deriv_check_result_;
}; // NormalModeMinimizer


} // namespace optimization
} // namespace core


#endif // INCLUDED_protocols_normalmode_NormalModeMinimizer_HH
