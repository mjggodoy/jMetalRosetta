// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief meta constraint where N constraints declared
/// @brief only the lowest K are evaluated
/// @author


#ifndef INCLUDED_core_scoring_constraints_KofNConstraint_hh
#define INCLUDED_core_scoring_constraints_KofNConstraint_hh

// Unit header
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/KofNConstraint.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>

//Utility Headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

class KofNConstraint : public MultiConstraint {
public:

	/// @brief Constructor
	KofNConstraint( core::Size K=0 );

	/// @brief Constructor
	KofNConstraint( ConstraintCOPs const & cst_in, core::Size K=0 );

	void
	setK( core::Size K );

	core::Size
	getK() const;

	virtual
	ConstraintOP clone() const;

	virtual
	MultiConstraintOP empty_clone() const;

	/// @brief
	void
	init_cst_score_types();

	std::string type() const;

	/// @brief read in constraint defiinition
	void
	read_def( std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory );

	/// @brief compute score
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	core::Real
	calculate_total_cst_score( EnergyMap const & weights, EnergyMap & emap) const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	void show( std::ostream& out) const;

	virtual
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	utility::vector1<ConstraintCOP> active_constraints() const;

	// void read_def( std::istream& in, pose::Pose const& pose, func::FuncFactory const& func_factory );
	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1.0 ) const;

	virtual bool operator == ( Constraint const & rhs ) const;
	virtual bool same_type_as_me( Constraint const & rhs ) const;

private:
	mutable utility::vector1<ConstraintCOP> active_constraints_;
	mutable core::Real cutoff_cst_score_;
	ScoreTypes cst_score_types_;
	core::Size K_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //KofNConstraint

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_KofNConstraint )
#endif // SERIALIZATION


#endif
