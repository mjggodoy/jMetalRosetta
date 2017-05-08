// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_core_scoring_constraints_AngleConstraint_hh
#define INCLUDED_core_scoring_constraints_AngleConstraint_hh

#include <core/scoring/constraints/AngleConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>

#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief An Angular Constraint.
class AngleConstraint : public Constraint {
public:

	/// @brief Constructor
	AngleConstraint(
		AtomID const & a1,
		AtomID const & a2,
		AtomID const & a3,
		core::scoring::func::FuncOP func_in, // we take ownership of this guy
		ScoreType scotype = angle_constraint
	);

	/// @brief Constructor without atom IDs -- if you create an AngleConstraint with
	/// this constructor, you must never call its score( XYZFunc ) method!  Dangerous and stupid!
	AngleConstraint(
		core::scoring::func::FuncOP func_in,
		ScoreType scoretype = angle_constraint
	);

	virtual std::string type() const;

	virtual ConstraintOP clone() const;

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( Constraint const & other ) const;

	virtual bool same_type_as_me( Constraint const & other ) const;

	/// @brief read in constraint defiinition
	void read_def( std::istream & data, pose::Pose const & pose, core::scoring::func::FuncFactory const & func_factory );

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	/// @brief compute score
	Real
	score(
		Vector const & xyz1,
		Vector const & xyz2,
		Vector const & xyz3
	) const;

	/// @brief compute score
	virtual
	void
	score( core::scoring::func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	virtual
	core::Real
	score( conformation::Conformation const & ) const;

	/// @brief compute
	Real
	angle(
		Vector const & xyz1,
		Vector const & xyz2,
		Vector const & xyz3
	) const;

	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual void setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	/// @brief number of atoms --- always 3 for angles
	Size
	natoms() const
	{
		return 3;
	}

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief return AtomID for atom 1,2,3
	AtomID const &
	atom( Size const n ) const;

	/// @brief output violation of constraint to out - returns 1 if violated ( i.e., func.show_violations() > 0 )
	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, core::Real threshold = 1  ) const;

	virtual void show(std::ostream& out ) const;

	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	virtual core::scoring::func::Func const& get_func() const {
		return *func_;
	}


public:

	/// Previously private member functions made public so that in the absence of atom_ids,
	/// these functions could still be called externally.

	/// @brief evaluate func at theta
	Real
	func( Real const theta ) const;

	/// @brief evaluate dfunc at theta
	Real
	dfunc( Real const theta ) const;

	static
	void
	p1_theta_deriv(
		Vector const & p1,
		Vector const & p2,
		Vector const & p3,
		Vector & f1,
		Vector & f2
	);


	static
	void
	helper(
		Vector const & M,
		Vector const & w,
		Vector & F1,
		Vector & F2
	);


	void
	p1_deriv(
		Vector const & p1,
		Vector const & p2,
		Vector const & p3,
		Vector & F1,
		Vector & F2
	) const;


	void
	p2_deriv(
		Vector const & p1,
		Vector const & p2,
		Vector const & p3,
		Vector & F1,
		Vector & F2
	) const;


	AtomID const & atom1() const { return atom1_; }
	AtomID const & atom2() const { return atom2_; }
	AtomID const & atom3() const { return atom3_; }

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	AngleConstraint( AngleConstraint const & src );

	/// @brief const access to func
	func::FuncCOP func() const;

	/// @brief set func
	void set_func( func::FuncOP f );

	void atom1( AtomID newid ) const;
	void atom2( AtomID newid ) const;
	void atom3( AtomID newid ) const;

private:
	// data
	mutable AtomID atom1_, atom2_, atom3_;
	core::scoring::func::FuncOP func_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AngleConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_AngleConstraint )
#endif // SERIALIZATION


#endif
