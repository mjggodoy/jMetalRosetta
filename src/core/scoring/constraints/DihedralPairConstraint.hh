// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/DihedralPairConstraint.hh
/// @brief Restrain a pair of residues to take the same torsion angles
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_constraints_DihedralPairConstraint_hh
#define INCLUDED_core_scoring_constraints_DihedralPairConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>

#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


/// constraint on dihedral angle formed by 4 points

class DihedralPairConstraint : public Constraint {

public:

	virtual std::string type() const {
		return "DihedralPair";
	}

	virtual ConstraintOP clone() const;

	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const;

	void read_def(
		std::istream & in,
		pose::Pose const & pose,
		func::FuncFactory const & func_factory
	);

	virtual bool operator == ( Constraint const & rhs ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	Real
	score(
		Vector const & X1, Vector const & X2, Vector const & X3, Vector const & X4,
		Vector const & Y1, Vector const & Y2, Vector const & Y3, Vector const & Y4
	) const;

	void
	score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	Real
	score( conformation::Conformation const & conformation  ) const;

	Real
	dihedral_diff(
		Vector const & X1, Vector const & X2, Vector const & X3, Vector const & X4,
		Vector const & Y1, Vector const & Y2, Vector const & Y3, Vector const & Y4
	) const;

	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1, Vector & F2,
		EnergyMap const & weights
	) const;

	///c-tor
	DihedralPairConstraint(
		AtomID const & a1, AtomID const & a2, AtomID const & a3, AtomID const & a4,
		AtomID const & b1, AtomID const & b2, AtomID const & b3, AtomID const & b4,
		func::FuncOP func,
		ScoreType scotype = dihedral_constraint
	):
		Constraint( scotype ),
		atomA1_(a1), atomA2_(a2), atomA3_(a3), atomA4_(a4),
		atomB1_(b1), atomB2_(b2), atomB3_(b3), atomB4_(b4),
		func_( func )
	{}

	Size
	natoms() const {
		return 8;
	}

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const & seqmap ) const;

	AtomID const &
	atom( Size const n ) const;

	virtual void show_def( std::ostream &, pose::Pose const & ) const;

	virtual void show( std::ostream & out ) const;

	// Coppied Remapped_clone from AtomPairConstraint
	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

private:
	Real
	func( Real const theta ) const;

	Real
	dfunc( Real const theta ) const;

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	DihedralPairConstraint( DihedralPairConstraint const & src );

private:
	// data
	AtomID atomA1_, atomA2_, atomA3_, atomA4_;
	AtomID atomB1_, atomB2_, atomB3_, atomB4_;
	func::FuncOP func_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DihedralPairConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_DihedralPairConstraint )
#endif // SERIALIZATION


#endif
