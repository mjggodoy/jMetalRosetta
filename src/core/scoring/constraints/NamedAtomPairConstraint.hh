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

#ifndef INCLUDED_core_scoring_constraints_NamedAtomPairConstraint_hh
#define INCLUDED_core_scoring_constraints_NamedAtomPairConstraint_hh

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <utility/vector1.hh>


// C++ Headers
//#include <cstdlib>
//#include <iostream>
//#include <map>
//#include <utility>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/* this helper class may have become obsoleted by the remapped_clone() method of ConstraintSet
but since it is still in use by abinitio::Template I haven't removed it.
*/
class Obsolet_NamedAtomPairConstraint : public utility::pointer::ReferenceCount {
	typedef core::scoring::constraints::AtomPairConstraintOP AtomPairConstraintOP;
	typedef core::id::NamedAtomID NamedAtomID;
	typedef core::id::SequenceMapping SequenceMapping;
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Obsolet_NamedAtomPairConstraint();
	Obsolet_NamedAtomPairConstraint( AtomPairConstraintOP, core::pose::Pose const& );
	Obsolet_NamedAtomPairConstraint( NamedAtomID const& atom1, NamedAtomID const& atom2, AtomPairConstraintOP cst );
	AtomPairConstraintOP mapto( SequenceMapping const&, core::pose::Pose const& ) const;
	Obsolet_NamedAtomPairConstraintOP mapto( SequenceMapping const& ) const;
	AtomPairConstraintOP mapto( core::pose::Pose const& ) const;

	friend std::ostream& operator<< ( std::ostream& out, Obsolet_NamedAtomPairConstraint const& cst );

	id::NamedAtomID const& atom1() {
		return atom1_;
	}

	id::NamedAtomID const& atom2() {
		return atom2_;
	}

private:

	id::NamedAtomID atom1_;
	id::NamedAtomID atom2_;
	AtomPairConstraintOP cst_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Obsolet_NamedAtomPairConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class NamedAtomPairConstraint : public AtomPairConstraint {
public:
	NamedAtomPairConstraint(
		id::NamedAtomID const& a1,
		id::NamedAtomID const& a2,
		func::FuncOP func,
		ScoreType scoretype = atom_pair_constraint
	);

	virtual ConstraintOP clone() const;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=NULL ) const;

	virtual bool operator == ( Constraint const & rhs ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	//@brief translates the atom-names into numbers
	virtual void setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const;

	virtual void show_def( std::ostream& out, pose::Pose const& pose ) const;
	void show_def_nopose( std::ostream& out ) const;

	virtual void read_def( std::istream& in, pose::Pose const& pose,func::FuncFactory const& func_factory );
	// //@brief set constraint such that the pose doesn't violate it.
	// virtual void steal( pose::Pose& );

	/// @brief This overrride updates the sequence numbering but not the atom names.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual ConstraintOP remap_resid( core::id::SequenceMapping const &seqmap ) const;

private:
	id::NamedAtomID named_atom1_;
	id::NamedAtomID named_atom2_;
	core::Size type1_id_;
	core::Size type2_id_;
	//@brief this could contain a checksum made from the "annotated-sequence"
	// unused Size pose_chemical_checksum_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NamedAtomPairConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_NamedAtomPairConstraint )
#endif // SERIALIZATION


#endif
