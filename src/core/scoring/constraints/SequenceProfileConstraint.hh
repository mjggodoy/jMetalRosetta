// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/SequenceProfileConstraint.hh
/// @brief  This is a constraint that refers to a core::sequence::SequenceProfile? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design).
/// @author ashworth

#ifndef INCLUDED_core_scoring_constraints_SequenceProfileConstraint_hh
#define INCLUDED_core_scoring_constraints_SequenceProfileConstraint_hh

#include <core/scoring/constraints/SequenceProfileConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/chemical/AA.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief
class SequenceProfileConstraint : public core::scoring::constraints::Constraint {
public:
	typedef core::sequence::SequenceProfile SequenceProfile;
	typedef core::sequence::SequenceProfileCOP SequenceProfileCOP;
	typedef core::id::SequenceMapping SequenceMapping;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Conformation Conformation;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::func::FuncFactory FuncFactory;
	typedef core::scoring::func::XYZ_Func XYZ_Func;
	typedef core::scoring::constraints::ConstraintOP ConstraintOP;

public:

	SequenceProfileConstraint();

	SequenceProfileConstraint(
		Pose const &,
		core::Size,
		SequenceProfileCOP profile = NULL,
		core::id::SequenceMappingCOP mapping = NULL // current pose numbers onto profile numbers.
	);

	SequenceProfileConstraint(
		core::Size,
		SequenceProfileCOP profile = NULL,
		core::id::SequenceMappingCOP mapping = NULL // current pose numbers onto profile numbers.
	);

	virtual ~SequenceProfileConstraint();

	virtual ConstraintOP clone() const;

	virtual bool operator == ( Constraint const & ) const;
	virtual bool same_type_as_me( Constraint const & ) const;

	virtual std::string type() const { return "SequenceProfile"; }

	/// @brief used by ConstraintIO and ConstraintFactory to construct this constraint from a input file stream (constraint file)
	virtual void
	read_def(
		std::istream &,
		Pose const &,
		func::FuncFactory const &
	);

	virtual void show_def( std::ostream &, Pose const & ) const;

	virtual void show( std::ostream & out ) const;

	core::Size seqpos() const { return seqpos_; }

	//hk
	void seqpos(core::Size seqpos) {
		seqpos_=seqpos;
	}
	//hk end

	void set_sequence_profile( SequenceProfileCOP profile, core::id::SequenceMappingCOP mapping = NULL);
	SequenceProfileCOP sequence_profile() const;
	core::id::SequenceMappingCOP profile_mapping() const;

	virtual core::Size natoms() const { return 0; }
	virtual AtomID const & atom( core::Size const ) const {
		utility_exit_with_message("SequenceProfileConstraint is not atom-based!.");
		return core::id::BOGUS_ATOM_ID; // required for compilation on Windows
	};

	virtual utility::vector1< core::Size > residues() const;

	virtual ConstraintOP remap_resid( SequenceMapping const & ) const;
	virtual ConstraintOP remapped_clone(pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual void
	score(
		XYZ_Func const &,
		EnergyMap const &,
		EnergyMap &
	) const;

	core::Real
	dist( core::chemical::AA aa ) const;

	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual void
	fill_f1_f2(
		AtomID const &,
		XYZ_Func const &,
		core::Vector &,
		core::Vector &,
		EnergyMap const &
	) const;

	void weight( core::Real const w );
	core::Real weight() const;

protected:
	SequenceProfileConstraint( SequenceProfileConstraint const & src );

private:
	core::Size seqpos_;
	// COP is deliberate, as sequence profiles are typically shared amoung many constraints, and as such probably shouldn't be fiddled with
	SequenceProfileCOP sequence_profile_;
	/// @brief A mapping of current pose numbers onto profile numbers.
	core::id::SequenceMappingCOP mapping_;
	core::Real weight_; //dflt 1; allows weighting the profile's contribution according to different criteria, for instance according to whether the position is close or far from the active site
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace constraints
} // namespace scoring
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_SequenceProfileConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_constraints_SequenceProfileConstraint_HH
