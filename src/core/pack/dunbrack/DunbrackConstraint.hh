// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackConstraint.hh
/// @author James Thompson

#ifndef INCLUDED_core_pack_dunbrack_DunbrackConstraint_hh
#define INCLUDED_core_pack_dunbrack_DunbrackConstraint_hh

#include <core/pack/dunbrack/DunbrackConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace dunbrack {

class DunbrackConstraint : public scoring::constraints::Constraint {
public:
	DunbrackConstraint();

	virtual ~DunbrackConstraint();

	virtual std::string type() const;

	virtual
	scoring::constraints::ConstraintOP
	clone() const;

	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	virtual
	Size
	natoms() const;

	virtual
	AtomID const &
	atom( Size const index ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual void score(
		scoring::func::XYZ_Func const & xyz_func,
		scoring::EnergyMap const & weights,
		scoring::EnergyMap & emap
	) const;

	bool
	bin_matches( core::conformation::Residue const & resi ) const;

	/// @details Will return 1.0 if the bonus will be added, 0.0 otherwise
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		scoring::EnergyMap const & weights
	) const;

	/// @brief This gets used to compare one constraint to
	/// another, so it should uniquely reflect all the
	/// parameters.
	virtual void show( std::ostream & out ) const;

	virtual void read_def(
		std::istream & in,
		pose::Pose const & pose,
		scoring::func::FuncFactory const & func_factory
	);

private:
	core::Real bonus_;
	core::Size seqpos_;
	core::Size rot_vec_pos_; // position in RotVector
	core::Size rot_bin_;     // desired bin at this position
	utility::vector1< AtomID > atom_ids_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // DunbrackConstraint

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_dunbrack_DunbrackConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_dunbrack_DunbrackConstraint_HH
