// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackEnergy.cc
/// @brief  Dunbrack energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/pack/dunbrack/DunbrackEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/id/TorsionID.hh>

// Utility headers
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

using namespace scoring;
using namespace scoring::methods;

/// @details This must return a fresh instance of the DunbrackEnergy class,
/// never an instance already in use
scoring::methods::EnergyMethodOP
DunbrackEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const &
) const {
	return scoring::methods::EnergyMethodOP( new DunbrackEnergy );
}

scoring::ScoreTypes
DunbrackEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_dun );
	sts.push_back( fa_dun_dev );
	sts.push_back( fa_dun_rot );
	sts.push_back( fa_dun_semi );
	return sts;
}


/// ctor
DunbrackEnergy::DunbrackEnergy() :
	parent( EnergyMethodCreatorOP( new DunbrackEnergyCreator ) )
{}

DunbrackEnergy::~DunbrackEnergy() {}

/// clone
EnergyMethodOP
DunbrackEnergy::clone() const
{
	return EnergyMethodOP( new DunbrackEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
DunbrackEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return;

	if ( rsd.is_virtual_residue() ) return;

	//Returns the equivalent L-amino acid library if a D-amino acid is provided
	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib = rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );

	if ( ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) ) return;

	dunbrack::RotamerLibraryScratchSpace scratch;
	emap[ fa_dun ] += rotlib->rotamer_energy( rsd, scratch );
	emap[ fa_dun_rot ] += scratch.fa_dun_rot();
	emap[ fa_dun_semi ] += scratch.fa_dun_semi();
	emap[ fa_dun_dev  ] += scratch.fa_dun_dev();
}

bool DunbrackEnergy::defines_dof_derivatives( pose::Pose const & ) const { return true; }

Real
DunbrackEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & ,//min_data,
	id::DOF_ID const & ,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & ,//pose,
	scoring::ScoreFunction const & ,//sfxn,
	scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return 0.0;

	Real deriv( 0.0 );
	Real deriv_dev( 0.0 );
	Real deriv_rot( 0.0 );
	Real deriv_semi( 0.0 );
	if ( ! tor_id.valid() ) return 0.0;

	debug_assert( rsd.seqpos() == tor_id.rsd() );

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );
	if ( ! rsd.is_protein() || ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) )  { return 0.0; }

	dunbrack::RotamerLibraryScratchSpace scratch;
	rotlib->rotamer_energy_deriv( rsd, scratch );
	if ( tor_id.type() == id::BB && tor_id.torsion() <= DUNBRACK_MAX_BBTOR ) {
		deriv      = scratch.dE_dbb()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dbb_dev()[ tor_id.torsion() ];
		deriv_rot  = scratch.dE_dbb_rot()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dbb_semi()[ tor_id.torsion() ];
	} else if ( tor_id.type() == id::CHI && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_SCTOR ) {
		deriv      = scratch.dE_dchi()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dchi_dev()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dchi_semi()[ tor_id.torsion() ];
	}

	return numeric::conversions::degrees( weights[ fa_dun ] * deriv + weights[ fa_dun_dev ] * deriv_dev + weights[ fa_dun_rot ] * deriv_rot + weights[ fa_dun_semi ] * deriv_semi);
}


Real
DunbrackEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	Real deriv( 0.0 );
	Real deriv_dev( 0.0 );
	Real deriv_rot( 0.0 );
	Real deriv_semi( 0.0 );
	if ( !tor_id.valid() )  return 0.0;

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( tor_id.rsd() ).type() );

	if ( pose.residue( tor_id.rsd() ).is_virtual_residue() ) return 0.0;

	/// ASSUMPTION: Derivatives for amino acids only!
	if ( ! rotlib || ! pose.residue_type( tor_id.rsd() ).is_protein() || pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::SC_BRANCH_POINT) )  { return 0.0; }

	dunbrack::RotamerLibraryScratchSpace scratch;
	rotlib->rotamer_energy_deriv( pose.residue( tor_id.rsd() ), scratch );

	if ( tor_id.type() == id::BB  && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_BBTOR ) {
		deriv      = scratch.dE_dbb()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dbb_dev()[ tor_id.torsion() ];
		deriv_rot  = scratch.dE_dbb_rot()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dbb_semi()[ tor_id.torsion() ];
	} else if ( tor_id.type() == id::CHI && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_SCTOR ) {
		deriv      = scratch.dE_dchi()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dchi_dev()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dchi_semi()[ tor_id.torsion() ];
	}

	return numeric::conversions::degrees( weights[ fa_dun ] * deriv + weights[ fa_dun_dev ] * deriv_dev + weights[ fa_dun_rot ] * deriv_rot + weights[ fa_dun_semi ] * deriv_semi );
}

/// @brief DunbrackEnergy is context independent; indicates that no context graphs are required
void
DunbrackEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}
core::Size
DunbrackEnergy::version() const
{
	return 1; // Initial versioning
}


} // dunbrack
} // pack
} // core

