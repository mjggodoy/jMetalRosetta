// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/SmoothEnvEnergy.cc
/// @brief  Smooth, differentiable, version of centroid env
/// @author Frank DiMaio


// Unit headers
#include <core/scoring/methods/SmoothEnvEnergy.hh>
#include <core/scoring/methods/SmoothEnvEnergyCreator.hh>

// Package headers
#include <core/scoring/SmoothEnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SmoothEnvEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SmoothEnvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SmoothEnvEnergy );
}

ScoreTypes
SmoothEnvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_env_smooth );
	sts.push_back( cbeta_smooth );
	return sts;
}


/// c-tor
SmoothEnvEnergy::SmoothEnvEnergy() :
	parent( methods::EnergyMethodCreatorOP( new SmoothEnvEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_SmoothEnvPairPotential() )
{}


/// clone
EnergyMethodOP
SmoothEnvEnergy::clone() const {
	return EnergyMethodOP( new SmoothEnvEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
SmoothEnvEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
}


void
SmoothEnvEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	// symmetry-specific code
	// since it is a whole-structure energy, special treatment is needed to make sure this is computed correctly
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation &symmconf =
			dynamic_cast<core::conformation::symmetry::SymmetricConformation & >( pose.conformation());
		symmconf.recalculate_transforms(); // this is needed by deriv calcs
	}

	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
	potential_.compute_dcentroid_environment( pose );
}

///////////////////////////////////////
//
// ENV SCORE AND CBETA SCORE
void
SmoothEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;

	Real env_score( 0.0 ), cb_score6( 0.0 ), cb_score12( 0.0 ), cb_score( 0.0 );

	potential_.evaluate_env_and_cbeta_scores( pose, rsd, env_score, cb_score6, cb_score12 );

	//fpd  constants match nonsmooth variants
	env_score *= 2.019;
	cb_score = 2.667 * ( cb_score6 + cb_score12 ) * 0.3;

	emap[ cen_env_smooth ] += env_score;
	emap[ cbeta_smooth ] += cb_score;
} // residue_energy


void
SmoothEnvEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;

	Real weight_env = weights[ cen_env_smooth ];
	Real weight_cbeta = weights[ cbeta_smooth ];

	numeric::xyzVector<Real> d_env_score, d_cb_score6, d_cb_score12, d_cb_score;
	potential_.evaluate_env_and_cbeta_deriv( pose, rsd, d_env_score, d_cb_score6, d_cb_score12);

	// again, multiply by magic constants
	d_env_score *= 2.019;
	d_cb_score = 2.667 * ( d_cb_score6 + d_cb_score12 ) * 0.3;

	Vector atom_x = rsd.atom( rsd.nbr_atom() ).xyz();

	Vector const f2_env( d_env_score );
	Vector atom_y = -f2_env + atom_x;
	Vector const f1_env( atom_x.cross( atom_y ) );

	Vector const f2_cbeta( d_cb_score );
	atom_y = -f2_cbeta + atom_x;
	Vector const f1_cbeta( atom_x.cross( atom_y ) );

	atom_derivs[ rsd.nbr_atom() ].f1() += weight_env * f1_env + weight_cbeta * f1_cbeta;
	atom_derivs[ rsd.nbr_atom() ].f2() += weight_env * f2_env + weight_cbeta * f2_cbeta;
}


void
SmoothEnvEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const {
	potential_.finalize( pose );
}

core::Size
SmoothEnvEnergy::version() const {
	return 1; // Initial versioning
}

}
}
}
