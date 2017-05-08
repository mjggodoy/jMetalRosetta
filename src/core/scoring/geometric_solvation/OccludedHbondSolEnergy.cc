// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/geometric_solvation/OccludedHbondSolEnergy.cc
/// @brief  Pairwise-decomposable, continuous solvation model based on penalizing potential for Hbonding to solvent.
///  This is the pwSHO model.
/// @author John Karanicolas
/// @author Andrea Bazzoli


// NOTES FOR IMPROVED PERFORMANCE / POTENTIAL IMPROVEMENTS....

// The GeometricSolvation implementation was context-dependent,
// because backbone groups participating in secondary structure
// were considered exempt. Here, though, we'll compute their solvation
// energy as for any other group (partly since it's not obvious how
// else they *should* be treated). This in turn allows this energy
// term to be context-independent. The real question is....
// what should be the solvation energy for CO in secondary structure??

// Probably the best alternative would be to *NOT* compute solvation energies
// here for backbone groups in secondary structure, and instead assign
// them a fixed solvation penalty. Not clear what this penalty should be though...

// It might make sense to precompute and cache scores and derivatives
// in eg. a "geometric solvation potential" object,
// so that they don't need to be computed over and over again

#include <core/scoring/geometric_solvation/OccludedHbondSolEnergy.hh>
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergyCreator.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.geometric_solvation.OccludedHbondSolEnergy" );

namespace core {
namespace scoring {
namespace geometric_solvation {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the OccludedHbondSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
OccludedHbondSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new geometric_solvation::OccludedHbondSolEnergy( options ) );
}

ScoreTypes
OccludedHbondSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( occ_sol_fitted );
	return sts;
}


Vector dummy_deriv_vector_;

// jumpouts will apply if this is the best possible energy.
// this value corresponds to the discontinuity we'll deem acceptable.
// deriv_check starts to give bad results with a value of 0.05 (dist=6.6), but is mostly acceptable with 0.01 (dist=7.5)
core::Real const MIN_OCC_ENERGY = { 0.01 };


OccludedHbondSolEnergy::OccludedHbondSolEnergy(
	methods::EnergyMethodOptions const & options,
	bool const verbose )
:
	parent( methods::EnergyMethodCreatorOP( new OccludedHbondSolEnergyCreator ) ),
	atom_type_set_ptr_( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ) ),
	amp_scaling_factors_(atom_type_set_ptr_->n_atomtypes(), 0),
	occ_hbond_sol_database_( ScoringManager::get_instance()->get_DatabaseOccSolEne( options.etable_type(), MIN_OCC_ENERGY ) ),
	verbose_( verbose )
{
	init_amp_scaling_factors();
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy constructor" << std::endl;
}

OccludedHbondSolEnergy::OccludedHbondSolEnergy( OccludedHbondSolEnergy const & src ):
	parent( src ),
	atom_type_set_ptr_(src.atom_type_set_ptr_),
	amp_scaling_factors_(src.amp_scaling_factors_),
	occ_hbond_sol_database_( src.occ_hbond_sol_database_ ),
	verbose_( src.verbose_ )
{
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy constructor" << std::endl;
}

methods::EnergyMethodOP
OccludedHbondSolEnergy::clone() const
{
	return methods::EnergyMethodOP( new OccludedHbondSolEnergy( *this ) );
}

///
/// @brief initializes amplitude scaling factors from command-line
///
void OccludedHbondSolEnergy::init_amp_scaling_factors() {

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Ntrp")] = option[score::occ_sol_fitted::Ntrp_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("NH2O")] = option[score::occ_sol_fitted::NH2O_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Nlys")] = option[score::occ_sol_fitted::Nlys_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Narg")] = option[score::occ_sol_fitted::Narg_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Nbb")] = option[score::occ_sol_fitted::Nbb_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Nhis")] = option[score::occ_sol_fitted::Nhis_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OH")] = option[score::occ_sol_fitted::OH_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("ONH2")] = option[score::occ_sol_fitted::ONH2_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OOC")] = option[score::occ_sol_fitted::OOC_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Oaro")] = option[score::occ_sol_fitted::Oaro_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Oet2")] = option[score::occ_sol_fitted::Oet2_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Oet3")] = option[score::occ_sol_fitted::Oet3_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OCbb")] = option[score::occ_sol_fitted::OCbb_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("HOH")] = option[score::occ_sol_fitted::HOH_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OPha")] = option[score::occ_sol_fitted::OPha_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OHha")] = option[score::occ_sol_fitted::OHha_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OC3")] = option[score::occ_sol_fitted::OC3_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("OSi")] = option[score::occ_sol_fitted::OSi_amp_scaling];
	amp_scaling_factors_[atom_type_set_ptr_->atom_type_index("Oice")] = option[score::occ_sol_fitted::Oice_amp_scaling];
}


void
OccludedHbondSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
OccludedHbondSolEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
OccludedHbondSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( verbose_ ) tr << "jk evaluating residue pair energy" << std::endl;

	// Note: no count-pair stuff, these will just be computed normally
	// jk is there double-counting with other stuff, eg. backbone-dependent Dunbrack when we include self-terms like this?

	debug_assert ( rsd1.seqpos() != rsd2.seqpos() ); // this should be computed via eval_intrares_energy

	Real occ_solE =
		res_res_occ_sol_one_way( rsd1, rsd2 ) +
		res_res_occ_sol_one_way( rsd2, rsd1 ) ;

	// store the energies
	emap[ occ_sol_fitted ] += occ_solE;
}

void
OccludedHbondSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap ) const {

	if ( verbose_ ) tr << "jk evaluating intrares energy" << std::endl;

	// behaves as not-same-residue, except that we only do the calculation once
	Real occ_solE = res_res_occ_sol_one_way( rsd, rsd );
	emap[ occ_sol_fitted ] += occ_solE;
}

/// @details return true if the two residues are moving with respect to each other.
bool
OccludedHbondSolEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool res_moving_wrt_eachother
) const
{
	return res_moving_wrt_eachother;
}

void
OccludedHbondSolEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	eval_residue_pair_derivatives_one_way( rsd1, rsd2, weights, r1_atom_derivs, r2_atom_derivs );
	eval_residue_pair_derivatives_one_way( rsd2, rsd1, weights, r2_atom_derivs, r1_atom_derivs );
}

void
OccludedHbondSolEnergy::eval_residue_pair_derivatives_one_way(
	conformation::Residue const & rsd1, // polar residue
	conformation::Residue const & rsd2, // occluding residue
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( rsd1.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return;
	if ( rsd2.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return;
	if ( rsd1.has_variant_type( chemical::REPLONLY ) ) return;
	if ( rsd2.has_variant_type( chemical::REPLONLY ) ) return;

	Real energy(0); // dummy variable
	Real const atom_pair_cutoff( occ_hbond_sol_database_.atomic_interaction_cutoff() );
	Real const atom_pair_cutoff2( atom_pair_cutoff*atom_pair_cutoff );
	Real dcut = ( atom_pair_cutoff + rsd2.nbr_radius() );
	Real d2cut = dcut * dcut;
	core::Vector r2_nb_xyz( rsd2.nbr_atom_xyz());
	for ( Size const don_h_atom : rsd1.Hpos_polar() ) {
		Size const don_base_atom( rsd1.atom_base( don_h_atom ) );
		if ( rsd1.is_virtual( don_base_atom ) ) continue;
		if ( rsd1.is_repulsive( don_base_atom ) ) continue;
		if ( rsd1.is_virtual( don_h_atom ) ) continue;
		if ( rsd1.is_repulsive( don_h_atom ) ) continue;

		if ( rsd1.xyz( don_h_atom ).distance_squared( r2_nb_xyz ) > d2cut ) continue;
		for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
			if ( rsd2.is_virtual( jj ) ) continue;
			if ( rsd2.is_repulsive( jj ) ) continue;
			if ( rsd1.xyz( don_h_atom ).distance_squared( rsd2.xyz(jj) ) > atom_pair_cutoff2 ) continue;
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, rsd1, jj, rsd2, energy, true, weights[ occ_sol_fitted ],
				r1_atom_derivs[ don_base_atom ].f1(), r1_atom_derivs[ don_base_atom ].f2(),
				r1_atom_derivs[ don_h_atom ].f1(), r1_atom_derivs[ don_h_atom ].f2(),
				r2_atom_derivs[ jj ].f1(), r2_atom_derivs[ jj ].f2() );
		}
	}

	for ( Size const acc_atom : rsd1.accpt_pos() ) {
		Size const base_atom ( rsd1.atom_base( acc_atom ) );
		if ( rsd1.is_virtual( base_atom ) ) continue;
		if ( rsd1.is_repulsive( base_atom ) ) continue;
		if ( rsd1.is_virtual( acc_atom ) ) continue;
		if ( rsd1.is_repulsive( acc_atom ) ) continue;
		if ( rsd1.xyz( acc_atom ).distance_squared( r2_nb_xyz ) > d2cut ) continue;
		for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
			if ( rsd2.is_virtual( jj ) ) continue;
			if ( rsd2.is_repulsive( jj ) ) continue;

			if ( rsd1.xyz( acc_atom ).distance_squared( rsd2.xyz(jj) ) > atom_pair_cutoff2 ) continue;
			get_atom_atom_occ_solvation( acc_atom, base_atom, rsd1, jj, rsd2, energy, true, weights[ occ_sol_fitted ],
				r1_atom_derivs[ base_atom ].f1(), r1_atom_derivs[ base_atom ].f2(),
				r1_atom_derivs[ acc_atom ].f1(), r1_atom_derivs[ acc_atom ].f2(),
				r2_atom_derivs[ jj ].f1(), r2_atom_derivs[ jj ].f2() );
		}
	}
}


void
OccludedHbondSolEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	eval_residue_pair_derivatives_one_way( rsd, rsd, weights, atom_derivs, atom_derivs );
}


Real
OccludedHbondSolEnergy::res_res_occ_sol_one_way(
	conformation::Residue const & polar_rsd,
	conformation::Residue const & occ_rsd
) const
{
	if ( polar_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return 0.0;
	if ( occ_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return 0.0;
	if ( polar_rsd.has_variant_type( chemical::REPLONLY ) ) return 0.0;
	if ( occ_rsd.has_variant_type( chemical::REPLONLY ) ) return 0.0;



	// Rhiju importantly notes: for GeometricSolvation he originally had the code in
	// the following functions written out inside these loop -- and packing was faster.
	// Perhaps something to do with inlining or compiler optimization.
	// I've left it this way for now, because it helps prevent copying too
	// much of the code shared between residue pair scoring and for the derivatives.
	// However, if speed becomes important, here's a place to start.

	// jk note: moved the loop over occluding atoms into the next fxn, this could be the speed diff...

	Real geo_solE(0.), energy(0.);

	Real dcut = ( occ_hbond_sol_database_.atomic_interaction_cutoff() + occ_rsd.nbr_radius() );
	Real d2cut = dcut * dcut;
	core::Vector occ_nb_xyz( occ_rsd.nbr_atom_xyz());

	// cycle through donors in polar_rsd
	for ( Size const don_h_atom : polar_rsd.Hpos_polar() ) {
		Size const don_base_atom( polar_rsd.atom_base( don_h_atom ) );
		if ( polar_rsd.is_virtual( don_h_atom ) ) continue;
		if ( polar_rsd.is_repulsive( don_h_atom ) ) continue;
		if ( polar_rsd.is_virtual( don_base_atom ) ) continue;
		if ( polar_rsd.is_repulsive( don_base_atom ) ) continue;

		if ( polar_rsd.xyz( don_h_atom ).distance_squared( occ_nb_xyz ) > d2cut ) continue;
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	// cycle through acceptors in polar_rsd
	for ( Size const acc_atom : polar_rsd.accpt_pos() ) {
		Size const base_atom ( polar_rsd.atom_base( acc_atom ) );
		if ( polar_rsd.is_virtual( acc_atom ) ) continue;
		if ( polar_rsd.is_repulsive( acc_atom ) ) continue;
		if ( polar_rsd.is_virtual( base_atom ) ) continue;
		if ( polar_rsd.is_repulsive( base_atom ) ) continue;

		if ( polar_rsd.xyz( acc_atom ).distance_squared( occ_nb_xyz ) > d2cut ) continue;
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( acc_atom, base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	return geo_solE;
}


void
OccludedHbondSolEnergy::get_atom_atom_occ_solvation(
	Size const polar_atom,
	Size const base_atom,
	conformation::Residue const & polar_rsd,
	Size const occ_atom,
	conformation::Residue const & occ_rsd,
	Real & energy,
	bool const update_deriv,  // = false
	Real const occ_sol_fitted_weight, // = 0.0
	//bool const update_deriv_base,  // = false
	//bool const update_deriv_occ,  // = false
	Vector & f1_base, // = dummy vector
	Vector & f2_base, // = dummy vector
	Vector & f1_polar, // = dummy vector
	Vector & f2_polar, // = dummy vector
	Vector & f1_occ, // = dummy vector
	Vector & f2_occ // = dummy vector
) const
{
	// In case of early return, initialize. Note that energy does NOT accumulate, but f1/f2 do.
	// Also note that f1 and f2 are returned unweighted.
	energy = 0.;

	if ( occ_rsd.has_variant_type( chemical::REPLONLY ) ) return;
	if ( occ_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return;
	if ( polar_rsd.has_variant_type( chemical::REPLONLY ) ) return;
	if ( polar_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return;

	// note: after testing, hydrogens need not occlude
	if ( occ_rsd.atom_is_hydrogen(occ_atom) ) return;

	if ( occ_rsd.is_repulsive(occ_atom) ) return;
	if ( occ_rsd.is_virtual(occ_atom) ) return;
	if ( polar_rsd.is_repulsive(polar_atom) ) return;
	if ( polar_rsd.is_virtual(polar_atom) ) return;
	if ( polar_rsd.is_repulsive(base_atom) ) return;
	if ( polar_rsd.is_virtual(base_atom) ) return;

	// note: the lines above don't exclude Proline NV...
	// catch proline NV here (and other virtual atoms, etc.)
	if ( occ_rsd.atom_type(occ_atom).lj_radius() < 0.1 ) return;

	// can be occluded by atoms directly bonded to this group, but not by self
	if ( polar_rsd.seqpos() == occ_rsd.seqpos() ) {
		if ( polar_atom == occ_atom ) return;
		if ( base_atom == occ_atom ) return;
	}

	bool polar_atom_donates = false;
	if ( polar_rsd.atom_is_hydrogen(polar_atom) ) polar_atom_donates = true;

	if ( polar_atom_donates ) {
		// polar donor cannot be occluded by an acceptor (analogous to exact_occ_skip_Hbonders in exact model, but not quite the same)
		if ( occ_rsd.heavyatom_is_an_acceptor( occ_atom ) ) return;
	} else {
		// polar acceptor cannot be occluded by a donor base (analogous to exact_occ_skip_Hbonders in exact model, but not quite the same)
		if ( occ_rsd.heavyatom_has_polar_hydrogens( occ_atom ) ) return;
	}

	debug_assert( ( polar_atom_donates && atom_is_donor_h( polar_rsd, polar_atom ) ) ||
		( ( ! polar_atom_donates ) && atom_is_acceptor( polar_rsd, polar_atom ) ) );

	// get atom type lookup index for polar and occluding atoms; map indexes to those of known atom types, if necessary
	Size polar_atom_type_lookup_index = polar_rsd.atom_type_index( polar_atom );
	Size occ_atom_type_index = occ_rsd.atom_type_index( occ_atom );
	if ( polar_atom_donates ) {

		// do donor lookup with base atom
		polar_atom_type_lookup_index = occ_hbond_sol_database_.don_type_mapping( polar_rsd.atom_type_index( base_atom ) );
		occ_atom_type_index = occ_hbond_sol_database_.don_occ_type_mapping( occ_atom_type_index );
	} else {
		polar_atom_type_lookup_index = occ_hbond_sol_database_.acc_type_mapping( polar_atom_type_lookup_index );
		occ_atom_type_index = occ_hbond_sol_database_.acc_occ_type_mapping( occ_atom_type_index );
	}

	Vector const & polar_atom_xyz( polar_rsd.atom( polar_atom ).xyz() );
	Vector const & base_atom_xyz( polar_rsd.atom( base_atom ).xyz() );
	Vector const & occ_atom_xyz( occ_rsd.atom( occ_atom ).xyz() );

	// jumpout with no calculations if easy tests are violated, ie. no contribution to solvation energy
	Real const dist_sq = ( occ_atom_xyz - polar_atom_xyz).length_squared();
	if ( dist_sq > occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_max_sq_dist ) ) return;
	Real const curr_cos_angle = get_cos_angle( base_atom_xyz, polar_atom_xyz, occ_atom_xyz );
	if ( curr_cos_angle < occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_min_cos_angle ) ) return;

	// geometric filters are met, compute energy (and derivatives, if desired)
	// get the appropriate parameters
	Real sf = amp_scaling_factors_[polar_atom_type_lookup_index];
	if ( !sf ) {
		tr << "Unsupported atom type index: " << polar_atom_type_lookup_index << std::endl;
		exit(0);
	}
	Real const amp = sf*occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_amp );
	Real const dist_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_dist_mu );
	Real const twice_dist_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_dist_sigma_sq );
	Real const cos_angle_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_cos_angle_mu );
	Real const twice_cos_angle_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_cos_angle_sigma_sq );

	// Note: differences are in different order. Doesn't matter for scores, does for derivatives (and these make derivatives work).
	// Briefly, we're in the regime where dist energy contribution gets small as we get big values,
	// while cos_angle contribution gets small as we get smaller values
	Real const dist_diff = sqrt(dist_sq) - dist_mu;
	Real const cos_angle_diff = cos_angle_mu - curr_cos_angle;

	Real dist_exp, cos_angle_exp;
	if ( update_deriv ) {
		dist_exp = exp( - ( dist_diff * dist_diff / twice_dist_sigma_sq ) );
		cos_angle_exp = exp( - ( cos_angle_diff * cos_angle_diff / twice_cos_angle_sigma_sq ) );
		energy = amp * dist_exp * cos_angle_exp;
	} else {
		// do the calculation with a single exp
		energy = amp *
			exp( - ( ( dist_diff * dist_diff / twice_dist_sigma_sq ) + ( cos_angle_diff * cos_angle_diff / twice_cos_angle_sigma_sq ) ) );
	}

	if ( verbose_ && ( energy > 0. ) ) {
		tr <<"jk res "<< polar_rsd.name1() << I(3,polar_rsd.seqpos()) <<
			" atom " << polar_rsd.atom_name( polar_atom ) << " is occluded by occ_res " <<
			occ_rsd.name1() << I(3, occ_rsd.seqpos()) <<
			" atom " << occ_rsd.atom_name( occ_atom ) <<
			" with energy " << F(8,3,energy) << std::endl;
	}

	if ( ! update_deriv ) return;

	// compute angle f1/f2
	// note: energy is what was computed above, since it does NOT accumulate
	Real const curr_angle = numeric::arccos(curr_cos_angle); // note: radians
	Real const cos_angle_dfunc = -2. * cos_angle_diff * energy / twice_cos_angle_sigma_sq;
	Real const angle_dfunc = -std::sin(curr_angle) * cos_angle_dfunc * occ_sol_fitted_weight;
	Real theta(0.);
	Vector angle_f1_p1(0.), angle_f2_p1(0.);
	Vector angle_f1_p2(0.), angle_f2_p2(0.);
	Vector angle_f1_p3(0.), angle_f2_p3(0.);

	numeric::deriv::angle_p1_p2_p3_deriv( base_atom_xyz, polar_atom_xyz, occ_atom_xyz, theta,
		angle_f1_p1, angle_f2_p1, angle_f1_p2, angle_f2_p2, angle_f1_p3, angle_f2_p3 );
	f1_polar += angle_dfunc * angle_f1_p2;
	f2_polar += angle_dfunc * angle_f2_p2;
	f1_base += angle_dfunc * angle_f1_p1;
	f2_base += angle_dfunc * angle_f2_p1;
	f1_occ += angle_dfunc * angle_f1_p3;
	f2_occ += angle_dfunc * angle_f2_p3;

	// compute distance f1/f2
	// note: energy is what was computed above, since it does NOT accumulate
	Real const dist_dfunc = -2. * dist_diff * energy * occ_sol_fitted_weight / twice_dist_sigma_sq;
	Real dist(0.);
	Vector dist_f1(0.), dist_f2(0.);

	numeric::deriv::distance_f1_f2_deriv( polar_atom_xyz, occ_atom_xyz, dist, dist_f1, dist_f2 );
	f1_polar += dist_dfunc * dist_f1;
	f2_polar += dist_dfunc * dist_f2;

	f1_occ -= dist_dfunc * dist_f1;
	f2_occ -= dist_dfunc * dist_f2;
}

Distance
OccludedHbondSolEnergy::atomic_interaction_cutoff() const
{
	// tr << "atomic_interaction_cutoff is:  " << occ_hbond_sol_database_.atomic_interaction_cutoff() << std::endl;
	// jk max interaction distance is computed using the hydrogen for donors - is this okay? or should we add one to get a heavyatom distance?
	// probably is doesn't matter, since at worst we'll just end up using an acceptor-based distance, which is fine...
	return occ_hbond_sol_database_.atomic_interaction_cutoff();
}

Real
OccludedHbondSolEnergy::get_cos_angle(
	Vector const & base_atom_xyz,
	Vector const & polar_atom_xyz,
	Vector const & occluding_atom_xyz ) const
{
	return dot( (polar_atom_xyz - base_atom_xyz).normalize(), (occluding_atom_xyz - polar_atom_xyz).normalize() );
}


// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_donor_h( conformation::Residue const & rsd, Size const atom ) const {
	for ( Size const don_h_atom : rsd.Hpos_polar() ) {
		if ( don_h_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_acceptor( conformation::Residue const & rsd, Size const atom ) const {
	for ( Size const acc_atom : rsd.accpt_pos() ) {
		if ( acc_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_valid_base( conformation::Residue const & rsd, Size const atom ) const {
	for ( Size const don_h_atom : rsd.Hpos_polar() ) {
		Size const base_atom ( rsd.atom_base( don_h_atom ) );
		if ( base_atom == atom ) return true;
	}
	for ( Size const acc_atom : rsd.accpt_pos() ) {
		Size const base_atom ( rsd.atom_base( acc_atom ) );
		if ( base_atom == atom ) return true;
	}
	return false;
}
core::Size
OccludedHbondSolEnergy::version() const
{
	return 1; // Initial versioning
}


} // geometric_solvation
} // scoring
} // core

