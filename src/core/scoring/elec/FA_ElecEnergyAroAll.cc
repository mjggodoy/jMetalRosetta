// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/FA_ElecEnergyAroAll.cc
/// @brief  Electrostatics energy method for aromatic side chain (stand-in for pi/pi interactions).
/// @author Rhiju Das


// Unit headers
#include <core/scoring/elec/FA_ElecEnergyAroAll.hh>
#include <core/scoring/elec/FA_ElecEnergyAroAllCreator.hh>

// Package headers
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>

// ObjexxFCL headers

// C++


/////////////////////////////////////////////////////////////////////////////////////////
// Quick version of fa_elec for just aromatic residues in proteins.
//  This could be made much faster (for packing) by using Andrew Leaver-Fay's trie stuff,
//   copying/pasting from FA_ElecEnergy.cc. For now, I just went with what I knew.
//         -- Rhiju, Nov. 2009
//
// Modified aroaro to calculate electrostatic interaction between an aromatic residue
// and all surrounding atoms from other residues. Derivitives are not implemented.
//     -- James Gleixner, April 2011

/////////////////////////////////////////////////////////////////////////////////////////
///
/// Hacky (hence the name) implementation of 10r dielectric model, cutoff at 5.5A
///

//     alternatives: WARSHEL (from ligand.cc)
//     E = 322.0637*q1*q2/r/e(r)
//     if ( r < 3 ) e(r) = 16.55
//     else         e(r) = 1 + 60*(1-exp(-0.1*r))
//     Warshel, A. Russell, S. T., Q. Rev. Biophys., 1984, 17, 283-422
//
//
//     sigmoidal dielectric: (hingerty 1985)
//
//     e(r) = D - (D-D0)/2 [ (rS)**2 + 2rS + 2]exp(-rS)
//     with eg:
//     D = 78, D0 = 1, S = 0.3565 (rouzina&bloomfield)
//     D = 80, D0 = 4, S = 0.4 (rohs)

namespace core {
namespace scoring {
namespace elec {


/// @details This must return a fresh instance of the FA_ElecEnergyAroAll class,
/// never an instance already in use
methods::EnergyMethodOP
FA_ElecEnergyAroAllCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new FA_ElecEnergyAroAll( options ) );
}

ScoreTypes
FA_ElecEnergyAroAllCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_elec_aro_all );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
FA_ElecEnergyAroAll::FA_ElecEnergyAroAll(
	methods::EnergyMethodOptions const & options
):
	parent( options )
{
	set_score_types( methods::EnergyMethodCreatorOP( new FA_ElecEnergyAroAllCreator ) );
}


////////////////////////////////////////////////////////////////////////////
FA_ElecEnergyAroAll::FA_ElecEnergyAroAll( FA_ElecEnergyAroAll const & src ):
	parent( src )
{
	set_score_types( methods::EnergyMethodCreatorOP( new FA_ElecEnergyAroAllCreator ) );
}

/// clone
methods::EnergyMethodOP
FA_ElecEnergyAroAll::clone() const
{
	return methods::EnergyMethodOP( new FA_ElecEnergyAroAll( *this ) );
}


void
FA_ElecEnergyAroAll::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
FA_ElecEnergyAroAll::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


// The FA_ElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
FA_ElecEnergyAroAll::setup_for_packing(
	pose::Pose &,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
FA_ElecEnergyAroAll::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
FA_ElecEnergyAroAll::update_residue_for_packing(
	pose::Pose &,
	Size
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
FA_ElecEnergyAroAll::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// Aromatic means: PHE, TRP, TYR (not HIS, currently).
	if ( rsd1.is_aromatic() ) {
		residue_pair_energy_aro_aro( rsd1, rsd2, emap ); // Only checks to make sure one of the residues is aromatic
	}
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}


//////////////////////////////////////
bool atom_is_aro2( conformation::Residue const & rsd, Size const i )
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// String lookup is slow, but overall won't be too bad here because there won't be many aro/aro pairs.
	// However, we could be fancy and do a one-pass lookup at the beginning of the function...
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::string const atom_type_name( rsd.atom_type( i ).name() );
	if ( atom_type_name == "aroC" ||
			atom_type_name == "Ntrp" ||
			atom_type_name == "Nhis" ||
			atom_type_name == "Oaro" ||
			atom_type_name == "Haro" ) return true;

	return false;
}

//////////////////////////////////////////////////////////////////////////////////
// Following copies a little code, but at least separates out this inelegant and
// probably useless RNA stuff from the usual stuff.
Real
FA_ElecEnergyAroAll::residue_pair_energy_aro_aro(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;

	Real total_score( 0.0 );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size i=1, i_end = rsd1.natoms(); i<= i_end; ++i ) {
		Vector const & i_xyz( rsd1.xyz(i) );
		Real const i_charge( rsd1.atomic_charge(i) );

		if ( !atom_is_aro2( rsd1, i ) ) continue;

		//debug output
		//std::cerr << "residue i " << rsd1.seqpos() << " is aromatic: " << rsd1.is_aromatic() << std::endl;
		//std::cerr << "residue j " << rsd2.seqpos() << " is aromatic: " << rsd2.is_aromatic() << std::endl;

		if ( i_charge == 0.0 ) continue;

		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

			Real const j_charge( rsd2.atomic_charge(j) );
			if ( j_charge == 0.0 ) continue;
			Real weight(1.0);
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, weight, path_dist ) ) {
				Real score = weight *
					coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge);

				total_score += score;
				emap[ fa_elec_aro_all ] += score;
				//debug output
				//std::cout << "elecenergy: " << rsd1.seqpos() << "  " << rsd1.atom_name(i) << "  " << rsd2.seqpos() << " " << rsd2.atom_name(j) << " is: " << score << std::endl;
			}
		}
	}
	return total_score;
}

void
FA_ElecEnergyAroAll::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	if ( set1.num_rotamers() >= 1 && set2.num_rotamers() >= 1 &&
			set1.rotamer(1)->is_aromatic() && set2.rotamer(1)->is_aromatic() ) {
		grandparent::evaluate_rotamer_pair_energies( set1, set2, pose, sfxn, weights, energy_table );
	} // else, non aromatic/aromatic interaction; early return
}

void
FA_ElecEnergyAroAll::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	if ( set.num_rotamers() >= 1 && set.rotamer(1)->is_aromatic() && residue.is_aromatic() ) {
		grandparent::evaluate_rotamer_background_energies( set, residue, pose, sfxn, weights, energy_vector );
	} // else, non aromatic/aromatic  interaction; early return
}


void
FA_ElecEnergyAroAll::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace etable::count_pair;

	// what is my charge?
	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	if ( ! rsd1.is_aromatic() ) return;

	Real const i_charge( rsd1.atomic_charge( i ) );
	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	if ( i_charge == 0.0 ) return;

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	//  kinematics::DomainMap const & domain_map( energies.domain_map() );
	//  bool const pos1_fixed( !energies.res_moved( pos1 ) );
	// debug_assert( pos1_fixed == ( domain_map(pos1) != 0 ) ); // this is probably not generally true but I'm curious

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( utility::graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		//EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		//Size const pos2( edge->get_second_node_ind() );
		Size const pos2( (*iru)->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );

		debug_assert( pos2 != pos1 );

		if ( rsd2.is_aromatic() ) {
			eval_atom_derivative_aro_aro( rsd1, i, rsd2, weights, F1, F2 );
		}

	} // loop over nbrs of rsd1

}


//////////////////////////////////////////////////
void
FA_ElecEnergyAroAll::eval_atom_derivative_aro_aro(
	conformation::Residue const & rsd1,
	Size const & i,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real const i_charge( rsd1.atomic_charge( i ) );
	Vector const & i_xyz( rsd1.xyz(i) );

	// NOTE THAT is_base, is_sugar, and is_phosphate contains MAGIC NUMBERS (for speed!)
	// Might be better to make it precomputed as part of the residue_type definition?
	//This repeats some unnecessary stuff, but hey, derivatives don't have to be that fast.
	if ( !atom_is_aro2( rsd1, i ) ) return;

	for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

		if ( !atom_is_aro2( rsd2, j ) ) continue;

		Real const j_charge( rsd2.atomic_charge(j) );
		if ( j_charge == 0.0 ) continue;

		Real weight(1.0);
		Size path_dist( 0 );
		if ( cpfxn->count( i, j, weight, path_dist ) ) {
			Vector const & j_xyz( rsd2.xyz(j) );
			Vector const f2( i_xyz - j_xyz );
			Real const dis2( f2.length_squared() );
			Real const dE_dr_over_r = weight *
				coulomb().eval_dfa_elecE_dr_over_r( dis2, i_charge, j_charge );

			if ( dE_dr_over_r == 0.0 ) continue;

			Vector const f1( i_xyz.cross( j_xyz ) );

			F1 += weights[ fa_elec_aro_all ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_aro_all ] * dE_dr_over_r * f2;
		}
	}
}


/// @brief FA_ElecEnergyAroAll is context independent; no context graphs required
void
FA_ElecEnergyAroAll::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
FA_ElecEnergyAroAll::version() const
{
	return 1; // Initial versioning
}


}
}
}
