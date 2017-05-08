// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/StackElecEnergy.cc
/// @brief  FA_Elec, but just 'perpendicular' to bases... trying to separate from H-bonds & geom_sol.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/StackElecEnergy.hh>
#include <core/scoring/rna/StackElecEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


////////////////////////////////////////////////////////////////////////////////////
// StackElecEnergy: Mash-up of FA_ElecEnergy & RNA_FullAtomStackingEnergy.
//
//  We need an approximate treatment of electrostatics, to handle several
//   puzzles in RNA modeling, including penalties  for 5'-CG-3' turner rule,
//   favored structures in the tandem GA motif, etc.
//
//  Unfortunately, just using fa_elec doesn't solve this -- severe convolution
//   with hydrogen bonds. This stack_elec function separates electrostatics
//   'perpendicular' to nucleobases away from hydrogen bonds & geom_sol which
//   have been balanced to account for in-plane interactions (base pairing).
//
//  Inspired by discussions with Fang-Chieh Chou.
//
//     -- Rhiju, July 2013
//
////////////////////////////////////////////////////////////////////////////////////

using namespace core::chemical;

// C++
static THREAD_LOCAL basic::Tracer tr( "core.scoring.rna.StackElecEnergy" );

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the StackElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
StackElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new StackElecEnergy( options ) );
}

ScoreTypes
StackElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( stack_elec );
	sts.push_back( stack_elec_base_base );
	sts.push_back( stack_elec_base_bb   );
	return sts;
}

/// c-tor
StackElecEnergy::StackElecEnergy( methods::EnergyMethodOptions const & options ) :
	parent( methods::EnergyMethodCreatorOP( new StackElecEnergyCreator ) ),
	coulomb_( options ),
	base_base_only_( false ), //true will be faster computation but appears less accurate (and less physically consistent)
	verbose_( false ),
	might_be_designing_( false ),
	using_extended_method_( false )
{
	coulomb_.initialize();
}

////////////////////////////////////////////////////////////////////////////
StackElecEnergy::StackElecEnergy( StackElecEnergy const & src ):
	parent( src ),
	coulomb_( src.coulomb() ),
	base_base_only_( src.base_base_only_ ),
	verbose_( src.verbose_ ),
	might_be_designing_( src.might_be_designing_ ),
	using_extended_method_( src.using_extended_method_ )
{
	coulomb_.initialize();
}


//clone
methods::EnergyMethodOP
StackElecEnergy::clone() const
{
	return methods::EnergyMethodOP( new StackElecEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
StackElecEnergy::setup_for_packing(
	pose::Pose  & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & designing_residues
) const {

	for ( bool const designing_this : designing_residues ) {
		if ( designing_this ) {
			might_be_designing_ = true;
			break;
		}
	}

	pose.update_residue_neighbors();

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );
}


void
StackElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	pose.update_residue_neighbors();

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );

	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::STACK_ELEC_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );
}

////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !pose.energies().use_nblist() ) return;

	// stash our nblist inside the pose's energies object
	Energies & energies( pose.energies() );

	// setup the atom-atom nblist
	NeighborListOP nblist;
	Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
	Real const XX = coulomb().max_dis() + 2 * tolerated_motion;
	nblist = NeighborListOP( new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX ) );
	if ( pose.energies().use_nblist_auto_update() ) {
		nblist->set_auto_update( tolerated_motion );
	}
	// this partially becomes the EtableEnergy classes's responsibility
	nblist->setup( pose, sfxn, *this );
	energies.set_nblist( EnergiesCacheableDataType::STACK_ELEC_NBLIST, nblist );
}

////////
bool
StackElecEnergy::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}


///////////////////////////////////////////////////////////////////////////////
bool
StackElecEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const {
	if ( rsd1.has_variant_type( REPLONLY ) ) return false;
	if ( rsd2.has_variant_type( REPLONLY ) ) return false;

	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		return false;
	}
	return res_moving_wrt_eachother;
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
StackElecEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const {
	using namespace etable::count_pair;
	if ( res1 == res2 ) {
		return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );
	}

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
StackElecEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair( rsd1, rsd2, true ) ) return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairAll ) );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
StackElecEnergy::get_intrares_countpair(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &
) const {
	utility_exit_with_message( "StackElecEnergy does not define intra - residue pair energies; do not call get_intrares_countpair()" );
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
bool
StackElecEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( pose.energies().use_nblist_auto_update() ) return;

	etable::count_pair::CountPairFunctionCOP count_pair =
		get_count_pair_function( rsd1, rsd2 );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( elec_pair_nblist ) ) );
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( coulomb().max_dis() + 2*tolerated_narrow_nblist_motion, 2 );

	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );

	pair_data.set_data( elec_pair_nblist, nblist );
}

////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using_extended_method_ = true;

	if ( pose.energies().use_nblist_auto_update() ) return;
	Real score( 0.0 ), score_base_base( 0.0 ), score_base_bb( 0.0 );

	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd2.has_variant_type( REPLONLY ) ) return;

	if ( !rsd1.is_RNA() || !rsd2.is_RNA() ) return;

	//debug_assert( dynamic_cast< ResiduePairNeighborList const * > (min_data.get_data( elec_pair_nblist )() ));
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( elec_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	if ( neighbs.size() <= 0 ) return;

	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );
	Size const i( rsd1.seqpos() );
	Size const j( rsd2.seqpos() );

	kinematics::Stub stub_i = base_stubs[i];
	kinematics::Stub stub_j = base_stubs[j];

	Matrix const M_i ( stub_i.M );
	Matrix const M_j ( stub_j.M );

	for ( auto const & neighb : neighbs ) {
		Size const m = neighb.atomno1();
		if ( rsd1.is_virtual( m ) ) continue;
		if ( rsd1.is_repulsive( m ) ) continue;
		if ( base_base_only_ && !is_rna_base( rsd1, m ) ) continue;
		Real const m_charge( rsd1.atomic_charge( m ) );
		if ( m_charge == 0.0 ) continue;

		Size const n = neighb.atomno2();
		if ( rsd2.is_virtual( n ) ) continue;
		if ( rsd2.is_repulsive( n ) ) continue;
		if ( base_base_only_ && !is_rna_base( rsd2, n ) ) continue;
		Real const n_charge( rsd2.atomic_charge( n ) );
		if ( n_charge == 0.0 ) continue;
		Vector const atom_m( rsd1.xyz( m ) );
		Vector const atom_n( rsd2.xyz( n ) );

		if ( is_rna_base( rsd1, m ) ) {
			Real cos_kappa2( 0.0 ); // useful for output...
			Real const stack_elec_score = get_stack_elec_score( atom_m, atom_n, m_charge, n_charge, M_i, cos_kappa2 );
			score += stack_elec_score;
			if ( is_rna_base( rsd2, n ) ) {
				score_base_base += stack_elec_score;
			} else {
				score_base_bb += stack_elec_score;
			}
		}
		if ( is_rna_base( rsd2, n ) ) {
			Real cos_kappa2( 0.0 ); // useful for output...
			Real const stack_elec_score = get_stack_elec_score( atom_n, atom_m, n_charge, m_charge, M_j, cos_kappa2 );
			score += stack_elec_score;
			if ( is_rna_base ( rsd1, m ) ) {
				score_base_base += stack_elec_score;
			} else {
				score_base_bb += stack_elec_score;
			}
		}
	}

	emap[ stack_elec ]           += score;
	emap[ stack_elec_base_base ] += score_base_base;
	emap[ stack_elec_base_bb ]   += score_base_bb;
}

//////////////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using_extended_method_ = false;
	if ( pose.energies().use_nblist() ) return;
	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd2.has_variant_type( REPLONLY ) ) return;

	Real score_base_base1( 0.0 ), score_base_base2( 0.0 );
	Real score_base_bb1( 0.0 ), score_base_bb2( 0.0 );

	Real const score = residue_pair_energy_one_way( rsd1, rsd2, pose, score_base_base1, score_base_bb1 ) +
		residue_pair_energy_one_way( rsd2, rsd1, pose, score_base_base2, score_base_bb2 ) ;

	emap[ stack_elec ]           += score;
	emap[ stack_elec_base_base ] += score_base_base1 + score_base_base2;
	emap[ stack_elec_base_bb ]   += score_base_bb1   + score_base_bb2;
}


//////////////////////////////////////////////////////////////////////////////////////////
Real
StackElecEnergy::residue_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & score_base_base,
	Real & score_base_bb
) const {
	score_base_base = 0.0;
	score_base_bb   = 0.0;

	// currently defined only for RNA -- could easily generalize to DNA or proteins,
	//  as long as we precompute centroid & base-normals for rings.
	if ( !rsd1.is_RNA() ) return 0.0;
	if ( rsd1.has_variant_type( REPLONLY ) ) return 0.0;
	if ( rsd2.has_variant_type( REPLONLY ) ) return 0.0;

	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );
	Size const i( rsd1.seqpos() );

	kinematics::Stub stub_i = base_stubs[i];

	if ( might_be_designing_ ) {
		Vector centroid1 = rna_centroid_info.get_base_centroid( rsd1 );
		stub_i = rna_centroid_info.get_base_coordinate_system( rsd1, centroid1 );
	}

	Real score( 0.0 );

	Matrix const M_i ( stub_i.M );

	// Loop over base heavy atoms.
	// If I want to generalize this to proteins, maybe could loop over "aromatic" atoms.
	for ( Size m = 1; m <= rsd1.natoms(); ++m )  {

		// following contains virtual check.
		if ( !is_rna_base( rsd1, m ) ) continue;
		if ( rsd1.is_repulsive( m ) ) continue;

		Real const i_charge( rsd1.atomic_charge( m ) );
		if ( i_charge == 0.0 ) continue;

		Vector const atom_i( rsd1.xyz( m ) );

		for ( Size n = 1; n <= rsd2.natoms(); ++n ) {

			if ( rsd2.is_virtual( n ) ) continue;
			if ( rsd2.is_repulsive( n ) ) continue;
			if ( base_base_only_ && !is_rna_base( rsd2, n ) ) continue;

			Real const j_charge( rsd2.atomic_charge( n ) );
			if ( j_charge == 0.0 ) continue;

			Vector const atom_j( rsd2.xyz( n ) );

			Real cos_kappa2( 0.0 ); // useful for output...
			Real const stack_elec_score = get_stack_elec_score( atom_i, atom_j, i_charge, j_charge, M_i, cos_kappa2 );

			score += stack_elec_score;

			if ( is_rna_base( rsd2, n ) ) {
				score_base_base += stack_elec_score;
			} else {
				score_base_bb   += stack_elec_score;
			}
		}
	}

	return score;
}


//////////////////////////////////////////////////////////////////////////////
bool
StackElecEnergy::is_rna_base(
	conformation::Residue const & rsd,
	Size const & m
) const {

	if ( !rsd.type().is_RNA() ) return false;

	//Need to be careful nheavyatoms count includes hydrogen! when the hydrogen is made virtual..
	if ( rsd.is_virtual( m ) ) return false;

	if ( rsd.atom_is_hydrogen( m ) ) return is_rna_base( rsd, rsd.atom_base( m ) );

	// Following is really really specific to RNA.
	return ( m > rsd.first_sidechain_atom() && m <= rsd.nheavyatoms() );
}

////////////////////////////////////////////////////////////////////////////// (Need to condition this? Parin Sep 2, 2009)
void
StackElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	if ( ! pose.energies().use_nblist_auto_update() ) return;

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd1.is_repulsive( m ) ) return;

	if ( rsd1.is_virtual( m ) ) return;

	Real const i_charge( rsd1.atomic_charge( m ) );
	if ( i_charge == 0.0 ) return;
	bool res1_is_base = is_rna_base( rsd1, m );
	if ( base_base_only_ && !res1_is_base ) return;

	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	Vector const atom_i( rsd1.xyz( m ) );
	kinematics::Stub stub_i = base_stubs[i];
	Matrix const M_i ( stub_i.M );

	//debug_assert( pose.energies().use_nblist() );
	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::STACK_ELEC_NBLIST ) );
	AtomNeighbors const & nbrs( nblist.atom_neighbors( i, m ) );

	for ( scoring::AtomNeighbors::const_iterator it2 = nbrs.begin(),
			it2e = nbrs.end(); it2 != it2e; ++it2 ) {
		scoring::AtomNeighbor const & nbr( *it2 );
		Size const j( nbr.rsd() );
		if ( i == j ) continue;
		Size const n( nbr.atomno() );
		conformation::Residue const & rsd2( pose.residue( j ) );

		if ( rsd2.is_virtual( n ) ) continue;
		if ( rsd2.is_repulsive( n ) ) continue;

		Real const j_charge( rsd2.atomic_charge( n ) );
		if ( j_charge == 0.0 ) continue; /// should prune out such atoms when constructing the neighborlist!
		bool res2_is_base = is_rna_base( rsd2, n );
		if ( base_base_only_ && !res2_is_base ) continue;

		Vector const & atom_j( rsd2.xyz( n ) );
		kinematics::Stub stub_j = base_stubs[j];
		Matrix const M_j ( stub_j.M );

		if ( res1_is_base ) {
			Vector const & deriv_vector_i = get_stack_elec_deriv( atom_i, atom_j, i_charge, j_charge, M_i );
			Vector force_vector_i = weights[ stack_elec ] * deriv_vector_i;

			if ( weights[ stack_elec_base_base ] != 0.0 &&  res2_is_base ) force_vector_i += weights[ stack_elec_base_base ] * deriv_vector_i;
			if ( weights[ stack_elec_base_bb   ] != 0.0 && !res2_is_base ) force_vector_i += weights[ stack_elec_base_bb   ] * deriv_vector_i;

			//Force/torque with which occluding atom j acts on "dipole" i.
			F1 += -1.0 * cross( force_vector_i, atom_j );
			F2 += -1.0 * force_vector_i;

		}

		if ( res2_is_base ) {
			Vector const & deriv_vector_j = get_stack_elec_deriv( atom_j, atom_i, j_charge, i_charge, M_j );
			Vector force_vector_j = weights[ stack_elec ] * deriv_vector_j;

			if ( weights[ stack_elec_base_base ] != 0.0  &&  res1_is_base )  force_vector_j += weights[ stack_elec_base_base ] * deriv_vector_j;
			if ( weights[ stack_elec_base_bb   ] != 0.0  && !res1_is_base )  force_vector_j += weights[ stack_elec_base_bb   ] * deriv_vector_j;

			F1 += cross( force_vector_j, atom_i );
			F2 += force_vector_j;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//get_stack_elec_deriv evaluates the stack_elec_score between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by residue_pair_energy_one_way
Real
StackElecEnergy::get_stack_elec_score( Vector const & r_i,
	Vector const & r_j,
	Real const & i_charge,
	Real const & j_charge,
	Matrix const & M_i,
	Real & cos_kappa2
) const {

	Vector const z_i = M_i.col_z();

	Vector const r_vec = r_j - r_i;
	Real const r = r_vec.length();
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	Real score  = coulomb().eval_atom_atom_fa_elecE( r_i, i_charge, r_j, j_charge );

	//Orientation dependence
	cos_kappa2 = cos_kappa * cos_kappa;

	score *= cos_kappa2;

	return score;
}


////////////////////////////////////////////////////////////////////////////////
//get_stack_elec_deriv evaluates the stack_elec_deriv between a pair of atoms.
//This function is called by eval_atom_derivative
Vector
StackElecEnergy::get_stack_elec_deriv( Vector const & r_i,
	Vector const & r_j,
	Real const & i_charge,
	Real const & j_charge,
	Matrix const & M_i
) const {

	//  Well, the energy is a function of the inter-atom distance and a special cos(angle)
	//            E = E ( r, cos(theta ) ).
	//  so
	//      dE/dx = dE/dr (x/r)   +   ( - x * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dy = dE/dr (y/r)   +   ( - y * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dz = dE/dr (z/r)   +   ( (r^2 - z^2) / r^3 ) ( dE/dcos(theta) )

	Vector const x_i = M_i.col_x();
	Vector const y_i = M_i.col_y();
	Vector const z_i = M_i.col_z();

	Vector const r_vec = r_j - r_i;
	Real const r = r_vec.length();
	Real const x = dot( x_i, r_vec );
	Real const y = dot( y_i, r_vec );
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	/////////////////////////////////
	//dE_dcoskappa
	/////////////////////////////////
	Real dE_dcoskappa = coulomb().eval_atom_atom_fa_elecE( r_i, i_charge, r_j, j_charge );
	dE_dcoskappa  *= 2 * cos_kappa;


	/////////////////////////////////
	//dE_dr
	/////////////////////////////////
	Real const r2 = r_vec.length_squared();
	Real dE_dr_over_r  = coulomb().eval_dfa_elecE_dr_over_r( r2, i_charge, j_charge );
	//Orientation dependence
	dE_dr_over_r *= cos_kappa * cos_kappa ;

	Real const dE_dx = ( dE_dr_over_r ) * x - ( dE_dcoskappa ) * ( x * z )/ ( r * r * r ) ;
	Real const dE_dy = ( dE_dr_over_r ) * y - ( dE_dcoskappa ) * ( y * z )/ ( r * r * r ) ;
	Real const dE_dz = ( dE_dr_over_r ) * z + ( dE_dcoskappa ) * ( x * x  +  y * y )/ ( r * r * r );

	return ( dE_dx * x_i + dE_dy * y_i + dE_dz * z_i );
}

/////////////////////////////
void
StackElecEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.calculated() = false;

	if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) return;

	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	// add in contributions from the nblist atom-pairs
	NeighborList const & nblist
		( pose.energies().nblist( EnergiesCacheableDataType::STACK_ELEC_NBLIST ) );

	nblist.check_domain_map( pose.energies().domain_map() );
	utility::vector1< conformation::Residue const * > resvect;
	resvect.reserve( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		resvect.push_back( & pose.residue( ii ) );
	}

	Real score( 0.0 ), score_base_base( 0.0 ), score_base_bb( 0.0 );

	for ( Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		conformation::Residue const & ires( *resvect[i] );
		for ( Size ii = 1, ii_end = ires.natoms(); ii <= ii_end; ++ii ) {
			if ( ires.is_virtual( ii ) ) continue;
			if ( ires.is_repulsive( ii ) ) continue;
			bool res1_is_base = is_rna_base( ires, ii );
			if ( base_base_only_ && !res1_is_base ) continue;
			Real const m_charge( ires.atomic_charge( ii ) );
			if ( m_charge == 0.0 ) continue;
			Vector const atom_m( ires.xyz( ii ) );
			kinematics::Stub stub_i = base_stubs[i];
			Matrix const M_i ( stub_i.M );
			AtomNeighbors const & nbrs( nblist.upper_atom_neighbors( i, ii ) );
			for ( AtomNeighbors::const_iterator nbr_iter = nbrs.begin(),
					nbr_end = nbrs.end(); nbr_iter != nbr_end; ++nbr_iter ) {
				AtomNeighbor const & nbr( *nbr_iter );

				Size const  j( nbr.rsd() );
				if ( i == j ) continue;
				Size const jj( nbr.atomno() );
				// could reorder the nbr lists so that we dont need this check:

				conformation::Residue const & jres( *resvect[j] );
				if ( jres.is_virtual ( jj ) ) continue;
				if ( jres.is_repulsive( jj ) ) continue;
				bool res2_is_base = is_rna_base( jres, jj );
				if ( base_base_only_ && !res2_is_base ) continue;
				Real const n_charge( jres.atomic_charge( jj ) );
				if ( n_charge == 0.0 ) continue;

				Vector const atom_n( jres.xyz( jj ) );
				kinematics::Stub stub_j = base_stubs[j];

				Matrix const M_j ( stub_j.M );

				if ( res1_is_base ) {
					Real cos_kappa2( 0.0 ); // useful for output...
					Real const stack_elec_score = get_stack_elec_score( atom_m, atom_n, m_charge, n_charge, M_i, cos_kappa2 );
					score += stack_elec_score;
					if ( res2_is_base ) {
						score_base_base += stack_elec_score;
					} else {
						score_base_bb += stack_elec_score;
					}
				}
				if ( res2_is_base ) {
					Real cos_kappa2( 0.0 ); // useful for output...
					Real const stack_elec_score = get_stack_elec_score( atom_n, atom_m, n_charge, m_charge, M_j, cos_kappa2 );
					score += stack_elec_score;
					if ( res1_is_base ) {
						score_base_base += stack_elec_score;
					} else {
						score_base_bb += stack_elec_score;
					}
				}
			}
		}
	}
	//std::cout << score;
	//std::cout << "\n";
	totals[ stack_elec ]           += score;
	totals[ stack_elec_base_base ] += score_base_base;
	totals[ stack_elec_base_bb ]   += score_base_bb;
}

/// @brief StackElecEnergy distance cutoff
Distance
StackElecEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}

core::Size
StackElecEnergy::version() const
{
	return 1; // Initial versioning
}


} //rna
} //scoring
} //core
