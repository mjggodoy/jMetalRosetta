// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMover.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.phosphate.MultiPhosphateSampler" );
using utility::tools::make_vector1;

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace phosphate {

MultiPhosphateSampler::MultiPhosphateSampler( pose::Pose const & reference_pose ):
	prepacked_( false )
{
	reset( reference_pose );
	initialize_parameters();
}

//Constructor
MultiPhosphateSampler::MultiPhosphateSampler( pose::Pose & pose_to_prepack,
	Size const moving_res /*sets partition*/ )
{
	initialize_parameters();
	initialize_by_prepack( pose_to_prepack, moving_res ); // prepacked --> true.
}


//Destructor
MultiPhosphateSampler::~MultiPhosphateSampler()
{}

//////////////////////////////
MultiPhosphateSamplerOP
MultiPhosphateSampler::clone_sampler() const {
	MultiPhosphateSamplerOP multi_phosphate_sampler( new MultiPhosphateSampler( *pose_with_original_phosphates_ ) );
	multi_phosphate_sampler->set_screen_all( screen_all_ );
	multi_phosphate_sampler->set_force_phosphate_instantiation( force_phosphate_instantiation_ );
	multi_phosphate_sampler->set_phosphate_move_list( phosphate_move_list_ );
	return multi_phosphate_sampler;
}

//////////////////////////////
void
MultiPhosphateSampler::initialize_parameters(){
	screen_all_ = true ;
	force_phosphate_instantiation_ = false;
	phosphate_takeoff_donor_distance_cutoff2_ = 8.0 * 8.0;
}

////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::sample_phosphates(){

	if ( scorefxn_ == 0 ) scorefxn_ = get_phosphate_scorefxn();

	phosphate_move_list_ = initialize_phosphate_move_list( *phosphate_sample_pose_ );
	copy_over_phosphate_variants( *phosphate_sample_pose_, *pose_with_original_phosphates_, phosphate_move_list_ );

	actual_phosphate_move_list_ = check_moved( phosphate_move_list_, *phosphate_sample_pose_ );

	// note that the mover will instantiate or virtualize phosphates as it sees fit.
	instantiated_some_phosphate_ = false;
	for ( Size i = 1; i <= actual_phosphate_move_list_.size(); i++ ) {
		PhosphateMover phosphate_mover( actual_phosphate_move_list_[i], scorefxn_ );
		phosphate_mover.set_force_phosphate_instantiation( force_phosphate_instantiation_ );
		phosphate_mover.apply( *phosphate_sample_pose_ );
		if ( phosphate_mover.instantiated_phosphate() ) instantiated_some_phosphate_ = true;
	}
}

////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::sample_phosphates( pose::PoseOP & viewer_pose_op ){
	phosphate_sample_pose_ = viewer_pose_op;
	sample_phosphates();
}

////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::reset_to_original_pose(){
	phosphate_sample_pose_ = pose_with_original_phosphates_->clone();
}

////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::copy_phosphates( pose::Pose & mod_pose ) const {
	copy_over_phosphate_variants( mod_pose, *phosphate_sample_pose_, phosphate_move_list_ );
}

////////////////////////////////////////////////////////////////////
utility::vector1< PhosphateMove >
MultiPhosphateSampler::initialize_phosphate_move_list( pose::Pose & pose ){

	using namespace pose::full_model_info;
	utility::vector1< PhosphateMove > phosphate_move_list;

	if ( screen_all_ ) {
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();

		Size const nres_full_model = full_model_info.size();

		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( !pose.residue( n ).is_RNA() ) continue;
			Size const seqpos_in_full_model = res_list[ n ];
			if ( pose.residue_type( n ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
			if ( ( n == 1 || pose.fold_tree().is_cutpoint( n - 1 ) ) &&
					seqpos_in_full_model > 1 &&
					!cutpoint_open_in_full_model.has_value( seqpos_in_full_model - 1 ) &&
					!pose.residue_type( n ).has_variant_type( "CUTPOINT_UPPER" ) &&
					!pose.residue_type( n ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ) {
				if ( static_cast<int>( n ) == pose.fold_tree().root() ) {
					// Note -- how to fix this? The only way would be to have the fold tree start inside the residue, e.g.
					// at the base, rather than at phosphate. In current Rosetta, the only way to do that is to create a
					// separate virtual residue, and have a jump into this 'root' residue. Sigh.
					TR.Debug << "skipping pack of root phosphate at " << n << ". In the future, fix this, rhiju!" << std::endl;
					continue;
				}
				phosphate_move_list.push_back( PhosphateMove( n, FIVE_PRIME_PHOSPHATE ) );
			}
			if ( ( n == pose.size() || pose.fold_tree().is_cutpoint( n ) ) &&
					seqpos_in_full_model < nres_full_model &&
					!cutpoint_open_in_full_model.has_value( seqpos_in_full_model ) &&
					!pose.residue_type( n ).has_variant_type( "CUTPOINT_LOWER" ) ) {
				phosphate_move_list.push_back( PhosphateMove( n, THREE_PRIME_PHOSPHATE ) );
			}
		}

	} else {
		for ( Size i = 1; i <= five_prime_phosphate_res_input_.size(); i++ ) {
			phosphate_move_list.push_back(  PhosphateMove( five_prime_phosphate_res_input_[i], FIVE_PRIME_PHOSPHATE ) );
		}
		for ( Size i = 1; i <= three_prime_phosphate_res_input_.size(); i++ ) {
			phosphate_move_list.push_back(  PhosphateMove( three_prime_phosphate_res_input_[i], THREE_PRIME_PHOSPHATE ) );
		}
	}

	return phosphate_move_list;
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PhosphateMove >
MultiPhosphateSampler::check_moved( utility::vector1< PhosphateMove > phosphate_move_list,
	pose::Pose const & pose ) const {

	utility::vector1< PhosphateMove > actual_phosphate_move_list;

	if ( moving_partition_res_.size() == 0 ) return phosphate_move_list;

	utility::vector1< Size > partition_res1, partition_res2;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !pose.residue( n ).is_RNA() ) continue;
		if ( moving_partition_res_.has_value( n ) ) {
			partition_res1.push_back( n );
		} else {
			partition_res2.push_back( n );
		}
	}

	if ( force_phosphate_instantiation_ ) find_uninstantiated_phosphates( pose, phosphate_move_list_, actual_phosphate_move_list );

	// check phosphates that make cross-partition contacts.
	find_phosphate_contacts_other_partition( partition_res1, partition_res2, pose,
		phosphate_move_list, actual_phosphate_move_list );
	find_phosphate_contacts_other_partition( partition_res2, partition_res1, pose,
		phosphate_move_list, actual_phosphate_move_list );

	if ( !prepacked_ ) {
		// need to also add in any phosphates that made cross-contact partitions in *original* pose.
		// those need to be re-evaluated.
		find_phosphate_contacts_other_partition( partition_res1, partition_res2, *pose_with_original_phosphates_,
			phosphate_move_list, actual_phosphate_move_list );
		find_phosphate_contacts_other_partition( partition_res2, partition_res1, *pose_with_original_phosphates_,
			phosphate_move_list, actual_phosphate_move_list );
	}

	return actual_phosphate_move_list;
}

///////////////////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::find_uninstantiated_phosphates( pose::Pose const & pose,
	utility::vector1< PhosphateMove > const & phosphate_move_list,
	utility::vector1< PhosphateMove > & actual_phosphate_move_list ) const {

	for ( PhosphateMove const & phosphate_move : phosphate_move_list ) {
		if ( actual_phosphate_move_list.has_value( phosphate_move ) ) continue;
		Size const n = phosphate_move.rsd();
		if ( phosphate_move.terminus() == FIVE_PRIME_PHOSPHATE &&
				pose.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ) continue;
		if ( phosphate_move.terminus() == THREE_PRIME_PHOSPHATE &&
				pose.residue( n ).has_variant_type( core::chemical::THREE_PRIME_PHOSPHATE ) ) continue;
		actual_phosphate_move_list.push_back( phosphate_move );
	}
}

///////////////////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::find_phosphate_contacts_other_partition( utility::vector1< Size > const & partition_res1,
	utility::vector1< Size > const & partition_res2,
	pose::Pose const & pose,
	utility::vector1< PhosphateMove > const & phosphate_move_list,
	utility::vector1< PhosphateMove > & actual_phosphate_move_list ) const {

	for ( PhosphateMove const & phosphate_move : phosphate_move_list ) {
		if ( actual_phosphate_move_list.has_value( phosphate_move ) ) continue;
		Size const n = phosphate_move.rsd();
		if ( !partition_res1.has_value( n ) ) continue;
		Vector const phosphate_takeoff_xyz = (phosphate_move.terminus() == FIVE_PRIME_PHOSPHATE) ?
			pose.residue( n ).xyz( " C5'" ) :
			pose.residue( n ).xyz( " O3'" );
		if ( !check_other_partition_for_contact( pose, partition_res2, phosphate_takeoff_xyz ) ) continue;
		actual_phosphate_move_list.push_back( phosphate_move );
	}
}

////////////////////////////////////////////////////////////////////
bool
MultiPhosphateSampler::check_other_partition_for_contact( pose::Pose const & pose,
	utility::vector1< Size > const & other_partition_res,
	Vector const & takeoff_xyz ) const {

	for ( Size const m : other_partition_res ) {
		core::chemical::AtomIndices Hpos_polar = pose.residue_type( m ).Hpos_polar();
		for ( Size const q : Hpos_polar ) {
			if ( ( pose.residue( m ).xyz( q ) - takeoff_xyz ).length_squared() < phosphate_takeoff_donor_distance_cutoff2_ ) {
				return true;
			}
		}
	}

	return false;
}


////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::initialize_by_prepack( pose::Pose & pose,
	Size const moving_res ){
	initialize_by_prepack( pose, make_vector1( moving_res ) );
}

////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::initialize_by_prepack( pose::Pose & pose,
	utility::vector1< Size > const & moving_res_list ){
	moving_partition_res_ = figure_out_moving_partition_res( pose, moving_res_list );
	do_prepack( pose, moving_res_list );
}

////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::do_prepack( pose::Pose & pose,
	utility::vector1< Size > const & moving_res_list ) {
	// sanity check:
	runtime_assert( moving_partition_res_ == figure_out_moving_partition_res( pose, moving_res_list ) );

	PoseOP pose_to_split = pose.clone();
	pose_to_split->remove_constraints(); // floating point errors if coordinate constraints are in there.
	split_pose( *pose_to_split, moving_res_list );

	utility::vector1< Size > moving_partition_res_save = moving_partition_res_;
	moving_partition_res_.clear(); // ensures that all phosphates will be packed -- could replace with inputted obligate_pack_res.

	pose_with_original_phosphates_ = pose_to_split->clone(); // dummy, will be replaced...
	sample_phosphates( pose_to_split );

	moving_partition_res_ = moving_partition_res_save;
	copy_over_phosphate_variants( pose, *pose_to_split, phosphate_move_list_ );

	reset( pose );
	prepacked_ = true;
}

////////////////////////////////////////////////////////////////////////
pose::Pose &
MultiPhosphateSampler::pose(){ return *phosphate_sample_pose_; }

////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::reset( pose::Pose const & pose ){
	pose_with_original_phosphates_ = pose.clone();
	phosphate_sample_pose_ = pose.clone();
}

////////////////////////////////////////////////////////////////////////
void
MultiPhosphateSampler::set_scorefxn( scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn;
}


} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols
