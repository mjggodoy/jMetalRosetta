// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/align/StepWiseClusterer.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/align/StepWiseClusterer.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/monte_carlo/util.hh> // for output_to_silent_file
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.align.StepWiseClusterer" );

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Very basic and align clusterer for StepWise Assembly & Monte Carlo.
//
//  Two interfaces ...
//
//     apply( pose )  - 'real-time' clustering, useful in SWA for preventing memory explosion.
//
//     set_pose_list( pose_list ), and then cluster()
//
// Takes advantage of StepWisePoseAligner, which holds all rules for defining superposition &
//  calc rmsd atoms (e.g., all heavy atoms for RNA, just backbone by default for protein,
//  no virtual atoms, and smart selection of superposition atoms based on fixed domain information
//  in full_model_info. )
//
//    -- rhiju, 2014
//
// Did not carry over:
//   'auto-tune' of RMSD (see StepWiseLegacyClusterer), explicit encoding of phosphate-base-phosphate
//   RMSDs (was not default anyway).
//
//////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

//Constructor
StepWiseClusterer::StepWiseClusterer():
	max_decoys_( 0 ),  // will be reset below.
	rmsd_( 0.0 ),
	cluster_rmsd_( 0.0 ), // will be reset below.
	score_diff_cut_( 0.0 ),
	do_checks_( true ),
	assume_atom_ids_invariant_( false ),
	initialized_( false ),
	count_( 0 ),
	cluster_size_( 0 ),
	output_cluster_size_( false )
{
}

StepWiseClusterer::StepWiseClusterer( options::StepWiseModelerOptionsCOP options ):
	max_decoys_( options->sampler_num_pose_kept() ),
	rmsd_( 0.0 ),
	cluster_rmsd_( options->cluster_rmsd() ),
	score_diff_cut_( 0.0 ),
	do_checks_( true ),
	assume_atom_ids_invariant_( false ),
	initialized_( false ),
	count_( 0 ),
	cluster_size_( 0 ),
	output_cluster_size_( options->output_cluster_size() ),
	silent_file_( options->sampler_silent_file() )
{
	runtime_assert( !options->cluster_by_all_atom_rmsd() ); // for now.
}

//Destructor
StepWiseClusterer::~StepWiseClusterer()
{}

//////////////////////////////////////////////////////////////////////////
void
StepWiseClusterer::apply( pose::Pose const & pose ) {
	initialize_parameters( pose );

	if ( silent_file_.size() > 0 )  monte_carlo::output_to_silent_file( "X_" + ObjexxFCL::lead_zero_string_of( count_++, 6 ), silent_file_, pose );

	if ( check_screen_and_kick_out_displaced_model( pose ) ) {
		//TR << "Inserting into list. Score:  " << total_energy_from_pose( pose ) << "  rmsd " << rmsd_ << "  cutoff: " << cluster_rmsd_ << "  list size " << pose_list_.size() << "  max_decoys " << max_decoys_ << std::endl;
		pose_list_.push_back( pose.clone() );
		if ( output_cluster_size_ ) setPoseExtraScore( *pose_list_[ pose_list_.size() ], "cluster_size", cluster_size_ );
		sort_pose_list();
	} else {
		//   TR << TR.Red << "not including in pose list. Score:  " << total_energy_from_pose( pose ) << "  rmsd " << rmsd_ << "  cutoff: " << cluster_rmsd_ << "  list size " << pose_list_.size() << "  max_decoys " << max_decoys_ << TR.Reset << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseClusterer::sort_pose_list() {
	sort( pose_list_.begin(), pose_list_.end(), core::pose::sort_pose_by_score );
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseClusterer::cluster()
{
	sort_pose_list();
	utility::vector1< pose::PoseOP > starting_pose_list = pose_list_;
	pose_list_.clear();

	for ( pose::PoseOP const & pose : starting_pose_list ) {
		runtime_assert( pose != NULL );
		if ( check_screen_and_kick_out_displaced_model( *pose ) ) pose_list_.push_back( pose );
		sort_pose_list();
	}
}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseClusterer::check_screen_and_kick_out_displaced_model( pose::Pose const & pose ) {

	Real const & pose_score = total_energy_from_pose( pose );

	// if we've got a full list, don't even consider the pose unless its energy is reasonable.
	if ( size() == max_decoys_ && size() > 0 ) {
		Real const & worst_score_in_list = total_energy_from_pose( *pose_list_[ size() ] );
		if ( pose_score > worst_score_in_list ) return false;
	}

	if ( score_diff_cut_ > 0.0 && size() > 0 ) {
		Real const & best_score = total_energy_from_pose( *pose_list_[1] );
		if ( pose_score - best_score > score_diff_cut_ ) return false;
	}

	// is the new pose a neighbor of an old pose? In that case should it replace the old one?
	for ( core::Size n = 1; n <= pose_list_.size(); n++ ) {
		if ( check_for_closeness( pose, *pose_list_[n] ) ) {
			if ( hasPoseExtraScore( *pose_list_[ n ], "cluster_size" ) ) getPoseExtraScore( *pose_list_[ n ], "cluster_size", cluster_size_ );
			cluster_size_ += 1;
			if ( pose_score < total_energy_from_pose( *pose_list_[n] ) ) {
				kick_out_pose_at_idx( n );
				return true;
			} else {
				if ( output_cluster_size_ ) setPoseExtraScore( *pose_list_[ n ], "cluster_size", cluster_size_ );
				return false;
			}
		}
	}

	// made it here -- pose has low enough energy to merit going on the list.
	// if we've got a full list,  kick out the worse scoring pose to make room.
	if ( size() > max_decoys_ ) kick_out_pose_at_idx( size() );
	cluster_size_ = 1;
	return true;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseClusterer::kick_out_pose_at_idx( core::Size const n ){

	runtime_assert( n >= 1 && n <= size() );
	utility::vector1< pose::PoseOP > const starting_pose_list = pose_list_;

	pose_list_.clear();
	for ( core::Size k = 1; k <= starting_pose_list.size(); k++ ) {
		if ( k != n ) pose_list_.push_back( starting_pose_list[k] );
	}
}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseClusterer::check_for_closeness( pose::Pose const & pose1, pose::Pose const & pose2 ) {
	bool set_rmsd( false );
	if ( assume_atom_ids_invariant_ ) {
		if ( pose_aligner_ == 0 ) {
			pose_aligner_ = StepWisePoseAlignerOP( new StepWisePoseAligner( pose2 ) );
			pose_aligner_->set_user_defined_calc_rms_res( calc_rms_res_ );
			pose_aligner_->set_root_partition_res( figure_out_root_partition_res( pose2, calc_rms_res_ ) );
			pose_aligner_->initialize( pose1 );
		}
		if ( pose_aligner_->check_matching_atom_names( pose1, pose2, false /*verbose*/ ) ) { // terminal phosphate variants can reorder atom order...
			rmsd_ = pose_aligner_->get_rmsd_no_superimpose( pose1, pose2, do_checks_ /*check_align_at_superimpose_res*/ );
			set_rmsd = true;
		}
	}

	if ( !set_rmsd ) {
		// this is the super-robust way to calculate RMSD.
		rmsd_ = get_rmsd( pose1, pose2, calc_rms_res_,
			do_checks_, //false /*check align at superimpose res*/,
			do_checks_ /*check switch*/ );
	}

	if ( rmsd_ <= cluster_rmsd_ ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseClusterer::initialize_parameters( pose::Pose const & pose ) {

	if ( initialized_ ) return;

	// default parameter setup, using 'historical' values tuned for protein/RNA applications.
	if ( calc_rms_res_.size() > 0 ) {
		if ( pose.residue_type( calc_rms_res_[1] ).is_protein() ) {
			if ( max_decoys_ == 0   ) max_decoys_ = 400;
			if ( cluster_rmsd_ == 0 ) cluster_rmsd_ = 0.1;
		} else if ( pose.residue_type( calc_rms_res_[1] ).is_RNA() ) {
			if ( max_decoys_ == 0   ) max_decoys_ = 108;
			if ( cluster_rmsd_ == 0 ) cluster_rmsd_ = 0.5;
		} else {
			utility_exit_with_message( "Not sure how to calculate RMSD of something that is not RNA or protein" );
			runtime_assert( max_decoys_ > 0 );
		}
	} else {
		runtime_assert(  max_decoys_ > 0 );
	}

	initialized_ = true;
}

} //align
} //modeler
} //stepwise
} //protocols
