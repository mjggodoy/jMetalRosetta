// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/RNA_KIC_Sampler.cc
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh>
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
using namespace protocols::stepwise::sampler::screener;
using namespace protocols::stepwise::sampler;

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.rna.RNA_KIC_Sampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

RNA_KIC_Sampler::RNA_KIC_Sampler(
	core::pose::PoseOP const & ref_pose,
	core::Size const moving_suite,
	core::Size const chainbreak_suite
):
	StepWiseSampler(),
	ref_pose_( ref_pose ),
	moving_suite_( moving_suite ),
	chainbreak_suite_( chainbreak_suite ),
	pucker_state_( ANY_PUCKER ), // ANY_PUCKER, NORTH, SOUTH, NO_PUCKER
	base_state_( ANY_CHI ), // ANY_CHI, ANTI, SYN, NO_CHI
	sample_nucleoside_( moving_suite + 1 ), // default, may be replaced.
	bin_size_( 20 ),
	max_tries_( 100 ),
	verbose_( false ),
	extra_epsilon_( false ),
	extra_chi_( false ),
	skip_same_pucker_( false ),
	idealize_coord_( false ),
	torsion_screen_( true ),
	random_chain_closed_( true ),
	screener_( screener::RNA_TorsionScreenerOP( new RNA_TorsionScreener ) )
{
	StepWiseSampler::set_random( false );
}

RNA_KIC_Sampler::~RNA_KIC_Sampler() {}

//////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::init() {
	using namespace core::id;

	////////// Backbone StepWiseSampler //////////
	bb_rotamer_ = StepWiseSamplerSizedCombOP( new StepWiseSamplerSizedComb );

	RNA_FittedTorsionInfo const torsion_info;

	TorsionList full_torsions;
	add_values_from_center( full_torsions, 0, 180, bin_size_ );
	//Avoid modeler both -180 and 180 deg
	Real const epsil =  0.001; //Arbitary small number
	if ( 360 - ( full_torsions.back() - full_torsions.front() ) < epsil ) {
		full_torsions.pop_back();
	}

	/////Epsilon1 rotamers/////
	//default: center +- 20 deg
	//extra_epsilon: center +- 60 deg
	TorsionList epsilon_torsions;
	Size const pucker1( assign_pucker( *ref_pose_, moving_suite_ ) );
	Real center = ( pucker1 == NORTH ) ? torsion_info.epsilon_north() :
		torsion_info.epsilon_south();
	Real max_range = 20;
	if ( extra_epsilon_ ) {
		max_range = 60;
		//Choice made by Parin, to cover the uneven tails of
		//epsilons distributions by extra_epsilon mode.
		center += ( pucker1 == NORTH ) ? -20 : 20;
	}
	add_values_from_center( epsilon_torsions, center, max_range, bin_size_ );
	StepWiseSamplerOneTorsionOP epsilon_rotamer( new StepWiseSamplerOneTorsion(
		TorsionID( moving_suite_, BB, EPSILON ), epsilon_torsions ) );
	bb_rotamer_->add_external_loop_rotamer( epsilon_rotamer );
	/////Zeta1 rotamers/////
	StepWiseSamplerOneTorsionOP zeta_rotamer( new StepWiseSamplerOneTorsion(
		TorsionID( moving_suite_, BB, ZETA ), full_torsions ) );
	bb_rotamer_->add_external_loop_rotamer( zeta_rotamer );
	/////Alpha1 rotamers/////
	StepWiseSamplerOneTorsionOP alpha1_rotamer( new StepWiseSamplerOneTorsion(
		TorsionID( moving_suite_ + 1, BB, ALPHA ), full_torsions ) );
	bb_rotamer_->add_external_loop_rotamer( alpha1_rotamer );
	/////Alpha2 rotamers/////
	StepWiseSamplerOneTorsionOP alpha2_rotamer( new StepWiseSamplerOneTorsion(
		TorsionID( chainbreak_suite_ + 1, BB, ALPHA ), full_torsions ) );
	bb_rotamer_->add_external_loop_rotamer( alpha2_rotamer );
	/////Pucker rotamers/////
	if ( pucker_state_ != NO_PUCKER ) {
		RNA_SugarStepWiseSamplerOP pucker_rotamer( new RNA_SugarStepWiseSampler(
			sample_nucleoside_, pucker_state_ ) );
		pucker_rotamer->set_skip_same_pucker( skip_same_pucker_ );
		pucker_rotamer->set_idealize_coord( idealize_coord_ );
		bb_rotamer_->add_external_loop_rotamer( pucker_rotamer );
	}
	bb_rotamer_->init();

	////////// Loop Closer //////////
	loop_closer_ = RNA_KinematicCloserOP( new RNA_KinematicCloser(
		*ref_pose_, moving_suite_, chainbreak_suite_ ) );
	loop_closer_->set_verbose( verbose_ );

	////////// Chi StepWiseSampler //////////
	if ( base_state_ != NO_CHI ) {
		chi_rotamer_ = RNA_ChiStepWiseSamplerOP( new RNA_ChiStepWiseSampler(
			sample_nucleoside_, NORTH/*arbitary*/, base_state_ ) );
		chi_rotamer_->set_extra_chi( extra_chi_ );
		chi_rotamer_->set_bin_size( bin_size_ );
		chi_rotamer_->init();
	}

	set_init( true );
	set_random( random() ); //update chi rotamer and init loop closer
}
//////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::reset() {
	runtime_assert( is_init() );
	if ( random() && not_end() ) {
		++( *this );
	} else {
		bb_rotamer_->reset();
		bb_rotamer_->apply( *ref_pose_ );
		loop_closer_->init();
		get_next_valid_bb();
	}
}
//////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		random_chain_closed_ = false;
		for ( Size i = 1; i <= max_tries_; ++i ) {
			++( *bb_rotamer_ );
			bb_rotamer_->apply( *ref_pose_ );

			loop_closer_->init();
			if ( !loop_closer_->not_end() ) continue;
			++( *loop_closer_ );

			if ( torsion_screen_ ) {
				loop_closer_->apply( *ref_pose_);
				if ( !screener_->screen( *ref_pose_, chainbreak_suite_ ) ||
						!screener_->screen( *ref_pose_, moving_suite_ ) ) {
					continue;
				}
			}

			if ( base_state_ != NO_CHI ) ++( *chi_rotamer_ );
			random_chain_closed_ = true;
			return;
		}
		TR.Debug << "Chain not closable after " << max_tries_ << " tries!" << std::endl;

	} else {
		if ( base_state_ != NO_CHI ) {
			++( *chi_rotamer_ );
			if ( chi_rotamer_->not_end() ) return;
		}
		++( *loop_closer_ );
		get_next_valid_bb();
	}
}
///////////////////////////////////////////////////////////////////////////
bool RNA_KIC_Sampler::not_end() const {
	runtime_assert( is_init() );
	// if ( random() ) return random_chain_closed_;
	return bb_rotamer_->not_end();
}
///////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::apply() {
	apply( *ref_pose_ );
}
///////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::apply( pose::Pose & pose ) {
	runtime_assert( is_init() );
	if ( random() && !random_chain_closed_ ) {
		TR.Debug << "Warning: Chain was not closable! Not doing anything to pose." << std::endl;
		return;
	}
	if ( base_state_ != NO_CHI ) chi_rotamer_->apply( pose );
	loop_closer_->apply( pose );
	bb_rotamer_->apply( pose );
}
///////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::set_random( bool const setting ) {
	StepWiseSampler::set_random( setting );
	if ( is_init() ) {
		bb_rotamer_->set_random( setting );
		loop_closer_->set_random( setting );
		if ( base_state_ != NO_CHI ) chi_rotamer_->set_random( setting );
		reset();
	}
}
///////////////////////////////////////////////////////////////////////////
void RNA_KIC_Sampler::get_next_valid_bb() {
	while ( true ) {
		for ( ; loop_closer_->not_end(); ++( *loop_closer_ ) ) {
			if ( torsion_screen_ ) {
				loop_closer_->apply( *ref_pose_);
				if ( !screener_->screen( *ref_pose_, chainbreak_suite_ ) ||
						!screener_->screen( *ref_pose_, moving_suite_ ) ) {
					continue;
				}
			}

			if ( base_state_ != NO_CHI ) {
				PuckerState const pucker_state( assign_pucker( *ref_pose_, moving_suite_ + 1 ) );
				if ( chi_rotamer_->pucker_state() != pucker_state ) {
					chi_rotamer_->set_pucker_state( pucker_state );
					chi_rotamer_->init();
				} else {
					chi_rotamer_->reset();
				}
			}
			return;
		}

		++( *bb_rotamer_ );
		if ( !bb_rotamer_->not_end() ) break;
		bb_rotamer_->apply( *ref_pose_ );
		loop_closer_->init();
	}
}
///////////////////////////////////////////////////////////////////////////
/// @brief Name of the class
std::string
RNA_KIC_Sampler::get_name() const {
	return "RNA_KIC_Sampler moving_suite:" + utility::to_string(moving_suite_) + " chainbreak_suite:" +
		utility::to_string(chainbreak_suite_);
}

} //rna
} //sampler
} //stepwise
} //protocols
