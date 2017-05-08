// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/RNA_NucleosideStepWiseSampler.cc
/// @brief Generate rotamers for one RNA nucleoside (pucker + glycosidic chi).
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_NucleosideStepWiseSampler.hh>

// Package headers
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/util.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.rna.RNA_NucleosideStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_NucleosideStepWiseSampler::RNA_NucleosideStepWiseSampler(
	core::Size const rsd_id,
	PuckerState const pucker_state, //ANY_PUCKER, NORTH, SOUTH, NO_PUCKER
	ChiState const base_state //ANY_CHI, ANTI, SYN, NO_CHI
):
	StepWiseSamplerSizedAny(),
	rsd_id_( rsd_id ),
	base_state_( base_state ),
	extra_chi_( false ),
	skip_same_pucker_( true ),
	idealize_coord_( true ),
	fast_( false ),
	bin_size_( 20 )
{
	runtime_assert( pucker_state <= 2 );
	runtime_assert( base_state <= 3 );

	if ( pucker_state == ANY_PUCKER ) {
		pucker_states_.push_back( NORTH );
		pucker_states_.push_back( SOUTH );
	} else {
		pucker_states_.push_back( pucker_state );
	}
}
////////////////////////////////////////////////////////////////////////////
void RNA_NucleosideStepWiseSampler::init() {
	clear_rotamer();

	//Setup the rotamer samplers
	for ( auto const & pucker_state : pucker_states_ ) {
		StepWiseSamplerSizedCombOP new_rotamer_agg( new StepWiseSamplerSizedComb );

		/////Chi rotamers/////
		if ( base_state_ != NO_CHI ) {
			RNA_ChiStepWiseSamplerOP chi_rotamer( new RNA_ChiStepWiseSampler(
				rsd_id_, pucker_state, base_state_ ) );
			chi_rotamer->set_bin_size( bin_size_ );
			chi_rotamer->set_extra_chi( extra_chi_);
			if ( fast_ ) chi_rotamer->set_max_range( 1 );
			new_rotamer_agg->add_external_loop_rotamer( chi_rotamer );
		}

		/////Pucker rotamers/////
		RNA_SugarStepWiseSamplerOP pucker_rotamer( new RNA_SugarStepWiseSampler(
			rsd_id_, pucker_state ) );
		pucker_rotamer->set_skip_same_pucker( skip_same_pucker_ );
		pucker_rotamer->set_idealize_coord( idealize_coord_ );
		new_rotamer_agg->add_external_loop_rotamer( pucker_rotamer );

		/////Add to the this sampler/////
		add_external_loop_rotamer( new_rotamer_agg );
	}

	StepWiseSamplerSizedAny::init();
}
//////////////////////////////////////////////////////////////////////////
/// @brief Name of the class
std::string
RNA_NucleosideStepWiseSampler::get_name() const {
	return "RNA_NucleosideStepWiseSampler residue:" + utility::to_string(rsd_id_);
}
//////////////////////////////////////////////////////////////////////////
} //rna
} //sampler
} //stepwise
} //protocols
