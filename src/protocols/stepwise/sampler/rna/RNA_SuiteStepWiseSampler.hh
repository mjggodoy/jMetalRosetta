// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.hh
/// @brief Generate rotamers for one RNA suite (from residue i to i+1).
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_RNA_SuiteStepWiseSampler_HH
#define INCLUDED_protocols_sampler_rna_RNA_SuiteStepWiseSampler_HH

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedAny.hh>

// Project headers
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

class RNA_SuiteStepWiseSampler : public StepWiseSamplerSizedAny {
public:
	using StepWiseSampler::TorsionList;

	RNA_SuiteStepWiseSampler(
		core::Size const rsd_id,
		core::chemical::rna::PuckerState const pucker_state_lower, //WHATEVER, NORTH, SOUTH
		core::chemical::rna::PuckerState const pucker_state_upper,
		core::chemical::rna::ChiState const base_state_lower, //WHATEVER, ANTI, SYN, NONE
		core::chemical::rna::ChiState const base_state_upper
	);

	/// @brief Initialization wrapper
	void init();

	/// @brief Initialization for standard modeler
	void init_standard();

	/// @brief Initialization for fast mode
	void init_fast();

	/// Set functions
	void set_sample_nucleoside_lower( bool const setting ) {
		set_and_reinit( sample_nucleoside_lower_, setting );
	}

	void set_sample_nucleoside_upper( bool const setting ) {
		set_and_reinit( sample_nucleoside_upper_, setting );
	}

	void set_extra_epsilon( bool const setting ) {
		set_and_reinit( extra_epsilon_, setting );
	}

	void set_extra_beta( bool const setting ) {
		set_and_reinit( extra_beta_, setting );
	}

	void set_extra_chi( bool const setting ) {
		set_and_reinit( extra_chi_, setting );
	}

	void set_skip_same_pucker( bool const setting ) {
		set_and_reinit( skip_same_pucker_, setting );
	}

	void set_idealize_coord( bool const setting ) {
		set_and_reinit( idealize_coord_, setting );
	}

	/// @brief Fast mode: just sample popular center torsions
	void set_fast( bool const setting ) {
		set_and_reinit( fast_, setting );
	}

	void set_bin_size( core::Real const setting ) {
		set_and_reinit( bin_size_, setting );
	}

	/// @brief Name of the class
	std::string get_name() const;

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::RNA_SUITE; }

private:
	TorsionList fast_sample_torsions_from_info(
		utility::vector1<core::chemical::rna::GaussianParameter>
		const & params );

	core::Size const rsd_id_;
	core::chemical::rna::ChiState base_state_lower_, base_state_upper_;

	bool sample_nucleoside_lower_, sample_nucleoside_upper_,
		extra_epsilon_, extra_beta_, extra_chi_,
		skip_same_pucker_, idealize_coord_, fast_;

	core::Real bin_size_;

	utility::vector1<core::chemical::rna::PuckerState> pucker_states_lower_, pucker_states_upper_;
};

} //rna
} //sampler
} //stepwise
} //protocols

#endif
