// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh
/// @brief Generate glycosidic chi rotamers for RNA.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_RNA_ChiStepWiseSampler_HH
#define INCLUDED_protocols_sampler_rna_RNA_ChiStepWiseSampler_HH

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneTorsion.hh>

// Project headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/util.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

class RNA_ChiStepWiseSampler : public StepWiseSamplerOneTorsion {
public:
	RNA_ChiStepWiseSampler(
		core::Size const rsd_id,
		core::chemical::rna::PuckerState const pucker_state,
		core::chemical::rna::ChiState const base_state
	);

	/// @brief Initialization
	void init();

	/// @brief Set the residue id
	void set_rsd_id( core::Size const setting ) {
		set_and_reinit( rsd_id_, setting );
	}

	/// @brief Get the pucker state (NORTH / SOUTH)
	core::chemical::rna::PuckerState pucker_state() { return pucker_state_; }

	/// @brief Set the pucker state (NORTH / SOUTH)
	void set_pucker_state( core::chemical::rna::PuckerState const setting ) {
		set_and_reinit( pucker_state_, setting );
	}

	/// @brief Get the base state (WHATEVER / ANTI / SYN)
	core::chemical::rna::ChiState base_state() { return base_state_; }

	/// @brief Set the base state (WHATEVER / ANTI / SYN)
	void set_base_state( core::chemical::rna::ChiState const setting ) {
		set_and_reinit( base_state_, setting );
	}

	/// @brief Set the bin_size (default: 20)
	void set_bin_size( core::Real const setting ) {
		set_and_reinit( bin_size_, setting );
	}

	/// @brief Set the max_range of modeler (default: +-20)
	void set_max_range( core::Real const setting ) {
		set_and_reinit( max_range_, setting );
	}

	/// @brief Set use extra chi (+-60)
	void set_extra_chi( bool const setting ) {
		if ( setting ) {
			set_max_range( 60 );
		} else {
			set_max_range( 20 );
		}
	}

	/// @brief Name of the class
	std::string get_name() const { return "RNA_ChiStepWiseSampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::RNA_CHI; }

private:
	core::Size rsd_id_;
	core::chemical::rna::ChiState base_state_;
	core::chemical::rna::PuckerState pucker_state_;
	core::Real bin_size_, max_range_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool extra_chi_;

	core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;
};

} //rna
} //sampler
} //stepwise
} //protocols

#endif
