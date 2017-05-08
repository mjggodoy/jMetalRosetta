// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh
/// @brief  Statistically derived chemical mapping score
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_ChemicalMappingEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_ChemicalMappingEnergy_hh

// Unit Headers
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

namespace core {
namespace scoring {
namespace rna {
namespace data {


class RNA_ChemicalMappingEnergy : public methods::WholeStructureEnergy {
public:

	typedef methods::WholeStructureEnergy parent;


	RNA_ChemicalMappingEnergy();

	~RNA_ChemicalMappingEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	Real
	calculate_energy( pose::Pose & pose,
		bool const use_low_res = false,
		bool const rna_base_pair_computed = false ) const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	virtual Size version() const { return 0; }

private:

	RNA_DMS_Potential & DMS_potential_;
	RNA_DMS_LowResolutionPotential & DMS_low_resolution_potential_;

};


} //data
} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
