// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_ContextIndependentOneBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextIndependentOneBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/EnergyMap.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class ContextIndependentOneBodyEnergy : public OneBodyEnergy {
public:
	typedef OneBodyEnergy parent;

public:
	/// @brief Constructor with an EnergyMethodCreator to inform the parent class
	/// which ScoreTypes this EnergyMethod is responsible for computing.
	ContextIndependentOneBodyEnergy( EnergyMethodCreatorOP );

	/// @brief Returns the ci_1b element of the EnergyMethodType enumeration; this
	/// method should NOT be overridden by derived classes.
	virtual
	EnergyMethodType
	method_type() const;


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap & emap
	) const = 0;


};

} // methods
} // scoring
} // core


#endif
