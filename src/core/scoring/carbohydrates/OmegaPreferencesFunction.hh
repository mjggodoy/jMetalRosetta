// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/OmegaPreferencesFunction.hh
/// @brief   Method declarations for OmegaPreferencesFunction.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_HH
#define INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_HH

// Unit Header
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.fwd.hh>
#include <core/scoring/carbohydrates/OmegaPreferenceType.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>

namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief  Struct for Omega-Preference dihedral sampling.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
struct OmegaPreferenceSamplingData {
	OmegaPreferenceType preference_type;
	Real step_size;

	utility::vector1< Real > probabilities;
	utility::vector1< Angle > angles;
};


/// @details  TBD
class OmegaPreferencesFunction : public utility::pointer::ReferenceCount {
public:  // Standard Methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	OmegaPreferencesFunction();

	virtual ~OmegaPreferencesFunction();


public:  // Other Public Methods //////////////////////////////////////////////
	Energy operator()( OmegaPreferenceType preference, Angle x ) const;

	Real evaluate_derivative( OmegaPreferenceType preference, Angle x ) const;

public: // Dihedral Sampling Methods //////////////////////////////////////////
	OmegaPreferenceSamplingData const & get_sampling_data( OmegaPreferenceType const type ) const;

	///@brief Set up CHI sampling data.
	void setup_for_sampling( core::Real step_size = 0.1 );

	bool sampling_data_setup() const;

	bool sampling_data_setup( OmegaPreferenceType const linkage_type ) const;

private:  // Private methods //////////////////////////////////////////////////
	void set_parameters( OmegaPreferenceType preference, Angle x ) const;

	Energy evaluate_function( OmegaPreferenceType preference, Angle x ) const;

	void init();


private:  // Private Data /////////////////////////////////////////////////////
	// Parabolic function parameters as defined in:
	// A.K. Nivedha et al. JCTC 2016, 12, 892-901.
	// For now, we will just use their k value, but I feel like we could do better.  ~Labonte
	core::Real const k_ = 0.0025;  // constant to determine parabola "width"
	mutable core::Angle theta_;  // location of the vertex of the parabola
	mutable core::Real b_;  // relative energy difference of minima from minimum

private:  // Formatted Sampling Data //////////////////////////////////////////
	std::map< OmegaPreferenceType, OmegaPreferenceSamplingData > dihedral_sampling_data_;


};  // class OmegaPreferencesFunction

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif // INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_HH
