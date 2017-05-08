// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/OmegaPreferencesFunction.fwd.hh
/// @brief   Forward declarations for OmegaPreferencesFunction.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_FWD_HH
#define INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief  A function to penalize glycosidic omega angles that do not demonstrate the gauche effect.
class OmegaPreferencesFunction;

typedef utility::pointer::shared_ptr< OmegaPreferencesFunction > OmegaPreferencesFunctionOP;
typedef utility::pointer::shared_ptr< OmegaPreferencesFunction const> OmegaPreferencesFunctionCOP;

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif // INCLUDED_core_scoring_carbohydrates_OmegaPreferencesFunction_FWD_HH
