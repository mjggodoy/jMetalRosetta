// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinRMSDNoSuperpositionFeatures.fwd.hh
/// @brief  report (no superposition) RMSD similarity of a protein structure against supplied reference structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_ProteinRMSDNoSuperpositionFeatures_fwd_hh
#define INCLUDED_protocols_features_ProteinRMSDNoSuperpositionFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class ProteinRMSDNoSuperpositionFeatures;
typedef utility::pointer::shared_ptr< ProteinRMSDNoSuperpositionFeatures > ProteinRMSDNoSuperpositionFeaturesOP;
typedef utility::pointer::shared_ptr< ProteinRMSDNoSuperpositionFeatures const > ProteinRMSDNoSuperpositionFeaturesCOP;

}//features
}//protocols

#endif // include guard
