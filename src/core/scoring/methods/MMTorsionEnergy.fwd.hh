// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMTorsionEnergy.fwd.hh
/// @brief  molecular mechanics torsion energy forward declaration
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_methods_MMTorsionEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_MMTorsionEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class MMTorsionEnergy;

typedef  utility::pointer::weak_ptr< MMTorsionEnergy > MMTorsionEnergyAP;
typedef  utility::pointer::weak_ptr< MMTorsionEnergy const > MMTorsionEnergyCAP;
typedef  utility::pointer::shared_ptr< MMTorsionEnergy > MMTorsionEnergyOP;
typedef  utility::pointer::shared_ptr< MMTorsionEnergy const > MMTorsionEnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_MMTorsionEnergy_HH
