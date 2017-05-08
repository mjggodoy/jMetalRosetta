// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/FullatomDisulfideEnergy.fwd.hh
/// @brief  Disulfide Energy class forward declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergy_fwd_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace disulfides {

class FullatomDisulfideEnergy;

typedef utility::pointer::shared_ptr< FullatomDisulfideEnergy > FullatomDisulfideEnergyOP;
typedef utility::pointer::shared_ptr< FullatomDisulfideEnergy const > FullatomDisulfideEnergyCOP;


} // namespace disulfides
} // namespace scoring
} // namespace core

#endif
