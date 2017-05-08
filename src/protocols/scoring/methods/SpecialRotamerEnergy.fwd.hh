// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/scoring/methods/SpecialRotamerEnergy.hh
/// @brief  Adds a bonus to any rotamer that is flagged
/// @author sthyme, sthyme@gmail.com, Feb 2010


#ifndef INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_fwd_hh
#define INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace scoring {
namespace methods {

class SpecialRotamerEnergy;

typedef utility::pointer::shared_ptr< SpecialRotamerEnergy > SpecialRotamerEnergyOP;

} // methods
} // scoring
} // protocols


#endif // INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_FWD_HH
