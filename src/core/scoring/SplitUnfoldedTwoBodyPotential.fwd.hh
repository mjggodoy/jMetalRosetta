// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/SplitUnfoldedTwoBodyPotential.cc
/// @brief  Reads in and stores the two body energies for each residue type for the split unfolded energy.
/// @author Riley Simmons-Edler (rse231@nyu.edu

#ifndef INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_fwd_hh
#define INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class SplitUnfoldedTwoBodyPotential;

typedef  utility::pointer::shared_ptr< SplitUnfoldedTwoBodyPotential > SplitUnfoldedTwoBodyPotentialOP;

} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_FWD_HH
