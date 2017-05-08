// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/TwelveANeighborGraph.fwd.hh
/// @brief  Twelve Angstrom Neighbor graph class forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_TwelveANeighborGraph_fwd_hh
#define INCLUDED_core_scoring_TwelveANeighborGraph_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class TwelveANeighborNode;
class TwelveANeighborEdge;
class TwelveANeighborGraph;

typedef utility::pointer::shared_ptr< TwelveANeighborGraph > TwelveANeighborGraphOP;
typedef utility::pointer::shared_ptr< TwelveANeighborGraph const > TwelveANeighborGraphCOP;

} // scoring
} // core

#endif
