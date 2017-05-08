// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_BetaAlphaBetaMotif_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_BetaAlphaBetaMotif_fwd_hh

// utitlity headers
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

class BetaAlphaBetaMotif;
class BetaAlphaBetaMotifSet;

typedef utility::pointer::shared_ptr< BetaAlphaBetaMotif > BetaAlphaBetaMotifOP;
typedef utility::pointer::shared_ptr< BetaAlphaBetaMotif const > BetaAlphaBetaMotifCOP;

typedef utility::vector1< BetaAlphaBetaMotifOP > BetaAlphaBetaMotifs;
typedef utility::pointer::shared_ptr< BetaAlphaBetaMotifSet > BetaAlphaBetaMotifSetOP;
typedef utility::pointer::shared_ptr< BetaAlphaBetaMotifSet const > BetaAlphaBetaMotifSetCOP;


} // topology
} // fldsgn
} // protocols

#endif
