// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DesignRepackMover.fwd.hh
/// @brief Forward declations for DesignRepackMover
/// @author Sarel Fleishman


#ifndef INCLUDED_protocols_simple_moves_DesignRepackMover_fwd_hh
#define INCLUDED_protocols_simple_moves_DesignRepackMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_moves {

class DesignRepackMover;

typedef utility::pointer::shared_ptr< DesignRepackMover > DesignRepackMoverOP;
typedef utility::pointer::shared_ptr< DesignRepackMover const > DesignRepackMoverCOP;

}
}

#endif /*INCLUDED_protocols_protein_interface_design_movers_DesignRepackMover_FWD_HH*/
