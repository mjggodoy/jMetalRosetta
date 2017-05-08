// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/docking/membrane/MPFindInterfaceMoverCreator.hh
/// @brief      Sample protein-protein interface in the membrane
/// @details The foldtree after the mover is reset to the original one as it
///    was before - so it should work with both fixed and movable membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_fwd_hh
#define INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace docking {
namespace membrane {

class MPFindInterfaceMover;
typedef utility::pointer::shared_ptr< MPFindInterfaceMover > MPFindInterfaceMoverOP;
typedef utility::pointer::shared_ptr< MPFindInterfaceMover const > MPFindInterfaceMoverCOP;

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_fwd_hh
