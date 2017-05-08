// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/membrane/AddMPLigandMover.fwd.hh
///
/// @brief  Add "single" ligand to to membrane pose
/// @details  Accommodate membrane protein ligand in the membrane framework by
///    reorganizing the current foldtree. Resulting foldtree will
///    keep the membrane attached to the COM and ligand to the closest
///    binding pocket residue, provided in the constructor.
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// #RosettaMPMover

#ifndef INCLUDED_protocols_membrane_AddMPLigandMover_fwd_hh
#define INCLUDED_protocols_membrane_AddMPLigandMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

/// @brief Add Membrane protein ligand mover
/// @details Accommodate membrane protein ligand in the membrane foldtree
class AddMPLigandMover;
typedef utility::pointer::shared_ptr< AddMPLigandMover > AddMPLigandMoverOP;
typedef utility::pointer::shared_ptr< AddMPLigandMover const > AddMPLigandMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMPLigandMover_fwd_hh

