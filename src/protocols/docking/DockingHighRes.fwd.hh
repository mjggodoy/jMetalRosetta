// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/DockinLowRes.fwd.hh
/// @brief forward declarations for low-res docking
/// @author Jeff Gray

#ifndef INCLUDED_protocols_docking_DockingHighRes_fwd_hh
#define INCLUDED_protocols_docking_DockingHighRes_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace docking {

/// @brief this mover does the high resolution refinement stage of the RosettaDock algorithm
class DockingHighRes;
typedef utility::pointer::shared_ptr< DockingHighRes > DockingHighResOP;
typedef utility::pointer::shared_ptr< DockingHighRes const > DockingHighResCOP;

} //docking
} //protocols

#endif
