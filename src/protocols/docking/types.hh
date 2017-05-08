// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/docking/types.hh
/// @brief  Additionl types for Docking protocol
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_docking_types_hh
#define INCLUDED_protocols_docking_types_hh

#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace docking {

typedef utility::vector1_int DockJumps;

} // docking
} // protocols

#endif


