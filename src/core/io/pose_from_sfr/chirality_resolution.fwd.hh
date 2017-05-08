// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pose_from_sfr/chirality_resolution.fwd.hh
/// @brief  Various utilities to accomodated PDB input issues.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_io_pose_from_sfr_chirality_resolution_FWD_HH
#define INCLUDED_core_io_pose_from_sfr_chirality_resolution_FWD_HH

//External headers
#include <boost/bimap.hpp>

// C++ headers
#include <string>

namespace core {
namespace io {
namespace pose_from_sfr {

// Note: boost::bimap has it's head so far up generalizability that it makes simple use cases
// horribly complex. (Primarily, support for operator[] is completely busted.)
// If you can replace this with a sane and usable bimap implementation, please do so.

typedef boost::bimap< std::string, std::string > NameBimap;

} // namespace pose_from_sfr
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pose_from_sfr_chirality_resolution_FWD_HH
