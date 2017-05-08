// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterCreator.cc
/// @brief  Base class for PosePropertyReporters for the load-time factory registration scheme
/// @author Luki Goldschmidt <lugo@uw.edu>

#include <protocols/rosetta_scripts/PosePropertyReporterCreator.hh>

namespace protocols {
namespace rosetta_scripts {

PosePropertyReporterCreator::PosePropertyReporterCreator() {}
PosePropertyReporterCreator::~PosePropertyReporterCreator() = default;

} //namespace rosetta_scripts
} //namespace protocols
