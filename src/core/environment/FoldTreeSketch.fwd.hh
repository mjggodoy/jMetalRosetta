// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FoldTreeSketch.fwd.hh
/// @brief definition of the FoldTreeSketch class
/// @author

#ifndef INCLUDED_core_environment_FoldTreeSketch_fwd_hh
#define INCLUDED_core_environment_FoldTreeSketch_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace core {
namespace environment {

class FoldTreeSketch;
typedef utility::pointer::shared_ptr< FoldTreeSketch > FoldTreeSketchOP;
typedef utility::pointer::shared_ptr< FoldTreeSketch const > FoldTreeSketchCOP;

} // core
} // environment

#endif //INCLUDED_core_moves_FoldTreeSketch_fwd_HH
