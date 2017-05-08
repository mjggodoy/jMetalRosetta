// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AbscriptLoopCloserCM.fwd.hh
/// @brief definition of the AbscriptLoopCloserCM class
/// @author

#ifndef INCLUDED_protocols_abinitio_abscript_AbscriptLoopCloserCM_fwd_hh
#define INCLUDED_protocols_abinitio_abscript_AbscriptLoopCloserCM_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace abinitio {
namespace abscript {

class AbscriptLoopCloserCM;
typedef utility::pointer::shared_ptr< AbscriptLoopCloserCM > AbscriptLoopCloserCMOP;
typedef utility::pointer::shared_ptr< AbscriptLoopCloserCM const > AbscriptLoopCloserCMCOP;

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_moves_AbscriptLoopCloserCM_fwd_HH
