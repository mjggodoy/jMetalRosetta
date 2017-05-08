// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/FlatHarmonicFunc.fwd.hh
/// @brief forward declaration for Mixture functions
/// @author Chris King modified from James Thompson


#ifndef INCLUDED_core_scoring_func_FlatHarmonicFunc_fwd_hh
#define INCLUDED_core_scoring_func_FlatHarmonicFunc_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

class FlatHarmonicFunc;
typedef utility::pointer::shared_ptr< FlatHarmonicFunc > FlatHarmonicFuncOP;
typedef utility::pointer::shared_ptr< FlatHarmonicFunc const > FlatHarmonicFuncCOP;

} // constraints
} // scoring
} // core

#endif
