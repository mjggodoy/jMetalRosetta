// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/MinimizerOptions.fwd.hh
/// @brief  core::optimization::MinimizerOptions forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_optimization_MinimizerOptions_fwd_hh
#define INCLUDED_core_optimization_MinimizerOptions_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace optimization {


// Forward
class MinimizerOptions;

typedef utility::pointer::shared_ptr< MinimizerOptions > MinimizerOptionsOP;
typedef utility::pointer::shared_ptr< MinimizerOptions const > MinimizerOptionsCOP;

} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_MinimizerOptions_FWD_HH
