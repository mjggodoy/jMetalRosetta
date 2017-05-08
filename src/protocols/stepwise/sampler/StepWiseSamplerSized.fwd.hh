// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSized.fwd.hh
/// @brief Abstract Base Class for Sampler generator with finite size.
/// @author Fang-Chieh Chou

#include <utility/pointer/owning_ptr.hh>

#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerSized_fwd_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerSized_fwd_HH

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerSized;
typedef utility::pointer::shared_ptr< StepWiseSamplerSized > StepWiseSamplerSizedOP;

} //sampler
} //stepwise
} //protocols
#endif
