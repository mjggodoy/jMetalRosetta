// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Sampler.fwd.hh
/// @brief Abstract Base Class for Markov chain rotamer sampler.
/// @author Fang-Chieh Chou

#include <utility/pointer/owning_ptr.hh>

#ifndef INCLUDED_protocols_sampler_MC_Sampler_fwd_HH
#define INCLUDED_protocols_sampler_MC_Sampler_fwd_HH

namespace protocols {
namespace recces {
namespace sampler {

class MC_Sampler;
typedef utility::pointer::shared_ptr< MC_Sampler > MC_SamplerOP;
typedef utility::pointer::shared_ptr< MC_Sampler const > MC_SamplerCOP;

} //sampler
} //recces
} //protocols
#endif
