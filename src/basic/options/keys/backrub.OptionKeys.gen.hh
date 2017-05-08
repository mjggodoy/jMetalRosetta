// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/backrub.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_backrub_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_backrub_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace backrub { extern BooleanOptionKey const backrub; }
namespace backrub { extern IntegerVectorOptionKey const pivot_residues; }
namespace backrub { extern StringVectorOptionKey const pivot_atoms; }
namespace backrub { extern IntegerOptionKey const min_atoms; }
namespace backrub { extern IntegerOptionKey const max_atoms; }
namespace backrub { extern IntegerOptionKey const ntrials; }
namespace backrub { extern RealOptionKey const sc_prob; }
namespace backrub { extern RealOptionKey const sm_prob; }
namespace backrub { extern RealOptionKey const sc_prob_uniform; }
namespace backrub { extern RealOptionKey const sc_prob_withinrot; }
namespace backrub { extern RealOptionKey const mc_kt; }
namespace backrub { extern RealOptionKey const mm_bend_weight; }
namespace backrub { extern BooleanOptionKey const initial_pack; }
namespace backrub { extern FileOptionKey const minimize_movemap; }
namespace backrub { extern BooleanOptionKey const trajectory; }
namespace backrub { extern BooleanOptionKey const trajectory_gz; }
namespace backrub { extern IntegerOptionKey const trajectory_stride; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
