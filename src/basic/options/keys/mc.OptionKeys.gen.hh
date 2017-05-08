// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/mc.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mc_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_mc_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mc { extern BooleanOptionKey const mc; }
namespace mc { extern BooleanOptionKey const log_scores_in_MC; }
namespace mc { extern StringOptionKey const hierarchical_pool; }
namespace mc { extern FileOptionKey const read_structures_into_pool; }
namespace mc { extern IntegerOptionKey const convergence_check_frequency; }
namespace mc { extern FileOptionKey const known_structures; }
namespace mc { extern RealOptionKey const max_rmsd_against_known_structures; }
namespace mc { extern IntegerVectorOptionKey const excluded_residues_from_rmsd; }
namespace mc { extern IntegerOptionKey const heat_convergence_check; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
