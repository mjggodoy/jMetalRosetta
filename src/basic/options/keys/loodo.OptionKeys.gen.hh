// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/loodo.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_loodo_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_loodo_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace loodo { extern BooleanOptionKey const loodo; }
namespace loodo { extern IntegerOptionKey const ins_begin; }
namespace loodo { extern StringOptionKey const cap; }
namespace loodo { extern StringOptionKey const bot; }
namespace loodo { extern IntegerVectorOptionKey const fragAlength; }
namespace loodo { extern IntegerVectorOptionKey const fragBlength; }
namespace loodo { extern IntegerOptionKey const known; }
namespace loodo { extern StringOptionKey const fragAnative; }
namespace loodo { extern StringOptionKey const fragBnative; }
namespace loodo { extern StringOptionKey const gridligpath; }
namespace loodo { extern BooleanOptionKey const debug; }
namespace loodo { extern RealOptionKey const ca_ratio; }
namespace loodo { extern RealOptionKey const distance_tolerance; }
namespace loodo { extern RealOptionKey const euler_tolerance; }
namespace loodo { extern IntegerOptionKey const num_frags; }
namespace loodo { extern StringOptionKey const use_fraglib; }
namespace loodo { extern StringOptionKey const use_fraglibsc; }
namespace loodo { extern BooleanOptionKey const com_in_grid; }
namespace loodo { extern BooleanOptionKey const loud; }
namespace loodo { extern BooleanOptionKey const dump_all_As; }
namespace loodo { extern BooleanOptionKey const dump_all_Bs; }
namespace loodo { extern StringOptionKey const caphit_rt_file; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
