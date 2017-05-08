// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/flexpack.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_flexpack_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_flexpack_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace flexpack { extern BooleanOptionKey const flexpack; }
namespace flexpack { namespace annealer { extern BooleanOptionKey const annealer; } }
namespace flexpack { namespace annealer { extern RealOptionKey const inner_iteration_scale; } }
namespace flexpack { namespace annealer { extern RealOptionKey const outer_iteration_scale; } }
namespace flexpack { namespace annealer { extern RealOptionKey const fixbb_substitutions_scale; } }
namespace flexpack { namespace annealer { extern RealOptionKey const pure_movebb_substitutions_scale; } }
namespace flexpack { namespace annealer { extern RealOptionKey const rotsub_movebb_substitutions_scale; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
