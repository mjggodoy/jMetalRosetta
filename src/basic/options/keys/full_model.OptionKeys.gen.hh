// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/full_model.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_full_model_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_full_model_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace full_model { extern BooleanOptionKey const full_model; }
namespace full_model { extern ResidueChainVectorOptionKey const cutpoint_open; }
namespace full_model { extern ResidueChainVectorOptionKey const cutpoint_closed; }
namespace full_model { extern ResidueChainVectorOptionKey const cyclize; }
namespace full_model { extern StringVectorOptionKey const other_poses; }
namespace full_model { extern ResidueChainVectorOptionKey const jump_res; }
namespace full_model { extern ResidueChainVectorOptionKey const extra_min_res; }
namespace full_model { extern ResidueChainVectorOptionKey const extra_min_jump_res; }
namespace full_model { extern ResidueChainVectorOptionKey const root_res; }
namespace full_model { extern ResidueChainVectorOptionKey const virtual_sugar_res; }
namespace full_model { extern ResidueChainVectorOptionKey const virtual_res; }
namespace full_model { extern ResidueChainVectorOptionKey const sample_res; }
namespace full_model { extern ResidueChainVectorOptionKey const calc_rms_res; }
namespace full_model { extern ResidueChainVectorOptionKey const working_res; }
namespace full_model { extern BooleanOptionKey const motif_mode; }
namespace full_model { extern BooleanOptionKey const allow_jump_in_numbering; }
namespace full_model { namespace rna { extern BooleanOptionKey const rna; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const terminal_res; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const block_stack_above_res; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const block_stack_below_res; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const force_syn_chi_res_list; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const force_anti_chi_res_list; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const force_north_sugar_list; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const force_south_sugar_list; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const bulge_res; } }
namespace full_model { namespace rna { extern ResidueChainVectorOptionKey const sample_sugar_res; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
