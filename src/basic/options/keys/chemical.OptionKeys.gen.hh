// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/chemical.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_chemical_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_chemical_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace chemical { extern BooleanOptionKey const chemical; }
namespace chemical { extern StringVectorOptionKey const exclude_patches; }
namespace chemical { extern StringVectorOptionKey const include_patches; }
namespace chemical { extern StringVectorOptionKey const add_atom_type_set_parameters; }
namespace chemical { extern StringVectorOptionKey const set_atom_properties; }
namespace chemical { extern StringVectorOptionKey const patch_selectors; }
namespace chemical { extern BooleanOptionKey const override_rsd_type_limit; }
namespace chemical { extern StringVectorOptionKey const clone_atom_types; }
namespace chemical { extern StringVectorOptionKey const reassign_atom_types; }
namespace chemical { extern StringVectorOptionKey const reassign_icoor; }
namespace chemical { extern StringVectorOptionKey const set_atomic_charge; }
namespace chemical { extern StringVectorOptionKey const set_patch_atomic_charge; }
namespace chemical { extern BooleanOptionKey const enlarge_H_lj; }
namespace chemical { extern BooleanOptionKey const no_hbonds_to_ether_oxygens; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
