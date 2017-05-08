// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/magnesium.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_magnesium_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_magnesium_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace magnesium { extern BooleanOptionKey const magnesium; }
namespace magnesium { extern BooleanOptionKey const scan; }
namespace magnesium { extern IntegerVectorOptionKey const mg_res; }
namespace magnesium { extern BooleanOptionKey const minimize_during_scoring; }
namespace magnesium { extern ResidueChainVectorOptionKey const ligand_res; }
namespace magnesium { extern IntegerVectorOptionKey const pose_ligand_res; }
namespace magnesium { extern BooleanOptionKey const lores_scan; }
namespace magnesium { extern RealOptionKey const xyz_step; }
namespace magnesium { extern RealOptionKey const score_cut; }
namespace magnesium { extern RealOptionKey const score_cut_PDB; }
namespace magnesium { extern BooleanOptionKey const integration_test; }
namespace magnesium { extern BooleanOptionKey const tether_to_closest_res; }
namespace magnesium { extern BooleanOptionKey const fixup; }
namespace magnesium { extern BooleanOptionKey const pack_water_hydrogens; }
namespace magnesium { extern BooleanOptionKey const hydrate; }
namespace magnesium { extern BooleanOptionKey const monte_carlo; }
namespace magnesium { extern BooleanOptionKey const scored_hydrogen_sampling; }
namespace magnesium { extern BooleanOptionKey const all_hydration_frames; }
namespace magnesium { extern BooleanOptionKey const leave_other_waters; }
namespace magnesium { extern BooleanOptionKey const minimize; }
namespace magnesium { extern RealOptionKey const minimize_mg_coord_constraint_distance; }
namespace magnesium { namespace montecarlo { extern BooleanOptionKey const montecarlo; } }
namespace magnesium { namespace montecarlo { extern RealOptionKey const temperature; } }
namespace magnesium { namespace montecarlo { extern IntegerOptionKey const cycles; } }
namespace magnesium { namespace montecarlo { extern BooleanOptionKey const dump; } }
namespace magnesium { namespace montecarlo { extern RealOptionKey const add_delete_frequency; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
