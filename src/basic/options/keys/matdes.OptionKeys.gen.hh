// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/matdes.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_matdes_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_matdes_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace matdes { extern BooleanOptionKey const matdes; }
namespace matdes { extern IntegerOptionKey const num_subs_building_block; }
namespace matdes { extern IntegerOptionKey const num_subs_total; }
namespace matdes { extern StringOptionKey const pdbID; }
namespace matdes { extern StringOptionKey const prefix; }
namespace matdes { extern RealVectorOptionKey const radial_disp; }
namespace matdes { extern RealVectorOptionKey const angle; }
namespace matdes { extern StringOptionKey const tag; }
namespace matdes { namespace dock { extern BooleanOptionKey const dock; } }
namespace matdes { namespace dock { extern RealOptionKey const neg_r; } }
namespace matdes { namespace dock { extern BooleanOptionKey const dump_pdb; } }
namespace matdes { namespace dock { extern BooleanOptionKey const dump_chainA_only; } }
namespace matdes { namespace design { extern BooleanOptionKey const design; } }
namespace matdes { namespace design { extern RealOptionKey const contact_dist; } }
namespace matdes { namespace design { extern RealOptionKey const grid_size_angle; } }
namespace matdes { namespace design { extern RealOptionKey const grid_size_radius; } }
namespace matdes { namespace design { extern IntegerOptionKey const grid_nsamp_angle; } }
namespace matdes { namespace design { extern IntegerOptionKey const grid_nsamp_radius; } }
namespace matdes { namespace design { extern RealOptionKey const fav_nat_bonus; } }
namespace matdes { namespace mutalyze { extern BooleanOptionKey const mutalyze; } }
namespace matdes { namespace mutalyze { extern BooleanOptionKey const calc_rot_boltz; } }
namespace matdes { namespace mutalyze { extern BooleanOptionKey const ala_scan; } }
namespace matdes { namespace mutalyze { extern BooleanOptionKey const revert_scan; } }
namespace matdes { namespace mutalyze { extern BooleanOptionKey const min_rb; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
