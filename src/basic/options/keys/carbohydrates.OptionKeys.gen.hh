// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/carbohydrates.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_carbohydrates_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_carbohydrates_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace carbohydrates { extern BooleanOptionKey const carbohydrates; }
namespace carbohydrates { extern BooleanOptionKey const glycam_pdb_format; }
namespace carbohydrates { extern StringOptionKey const linkage_conformer_data_file; }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const glycan_relax; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const glycan_relax_test; } }
namespace carbohydrates { namespace glycan_relax { extern IntegerOptionKey const glycan_relax_rounds; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const pack_glycans; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const final_min_glycans; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const glycan_relax_movie; } }
namespace carbohydrates { namespace glycan_relax { extern RealOptionKey const glycan_relax_kt; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const glycan_relax_refine; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const cartmin; } }
namespace carbohydrates { namespace glycan_relax { extern BooleanOptionKey const tree_based_min_pack; } }
namespace carbohydrates { namespace clash_check { extern BooleanOptionKey const clash_check; } }
namespace carbohydrates { namespace clash_check { extern StringVectorOptionKey const glycan_branches; } }
namespace carbohydrates { namespace clash_check { extern StringVectorOptionKey const check_chains; } }
namespace carbohydrates { namespace clash_check { extern RealOptionKey const soft_clash; } }
namespace carbohydrates { namespace clash_check { extern RealOptionKey const cb_clash_distance; } }
namespace carbohydrates { namespace clash_check { extern BooleanOptionKey const ignore_hydrogens; } }
namespace carbohydrates { namespace clash_check { extern BooleanOptionKey const ignore_full_res_output; } }
namespace carbohydrates { namespace clash_check { extern BooleanOptionKey const output_per_glycan_data; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
