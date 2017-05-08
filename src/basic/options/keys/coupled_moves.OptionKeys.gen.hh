// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/coupled_moves.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_coupled_moves_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_coupled_moves_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace coupled_moves { extern BooleanOptionKey const coupled_moves; }
namespace coupled_moves { extern IntegerOptionKey const ntrials; }
namespace coupled_moves { extern IntegerOptionKey const number_ligands; }
namespace coupled_moves { extern RealOptionKey const mc_kt; }
namespace coupled_moves { extern RealOptionKey const boltzmann_kt; }
namespace coupled_moves { extern RealOptionKey const mm_bend_weight; }
namespace coupled_moves { extern BooleanOptionKey const trajectory; }
namespace coupled_moves { extern BooleanOptionKey const trajectory_gz; }
namespace coupled_moves { extern IntegerOptionKey const trajectory_stride; }
namespace coupled_moves { extern StringOptionKey const trajectory_file; }
namespace coupled_moves { extern StringOptionKey const output_fasta; }
namespace coupled_moves { extern StringOptionKey const output_stats; }
namespace coupled_moves { extern BooleanOptionKey const ligand_mode; }
namespace coupled_moves { extern BooleanOptionKey const initial_repack; }
namespace coupled_moves { extern BooleanOptionKey const min_pack; }
namespace coupled_moves { extern BooleanOptionKey const save_sequences; }
namespace coupled_moves { extern BooleanOptionKey const save_structures; }
namespace coupled_moves { extern RealOptionKey const ligand_prob; }
namespace coupled_moves { extern BooleanOptionKey const fix_backbone; }
namespace coupled_moves { extern BooleanOptionKey const uniform_backrub; }
namespace coupled_moves { extern BooleanOptionKey const bias_sampling; }
namespace coupled_moves { extern BooleanOptionKey const bump_check; }
namespace coupled_moves { extern RealOptionKey const ligand_weight; }
namespace coupled_moves { extern StringOptionKey const output_prefix; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
