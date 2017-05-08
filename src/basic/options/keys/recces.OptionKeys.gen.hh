// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/recces.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_recces_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_recces_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace recces { extern StringOptionKey const seq1; }
namespace recces { extern StringOptionKey const seq2; }
namespace recces { extern IntegerOptionKey const n_cycle; }
namespace recces { extern RealOptionKey const a_form_range; }
namespace recces { extern BooleanOptionKey const dump_pdb; }
namespace recces { extern BooleanOptionKey const dump_silent; }
namespace recces { extern BooleanOptionKey const out_torsions; }
namespace recces { extern RealVectorOptionKey const temps; }
namespace recces { extern RealVectorOptionKey const st_weights; }
namespace recces { extern BooleanOptionKey const save_score_terms; }
namespace recces { extern StringOptionKey const out_prefix; }
namespace recces { extern IntegerOptionKey const dump_freq; }
namespace recces { extern IntegerOptionKey const n_intermediate_dump; }
namespace recces { extern BooleanOptionKey const output_min_pose; }
namespace recces { extern BooleanOptionKey const accept_no_op_moves; }
namespace recces { extern RealOptionKey const histogram_min; }
namespace recces { extern RealOptionKey const histogram_max; }
namespace recces { extern RealOptionKey const histogram_spacing; }
namespace recces { namespace base_pair { extern BooleanOptionKey const base_pair; } }
namespace recces { namespace base_pair { extern RealOptionKey const rmsd_cutoff; } }
namespace recces { namespace base_pair { extern RealOptionKey const translation_mag; } }
namespace recces { namespace base_pair { extern RealOptionKey const rotation_mag; } }
namespace recces { namespace base_pair { extern BooleanOptionKey const recces; } }
namespace recces { namespace base_pair { extern BooleanOptionKey const block_stack; } }
namespace recces { namespace base_pair { extern BooleanOptionKey const sample_jump; } }
namespace recces { namespace thermal_sampling { extern BooleanOptionKey const thermal_sampling; } }
namespace recces { namespace thermal_sampling { extern IntegerVectorOptionKey const sample_residues; } }
namespace recces { namespace thermal_sampling { extern IntegerVectorOptionKey const free_residues; } }
namespace recces { namespace thermal_sampling { extern IntegerVectorOptionKey const loop_residues; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const angle_range_bb; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const angle_range_free_bb; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const angle_range_chi; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const angle_range_free_chi; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const chi_stdev; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const bb_stdev; } }
namespace recces { namespace thermal_sampling { extern RealOptionKey const standard_bb_stdev; } }
namespace recces { namespace thermal_sampling { extern BooleanOptionKey const setup_base_pair_constraints; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
