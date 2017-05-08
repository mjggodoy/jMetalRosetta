// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/mp.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mp_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_mp_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mp { extern BooleanOptionKey const mp; }
namespace mp { extern RealOptionKey const thickness; }
namespace mp { extern RealOptionKey const steepness; }
namespace mp { extern RealVectorOptionKey const center_start; }
namespace mp { extern RealOptionKey const center_delta; }
namespace mp { extern RealOptionKey const center_search_cycles; }
namespace mp { extern RealVectorOptionKey const normal_start; }
namespace mp { extern RealOptionKey const normal_angle_start; }
namespace mp { extern RealOptionKey const normal_angle_delta; }
namespace mp { extern RealOptionKey const normal_search_cycles; }
namespace mp { extern RealOptionKey const chain_normal_angle_max; }
namespace mp { extern RealOptionKey const pose_normal_angle_max; }
namespace mp { extern BooleanOptionKey const no_interpolate_Mpair; }
namespace mp { extern BooleanOptionKey const Hbond_depth_correction; }
namespace mp { extern BooleanOptionKey const TMprojection; }
namespace mp { extern RealOptionKey const wt_TMprojection; }
namespace mp { extern BooleanOptionKey const non_helix; }
namespace mp { extern RealOptionKey const wt_non_helix; }
namespace mp { extern BooleanOptionKey const termini; }
namespace mp { extern RealOptionKey const wt_termini; }
namespace mp { extern BooleanOptionKey const secstruct; }
namespace mp { extern RealOptionKey const wt_secstruct; }
namespace mp { extern BooleanOptionKey const spanning; }
namespace mp { extern RealOptionKey const wt_spanning; }
namespace mp { namespace viewer { extern BooleanOptionKey const viewer; } }
namespace mp { namespace viewer { extern RealOptionKey const thickness; } }
namespace mp { namespace viewer { extern IntegerOptionKey const num_points; } }
namespace mp { namespace visualize { extern BooleanOptionKey const visualize; } }
namespace mp { namespace visualize { extern BooleanOptionKey const embedding; } }
namespace mp { namespace visualize { extern RealOptionKey const spacing; } }
namespace mp { namespace visualize { extern RealOptionKey const width; } }
namespace mp { namespace visualize { extern RealOptionKey const thickness; } }
namespace mp { namespace visualize { extern RealOptionKey const plane_radius; } }
namespace mp { namespace scoring { extern BooleanOptionKey const scoring; } }
namespace mp { namespace scoring { extern BooleanOptionKey const hbond; } }
namespace mp { namespace setup { extern BooleanOptionKey const setup; } }
namespace mp { namespace setup { extern StringVectorOptionKey const spanfiles; } }
namespace mp { namespace setup { extern BooleanOptionKey const spans_from_structure; } }
namespace mp { namespace setup { extern StringOptionKey const lipsfile; } }
namespace mp { namespace setup { extern RealVectorOptionKey const center; } }
namespace mp { namespace setup { extern RealVectorOptionKey const normal; } }
namespace mp { namespace setup { extern RealOptionKey const membrane_rsd; } }
namespace mp { namespace setup { extern BooleanOptionKey const transform_into_membrane; } }
namespace mp { namespace setup { extern BooleanOptionKey const position_from_topo; } }
namespace mp { namespace setup { extern BooleanOptionKey const optimize1; } }
namespace mp { namespace setup { extern BooleanOptionKey const optimize2; } }
namespace mp { namespace lipid_acc { extern BooleanOptionKey const lipid_acc; } }
namespace mp { namespace lipid_acc { extern RealOptionKey const angle_cutoff; } }
namespace mp { namespace lipid_acc { extern RealOptionKey const slice_width; } }
namespace mp { namespace lipid_acc { extern RealOptionKey const shell_radius; } }
namespace mp { namespace lipid_acc { extern RealOptionKey const dist_cutoff; } }
namespace mp { namespace lipid_acc { extern BooleanOptionKey const tm_alpha; } }
namespace mp { namespace transform { extern BooleanOptionKey const transform; } }
namespace mp { namespace transform { extern BooleanOptionKey const optimize_embedding; } }
namespace mp { namespace dock { extern BooleanOptionKey const dock; } }
namespace mp { namespace dock { extern StringOptionKey const weights_cen; } }
namespace mp { namespace dock { extern StringOptionKey const weights_fa; } }
namespace mp { namespace dock { extern BooleanOptionKey const lowres; } }
namespace mp { namespace dock { extern BooleanOptionKey const allow_flips; } }
namespace mp { namespace dock { extern BooleanOptionKey const flexible_bb; } }
namespace mp { namespace dock { extern BooleanOptionKey const flexible_sc; } }
namespace mp { namespace dock { extern RealOptionKey const slide_threshold; } }
namespace mp { namespace quickrelax { extern BooleanOptionKey const quickrelax; } }
namespace mp { namespace quickrelax { extern RealOptionKey const angle_max; } }
namespace mp { namespace quickrelax { extern StringOptionKey const nmoves; } }
namespace mp { namespace quickrelax { extern BooleanOptionKey const repack_again; } }
namespace mp { namespace mutate_relax { extern BooleanOptionKey const mutate_relax; } }
namespace mp { namespace mutate_relax { extern StringOptionKey const mutation; } }
namespace mp { namespace mutate_relax { extern StringOptionKey const mutant_file; } }
namespace mp { namespace mutate_relax { extern IntegerOptionKey const nmodels; } }
namespace mp { namespace mutate_relax { extern BooleanOptionKey const repack_mutation_only; } }
namespace mp { namespace mutate_relax { extern RealOptionKey const repack_radius; } }
namespace mp { namespace mutate_relax { extern BooleanOptionKey const relax; } }
namespace mp { namespace benchmark { extern BooleanOptionKey const benchmark; } }
namespace mp { namespace benchmark { namespace ideal_helix { extern BooleanOptionKey const ideal_helix; } } }
namespace mp { namespace benchmark { namespace ideal_helix { extern RealOptionKey const helix_start; } } }
namespace mp { namespace benchmark { namespace ideal_helix { extern RealOptionKey const helix_end; } } }
namespace mp { namespace benchmark { namespace tilt_angle { extern BooleanOptionKey const tilt_angle; } } }
namespace mp { namespace benchmark { namespace tilt_angle { extern StringOptionKey const output; } } }
namespace mp { namespace output { extern BooleanOptionKey const output; } }
namespace mp { namespace output { extern BooleanOptionKey const normalize_to_thk; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
