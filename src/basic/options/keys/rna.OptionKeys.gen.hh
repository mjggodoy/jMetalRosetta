// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/rna.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_rna_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_rna_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace rna { extern BooleanOptionKey const rna; }
namespace rna { extern BooleanOptionKey const corrected_geo; }
namespace rna { extern BooleanOptionKey const rna_prot_erraser; }
namespace rna { extern BooleanOptionKey const vary_geometry; }
namespace rna { extern StringOptionKey const data_file; }
namespace rna { extern BooleanOptionKey const cut_at_rna_chainbreak; }
namespace rna { extern BooleanOptionKey const evaluate_base_pairs; }
namespace rna { namespace farna { extern BooleanOptionKey const farna; } }
namespace rna { namespace farna { extern IntegerOptionKey const cycles; } }
namespace rna { namespace farna { extern IntegerOptionKey const rna_protein_docking_freq; } }
namespace rna { namespace farna { extern IntegerOptionKey const rounds; } }
namespace rna { namespace farna { extern RealOptionKey const temperature; } }
namespace rna { namespace farna { extern BooleanOptionKey const minimize_rna; } }
namespace rna { namespace farna { extern StringVectorOptionKey const sequence; } }
namespace rna { namespace farna { extern StringOptionKey const secstruct; } }
namespace rna { namespace farna { extern StringOptionKey const secstruct_general; } }
namespace rna { namespace farna { extern StringOptionKey const secstruct_file; } }
namespace rna { namespace farna { extern StringOptionKey const secstruct_general_file; } }
namespace rna { namespace farna { extern StringOptionKey const secstruct_legacy; } }
namespace rna { namespace farna { extern StringOptionKey const lores_scorefxn; } }
namespace rna { namespace farna { extern StringOptionKey const params_file; } }
namespace rna { namespace farna { extern BooleanOptionKey const filter_lores_base_pairs; } }
namespace rna { namespace farna { extern BooleanOptionKey const filter_lores_base_pairs_early; } }
namespace rna { namespace farna { extern BooleanOptionKey const filter_chain_closure; } }
namespace rna { namespace farna { extern BooleanOptionKey const filter_chain_closure_halfway; } }
namespace rna { namespace farna { extern RealOptionKey const filter_chain_closure_distance; } }
namespace rna { namespace farna { extern BooleanOptionKey const relax_rna; } }
namespace rna { namespace farna { extern BooleanOptionKey const simple_relax; } }
namespace rna { namespace farna { extern BooleanOptionKey const ignore_secstruct; } }
namespace rna { namespace farna { extern RealOptionKey const jump_change_frequency; } }
namespace rna { namespace farna { extern BooleanOptionKey const close_loops; } }
namespace rna { namespace farna { extern BooleanOptionKey const close_loops_after_each_move; } }
namespace rna { namespace farna { extern BooleanOptionKey const heat; } }
namespace rna { namespace farna { extern BooleanOptionKey const staged_constraints; } }
namespace rna { namespace farna { extern StringOptionKey const jump_library_file; } }
namespace rna { namespace farna { extern StringOptionKey const vall_torsions; } }
namespace rna { namespace farna { extern BooleanOptionKey const use_1jj2_torsions; } }
namespace rna { namespace farna { extern RealOptionKey const rna_lores_chainbreak_weight; } }
namespace rna { namespace farna { extern RealOptionKey const rna_lores_linear_chainbreak_weight; } }
namespace rna { namespace farna { extern BooleanOptionKey const fixed_stems; } }
namespace rna { namespace farna { extern BooleanOptionKey const allow_bulge; } }
namespace rna { namespace farna { extern IntegerVectorOptionKey const allowed_bulge_res; } }
namespace rna { namespace farna { extern BooleanOptionKey const allow_consecutive_bulges; } }
namespace rna { namespace farna { extern BooleanOptionKey const move_first_rigid_body; } }
namespace rna { namespace farna { extern BooleanOptionKey const root_at_first_rigid_body; } }
namespace rna { namespace farna { extern RealOptionKey const suppress_bp_constraint; } }
namespace rna { namespace farna { extern BooleanOptionKey const output_filters; } }
namespace rna { namespace farna { extern BooleanOptionKey const autofilter; } }
namespace rna { namespace farna { extern BooleanOptionKey const no_filters; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const output_res_num; } }
namespace rna { namespace farna { extern IntegerOptionKey const offset; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const input_silent_res; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const virtual_anchor; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const obligate_pair; } }
namespace rna { namespace farna { extern StringVectorOptionKey const obligate_pair_explicit; } }
namespace rna { namespace farna { extern StringVectorOptionKey const chain_connection; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const remove_pair; } }
namespace rna { namespace farna { extern ResidueChainVectorOptionKey const remove_obligate_pair; } }
namespace rna { namespace farna { extern StringOptionKey const refine_silent_file; } }
namespace rna { namespace farna { extern BooleanOptionKey const refine_native; } }
namespace rna { namespace farna { extern BooleanOptionKey const bps_moves; } }
namespace rna { namespace farna { extern BooleanOptionKey const disallow_bps_at_extra_min_res; } }
namespace rna { namespace farna { extern BooleanOptionKey const allow_fragment_moves_in_bps; } }
namespace rna { namespace farna { extern IntegerOptionKey const frag_size; } }
namespace rna { namespace farna { extern BooleanOptionKey const VDW_rep_screen_include_sidechains; } }
namespace rna { namespace farna { extern BooleanOptionKey const gradual_constraints; } }
namespace rna { namespace farna { extern RealOptionKey const grid_vdw_weight; } }
namespace rna { namespace farna { extern StringOptionKey const tag; } }
namespace rna { namespace farna { extern StringOptionKey const working_native; } }
namespace rna { namespace farna { extern BooleanOptionKey const use_legacy_setup; } }
namespace rna { namespace farna { extern BooleanOptionKey const cst_gap; } }
namespace rna { namespace farna { extern BooleanOptionKey const convert_protein_CEN; } }
namespace rna { namespace farna { extern BooleanOptionKey const rna_protein_docking; } }
namespace rna { namespace farna { extern IntegerVectorOptionKey const exclude_fragments; } }
namespace rna { namespace farna { extern StringOptionKey const exclusion_match_type; } }
namespace rna { namespace farna { extern RealOptionKey const fragment_homology_rmsd; } }
namespace rna { namespace farna { extern BooleanOptionKey const exclude_native_fragments; } }
namespace rna { namespace farna { extern StringVectorOptionKey const exclude_fragment_files; } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const out; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const output_lores_silent_file; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const dump; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const binary_output; } } }
namespace rna { namespace farna { namespace out { extern StringOptionKey const output_score_file; } } }
namespace rna { namespace farna { namespace out { extern IntegerOptionKey const output_score_frequency; } } }
namespace rna { namespace farna { namespace out { extern IntegerVectorOptionKey const output_jump_res; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const output_jump_o3p_to_o5p; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const output_rotation_vector; } } }
namespace rna { namespace farna { namespace out { extern RealVectorOptionKey const target_xyz; } } }
namespace rna { namespace farna { namespace out { extern BooleanOptionKey const save_jump_histogram; } } }
namespace rna { namespace farna { namespace out { extern StringOptionKey const output_histogram_file; } } }
namespace rna { namespace farna { namespace out { extern RealOptionKey const jump_histogram_boxsize; } } }
namespace rna { namespace farna { namespace out { extern RealOptionKey const jump_histogram_binwidth; } } }
namespace rna { namespace farna { namespace out { extern RealOptionKey const jump_histogram_binwidth_rotvector; } } }
namespace rna { namespace farna { namespace db { extern BooleanOptionKey const db; } } }
namespace rna { namespace farna { namespace db { extern BooleanOptionKey const jump_database; } } }
namespace rna { namespace farna { namespace db { extern BooleanOptionKey const bps_database; } } }
namespace rna { namespace farna { namespace erraser { extern BooleanOptionKey const erraser; } } }
namespace rna { namespace farna { namespace erraser { extern BooleanOptionKey const constrain_P; } } }
namespace rna { namespace farna { namespace erraser { extern IntegerVectorOptionKey const fixed_res; } } }
namespace rna { namespace farna { namespace erraser { extern BooleanOptionKey const ready_set_only; } } }
namespace rna { namespace farna { namespace erraser { extern BooleanOptionKey const skip_minimize; } } }
namespace rna { namespace farna { namespace erraser { extern BooleanOptionKey const attempt_pyrimidine_flip; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const minimize; } } }
namespace rna { namespace farna { namespace minimize { extern IntegerOptionKey const minimize_rounds; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const skip_coord_constraints; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const skip_o2prime_trials; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const deriv_check; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const minimizer_use_coordinate_constraints; } } }
namespace rna { namespace farna { namespace minimize { extern StringOptionKey const min_type; } } }
namespace rna { namespace farna { namespace minimize { extern BooleanOptionKey const minimize_bps; } } }
namespace rna { namespace farna { namespace minimize { extern ResidueChainVectorOptionKey const extra_minimize_res; } } }
namespace rna { namespace farna { namespace minimize { extern ResidueChainVectorOptionKey const extra_minimize_chi_res; } } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
