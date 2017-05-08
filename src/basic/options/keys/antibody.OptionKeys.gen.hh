// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/antibody.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_antibody_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_antibody_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace antibody { extern BooleanOptionKey const antibody; }
namespace antibody { extern StringOptionKey const input_ab_scheme; }
namespace antibody { extern StringOptionKey const output_ab_scheme; }
namespace antibody { extern StringOptionKey const cdr_definition; }
namespace antibody { extern StringOptionKey const light_chain; }
namespace antibody { extern BooleanOptionKey const check_cdr_chainbreaks; }
namespace antibody { extern BooleanOptionKey const check_cdr_pep_bond_geom; }
namespace antibody { extern StringOptionKey const numbering_scheme; }
namespace antibody { extern BooleanOptionKey const graft_l1; }
namespace antibody { extern StringOptionKey const l1_template; }
namespace antibody { extern BooleanOptionKey const graft_l2; }
namespace antibody { extern StringOptionKey const l2_template; }
namespace antibody { extern BooleanOptionKey const graft_l3; }
namespace antibody { extern StringOptionKey const l3_template; }
namespace antibody { extern BooleanOptionKey const graft_h1; }
namespace antibody { extern StringOptionKey const h1_template; }
namespace antibody { extern BooleanOptionKey const graft_h2; }
namespace antibody { extern StringOptionKey const h2_template; }
namespace antibody { extern BooleanOptionKey const graft_h3; }
namespace antibody { extern StringOptionKey const h3_template; }
namespace antibody { extern StringOptionKey const light_heavy_template; }
namespace antibody { extern StringOptionKey const frl_template; }
namespace antibody { extern StringOptionKey const frh_template; }
namespace antibody { extern BooleanOptionKey const h3_no_stem_graft; }
namespace antibody { extern BooleanOptionKey const packonly_after_graft; }
namespace antibody { extern BooleanOptionKey const stem_optimize; }
namespace antibody { extern IntegerOptionKey const stem_optimize_size; }
namespace antibody { extern StringOptionKey const preprocessing_script_version; }
namespace antibody { extern BooleanOptionKey const model_h3; }
namespace antibody { extern BooleanOptionKey const snugfit; }
namespace antibody { extern BooleanOptionKey const refine_h3; }
namespace antibody { extern BooleanOptionKey const h3_filter; }
namespace antibody { extern RealOptionKey const h3_filter_tolerance; }
namespace antibody { extern BooleanOptionKey const cter_insert; }
namespace antibody { extern BooleanOptionKey const flank_residue_min; }
namespace antibody { extern BooleanOptionKey const sc_min; }
namespace antibody { extern BooleanOptionKey const rt_min; }
namespace antibody { extern BooleanOptionKey const bad_nter; }
namespace antibody { extern BooleanOptionKey const extend_h3_before_modeling; }
namespace antibody { extern BooleanOptionKey const idealize_h3_stems_before_modeling; }
namespace antibody { extern StringOptionKey const remodel; }
namespace antibody { extern StringOptionKey const refine; }
namespace antibody { extern StringOptionKey const centroid_refine; }
namespace antibody { extern BooleanOptionKey const constrain_cter; }
namespace antibody { extern BooleanOptionKey const constrain_vlvh_qq; }
namespace antibody { extern BooleanOptionKey const auto_generate_kink_constraint; }
namespace antibody { extern BooleanOptionKey const all_atom_mode_kink_constraint; }
namespace antibody { extern BooleanOptionKey const snug_loops; }
namespace antibody { extern FileOptionKey const input_fv; }
namespace antibody { extern BooleanOptionKey const camelid; }
namespace antibody { extern BooleanOptionKey const camelid_constraints; }
namespace antibody { extern BooleanOptionKey const use_mean_cluster_cst_data; }
namespace antibody { extern BooleanOptionKey const force_use_of_cluster_csts_with_outliers; }
namespace antibody { extern IntegerOptionKey const cluster_csts_stats_cutoff; }
namespace antibody { extern RealOptionKey const general_dihedral_cst_phi_sd; }
namespace antibody { extern RealOptionKey const general_dihedral_cst_psi_sd; }
namespace antibody { extern BooleanOptionKey const allow_omega_mismatches_for_north_clusters; }
namespace antibody { extern StringOptionKey const prefix; }
namespace antibody { extern StringOptionKey const grafting_database; }
namespace antibody { extern StringOptionKey const blastp; }
namespace antibody { extern BooleanOptionKey const exclude_homologs; }
namespace antibody { extern RealOptionKey const exclude_homologs_cdr_cutoff; }
namespace antibody { extern RealOptionKey const exclude_homologs_fr_cutoff; }
namespace antibody { extern RealOptionKey const ocd_cutoff; }
namespace antibody { extern IntegerOptionKey const n_multi_templates; }
namespace antibody { namespace design { extern BooleanOptionKey const design; } }
namespace antibody { namespace design { extern StringOptionKey const base_cdr_instructions; } }
namespace antibody { namespace design { extern StringOptionKey const cdr_instructions; } }
namespace antibody { namespace design { extern StringOptionKey const antibody_database; } }
namespace antibody { namespace design { extern BooleanOptionKey const paper_ab_db; } }
namespace antibody { namespace design { extern StringOptionKey const paper_ab_db_path; } }
namespace antibody { namespace design { extern StringVectorOptionKey const design_cdrs; } }
namespace antibody { namespace design { extern IntegerOptionKey const top_designs; } }
namespace antibody { namespace design { extern StringOptionKey const design_protocol; } }
namespace antibody { namespace design { extern BooleanOptionKey const run_snugdock; } }
namespace antibody { namespace design { extern BooleanOptionKey const run_relax; } }
namespace antibody { namespace design { extern BooleanOptionKey const run_interface_analyzer; } }
namespace antibody { namespace design { extern StringVectorOptionKey const paratope; } }
namespace antibody { namespace design { extern StringVectorOptionKey const epitope; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_epitope_constraints; } }
namespace antibody { namespace design { extern RealOptionKey const dihedral_cst_weight; } }
namespace antibody { namespace design { extern RealOptionKey const atom_pair_cst_weight; } }
namespace antibody { namespace design { extern BooleanOptionKey const global_dihedral_cst_scoring; } }
namespace antibody { namespace design { extern BooleanOptionKey const global_atom_pair_cst_scoring; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_dock; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_rb_min; } }
namespace antibody { namespace design { extern BooleanOptionKey const dock_min_dock; } }
namespace antibody { namespace design { extern IntegerOptionKey const outer_cycle_rounds; } }
namespace antibody { namespace design { extern IntegerOptionKey const inner_cycle_rounds; } }
namespace antibody { namespace design { extern IntegerOptionKey const dock_cycle_rounds; } }
namespace antibody { namespace design { extern RealOptionKey const interface_dis; } }
namespace antibody { namespace design { extern RealOptionKey const neighbor_dis; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_outliers; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_H3_graft_outliers; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_only_H3_kinked; } }
namespace antibody { namespace design { extern BooleanOptionKey const design_antigen; } }
namespace antibody { namespace design { extern BooleanOptionKey const design_framework; } }
namespace antibody { namespace design { extern BooleanOptionKey const conservative_framework_design; } }
namespace antibody { namespace design { extern BooleanOptionKey const design_H3_stem; } }
namespace antibody { namespace design { extern BooleanOptionKey const design_proline; } }
namespace antibody { namespace design { extern RealOptionKey const sample_zero_probs_at; } }
namespace antibody { namespace design { extern BooleanOptionKey const force_mutate_framework_for_cluster; } }
namespace antibody { namespace design { extern IntegerOptionKey const seq_design_stats_cutoff; } }
namespace antibody { namespace design { extern IntegerOptionKey const seq_design_profile_samples; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_light_chain_type; } }
namespace antibody { namespace design { extern BooleanOptionKey const idealize_graft_cdrs; } }
namespace antibody { namespace design { extern StringVectorOptionKey const add_backrub_pivots; } }
namespace antibody { namespace design { extern RealOptionKey const inner_kt; } }
namespace antibody { namespace design { extern RealOptionKey const outer_kt; } }
namespace antibody { namespace design { extern BooleanOptionKey const random_start; } }
namespace antibody { namespace design { extern BooleanOptionKey const adapt_graft; } }
namespace antibody { namespace design { extern BooleanOptionKey const enable_adapt_graft_cartesian; } }
namespace antibody { namespace design { extern BooleanOptionKey const remove_antigen; } }
namespace antibody { namespace design { extern BooleanOptionKey const add_graft_log_to_pdb; } }
namespace antibody { namespace design { extern BooleanOptionKey const mutate_framework_for_cluster; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
