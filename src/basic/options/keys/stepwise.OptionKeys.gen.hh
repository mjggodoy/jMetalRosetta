// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/stepwise.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_stepwise_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_stepwise_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace stepwise { extern BooleanOptionKey const stepwise; }
namespace stepwise { extern StringVectorOptionKey const s1; }
namespace stepwise { extern StringVectorOptionKey const s2; }
namespace stepwise { extern StringVectorOptionKey const silent1; }
namespace stepwise { extern StringVectorOptionKey const silent2; }
namespace stepwise { extern StringVectorOptionKey const tags1; }
namespace stepwise { extern StringVectorOptionKey const tags2; }
namespace stepwise { extern IntegerVectorOptionKey const slice_res1; }
namespace stepwise { extern IntegerVectorOptionKey const slice_res2; }
namespace stepwise { extern IntegerVectorOptionKey const input_res1; }
namespace stepwise { extern IntegerVectorOptionKey const input_res2; }
namespace stepwise { extern BooleanOptionKey const backbone_only1; }
namespace stepwise { extern BooleanOptionKey const backbone_only2; }
namespace stepwise { extern IntegerVectorOptionKey const fixed_res; }
namespace stepwise { extern BooleanOptionKey const test_encapsulation; }
namespace stepwise { extern BooleanOptionKey const choose_random; }
namespace stepwise { extern IntegerOptionKey const num_random_samples; }
namespace stepwise { extern IntegerOptionKey const max_tries_multiplier_for_ccd; }
namespace stepwise { extern IntegerOptionKey const num_pose_minimize; }
namespace stepwise { extern BooleanOptionKey const atr_rep_screen; }
namespace stepwise { extern BooleanOptionKey const atr_rep_screen_for_docking; }
namespace stepwise { extern StringOptionKey const align_pdb; }
namespace stepwise { extern BooleanOptionKey const enumerate; }
namespace stepwise { extern BooleanOptionKey const preminimize; }
namespace stepwise { extern BooleanOptionKey const skip_preminimize; }
namespace stepwise { extern BooleanOptionKey const minimize_waters; }
namespace stepwise { extern BooleanOptionKey const test_all_moves; }
namespace stepwise { extern BooleanOptionKey const new_move_selector; }
namespace stepwise { extern BooleanOptionKey const dump; }
namespace stepwise { extern BooleanOptionKey const VERBOSE; }
namespace stepwise { extern BooleanOptionKey const use_green_packer; }
namespace stepwise { extern RealOptionKey const rmsd_screen; }
namespace stepwise { extern BooleanOptionKey const skip_minimize; }
namespace stepwise { extern BooleanOptionKey const virtualize_packable_moieties_in_screening_pose; }
namespace stepwise { extern StringOptionKey const minimizer_mode; }
namespace stepwise { extern StringOptionKey const sampler_silent_file; }
namespace stepwise { extern BooleanOptionKey const superimpose_over_all; }
namespace stepwise { extern StringVectorOptionKey const move; }
namespace stepwise { extern StringOptionKey const min_type; }
namespace stepwise { extern RealOptionKey const min_tolerance; }
namespace stepwise { extern BooleanOptionKey const vary_polar_hydrogen_geometry; }
namespace stepwise { extern BooleanOptionKey const output_minimized_pose_list; }
namespace stepwise { extern BooleanOptionKey const virtualize_free_moieties_in_native; }
namespace stepwise { extern BooleanOptionKey const output_cluster_size; }
namespace stepwise { extern BooleanOptionKey const lores; }
namespace stepwise { extern BooleanOptionKey const verbose_sampler; }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const monte_carlo; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const verbose_scores; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const skip_deletions; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_internal_hinge_moves; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_internal_local_moves; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_skip_bulge; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const skip_bulge_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const from_scratch_frequency; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_split_off; } }
namespace stepwise { namespace monte_carlo { extern IntegerOptionKey const cycles; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const temperature; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const add_proposal_density_factor; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const add_delete_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const docking_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const submotif_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const intermolecular_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const minimize_single_res_frequency; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_variable_bond_geometry; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const switch_focus_frequency; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const just_min_after_mutation_frequency; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const local_redock_only; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const make_movie; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const recover_low; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const use_precomputed_library; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const allow_submotif_split; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const force_submotif_without_intervening_bulge; } }
namespace stepwise { namespace monte_carlo { extern BooleanOptionKey const use_first_jump_for_submotif; } }
namespace stepwise { namespace monte_carlo { extern RealOptionKey const vary_loop_length_frequency; } }
namespace stepwise { namespace monte_carlo { extern IntegerOptionKey const checkpointing_frequency; } }
namespace stepwise { namespace monte_carlo { namespace csa { extern BooleanOptionKey const csa; } } }
namespace stepwise { namespace monte_carlo { namespace csa { extern IntegerOptionKey const csa_bank_size; } } }
namespace stepwise { namespace monte_carlo { namespace csa { extern RealOptionKey const csa_rmsd; } } }
namespace stepwise { namespace monte_carlo { namespace csa { extern BooleanOptionKey const csa_output_rounds; } } }
namespace stepwise { namespace rna { extern BooleanOptionKey const rna; } }
namespace stepwise { namespace rna { extern IntegerOptionKey const sampler_num_pose_kept; } }
namespace stepwise { namespace rna { extern RealOptionKey const native_edensity_score_cutoff; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const o2prime_legacy_mode; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const allow_virtual_o2prime_hydrogens; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_perform_phosphate_pack; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const force_phosphate_instantiation; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const distinguish_pucker; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const finer_sampling_at_chain_closure; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const PBP_clustering_at_chain_closure; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_allow_syn_pyrimidine; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_extra_chi_rotamer; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_extra_beta_rotamer; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_extra_epsilon_rotamer; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const force_centroid_interaction; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const virtual_sugar_legacy_mode; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const VDW_rep_optimize_memory_usage; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const erraser; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const centroid_screen; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const VDW_atr_rep_screen; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const minimize_and_score_native_pose; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const rm_virt_phosphate; } }
namespace stepwise { namespace rna { extern StringVectorOptionKey const VDW_rep_screen_info; } }
namespace stepwise { namespace rna { extern RealOptionKey const VDW_rep_alignment_RMSD_CUTOFF; } }
namespace stepwise { namespace rna { extern StringVectorOptionKey const VDW_rep_delete_matching_res; } }
namespace stepwise { namespace rna { extern RealOptionKey const VDW_rep_screen_physical_pose_clash_dist_cutoff; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const integration_test; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const allow_bulge_at_chainbreak; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const parin_favorite_output; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const reinitialize_CCD_torsions; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sample_both_sugar_base_rotamer; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_include_torsion_value_in_tag; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_assert_no_virt_sugar_sampling; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sampler_try_sugar_instantiation; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const do_not_sample_multiple_virtual_sugar; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const sample_ONLY_multiple_virtual_sugar; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const allow_base_pair_only_centroid_screen; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const minimizer_rename_tag; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const minimize_res; } }
namespace stepwise { namespace rna { extern StringVectorOptionKey const alignment_res; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const native_alignment_res; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const rmsd_res; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const missing_res; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const missing_res2; } }
namespace stepwise { namespace rna { extern IntegerOptionKey const job_queue_ID; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const minimize_and_score_sugar; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const global_sample_res_list; } }
namespace stepwise { namespace rna { extern FileOptionKey const filter_output_filename; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const combine_long_loop_mode; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const combine_helical_silent_file; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const output_extra_RMSDs; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const protonated_H1_adenosine_list; } }
namespace stepwise { namespace rna { extern IntegerVectorOptionKey const native_virtual_res; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const simple_append_map; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const allow_fixed_res_at_moving_res; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const force_user_defined_jumps; } }
namespace stepwise { namespace rna { extern StringVectorOptionKey const jump_point_pairs; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const add_virt_root; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const floating_base; } }
namespace stepwise { namespace rna { extern IntegerOptionKey const floating_base_anchor_res; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const allow_chain_boundary_jump_partner_right_at_fixed_BP; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const rebuild_bulge_mode; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const virtual_sugar_keep_base_fixed; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const virtual_sugar_do_minimize; } }
namespace stepwise { namespace rna { extern RealOptionKey const sampler_max_centroid_distance; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const filter_user_alignment_res; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const tether_jump; } }
namespace stepwise { namespace rna { extern BooleanOptionKey const turn_off_rna_chem_map_during_optimize; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const protein; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const global_optimize; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const disable_sampling_of_loop_takeoff; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const sample_beta; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const ghost_loops; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const centroid_screen; } }
namespace stepwise { namespace protein { extern RealOptionKey const centroid_score_diff_cut; } }
namespace stepwise { namespace protein { extern StringOptionKey const centroid_weights; } }
namespace stepwise { namespace protein { extern RealOptionKey const score_diff_cut; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const filter_native_big_bins; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const cluster_by_all_atom_rmsd; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const centroid_output; } }
namespace stepwise { namespace protein { extern IntegerOptionKey const n_sample; } }
namespace stepwise { namespace protein { extern IntegerOptionKey const nstruct_centroid; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const ccd_close; } }
namespace stepwise { namespace protein { extern IntegerVectorOptionKey const bridge_res; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const cart_min; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const move_jumps_between_chains; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const use_packer_instead_of_rotamer_trials; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const expand_loop_takeoff; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const skip_coord_constraints; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const allow_virtual_side_chains; } }
namespace stepwise { namespace protein { extern BooleanOptionKey const protein_prepack; } }
namespace stepwise { namespace protein { extern StringOptionKey const disulfide_file; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
