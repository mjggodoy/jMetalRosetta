// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/score.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_score_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_score_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace score { extern BooleanOptionKey const score_pose_cutpoint_variants; }
namespace score { extern BooleanOptionKey const score; }
namespace score { extern StringOptionKey const weights; }
namespace score { extern StringVectorOptionKey const set_weights; }
namespace score { extern StringOptionKey const pack_weights; }
namespace score { extern StringOptionKey const soft_wts; }
namespace score { extern BooleanOptionKey const docking_interface_score; }
namespace score { extern RealOptionKey const min_score_score; }
namespace score { extern StringOptionKey const custom_atom_pair; }
namespace score { extern FileVectorOptionKey const patch; }
namespace score { extern BooleanOptionKey const empty; }
namespace score { extern RealOptionKey const fa_max_dis; }
namespace score { extern BooleanOptionKey const fa_Hatr; }
namespace score { extern BooleanOptionKey const no_smooth_etables; }
namespace score { extern BooleanOptionKey const no_lk_polar_desolvation; }
namespace score { extern BooleanOptionKey const lk_polar_without_proline_N; }
namespace score { extern StringOptionKey const input_etables; }
namespace score { extern StringOptionKey const output_etables; }
namespace score { extern BooleanOptionKey const analytic_etable_evaluation; }
namespace score { extern BooleanOptionKey const put_intra_into_total; }
namespace score { extern BooleanOptionKey const include_intra_res_protein; }
namespace score { extern BooleanOptionKey const fa_stack_base_base_only; }
namespace score { extern RealOptionKey const fa_stack_sol_prefactor; }
namespace score { extern RealOptionKey const fa_stack_sol_stack_cutoff; }
namespace score { extern RealOptionKey const fa_stack_sol_dist_cutoff; }
namespace score { extern RealOptionKey const fa_stack_lr_prefactor; }
namespace score { extern RealOptionKey const fa_stack_lr_stack_cutoff; }
namespace score { extern RealOptionKey const fa_stack_lr_dist_cutoff; }
namespace score { extern IntegerOptionKey const geom_sol_interres_path_distance_cutoff; }
namespace score { extern IntegerOptionKey const geom_sol_intrares_path_distance_cutoff; }
namespace score { extern RealOptionKey const rms_target; }
namespace score { extern BooleanOptionKey const ramaneighbors; }
namespace score { extern StringOptionKey const optH_weights; }
namespace score { extern StringOptionKey const optH_patch; }
namespace score { extern StringVectorOptionKey const hb_don_strength; }
namespace score { extern StringVectorOptionKey const hb_acc_strength; }
namespace score { extern StringOptionKey const hbe_for_dH2O_aGEN_SP3SC_ssother; }
namespace score { extern StringOptionKey const hbond_params; }
namespace score { extern BooleanOptionKey const hbond_bb_per_residue_energy; }
namespace score { extern BooleanOptionKey const hbond_disable_bbsc_exclusion_rule; }
namespace score { extern IntegerOptionKey const symE_units; }
namespace score { extern RealOptionKey const symE_bonus; }
namespace score { extern BooleanOptionKey const symmetric_gly_tables; }
namespace score { extern RealOptionKey const NV_lbound; }
namespace score { extern RealOptionKey const NV_ubound; }
namespace score { extern StringOptionKey const NV_table; }
namespace score { extern BooleanOptionKey const disable_orientation_dependent_rna_ch_o_bonds; }
namespace score { extern StringOptionKey const rna_torsion_potential; }
namespace score { extern StringOptionKey const rna_suite_potential; }
namespace score { extern StringOptionKey const suiteness_bonus; }
namespace score { extern BooleanOptionKey const rna_torsion_skip_chainbreak; }
namespace score { extern BooleanOptionKey const rna_chemical_shift_verbose; }
namespace score { extern BooleanOptionKey const rna_chemical_shift_larmord; }
namespace score { extern StringOptionKey const rna_chemical_shift_exp_data; }
namespace score { extern StringOptionKey const rna_chemical_shift_larmord_par; }
namespace score { extern StringOptionKey const rna_chemical_shift_larmord_wt; }
namespace score { extern StringOptionKey const rna_chemical_shift_larmord_ref; }
namespace score { extern StringOptionKey const rna_chemical_shift_larmord_nei; }
namespace score { extern RealOptionKey const rna_chemical_shift_larmord_cut; }
namespace score { extern RealOptionKey const rna_chemical_shift_larmord_beta; }
namespace score { extern StringOptionKey const rna_chemical_shift_H5_prime_mode; }
namespace score { extern IntegerVectorOptionKey const rna_chemical_shift_include_res; }
namespace score { extern BooleanOptionKey const DMS_separate_features; }
namespace score { extern BooleanOptionKey const DMS_careful_base_pair_classifier; }
namespace score { extern RealOptionKey const rna_chem_map_lores_weight; }
namespace score { extern BooleanOptionKey const use_2prime_OH_potential; }
namespace score { extern BooleanOptionKey const include_neighbor_base_stacks; }
namespace score { extern BooleanOptionKey const FA_low_res_rnp_scoring; }
namespace score { extern BooleanOptionKey const find_neighbors_3dgrid; }
namespace score { extern BooleanOptionKey const find_neighbors_stripehash; }
namespace score { extern StringOptionKey const seqdep_refene_fname; }
namespace score { extern StringOptionKey const secondary_seqdep_refene_fname; }
namespace score { extern BooleanOptionKey const exact_occ_pairwise; }
namespace score { extern BooleanOptionKey const exact_occ_skip_Hbonders; }
namespace score { extern BooleanOptionKey const exact_occ_pairwise_by_res; }
namespace score { extern BooleanOptionKey const exact_occ_split_between_res; }
namespace score { extern BooleanOptionKey const exact_occ_self_res_no_occ; }
namespace score { extern RealOptionKey const exact_occ_radius_scaling; }
namespace score { extern StringVectorOptionKey const ref_offsets; }
namespace score { extern RealOptionKey const ref_offset; }
namespace score { extern BooleanOptionKey const output_residue_energies; }
namespace score { extern StringOptionKey const fa_custom_pair_distance_file; }
namespace score { extern RealOptionKey const disulf_matching_probe; }
namespace score { extern RealVectorOptionKey const bonded_params; }
namespace score { extern StringOptionKey const bonded_params_dir; }
namespace score { extern StringOptionKey const extra_improper_file; }
namespace score { extern RealOptionKey const pro_close_planar_constraint; }
namespace score { extern BooleanOptionKey const no_pro_close_ring_closure; }
namespace score { extern RealOptionKey const ring_close_shadow_constraint; }
namespace score { extern BooleanOptionKey const linear_bonded_potential; }
namespace score { extern RealOptionKey const free_suite_bonus; }
namespace score { extern RealOptionKey const free_sugar_bonus; }
namespace score { extern RealOptionKey const free_2HOprime_bonus; }
namespace score { extern RealOptionKey const syn_G_potential_bonus; }
namespace score { extern RealOptionKey const pack_phosphate_penalty; }
namespace score { extern RealOptionKey const free_side_chain_bonus; }
namespace score { extern RealOptionKey const bond_angle_sd_polar_hydrogen; }
namespace score { extern RealOptionKey const bond_torsion_sd_polar_hydrogen; }
namespace score { extern BooleanOptionKey const rna_bulge_bonus_once_per_loop; }
namespace score { extern BooleanOptionKey const compute_mg_sol_for_hydrogens; }
namespace score { extern IntegerVectorOptionKey const rg_local_span; }
namespace score { extern BooleanOptionKey const unmodifypot; }
namespace score { extern RealOptionKey const conc; }
namespace score { extern IntegerVectorOptionKey const sidechain_buried; }
namespace score { extern IntegerVectorOptionKey const sidechain_exposed; }
namespace score { extern StringVectorOptionKey const aa_composition_setup_file; }
namespace score { extern StringOptionKey const aa_repeat_energy_penalty_file; }
namespace score { extern RealOptionKey const aspartimide_penalty_value; }
namespace score { extern RealOptionKey const elec_min_dis; }
namespace score { extern RealOptionKey const elec_max_dis; }
namespace score { extern RealOptionKey const elec_die; }
namespace score { extern BooleanOptionKey const elec_r_option; }
namespace score { extern BooleanOptionKey const elec_sigmoidal_die; }
namespace score { extern RealOptionKey const elec_sigmoidal_die_D; }
namespace score { extern RealOptionKey const elec_sigmoidal_die_D0; }
namespace score { extern RealOptionKey const elec_sigmoidal_die_S; }
namespace score { extern BooleanOptionKey const elec_representative_cp; }
namespace score { extern BooleanOptionKey const elec_representative_cp_flip; }
namespace score { extern BooleanOptionKey const eval_intrares_elec_ST_only; }
namespace score { extern BooleanOptionKey const smooth_fa_elec; }
namespace score { extern StringOptionKey const grpelec_fade_type; }
namespace score { extern RealOptionKey const grpelec_fade_param1; }
namespace score { extern RealOptionKey const grpelec_fade_param2; }
namespace score { extern StringOptionKey const elec_group_file; }
namespace score { extern StringOptionKey const elec_group_extrafile; }
namespace score { extern BooleanOptionKey const grpelec_fade_hbond; }
namespace score { extern RealVectorOptionKey const grpelec_max_qeps; }
namespace score { extern BooleanOptionKey const grpelec_context_dependent; }
namespace score { extern BooleanOptionKey const grp_cpfxn; }
namespace score { extern RealVectorOptionKey const grpelec_cpfxn_weight; }
namespace score { extern RealOptionKey const elec_context_minstrength; }
namespace score { extern RealOptionKey const elec_context_minburial; }
namespace score { extern RealOptionKey const elec_context_maxburial; }
namespace score { extern BooleanOptionKey const use_polarization; }
namespace score { extern BooleanOptionKey const use_gen_kirkwood; }
namespace score { extern RealOptionKey const protein_dielectric; }
namespace score { extern RealOptionKey const water_dielectric; }
namespace score { extern RealOptionKey const facts_GBpair_cut; }
namespace score { extern RealOptionKey const facts_kappa; }
namespace score { extern IntegerOptionKey const facts_asp_patch; }
namespace score { extern BooleanOptionKey const facts_plane_to_self; }
namespace score { extern RealOptionKey const facts_saltbridge_correction; }
namespace score { extern RealVectorOptionKey const facts_dshift; }
namespace score { extern RealOptionKey const facts_die; }
namespace score { extern BooleanOptionKey const facts_binding_affinity; }
namespace score { extern BooleanOptionKey const facts_intrascale_by_level; }
namespace score { extern RealVectorOptionKey const facts_intbb_elec_scale; }
namespace score { extern RealVectorOptionKey const facts_intbb_solv_scale; }
namespace score { extern RealVectorOptionKey const facts_adjbb_elec_scale; }
namespace score { extern RealVectorOptionKey const facts_adjbb_solv_scale; }
namespace score { extern RealVectorOptionKey const facts_intbs_elec_scale; }
namespace score { extern RealVectorOptionKey const facts_intbs_solv_scale; }
namespace score { extern RealVectorOptionKey const facts_adjbs_elec_scale; }
namespace score { extern RealVectorOptionKey const facts_adjbs_solv_scale; }
namespace score { extern RealVectorOptionKey const facts_intsc_elec_scale; }
namespace score { extern RealVectorOptionKey const facts_intsc_solv_scale; }
namespace score { extern StringOptionKey const facts_charge_dir; }
namespace score { extern StringOptionKey const facts_eff_charge_dir; }
namespace score { extern StringVectorOptionKey const facts_plane_aa; }
namespace score { extern StringOptionKey const facts_eq_type; }
namespace score { extern IntegerOptionKey const ignore_terminal_ss_elements; }
namespace score { extern BooleanOptionKey const length_dep_srbb; }
namespace score { extern RealOptionKey const ldsrbb_low_scale; }
namespace score { extern RealOptionKey const ldsrbb_high_scale; }
namespace score { extern IntegerOptionKey const ldsrbb_minlength; }
namespace score { extern IntegerOptionKey const ldsrbb_maxlength; }
namespace score { extern RealOptionKey const max_motif_per_res; }
namespace score { extern IntegerOptionKey const max_contacting_ss; }
namespace score { extern BooleanOptionKey const motif_ignore_symmmetry; }
namespace score { extern IntegerVectorOptionKey const motif_residues; }
namespace score { extern StringOptionKey const nmer_ref_energies; }
namespace score { extern StringOptionKey const nmer_ref_energies_list; }
namespace score { extern StringOptionKey const nmer_pssm; }
namespace score { extern StringOptionKey const nmer_pssm_list; }
namespace score { extern RealOptionKey const nmer_pssm_scorecut; }
namespace score { extern StringOptionKey const nmer_svm; }
namespace score { extern StringOptionKey const nmer_svm_list; }
namespace score { extern RealOptionKey const nmer_svm_scorecut; }
namespace score { extern StringOptionKey const nmer_svm_aa_matrix; }
namespace score { extern IntegerOptionKey const nmer_svm_term_length; }
namespace score { extern BooleanOptionKey const nmer_svm_pssm_feat; }
namespace score { extern IntegerOptionKey const nmer_ref_seq_length; }
namespace score { extern BooleanOptionKey const just_calc_rmsd; }
namespace score { extern BooleanOptionKey const envsmooth_zero_negatives; }
namespace score { extern RealOptionKey const rama_power; }
namespace score { extern RealOptionKey const hbond_fade; }
namespace score { extern BooleanOptionKey const hbond_new_sp3_acc; }
namespace score { namespace saxs { extern BooleanOptionKey const saxs; } }
namespace score { namespace saxs { extern RealOptionKey const min_score; } }
namespace score { namespace saxs { extern StringOptionKey const custom_ff; } }
namespace score { namespace saxs { extern StringOptionKey const print_i_calc; } }
namespace score { namespace saxs { extern FileOptionKey const ref_fa_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_cen_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_pddf; } }
namespace score { namespace saxs { extern RealOptionKey const d_min; } }
namespace score { namespace saxs { extern RealOptionKey const d_max; } }
namespace score { namespace saxs { extern RealOptionKey const d_step; } }
namespace score { namespace saxs { extern RealOptionKey const q_min; } }
namespace score { namespace saxs { extern RealOptionKey const q_max; } }
namespace score { namespace saxs { extern RealOptionKey const q_step; } }
namespace score { namespace saxs { extern BooleanOptionKey const fit_pddf_area; } }
namespace score { namespace fiber_diffraction { extern BooleanOptionKey const fiber_diffraction; } }
namespace score { namespace fiber_diffraction { extern FileOptionKey const layer_lines; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const a; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const b; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const p; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const radius; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const resolution_cutoff_low; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const resolution_cutoff_high; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const max_bessel_order; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const cn_symmetry; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const b_factor; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const b_factor_solv; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const b_factor_solv_K; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const grid_reso; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const grid_r; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const grid_phi; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const grid_z; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const qfht_K1; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const qfht_K2; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const chi_iterations; } }
namespace score { namespace fiber_diffraction { extern BooleanOptionKey const rfactor_refinement; } }
namespace score { namespace fiber_diffraction { extern BooleanOptionKey const output_fiber_spectra; } }
namespace score { namespace fiber_diffraction { extern IntegerOptionKey const gpu_processor; } }
namespace score { namespace fiber_diffraction { extern RealOptionKey const centroid_density_mass; } }
namespace score { namespace occ_sol_fitted { extern BooleanOptionKey const occ_sol_fitted; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Ntrp_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const NH2O_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Nlys_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Narg_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Nbb_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Nhis_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OH_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const ONH2_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OOC_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Oaro_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Oet2_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Oet3_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OCbb_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const HOH_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OPha_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OHha_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OC3_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const OSi_amp_scaling; } }
namespace score { namespace occ_sol_fitted { extern RealOptionKey const Oice_amp_scaling; } }
namespace score { namespace loop_close { extern BooleanOptionKey const loop_close; } }
namespace score { namespace loop_close { extern RealOptionKey const loop_fixed_cost; } }
namespace score { namespace loop_close { extern BooleanOptionKey const allow_complex_loop_graph; } }
namespace score { namespace loop_close { extern BooleanOptionKey const use_6D_potential; } }
namespace score { namespace loop_close { extern StringOptionKey const force_6D_potential_file; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
