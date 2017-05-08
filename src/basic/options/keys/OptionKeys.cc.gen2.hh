namespace score { BooleanOptionKey const grpelec_fade_hbond( "score:grpelec_fade_hbond" );  }
namespace score { RealVectorOptionKey const grpelec_max_qeps( "score:grpelec_max_qeps" );  }
namespace score { BooleanOptionKey const grpelec_context_dependent( "score:grpelec_context_dependent" );  }
namespace score { BooleanOptionKey const grp_cpfxn( "score:grp_cpfxn" );  }
namespace score { RealVectorOptionKey const grpelec_cpfxn_weight( "score:grpelec_cpfxn_weight" );  }
namespace score { RealOptionKey const elec_context_minstrength( "score:elec_context_minstrength" );  }
namespace score { RealOptionKey const elec_context_minburial( "score:elec_context_minburial" );  }
namespace score { RealOptionKey const elec_context_maxburial( "score:elec_context_maxburial" );  }
namespace score { BooleanOptionKey const use_polarization( "score:use_polarization" );  }
namespace score { BooleanOptionKey const use_gen_kirkwood( "score:use_gen_kirkwood" );  }
namespace score { RealOptionKey const protein_dielectric( "score:protein_dielectric" );  }
namespace score { RealOptionKey const water_dielectric( "score:water_dielectric" );  }
namespace score { RealOptionKey const facts_GBpair_cut( "score:facts_GBpair_cut" );  }
namespace score { RealOptionKey const facts_kappa( "score:facts_kappa" );  }
namespace score { IntegerOptionKey const facts_asp_patch( "score:facts_asp_patch" );  }
namespace score { BooleanOptionKey const facts_plane_to_self( "score:facts_plane_to_self" );  }
namespace score { RealOptionKey const facts_saltbridge_correction( "score:facts_saltbridge_correction" );  }
namespace score { RealVectorOptionKey const facts_dshift( "score:facts_dshift" );  }
namespace score { RealOptionKey const facts_die( "score:facts_die" );  }
namespace score { BooleanOptionKey const facts_binding_affinity( "score:facts_binding_affinity" );  }
namespace score { BooleanOptionKey const facts_intrascale_by_level( "score:facts_intrascale_by_level" );  }
namespace score { RealVectorOptionKey const facts_intbb_elec_scale( "score:facts_intbb_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intbb_solv_scale( "score:facts_intbb_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbb_elec_scale( "score:facts_adjbb_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbb_solv_scale( "score:facts_adjbb_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_intbs_elec_scale( "score:facts_intbs_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intbs_solv_scale( "score:facts_intbs_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbs_elec_scale( "score:facts_adjbs_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbs_solv_scale( "score:facts_adjbs_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_intsc_elec_scale( "score:facts_intsc_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intsc_solv_scale( "score:facts_intsc_solv_scale" );  }
namespace score { StringOptionKey const facts_charge_dir( "score:facts_charge_dir" );  }
namespace score { StringOptionKey const facts_eff_charge_dir( "score:facts_eff_charge_dir" );  }
namespace score { StringVectorOptionKey const facts_plane_aa( "score:facts_plane_aa" );  }
namespace score { StringOptionKey const facts_eq_type( "score:facts_eq_type" );  }
namespace score { IntegerOptionKey const ignore_terminal_ss_elements( "score:ignore_terminal_ss_elements" );  }
namespace score { BooleanOptionKey const length_dep_srbb( "score:length_dep_srbb" );  }
namespace score { RealOptionKey const ldsrbb_low_scale( "score:ldsrbb_low_scale" );  }
namespace score { RealOptionKey const ldsrbb_high_scale( "score:ldsrbb_high_scale" );  }
namespace score { IntegerOptionKey const ldsrbb_minlength( "score:ldsrbb_minlength" );  }
namespace score { IntegerOptionKey const ldsrbb_maxlength( "score:ldsrbb_maxlength" );  }
namespace score { RealOptionKey const max_motif_per_res( "score:max_motif_per_res" );  }
namespace score { IntegerOptionKey const max_contacting_ss( "score:max_contacting_ss" );  }
namespace score { BooleanOptionKey const motif_ignore_symmmetry( "score:motif_ignore_symmmetry" );  }
namespace score { IntegerVectorOptionKey const motif_residues( "score:motif_residues" );  }
namespace score { StringOptionKey const nmer_ref_energies( "score:nmer_ref_energies" );  }
namespace score { StringOptionKey const nmer_ref_energies_list( "score:nmer_ref_energies_list" );  }
namespace score { StringOptionKey const nmer_pssm( "score:nmer_pssm" );  }
namespace score { StringOptionKey const nmer_pssm_list( "score:nmer_pssm_list" );  }
namespace score { RealOptionKey const nmer_pssm_scorecut( "score:nmer_pssm_scorecut" );  }
namespace score { StringOptionKey const nmer_svm( "score:nmer_svm" );  }
namespace score { StringOptionKey const nmer_svm_list( "score:nmer_svm_list" );  }
namespace score { RealOptionKey const nmer_svm_scorecut( "score:nmer_svm_scorecut" );  }
namespace score { StringOptionKey const nmer_svm_aa_matrix( "score:nmer_svm_aa_matrix" );  }
namespace score { IntegerOptionKey const nmer_svm_term_length( "score:nmer_svm_term_length" );  }
namespace score { BooleanOptionKey const nmer_svm_pssm_feat( "score:nmer_svm_pssm_feat" );  }
namespace score { IntegerOptionKey const nmer_ref_seq_length( "score:nmer_ref_seq_length" );  }
namespace score { BooleanOptionKey const just_calc_rmsd( "score:just_calc_rmsd" );  }
namespace score { BooleanOptionKey const envsmooth_zero_negatives( "score:envsmooth_zero_negatives" );  }
namespace score { RealOptionKey const rama_power( "score:rama_power" );  }
namespace score { RealOptionKey const hbond_fade( "score:hbond_fade" );  }
namespace score { BooleanOptionKey const hbond_new_sp3_acc( "score:hbond_new_sp3_acc" );  }
namespace score { namespace saxs { BooleanOptionKey const saxs( "score:saxs" );  } }
namespace score { namespace saxs { RealOptionKey const min_score( "score:saxs:min_score" );  } }
namespace score { namespace saxs { StringOptionKey const custom_ff( "score:saxs:custom_ff" );  } }
namespace score { namespace saxs { StringOptionKey const print_i_calc( "score:saxs:print_i_calc" );  } }
namespace score { namespace saxs { FileOptionKey const ref_fa_spectrum( "score:saxs:ref_fa_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_cen_spectrum( "score:saxs:ref_cen_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_spectrum( "score:saxs:ref_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_pddf( "score:saxs:ref_pddf" );  } }
namespace score { namespace saxs { RealOptionKey const d_min( "score:saxs:d_min" );  } }
namespace score { namespace saxs { RealOptionKey const d_max( "score:saxs:d_max" );  } }
namespace score { namespace saxs { RealOptionKey const d_step( "score:saxs:d_step" );  } }
namespace score { namespace saxs { RealOptionKey const q_min( "score:saxs:q_min" );  } }
namespace score { namespace saxs { RealOptionKey const q_max( "score:saxs:q_max" );  } }
namespace score { namespace saxs { RealOptionKey const q_step( "score:saxs:q_step" );  } }
namespace score { namespace saxs { BooleanOptionKey const fit_pddf_area( "score:saxs:fit_pddf_area" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const fiber_diffraction( "score:fiber_diffraction" );  } }
namespace score { namespace fiber_diffraction { FileOptionKey const layer_lines( "score:fiber_diffraction:layer_lines" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const a( "score:fiber_diffraction:a" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b( "score:fiber_diffraction:b" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const p( "score:fiber_diffraction:p" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const radius( "score:fiber_diffraction:radius" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const resolution_cutoff_low( "score:fiber_diffraction:resolution_cutoff_low" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const resolution_cutoff_high( "score:fiber_diffraction:resolution_cutoff_high" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const max_bessel_order( "score:fiber_diffraction:max_bessel_order" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const cn_symmetry( "score:fiber_diffraction:cn_symmetry" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor( "score:fiber_diffraction:b_factor" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor_solv( "score:fiber_diffraction:b_factor_solv" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor_solv_K( "score:fiber_diffraction:b_factor_solv_K" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const grid_reso( "score:fiber_diffraction:grid_reso" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_r( "score:fiber_diffraction:grid_r" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_phi( "score:fiber_diffraction:grid_phi" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_z( "score:fiber_diffraction:grid_z" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const qfht_K1( "score:fiber_diffraction:qfht_K1" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const qfht_K2( "score:fiber_diffraction:qfht_K2" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const chi_iterations( "score:fiber_diffraction:chi_iterations" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const rfactor_refinement( "score:fiber_diffraction:rfactor_refinement" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const output_fiber_spectra( "score:fiber_diffraction:output_fiber_spectra" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const gpu_processor( "score:fiber_diffraction:gpu_processor" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const centroid_density_mass( "score:fiber_diffraction:centroid_density_mass" );  } }
namespace score { namespace occ_sol_fitted { BooleanOptionKey const occ_sol_fitted( "score:occ_sol_fitted" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Ntrp_amp_scaling( "score:occ_sol_fitted:Ntrp_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const NH2O_amp_scaling( "score:occ_sol_fitted:NH2O_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Nlys_amp_scaling( "score:occ_sol_fitted:Nlys_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Narg_amp_scaling( "score:occ_sol_fitted:Narg_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Nbb_amp_scaling( "score:occ_sol_fitted:Nbb_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Nhis_amp_scaling( "score:occ_sol_fitted:Nhis_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OH_amp_scaling( "score:occ_sol_fitted:OH_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const ONH2_amp_scaling( "score:occ_sol_fitted:ONH2_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OOC_amp_scaling( "score:occ_sol_fitted:OOC_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Oaro_amp_scaling( "score:occ_sol_fitted:Oaro_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Oet2_amp_scaling( "score:occ_sol_fitted:Oet2_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Oet3_amp_scaling( "score:occ_sol_fitted:Oet3_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OCbb_amp_scaling( "score:occ_sol_fitted:OCbb_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const HOH_amp_scaling( "score:occ_sol_fitted:HOH_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OPha_amp_scaling( "score:occ_sol_fitted:OPha_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OHha_amp_scaling( "score:occ_sol_fitted:OHha_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OC3_amp_scaling( "score:occ_sol_fitted:OC3_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const OSi_amp_scaling( "score:occ_sol_fitted:OSi_amp_scaling" );  } }
namespace score { namespace occ_sol_fitted { RealOptionKey const Oice_amp_scaling( "score:occ_sol_fitted:Oice_amp_scaling" );  } }
namespace score { namespace loop_close { BooleanOptionKey const loop_close( "score:loop_close" );  } }
namespace score { namespace loop_close { RealOptionKey const loop_fixed_cost( "score:loop_close:loop_fixed_cost" );  } }
namespace score { namespace loop_close { BooleanOptionKey const allow_complex_loop_graph( "score:loop_close:allow_complex_loop_graph" );  } }
namespace score { namespace loop_close { BooleanOptionKey const use_6D_potential( "score:loop_close:use_6D_potential" );  } }
namespace score { namespace loop_close { StringOptionKey const force_6D_potential_file( "score:loop_close:force_6D_potential_file" );  } }
namespace packing { BooleanOptionKey const packing( "packing" );  }
namespace packing { BooleanOptionKey const repack_only( "packing:repack_only" );  }
namespace packing { BooleanOptionKey const prevent_repacking( "packing:prevent_repacking" );  }
namespace packing { RealOptionKey const cenrot_cutoff( "packing:cenrot_cutoff" );  }
namespace packing { BooleanOptionKey const ignore_ligand_chi( "packing:ignore_ligand_chi" );  }
namespace packing { BooleanOptionKey const quasisymmetry( "packing:quasisymmetry" );  }
namespace packing { IntegerOptionKey const ndruns( "packing:ndruns" );  }
namespace packing { BooleanOptionKey const soft_rep_design( "packing:soft_rep_design" );  }
namespace packing { RealOptionKey const mainchain_h_rebuild_threshold( "packing:mainchain_h_rebuild_threshold" );  }
namespace packing { BooleanOptionKey const use_electrostatic_repulsion( "packing:use_electrostatic_repulsion" );  }
namespace packing { BooleanOptionKey const dump_rotamer_sets( "packing:dump_rotamer_sets" );  }
namespace packing { RealOptionKey const dunbrack_prob_buried( "packing:dunbrack_prob_buried" );  }
namespace packing { RealOptionKey const dunbrack_prob_nonburied( "packing:dunbrack_prob_nonburied" );  }
namespace packing { BooleanOptionKey const no_optH( "packing:no_optH" );  }
namespace packing { BooleanOptionKey const optH_MCA( "packing:optH_MCA" );  }
namespace packing { BooleanOptionKey const pack_missing_sidechains( "packing:pack_missing_sidechains" );  }
namespace packing { BooleanOptionKey const preserve_c_beta( "packing:preserve_c_beta" );  }
namespace packing { BooleanOptionKey const flip_HNQ( "packing:flip_HNQ" );  }
namespace packing { IntegerVectorOptionKey const fix_his_tautomer( "packing:fix_his_tautomer" );  }
namespace packing { BooleanOptionKey const print_pymol_selection( "packing:print_pymol_selection" );  }
namespace packing { IntegerOptionKey const extrachi_cutoff( "packing:extrachi_cutoff" );  }
namespace packing { FileVectorOptionKey const resfile( "packing:resfile" );  }
namespace packing { RealOptionKey const outeriterations_scaling( "packing:outeriterations_scaling" );  }
namespace packing { RealOptionKey const inneriterations_scaling( "packing:inneriterations_scaling" );  }
namespace packing { StringVectorOptionKey const adducts( "packing:adducts" );  }
namespace packing { BooleanOptionKey const use_input_sc( "packing:use_input_sc" );  }
namespace packing { FileVectorOptionKey const unboundrot( "packing:unboundrot" );  }
namespace packing { RealOptionKey const max_rotbump_energy( "packing:max_rotbump_energy" );  }
namespace packing { BooleanOptionKey const lazy_ig( "packing:lazy_ig" );  }
namespace packing { BooleanOptionKey const double_lazy_ig( "packing:double_lazy_ig" );  }
namespace packing { IntegerOptionKey const linmem_ig( "packing:linmem_ig" );  }
namespace packing { IntegerOptionKey const multi_cool_annealer( "packing:multi_cool_annealer" );  }
namespace packing { RealVectorOptionKey const minpack_temp_schedule( "packing:minpack_temp_schedule" );  }
namespace packing { IntegerOptionKey const minpack_inner_iteration_scale( "packing:minpack_inner_iteration_scale" );  }
namespace packing { BooleanOptionKey const minpack_disable_bumpcheck( "packing:minpack_disable_bumpcheck" );  }
namespace packing { namespace ex1 { BooleanOptionKey const ex1( "packing:ex1" );  } }
namespace packing { namespace ex1 { IntegerOptionKey const level( "packing:ex1:level" );  } }
namespace packing { namespace ex1 { BooleanOptionKey const operate( "packing:ex1:operate" );  } }
namespace packing { namespace ex2 { BooleanOptionKey const ex2( "packing:ex2" );  } }
namespace packing { namespace ex2 { IntegerOptionKey const level( "packing:ex2:level" );  } }
namespace packing { namespace ex2 { BooleanOptionKey const operate( "packing:ex2:operate" );  } }
namespace packing { namespace ex3 { BooleanOptionKey const ex3( "packing:ex3" );  } }
namespace packing { namespace ex3 { IntegerOptionKey const level( "packing:ex3:level" );  } }
namespace packing { namespace ex3 { BooleanOptionKey const operate( "packing:ex3:operate" );  } }
namespace packing { namespace ex4 { BooleanOptionKey const ex4( "packing:ex4" );  } }
namespace packing { namespace ex4 { IntegerOptionKey const level( "packing:ex4:level" );  } }
namespace packing { namespace ex4 { BooleanOptionKey const operate( "packing:ex4:operate" );  } }
namespace packing { namespace ex1aro { BooleanOptionKey const ex1aro( "packing:ex1aro" );  } }
namespace packing { namespace ex1aro { IntegerOptionKey const level( "packing:ex1aro:level" );  } }
namespace packing { namespace ex2aro { BooleanOptionKey const ex2aro( "packing:ex2aro" );  } }
namespace packing { namespace ex2aro { IntegerOptionKey const level( "packing:ex2aro:level" );  } }
namespace packing { namespace ex1aro_exposed { BooleanOptionKey const ex1aro_exposed( "packing:ex1aro_exposed" );  } }
namespace packing { namespace ex1aro_exposed { IntegerOptionKey const level( "packing:ex1aro_exposed:level" );  } }
namespace packing { namespace ex2aro_exposed { BooleanOptionKey const ex2aro_exposed( "packing:ex2aro_exposed" );  } }
namespace packing { namespace ex2aro_exposed { IntegerOptionKey const level( "packing:ex2aro_exposed:level" );  } }
namespace packing { namespace exdna { BooleanOptionKey const exdna( "packing:exdna" );  } }
namespace packing { namespace exdna { IntegerOptionKey const level( "packing:exdna:level" );  } }
namespace archive { BooleanOptionKey const archive( "archive" );  }
namespace archive { BooleanOptionKey const reread_all_structures( "archive:reread_all_structures" );  }
namespace archive { IntegerOptionKey const completion_notify_frequency( "archive:completion_notify_frequency" );  }
namespace carbohydrates { BooleanOptionKey const carbohydrates( "carbohydrates" );  }
namespace carbohydrates { BooleanOptionKey const glycam_pdb_format( "carbohydrates:glycam_pdb_format" );  }
namespace carbohydrates { StringOptionKey const linkage_conformer_data_file( "carbohydrates:linkage_conformer_data_file" );  }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const glycan_relax( "carbohydrates:glycan_relax" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const glycan_relax_test( "carbohydrates:glycan_relax:glycan_relax_test" );  } }
namespace carbohydrates { namespace glycan_relax { IntegerOptionKey const glycan_relax_rounds( "carbohydrates:glycan_relax:glycan_relax_rounds" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const pack_glycans( "carbohydrates:glycan_relax:pack_glycans" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const final_min_glycans( "carbohydrates:glycan_relax:final_min_glycans" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const glycan_relax_movie( "carbohydrates:glycan_relax:glycan_relax_movie" );  } }
namespace carbohydrates { namespace glycan_relax { RealOptionKey const glycan_relax_kt( "carbohydrates:glycan_relax:glycan_relax_kt" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const glycan_relax_refine( "carbohydrates:glycan_relax:glycan_relax_refine" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const cartmin( "carbohydrates:glycan_relax:cartmin" );  } }
namespace carbohydrates { namespace glycan_relax { BooleanOptionKey const tree_based_min_pack( "carbohydrates:glycan_relax:tree_based_min_pack" );  } }
namespace carbohydrates { namespace clash_check { BooleanOptionKey const clash_check( "carbohydrates:clash_check" );  } }
namespace carbohydrates { namespace clash_check { StringVectorOptionKey const glycan_branches( "carbohydrates:clash_check:glycan_branches" );  } }
namespace carbohydrates { namespace clash_check { StringVectorOptionKey const check_chains( "carbohydrates:clash_check:check_chains" );  } }
namespace carbohydrates { namespace clash_check { RealOptionKey const soft_clash( "carbohydrates:clash_check:soft_clash" );  } }
namespace carbohydrates { namespace clash_check { RealOptionKey const cb_clash_distance( "carbohydrates:clash_check:cb_clash_distance" );  } }
namespace carbohydrates { namespace clash_check { BooleanOptionKey const ignore_hydrogens( "carbohydrates:clash_check:ignore_hydrogens" );  } }
namespace carbohydrates { namespace clash_check { BooleanOptionKey const ignore_full_res_output( "carbohydrates:clash_check:ignore_full_res_output" );  } }
namespace carbohydrates { namespace clash_check { BooleanOptionKey const output_per_glycan_data( "carbohydrates:clash_check:output_per_glycan_data" );  } }
namespace rings { BooleanOptionKey const rings( "rings" );  }
namespace rings { BooleanOptionKey const lock_rings( "rings:lock_rings" );  }
namespace rings { BooleanOptionKey const idealize_rings( "rings:idealize_rings" );  }
namespace rings { BooleanOptionKey const sample_high_energy_conformers( "rings:sample_high_energy_conformers" );  }
namespace chemical { BooleanOptionKey const chemical( "chemical" );  }
namespace chemical { StringVectorOptionKey const exclude_patches( "chemical:exclude_patches" );  }
namespace chemical { StringVectorOptionKey const include_patches( "chemical:include_patches" );  }
namespace chemical { StringVectorOptionKey const add_atom_type_set_parameters( "chemical:add_atom_type_set_parameters" );  }
namespace chemical { StringVectorOptionKey const set_atom_properties( "chemical:set_atom_properties" );  }
namespace chemical { StringVectorOptionKey const patch_selectors( "chemical:patch_selectors" );  }
namespace chemical { BooleanOptionKey const override_rsd_type_limit( "chemical:override_rsd_type_limit" );  }
namespace chemical { StringVectorOptionKey const clone_atom_types( "chemical:clone_atom_types" );  }
namespace chemical { StringVectorOptionKey const reassign_atom_types( "chemical:reassign_atom_types" );  }
namespace chemical { StringVectorOptionKey const reassign_icoor( "chemical:reassign_icoor" );  }
namespace chemical { StringVectorOptionKey const set_atomic_charge( "chemical:set_atomic_charge" );  }
namespace chemical { StringVectorOptionKey const set_patch_atomic_charge( "chemical:set_patch_atomic_charge" );  }
namespace chemical { BooleanOptionKey const enlarge_H_lj( "chemical:enlarge_H_lj" );  }
namespace chemical { BooleanOptionKey const no_hbonds_to_ether_oxygens( "chemical:no_hbonds_to_ether_oxygens" );  }
namespace constraints { BooleanOptionKey const constraints( "constraints" );  }
namespace constraints { BooleanOptionKey const exit_on_bad_read( "constraints:exit_on_bad_read" );  }
namespace constraints { StringVectorOptionKey const cst_file( "constraints:cst_file" );  }
namespace constraints { RealOptionKey const cst_weight( "constraints:cst_weight" );  }
namespace constraints { RealOptionKey const max_cst_dist( "constraints:max_cst_dist" );  }
namespace constraints { StringVectorOptionKey const cst_fa_file( "constraints:cst_fa_file" );  }
namespace constraints { RealOptionKey const cst_fa_weight( "constraints:cst_fa_weight" );  }
namespace constraints { BooleanOptionKey const normalize_mixture_func( "constraints:normalize_mixture_func" );  }
namespace constraints { BooleanOptionKey const penalize_mixture_func( "constraints:penalize_mixture_func" );  }
namespace constraints { FileOptionKey const forest_file( "constraints:forest_file" );  }
namespace constraints { BooleanOptionKey const compute_total_dist_cst( "constraints:compute_total_dist_cst" );  }
namespace constraints { IntegerOptionKey const cull_with_native( "constraints:cull_with_native" );  }
namespace constraints { FileOptionKey const dump_cst_set( "constraints:dump_cst_set" );  }
namespace constraints { IntegerVectorOptionKey const evaluate_max_seq_sep( "constraints:evaluate_max_seq_sep" );  }
namespace constraints { BooleanOptionKey const named( "constraints:named" );  }
namespace constraints { BooleanOptionKey const no_cst_in_relax( "constraints:no_cst_in_relax" );  }
namespace constraints { BooleanOptionKey const no_linearize_bounded( "constraints:no_linearize_bounded" );  }
namespace constraints { RealOptionKey const pocket_constraint_weight( "constraints:pocket_constraint_weight" );  }
namespace constraints { BooleanOptionKey const pocket_zero_derivatives( "constraints:pocket_zero_derivatives" );  }
namespace constraints { BooleanOptionKey const viol( "constraints:viol" );  }
namespace constraints { IntegerOptionKey const viol_level( "constraints:viol_level" );  }
