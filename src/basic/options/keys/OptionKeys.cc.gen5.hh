namespace abinitio { IntegerOptionKey const number_3mer_frags( "abinitio:number_3mer_frags" );  }
namespace abinitio { IntegerOptionKey const number_9mer_frags( "abinitio:number_9mer_frags" );  }
namespace abinitio { RealOptionKey const temperature( "abinitio:temperature" );  }
namespace abinitio { RealOptionKey const rg_reweight( "abinitio:rg_reweight" );  }
namespace abinitio { RealOptionKey const rsd_wt_helix( "abinitio:rsd_wt_helix" );  }
namespace abinitio { RealOptionKey const rsd_wt_strand( "abinitio:rsd_wt_strand" );  }
namespace abinitio { RealOptionKey const rsd_wt_loop( "abinitio:rsd_wt_loop" );  }
namespace abinitio { BooleanOptionKey const skip_convergence_check( "abinitio:skip_convergence_check" );  }
namespace abinitio { FileVectorOptionKey const stage1_patch( "abinitio:stage1_patch" );  }
namespace abinitio { FileVectorOptionKey const stage2_patch( "abinitio:stage2_patch" );  }
namespace abinitio { FileVectorOptionKey const stage3a_patch( "abinitio:stage3a_patch" );  }
namespace abinitio { FileVectorOptionKey const stage3b_patch( "abinitio:stage3b_patch" );  }
namespace abinitio { FileVectorOptionKey const stage4_patch( "abinitio:stage4_patch" );  }
namespace abinitio { FileVectorOptionKey const stage5_patch( "abinitio:stage5_patch" );  }
namespace abinitio { BooleanOptionKey const steal_3mers( "abinitio:steal_3mers" );  }
namespace abinitio { BooleanOptionKey const steal_9mers( "abinitio:steal_9mers" );  }
namespace abinitio { BooleanOptionKey const no_write_failures( "abinitio:no_write_failures" );  }
namespace abinitio { BooleanOptionKey const relax_failures( "abinitio:relax_failures" );  }
namespace abinitio { BooleanOptionKey const relax_with_jumps( "abinitio:relax_with_jumps" );  }
namespace abinitio { BooleanOptionKey const process_store( "abinitio:process_store" );  }
namespace abinitio { IntegerVectorOptionKey const fix_residues_to_native( "abinitio:fix_residues_to_native" );  }
namespace abinitio { BooleanOptionKey const return_full_atom( "abinitio:return_full_atom" );  }
namespace abinitio { BooleanOptionKey const detect_disulfide_before_relax( "abinitio:detect_disulfide_before_relax" );  }
namespace abinitio { BooleanOptionKey const close_loops( "abinitio:close_loops" );  }
namespace abinitio { BooleanOptionKey const bGDT( "abinitio:bGDT" );  }
namespace abinitio { BooleanOptionKey const dump_frags( "abinitio:dump_frags" );  }
namespace abinitio { BooleanOptionKey const jdist_rerun( "abinitio:jdist_rerun" );  }
namespace abinitio { RealOptionKey const perturb( "abinitio:perturb" );  }
namespace abinitio { BooleanOptionKey const rerun( "abinitio:rerun" );  }
namespace abinitio { IntegerVectorOptionKey const rmsd_residues( "abinitio:rmsd_residues" );  }
namespace abinitio { BooleanOptionKey const start_native( "abinitio:start_native" );  }
namespace abinitio { BooleanOptionKey const cyclic_peptide( "abinitio:cyclic_peptide" );  }
namespace abinitio { BooleanOptionKey const debug_structures( "abinitio:debug_structures" );  }
namespace abinitio { FileOptionKey const log_frags( "abinitio:log_frags" );  }
namespace abinitio { BooleanOptionKey const only_stage1( "abinitio:only_stage1" );  }
namespace abinitio { RealOptionKey const end_bias( "abinitio:end_bias" );  }
namespace abinitio { BooleanOptionKey const apply_end_bias( "abinitio:apply_end_bias" );  }
namespace abinitio { IntegerOptionKey const symmetry_residue( "abinitio:symmetry_residue" );  }
namespace abinitio { RealOptionKey const vdw_weight_stage1( "abinitio:vdw_weight_stage1" );  }
namespace abinitio { BooleanOptionKey const override_vdw_all_stages( "abinitio:override_vdw_all_stages" );  }
namespace abinitio { IntegerVectorOptionKey const recover_low_in_stages( "abinitio:recover_low_in_stages" );  }
namespace abinitio { IntegerVectorOptionKey const skip_stages( "abinitio:skip_stages" );  }
namespace abinitio { BooleanOptionKey const close_chbrk( "abinitio:close_chbrk" );  }
namespace abinitio { BooleanOptionKey const include_stage5( "abinitio:include_stage5" );  }
namespace abinitio { BooleanOptionKey const close_loops_by_idealizing( "abinitio:close_loops_by_idealizing" );  }
namespace abinitio { BooleanOptionKey const optimize_cutpoints_using_kic( "abinitio:optimize_cutpoints_using_kic" );  }
namespace abinitio { IntegerOptionKey const optimize_cutpoints_margin( "abinitio:optimize_cutpoints_margin" );  }
namespace abinitio { namespace star { BooleanOptionKey const star( "abinitio:star" );  } }
namespace abinitio { namespace star { RealOptionKey const initial_dist_cutoff( "abinitio:star:initial_dist_cutoff" );  } }
namespace abinitio { namespace star { IntegerOptionKey const min_unaligned_len( "abinitio:star:min_unaligned_len" );  } }
namespace abrelax { BooleanOptionKey const abrelax( "abrelax" );  }
namespace abrelax { BooleanOptionKey const fail_unclosed( "abrelax:fail_unclosed" );  }
namespace AnchoredDesign { BooleanOptionKey const AnchoredDesign( "AnchoredDesign" );  }
namespace AnchoredDesign { FileOptionKey const anchor( "AnchoredDesign:anchor" );  }
namespace AnchoredDesign { BooleanOptionKey const allow_anchor_repack( "AnchoredDesign:allow_anchor_repack" );  }
namespace AnchoredDesign { BooleanOptionKey const vary_cutpoints( "AnchoredDesign:vary_cutpoints" );  }
namespace AnchoredDesign { BooleanOptionKey const no_frags( "AnchoredDesign:no_frags" );  }
namespace AnchoredDesign { BooleanOptionKey const debug( "AnchoredDesign:debug" );  }
namespace AnchoredDesign { BooleanOptionKey const show_extended( "AnchoredDesign:show_extended" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_only( "AnchoredDesign:refine_only" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_show( "AnchoredDesign:perturb_show" );  }
namespace AnchoredDesign { IntegerOptionKey const perturb_cycles( "AnchoredDesign:perturb_cycles" );  }
namespace AnchoredDesign { RealOptionKey const perturb_temp( "AnchoredDesign:perturb_temp" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_CCD_off( "AnchoredDesign:perturb_CCD_off" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_KIC_off( "AnchoredDesign:perturb_KIC_off" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_CCD_off( "AnchoredDesign:refine_CCD_off" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_KIC_off( "AnchoredDesign:refine_KIC_off" );  }
namespace AnchoredDesign { IntegerOptionKey const refine_cycles( "AnchoredDesign:refine_cycles" );  }
namespace AnchoredDesign { RealOptionKey const refine_temp( "AnchoredDesign:refine_temp" );  }
namespace AnchoredDesign { IntegerOptionKey const refine_repack_cycles( "AnchoredDesign:refine_repack_cycles" );  }
namespace AnchoredDesign { BooleanOptionKey const rmsd( "AnchoredDesign:rmsd" );  }
namespace AnchoredDesign { BooleanOptionKey const unbound_mode( "AnchoredDesign:unbound_mode" );  }
namespace AnchoredDesign { RealOptionKey const chainbreak_weight( "AnchoredDesign:chainbreak_weight" );  }
namespace AnchoredDesign { namespace filters { BooleanOptionKey const filters( "AnchoredDesign:filters" );  } }
namespace AnchoredDesign { namespace filters { RealOptionKey const score( "AnchoredDesign:filters:score" );  } }
namespace AnchoredDesign { namespace filters { RealOptionKey const sasa( "AnchoredDesign:filters:sasa" );  } }
namespace AnchoredDesign { namespace filters { BooleanOptionKey const omega( "AnchoredDesign:filters:omega" );  } }
namespace AnchoredDesign { namespace akash { BooleanOptionKey const akash( "AnchoredDesign:akash" );  } }
namespace AnchoredDesign { namespace akash { IntegerOptionKey const dyepos( "AnchoredDesign:akash:dyepos" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const testing( "AnchoredDesign:testing" );  } }
namespace AnchoredDesign { namespace testing { RealOptionKey const VDW_weight( "AnchoredDesign:testing:VDW_weight" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const anchor_via_constraints( "AnchoredDesign:testing:anchor_via_constraints" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const delete_interface_native_sidechains( "AnchoredDesign:testing:delete_interface_native_sidechains" );  } }
namespace AnchoredDesign { namespace testing { FileOptionKey const RMSD_only_this( "AnchoredDesign:testing:RMSD_only_this" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const anchor_noise_constraints_mode( "AnchoredDesign:testing:anchor_noise_constraints_mode" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const super_secret_fixed_interface_mode( "AnchoredDesign:testing:super_secret_fixed_interface_mode" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const randomize_input_sequence( "AnchoredDesign:testing:randomize_input_sequence" );  } }
namespace antibody { BooleanOptionKey const antibody( "antibody" );  }
namespace antibody { StringOptionKey const input_ab_scheme( "antibody:input_ab_scheme" );  }
namespace antibody { StringOptionKey const output_ab_scheme( "antibody:output_ab_scheme" );  }
namespace antibody { StringOptionKey const cdr_definition( "antibody:cdr_definition" );  }
namespace antibody { StringOptionKey const light_chain( "antibody:light_chain" );  }
namespace antibody { BooleanOptionKey const check_cdr_chainbreaks( "antibody:check_cdr_chainbreaks" );  }
namespace antibody { BooleanOptionKey const check_cdr_pep_bond_geom( "antibody:check_cdr_pep_bond_geom" );  }
namespace antibody { StringOptionKey const numbering_scheme( "antibody:numbering_scheme" );  }
namespace antibody { BooleanOptionKey const graft_l1( "antibody:graft_l1" );  }
namespace antibody { StringOptionKey const l1_template( "antibody:l1_template" );  }
namespace antibody { BooleanOptionKey const graft_l2( "antibody:graft_l2" );  }
namespace antibody { StringOptionKey const l2_template( "antibody:l2_template" );  }
namespace antibody { BooleanOptionKey const graft_l3( "antibody:graft_l3" );  }
namespace antibody { StringOptionKey const l3_template( "antibody:l3_template" );  }
namespace antibody { BooleanOptionKey const graft_h1( "antibody:graft_h1" );  }
namespace antibody { StringOptionKey const h1_template( "antibody:h1_template" );  }
namespace antibody { BooleanOptionKey const graft_h2( "antibody:graft_h2" );  }
namespace antibody { StringOptionKey const h2_template( "antibody:h2_template" );  }
namespace antibody { BooleanOptionKey const graft_h3( "antibody:graft_h3" );  }
namespace antibody { StringOptionKey const h3_template( "antibody:h3_template" );  }
namespace antibody { StringOptionKey const light_heavy_template( "antibody:light_heavy_template" );  }
namespace antibody { StringOptionKey const frl_template( "antibody:frl_template" );  }
namespace antibody { StringOptionKey const frh_template( "antibody:frh_template" );  }
namespace antibody { BooleanOptionKey const h3_no_stem_graft( "antibody:h3_no_stem_graft" );  }
namespace antibody { BooleanOptionKey const packonly_after_graft( "antibody:packonly_after_graft" );  }
namespace antibody { BooleanOptionKey const stem_optimize( "antibody:stem_optimize" );  }
namespace antibody { IntegerOptionKey const stem_optimize_size( "antibody:stem_optimize_size" );  }
namespace antibody { StringOptionKey const preprocessing_script_version( "antibody:preprocessing_script_version" );  }
namespace antibody { BooleanOptionKey const model_h3( "antibody:model_h3" );  }
namespace antibody { BooleanOptionKey const snugfit( "antibody:snugfit" );  }
namespace antibody { BooleanOptionKey const refine_h3( "antibody:refine_h3" );  }
namespace antibody { BooleanOptionKey const h3_filter( "antibody:h3_filter" );  }
namespace antibody { RealOptionKey const h3_filter_tolerance( "antibody:h3_filter_tolerance" );  }
namespace antibody { BooleanOptionKey const cter_insert( "antibody:cter_insert" );  }
namespace antibody { BooleanOptionKey const flank_residue_min( "antibody:flank_residue_min" );  }
namespace antibody { BooleanOptionKey const sc_min( "antibody:sc_min" );  }
namespace antibody { BooleanOptionKey const rt_min( "antibody:rt_min" );  }
namespace antibody { BooleanOptionKey const bad_nter( "antibody:bad_nter" );  }
namespace antibody { BooleanOptionKey const extend_h3_before_modeling( "antibody:extend_h3_before_modeling" );  }
namespace antibody { BooleanOptionKey const idealize_h3_stems_before_modeling( "antibody:idealize_h3_stems_before_modeling" );  }
namespace antibody { StringOptionKey const remodel( "antibody:remodel" );  }
namespace antibody { StringOptionKey const refine( "antibody:refine" );  }
namespace antibody { StringOptionKey const centroid_refine( "antibody:centroid_refine" );  }
namespace antibody { BooleanOptionKey const constrain_cter( "antibody:constrain_cter" );  }
namespace antibody { BooleanOptionKey const constrain_vlvh_qq( "antibody:constrain_vlvh_qq" );  }
namespace antibody { BooleanOptionKey const auto_generate_kink_constraint( "antibody:auto_generate_kink_constraint" );  }
namespace antibody { BooleanOptionKey const all_atom_mode_kink_constraint( "antibody:all_atom_mode_kink_constraint" );  }
namespace antibody { BooleanOptionKey const snug_loops( "antibody:snug_loops" );  }
namespace antibody { FileOptionKey const input_fv( "antibody:input_fv" );  }
namespace antibody { BooleanOptionKey const camelid( "antibody:camelid" );  }
namespace antibody { BooleanOptionKey const camelid_constraints( "antibody:camelid_constraints" );  }
namespace antibody { BooleanOptionKey const use_mean_cluster_cst_data( "antibody:use_mean_cluster_cst_data" );  }
namespace antibody { BooleanOptionKey const force_use_of_cluster_csts_with_outliers( "antibody:force_use_of_cluster_csts_with_outliers" );  }
namespace antibody { IntegerOptionKey const cluster_csts_stats_cutoff( "antibody:cluster_csts_stats_cutoff" );  }
namespace antibody { RealOptionKey const general_dihedral_cst_phi_sd( "antibody:general_dihedral_cst_phi_sd" );  }
namespace antibody { RealOptionKey const general_dihedral_cst_psi_sd( "antibody:general_dihedral_cst_psi_sd" );  }
namespace antibody { BooleanOptionKey const allow_omega_mismatches_for_north_clusters( "antibody:allow_omega_mismatches_for_north_clusters" );  }
namespace antibody { StringOptionKey const prefix( "antibody:prefix" );  }
namespace antibody { StringOptionKey const grafting_database( "antibody:grafting_database" );  }
namespace antibody { StringOptionKey const blastp( "antibody:blastp" );  }
namespace antibody { BooleanOptionKey const exclude_homologs( "antibody:exclude_homologs" );  }
namespace antibody { RealOptionKey const exclude_homologs_cdr_cutoff( "antibody:exclude_homologs_cdr_cutoff" );  }
namespace antibody { RealOptionKey const exclude_homologs_fr_cutoff( "antibody:exclude_homologs_fr_cutoff" );  }
namespace antibody { RealOptionKey const ocd_cutoff( "antibody:ocd_cutoff" );  }
namespace antibody { IntegerOptionKey const n_multi_templates( "antibody:n_multi_templates" );  }
namespace antibody { namespace design { BooleanOptionKey const design( "antibody:design" );  } }
namespace antibody { namespace design { StringOptionKey const base_cdr_instructions( "antibody:design:base_cdr_instructions" );  } }
namespace antibody { namespace design { StringOptionKey const cdr_instructions( "antibody:design:cdr_instructions" );  } }
namespace antibody { namespace design { StringOptionKey const antibody_database( "antibody:design:antibody_database" );  } }
namespace antibody { namespace design { BooleanOptionKey const paper_ab_db( "antibody:design:paper_ab_db" );  } }
namespace antibody { namespace design { StringOptionKey const paper_ab_db_path( "antibody:design:paper_ab_db_path" );  } }
namespace antibody { namespace design { StringVectorOptionKey const design_cdrs( "antibody:design:design_cdrs" );  } }
namespace antibody { namespace design { IntegerOptionKey const top_designs( "antibody:design:top_designs" );  } }
namespace antibody { namespace design { StringOptionKey const design_protocol( "antibody:design:design_protocol" );  } }
namespace antibody { namespace design { BooleanOptionKey const run_snugdock( "antibody:design:run_snugdock" );  } }
namespace antibody { namespace design { BooleanOptionKey const run_relax( "antibody:design:run_relax" );  } }
namespace antibody { namespace design { BooleanOptionKey const run_interface_analyzer( "antibody:design:run_interface_analyzer" );  } }
namespace antibody { namespace design { StringVectorOptionKey const paratope( "antibody:design:paratope" );  } }
namespace antibody { namespace design { StringVectorOptionKey const epitope( "antibody:design:epitope" );  } }
namespace antibody { namespace design { BooleanOptionKey const use_epitope_constraints( "antibody:design:use_epitope_constraints" );  } }
namespace antibody { namespace design { RealOptionKey const dihedral_cst_weight( "antibody:design:dihedral_cst_weight" );  } }
namespace antibody { namespace design { RealOptionKey const atom_pair_cst_weight( "antibody:design:atom_pair_cst_weight" );  } }
namespace antibody { namespace design { BooleanOptionKey const global_dihedral_cst_scoring( "antibody:design:global_dihedral_cst_scoring" );  } }
namespace antibody { namespace design { BooleanOptionKey const global_atom_pair_cst_scoring( "antibody:design:global_atom_pair_cst_scoring" );  } }
namespace antibody { namespace design { BooleanOptionKey const do_dock( "antibody:design:do_dock" );  } }
namespace antibody { namespace design { BooleanOptionKey const do_rb_min( "antibody:design:do_rb_min" );  } }
namespace antibody { namespace design { BooleanOptionKey const dock_min_dock( "antibody:design:dock_min_dock" );  } }
namespace antibody { namespace design { IntegerOptionKey const outer_cycle_rounds( "antibody:design:outer_cycle_rounds" );  } }
namespace antibody { namespace design { IntegerOptionKey const inner_cycle_rounds( "antibody:design:inner_cycle_rounds" );  } }
namespace antibody { namespace design { IntegerOptionKey const dock_cycle_rounds( "antibody:design:dock_cycle_rounds" );  } }
namespace antibody { namespace design { RealOptionKey const interface_dis( "antibody:design:interface_dis" );  } }
namespace antibody { namespace design { RealOptionKey const neighbor_dis( "antibody:design:neighbor_dis" );  } }
namespace antibody { namespace design { BooleanOptionKey const use_outliers( "antibody:design:use_outliers" );  } }
namespace antibody { namespace design { BooleanOptionKey const use_H3_graft_outliers( "antibody:design:use_H3_graft_outliers" );  } }
namespace antibody { namespace design { BooleanOptionKey const use_only_H3_kinked( "antibody:design:use_only_H3_kinked" );  } }
namespace antibody { namespace design { BooleanOptionKey const design_antigen( "antibody:design:design_antigen" );  } }
namespace antibody { namespace design { BooleanOptionKey const design_framework( "antibody:design:design_framework" );  } }
namespace antibody { namespace design { BooleanOptionKey const conservative_framework_design( "antibody:design:conservative_framework_design" );  } }
namespace antibody { namespace design { BooleanOptionKey const design_H3_stem( "antibody:design:design_H3_stem" );  } }
namespace antibody { namespace design { BooleanOptionKey const design_proline( "antibody:design:design_proline" );  } }
namespace antibody { namespace design { RealOptionKey const sample_zero_probs_at( "antibody:design:sample_zero_probs_at" );  } }
namespace antibody { namespace design { BooleanOptionKey const force_mutate_framework_for_cluster( "antibody:design:force_mutate_framework_for_cluster" );  } }
namespace antibody { namespace design { IntegerOptionKey const seq_design_stats_cutoff( "antibody:design:seq_design_stats_cutoff" );  } }
namespace antibody { namespace design { IntegerOptionKey const seq_design_profile_samples( "antibody:design:seq_design_profile_samples" );  } }
namespace antibody { namespace design { BooleanOptionKey const use_light_chain_type( "antibody:design:use_light_chain_type" );  } }
namespace antibody { namespace design { BooleanOptionKey const idealize_graft_cdrs( "antibody:design:idealize_graft_cdrs" );  } }
namespace antibody { namespace design { StringVectorOptionKey const add_backrub_pivots( "antibody:design:add_backrub_pivots" );  } }
namespace antibody { namespace design { RealOptionKey const inner_kt( "antibody:design:inner_kt" );  } }
namespace antibody { namespace design { RealOptionKey const outer_kt( "antibody:design:outer_kt" );  } }
namespace antibody { namespace design { BooleanOptionKey const random_start( "antibody:design:random_start" );  } }
namespace antibody { namespace design { BooleanOptionKey const adapt_graft( "antibody:design:adapt_graft" );  } }
namespace antibody { namespace design { BooleanOptionKey const enable_adapt_graft_cartesian( "antibody:design:enable_adapt_graft_cartesian" );  } }
namespace antibody { namespace design { BooleanOptionKey const remove_antigen( "antibody:design:remove_antigen" );  } }
namespace antibody { namespace design { BooleanOptionKey const add_graft_log_to_pdb( "antibody:design:add_graft_log_to_pdb" );  } }
namespace antibody { namespace design { BooleanOptionKey const mutate_framework_for_cluster( "antibody:design:mutate_framework_for_cluster" );  } }
namespace task_operations { BooleanOptionKey const task_operations( "task_operations" );  }
namespace task_operations { StringOptionKey const cons_design_data_source( "task_operations:cons_design_data_source" );  }
namespace assembly { BooleanOptionKey const assembly( "assembly" );  }
namespace assembly { FileOptionKey const pdb1( "assembly:pdb1" );  }
namespace assembly { FileOptionKey const pdb2( "assembly:pdb2" );  }
namespace backrub { BooleanOptionKey const backrub( "backrub" );  }
namespace backrub { IntegerVectorOptionKey const pivot_residues( "backrub:pivot_residues" );  }
namespace backrub { StringVectorOptionKey const pivot_atoms( "backrub:pivot_atoms" );  }
namespace backrub { IntegerOptionKey const min_atoms( "backrub:min_atoms" );  }
namespace backrub { IntegerOptionKey const max_atoms( "backrub:max_atoms" );  }
namespace backrub { IntegerOptionKey const ntrials( "backrub:ntrials" );  }
namespace backrub { RealOptionKey const sc_prob( "backrub:sc_prob" );  }
namespace backrub { RealOptionKey const sm_prob( "backrub:sm_prob" );  }
namespace backrub { RealOptionKey const sc_prob_uniform( "backrub:sc_prob_uniform" );  }
namespace backrub { RealOptionKey const sc_prob_withinrot( "backrub:sc_prob_withinrot" );  }
namespace backrub { RealOptionKey const mc_kt( "backrub:mc_kt" );  }
namespace backrub { RealOptionKey const mm_bend_weight( "backrub:mm_bend_weight" );  }
namespace backrub { BooleanOptionKey const initial_pack( "backrub:initial_pack" );  }
namespace backrub { FileOptionKey const minimize_movemap( "backrub:minimize_movemap" );  }
namespace backrub { BooleanOptionKey const trajectory( "backrub:trajectory" );  }
namespace backrub { BooleanOptionKey const trajectory_gz( "backrub:trajectory_gz" );  }
namespace backrub { IntegerOptionKey const trajectory_stride( "backrub:trajectory_stride" );  }
namespace batch_relax { BooleanOptionKey const batch_relax( "batch_relax" );  }
namespace batch_relax { IntegerOptionKey const batch_size( "batch_relax:batch_size" );  }
namespace bbg { BooleanOptionKey const bbg( "bbg" );  }
namespace bbg { RealOptionKey const factorA( "bbg:factorA" );  }
namespace bbg { RealOptionKey const factorB( "bbg:factorB" );  }
namespace bbg { BooleanOptionKey const ignore_improper_res( "bbg:ignore_improper_res" );  }
namespace bbg { BooleanOptionKey const fix_short_segment( "bbg:fix_short_segment" );  }
namespace boinc { BooleanOptionKey const boinc( "boinc" );  }
namespace boinc { BooleanOptionKey const graphics( "boinc:graphics" );  }
namespace boinc { BooleanOptionKey const fullscreen( "boinc:fullscreen" );  }
namespace boinc { IntegerOptionKey const max_fps( "boinc:max_fps" );  }
namespace boinc { IntegerOptionKey const max_cpu( "boinc:max_cpu" );  }
namespace boinc { BooleanOptionKey const noshmem( "boinc:noshmem" );  }
namespace boinc { IntegerOptionKey const cpu_run_time( "boinc:cpu_run_time" );  }
namespace boinc { IntegerOptionKey const max_nstruct( "boinc:max_nstruct" );  }
namespace boinc { RealOptionKey const cpu_frac( "boinc:cpu_frac" );  }
namespace boinc { RealOptionKey const frame_rate( "boinc:frame_rate" );  }
namespace boinc { BooleanOptionKey const watchdog( "boinc:watchdog" );  }
namespace boinc { IntegerOptionKey const watchdog_time( "boinc:watchdog_time" );  }
namespace boinc { IntegerOptionKey const cpu_run_timeout( "boinc:cpu_run_timeout" );  }
namespace boinc { FileOptionKey const description_file( "boinc:description_file" );  }
namespace boinc { RealOptionKey const score_cut_pct( "boinc:score_cut_pct" );  }
