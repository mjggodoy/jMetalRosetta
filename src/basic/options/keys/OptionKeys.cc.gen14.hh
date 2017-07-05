namespace stepwise { BooleanOptionKey const skip_preminimize( "stepwise:skip_preminimize" );  }
namespace stepwise { BooleanOptionKey const minimize_waters( "stepwise:minimize_waters" );  }
namespace stepwise { BooleanOptionKey const test_all_moves( "stepwise:test_all_moves" );  }
namespace stepwise { BooleanOptionKey const new_move_selector( "stepwise:new_move_selector" );  }
namespace stepwise { BooleanOptionKey const dump( "stepwise:dump" );  }
namespace stepwise { BooleanOptionKey const VERBOSE( "stepwise:VERBOSE" );  }
namespace stepwise { BooleanOptionKey const use_green_packer( "stepwise:use_green_packer" );  }
namespace stepwise { RealOptionKey const rmsd_screen( "stepwise:rmsd_screen" );  }
namespace stepwise { BooleanOptionKey const skip_minimize( "stepwise:skip_minimize" );  }
namespace stepwise { BooleanOptionKey const virtualize_packable_moieties_in_screening_pose( "stepwise:virtualize_packable_moieties_in_screening_pose" );  }
namespace stepwise { StringOptionKey const minimizer_mode( "stepwise:minimizer_mode" );  }
namespace stepwise { StringOptionKey const sampler_silent_file( "stepwise:sampler_silent_file" );  }
namespace stepwise { BooleanOptionKey const superimpose_over_all( "stepwise:superimpose_over_all" );  }
namespace stepwise { StringVectorOptionKey const move( "stepwise:move" );  }
namespace stepwise { StringOptionKey const min_type( "stepwise:min_type" );  }
namespace stepwise { RealOptionKey const min_tolerance( "stepwise:min_tolerance" );  }
namespace stepwise { BooleanOptionKey const vary_polar_hydrogen_geometry( "stepwise:vary_polar_hydrogen_geometry" );  }
namespace stepwise { BooleanOptionKey const output_minimized_pose_list( "stepwise:output_minimized_pose_list" );  }
namespace stepwise { BooleanOptionKey const virtualize_free_moieties_in_native( "stepwise:virtualize_free_moieties_in_native" );  }
namespace stepwise { BooleanOptionKey const output_cluster_size( "stepwise:output_cluster_size" );  }
namespace stepwise { BooleanOptionKey const lores( "stepwise:lores" );  }
namespace stepwise { BooleanOptionKey const verbose_sampler( "stepwise:verbose_sampler" );  }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const monte_carlo( "stepwise:monte_carlo" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const verbose_scores( "stepwise:monte_carlo:verbose_scores" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const skip_deletions( "stepwise:monte_carlo:skip_deletions" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_internal_hinge_moves( "stepwise:monte_carlo:allow_internal_hinge_moves" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_internal_local_moves( "stepwise:monte_carlo:allow_internal_local_moves" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_skip_bulge( "stepwise:monte_carlo:allow_skip_bulge" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const skip_bulge_frequency( "stepwise:monte_carlo:skip_bulge_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const from_scratch_frequency( "stepwise:monte_carlo:from_scratch_frequency" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_split_off( "stepwise:monte_carlo:allow_split_off" );  } }
namespace stepwise { namespace monte_carlo { IntegerOptionKey const cycles( "stepwise:monte_carlo:cycles" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const temperature( "stepwise:monte_carlo:temperature" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const add_proposal_density_factor( "stepwise:monte_carlo:add_proposal_density_factor" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const add_delete_frequency( "stepwise:monte_carlo:add_delete_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const docking_frequency( "stepwise:monte_carlo:docking_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const submotif_frequency( "stepwise:monte_carlo:submotif_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const intermolecular_frequency( "stepwise:monte_carlo:intermolecular_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const minimize_single_res_frequency( "stepwise:monte_carlo:minimize_single_res_frequency" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_variable_bond_geometry( "stepwise:monte_carlo:allow_variable_bond_geometry" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const switch_focus_frequency( "stepwise:monte_carlo:switch_focus_frequency" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const just_min_after_mutation_frequency( "stepwise:monte_carlo:just_min_after_mutation_frequency" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const local_redock_only( "stepwise:monte_carlo:local_redock_only" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const make_movie( "stepwise:monte_carlo:make_movie" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const recover_low( "stepwise:monte_carlo:recover_low" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const use_precomputed_library( "stepwise:monte_carlo:use_precomputed_library" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const allow_submotif_split( "stepwise:monte_carlo:allow_submotif_split" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const force_submotif_without_intervening_bulge( "stepwise:monte_carlo:force_submotif_without_intervening_bulge" );  } }
namespace stepwise { namespace monte_carlo { BooleanOptionKey const use_first_jump_for_submotif( "stepwise:monte_carlo:use_first_jump_for_submotif" );  } }
namespace stepwise { namespace monte_carlo { RealOptionKey const vary_loop_length_frequency( "stepwise:monte_carlo:vary_loop_length_frequency" );  } }
namespace stepwise { namespace monte_carlo { IntegerOptionKey const checkpointing_frequency( "stepwise:monte_carlo:checkpointing_frequency" );  } }
namespace stepwise { namespace monte_carlo { namespace csa { BooleanOptionKey const csa( "stepwise:monte_carlo:csa" );  } } }
namespace stepwise { namespace monte_carlo { namespace csa { IntegerOptionKey const csa_bank_size( "stepwise:monte_carlo:csa:csa_bank_size" );  } } }
namespace stepwise { namespace monte_carlo { namespace csa { RealOptionKey const csa_rmsd( "stepwise:monte_carlo:csa:csa_rmsd" );  } } }
namespace stepwise { namespace monte_carlo { namespace csa { BooleanOptionKey const csa_output_rounds( "stepwise:monte_carlo:csa:csa_output_rounds" );  } } }
namespace stepwise { namespace rna { BooleanOptionKey const rna( "stepwise:rna" );  } }
namespace stepwise { namespace rna { IntegerOptionKey const sampler_num_pose_kept( "stepwise:rna:sampler_num_pose_kept" );  } }
namespace stepwise { namespace rna { RealOptionKey const native_edensity_score_cutoff( "stepwise:rna:native_edensity_score_cutoff" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const o2prime_legacy_mode( "stepwise:rna:o2prime_legacy_mode" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const allow_virtual_o2prime_hydrogens( "stepwise:rna:allow_virtual_o2prime_hydrogens" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_perform_phosphate_pack( "stepwise:rna:sampler_perform_phosphate_pack" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const force_phosphate_instantiation( "stepwise:rna:force_phosphate_instantiation" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const distinguish_pucker( "stepwise:rna:distinguish_pucker" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const finer_sampling_at_chain_closure( "stepwise:rna:finer_sampling_at_chain_closure" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const PBP_clustering_at_chain_closure( "stepwise:rna:PBP_clustering_at_chain_closure" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_allow_syn_pyrimidine( "stepwise:rna:sampler_allow_syn_pyrimidine" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_extra_chi_rotamer( "stepwise:rna:sampler_extra_chi_rotamer" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_extra_beta_rotamer( "stepwise:rna:sampler_extra_beta_rotamer" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_extra_epsilon_rotamer( "stepwise:rna:sampler_extra_epsilon_rotamer" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const force_centroid_interaction( "stepwise:rna:force_centroid_interaction" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const virtual_sugar_legacy_mode( "stepwise:rna:virtual_sugar_legacy_mode" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const VDW_rep_optimize_memory_usage( "stepwise:rna:VDW_rep_optimize_memory_usage" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const erraser( "stepwise:rna:erraser" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const centroid_screen( "stepwise:rna:centroid_screen" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const VDW_atr_rep_screen( "stepwise:rna:VDW_atr_rep_screen" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const minimize_and_score_native_pose( "stepwise:rna:minimize_and_score_native_pose" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const rm_virt_phosphate( "stepwise:rna:rm_virt_phosphate" );  } }
namespace stepwise { namespace rna { StringVectorOptionKey const VDW_rep_screen_info( "stepwise:rna:VDW_rep_screen_info" );  } }
namespace stepwise { namespace rna { RealOptionKey const VDW_rep_alignment_RMSD_CUTOFF( "stepwise:rna:VDW_rep_alignment_RMSD_CUTOFF" );  } }
namespace stepwise { namespace rna { StringVectorOptionKey const VDW_rep_delete_matching_res( "stepwise:rna:VDW_rep_delete_matching_res" );  } }
namespace stepwise { namespace rna { RealOptionKey const VDW_rep_screen_physical_pose_clash_dist_cutoff( "stepwise:rna:VDW_rep_screen_physical_pose_clash_dist_cutoff" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const integration_test( "stepwise:rna:integration_test" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const allow_bulge_at_chainbreak( "stepwise:rna:allow_bulge_at_chainbreak" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const parin_favorite_output( "stepwise:rna:parin_favorite_output" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const reinitialize_CCD_torsions( "stepwise:rna:reinitialize_CCD_torsions" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sample_both_sugar_base_rotamer( "stepwise:rna:sample_both_sugar_base_rotamer" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_include_torsion_value_in_tag( "stepwise:rna:sampler_include_torsion_value_in_tag" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_assert_no_virt_sugar_sampling( "stepwise:rna:sampler_assert_no_virt_sugar_sampling" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sampler_try_sugar_instantiation( "stepwise:rna:sampler_try_sugar_instantiation" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const do_not_sample_multiple_virtual_sugar( "stepwise:rna:do_not_sample_multiple_virtual_sugar" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const sample_ONLY_multiple_virtual_sugar( "stepwise:rna:sample_ONLY_multiple_virtual_sugar" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const allow_base_pair_only_centroid_screen( "stepwise:rna:allow_base_pair_only_centroid_screen" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const minimizer_rename_tag( "stepwise:rna:minimizer_rename_tag" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const minimize_res( "stepwise:rna:minimize_res" );  } }
namespace stepwise { namespace rna { StringVectorOptionKey const alignment_res( "stepwise:rna:alignment_res" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const native_alignment_res( "stepwise:rna:native_alignment_res" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const rmsd_res( "stepwise:rna:rmsd_res" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const missing_res( "stepwise:rna:missing_res" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const missing_res2( "stepwise:rna:missing_res2" );  } }
namespace stepwise { namespace rna { IntegerOptionKey const job_queue_ID( "stepwise:rna:job_queue_ID" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const minimize_and_score_sugar( "stepwise:rna:minimize_and_score_sugar" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const global_sample_res_list( "stepwise:rna:global_sample_res_list" );  } }
namespace stepwise { namespace rna { FileOptionKey const filter_output_filename( "stepwise:rna:filter_output_filename" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const combine_long_loop_mode( "stepwise:rna:combine_long_loop_mode" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const combine_helical_silent_file( "stepwise:rna:combine_helical_silent_file" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const output_extra_RMSDs( "stepwise:rna:output_extra_RMSDs" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const protonated_H1_adenosine_list( "stepwise:rna:protonated_H1_adenosine_list" );  } }
namespace stepwise { namespace rna { IntegerVectorOptionKey const native_virtual_res( "stepwise:rna:native_virtual_res" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const simple_append_map( "stepwise:rna:simple_append_map" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const allow_fixed_res_at_moving_res( "stepwise:rna:allow_fixed_res_at_moving_res" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const force_user_defined_jumps( "stepwise:rna:force_user_defined_jumps" );  } }
namespace stepwise { namespace rna { StringVectorOptionKey const jump_point_pairs( "stepwise:rna:jump_point_pairs" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const add_virt_root( "stepwise:rna:add_virt_root" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const floating_base( "stepwise:rna:floating_base" );  } }
namespace stepwise { namespace rna { IntegerOptionKey const floating_base_anchor_res( "stepwise:rna:floating_base_anchor_res" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const allow_chain_boundary_jump_partner_right_at_fixed_BP( "stepwise:rna:allow_chain_boundary_jump_partner_right_at_fixed_BP" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const rebuild_bulge_mode( "stepwise:rna:rebuild_bulge_mode" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const virtual_sugar_keep_base_fixed( "stepwise:rna:virtual_sugar_keep_base_fixed" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const virtual_sugar_do_minimize( "stepwise:rna:virtual_sugar_do_minimize" );  } }
namespace stepwise { namespace rna { RealOptionKey const sampler_max_centroid_distance( "stepwise:rna:sampler_max_centroid_distance" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const filter_user_alignment_res( "stepwise:rna:filter_user_alignment_res" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const tether_jump( "stepwise:rna:tether_jump" );  } }
namespace stepwise { namespace rna { BooleanOptionKey const turn_off_rna_chem_map_during_optimize( "stepwise:rna:turn_off_rna_chem_map_during_optimize" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const protein( "stepwise:protein" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const global_optimize( "stepwise:protein:global_optimize" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const disable_sampling_of_loop_takeoff( "stepwise:protein:disable_sampling_of_loop_takeoff" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const sample_beta( "stepwise:protein:sample_beta" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const ghost_loops( "stepwise:protein:ghost_loops" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const centroid_screen( "stepwise:protein:centroid_screen" );  } }
namespace stepwise { namespace protein { RealOptionKey const centroid_score_diff_cut( "stepwise:protein:centroid_score_diff_cut" );  } }
namespace stepwise { namespace protein { StringOptionKey const centroid_weights( "stepwise:protein:centroid_weights" );  } }
namespace stepwise { namespace protein { RealOptionKey const score_diff_cut( "stepwise:protein:score_diff_cut" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const filter_native_big_bins( "stepwise:protein:filter_native_big_bins" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const cluster_by_all_atom_rmsd( "stepwise:protein:cluster_by_all_atom_rmsd" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const centroid_output( "stepwise:protein:centroid_output" );  } }
namespace stepwise { namespace protein { IntegerOptionKey const n_sample( "stepwise:protein:n_sample" );  } }
namespace stepwise { namespace protein { IntegerOptionKey const nstruct_centroid( "stepwise:protein:nstruct_centroid" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const ccd_close( "stepwise:protein:ccd_close" );  } }
namespace stepwise { namespace protein { IntegerVectorOptionKey const bridge_res( "stepwise:protein:bridge_res" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const cart_min( "stepwise:protein:cart_min" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const move_jumps_between_chains( "stepwise:protein:move_jumps_between_chains" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const use_packer_instead_of_rotamer_trials( "stepwise:protein:use_packer_instead_of_rotamer_trials" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const expand_loop_takeoff( "stepwise:protein:expand_loop_takeoff" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const skip_coord_constraints( "stepwise:protein:skip_coord_constraints" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const allow_virtual_side_chains( "stepwise:protein:allow_virtual_side_chains" );  } }
namespace stepwise { namespace protein { BooleanOptionKey const protein_prepack( "stepwise:protein:protein_prepack" );  } }
namespace stepwise { namespace protein { StringOptionKey const disulfide_file( "stepwise:protein:disulfide_file" );  } }
namespace full_model { BooleanOptionKey const full_model( "full_model" );  }
namespace full_model { ResidueChainVectorOptionKey const cutpoint_open( "full_model:cutpoint_open" );  }
namespace full_model { ResidueChainVectorOptionKey const cutpoint_closed( "full_model:cutpoint_closed" );  }
namespace full_model { ResidueChainVectorOptionKey const cyclize( "full_model:cyclize" );  }
namespace full_model { StringVectorOptionKey const other_poses( "full_model:other_poses" );  }
namespace full_model { ResidueChainVectorOptionKey const jump_res( "full_model:jump_res" );  }
namespace full_model { ResidueChainVectorOptionKey const extra_min_res( "full_model:extra_min_res" );  }
namespace full_model { ResidueChainVectorOptionKey const extra_min_jump_res( "full_model:extra_min_jump_res" );  }
namespace full_model { ResidueChainVectorOptionKey const root_res( "full_model:root_res" );  }
namespace full_model { ResidueChainVectorOptionKey const virtual_sugar_res( "full_model:virtual_sugar_res" );  }
namespace full_model { ResidueChainVectorOptionKey const virtual_res( "full_model:virtual_res" );  }
namespace full_model { ResidueChainVectorOptionKey const sample_res( "full_model:sample_res" );  }
namespace full_model { ResidueChainVectorOptionKey const calc_rms_res( "full_model:calc_rms_res" );  }
namespace full_model { ResidueChainVectorOptionKey const working_res( "full_model:working_res" );  }
namespace full_model { BooleanOptionKey const motif_mode( "full_model:motif_mode" );  }
namespace full_model { BooleanOptionKey const allow_jump_in_numbering( "full_model:allow_jump_in_numbering" );  }
namespace full_model { namespace rna { BooleanOptionKey const rna( "full_model:rna" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const terminal_res( "full_model:rna:terminal_res" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const block_stack_above_res( "full_model:rna:block_stack_above_res" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const block_stack_below_res( "full_model:rna:block_stack_below_res" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const force_syn_chi_res_list( "full_model:rna:force_syn_chi_res_list" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const force_anti_chi_res_list( "full_model:rna:force_anti_chi_res_list" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const force_north_sugar_list( "full_model:rna:force_north_sugar_list" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const force_south_sugar_list( "full_model:rna:force_south_sugar_list" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const bulge_res( "full_model:rna:bulge_res" );  } }
namespace full_model { namespace rna { ResidueChainVectorOptionKey const sample_sugar_res( "full_model:rna:sample_sugar_res" );  } }
namespace recces { StringOptionKey const seq1( "recces:seq1" );  }
namespace recces { StringOptionKey const seq2( "recces:seq2" );  }
namespace recces { IntegerOptionKey const n_cycle( "recces:n_cycle" );  }
namespace recces { RealOptionKey const a_form_range( "recces:a_form_range" );  }
namespace recces { BooleanOptionKey const dump_pdb( "recces:dump_pdb" );  }
namespace recces { BooleanOptionKey const dump_silent( "recces:dump_silent" );  }
namespace recces { BooleanOptionKey const out_torsions( "recces:out_torsions" );  }
namespace recces { RealVectorOptionKey const temps( "recces:temps" );  }
namespace recces { RealVectorOptionKey const st_weights( "recces:st_weights" );  }
namespace recces { BooleanOptionKey const save_score_terms( "recces:save_score_terms" );  }
namespace recces { StringOptionKey const out_prefix( "recces:out_prefix" );  }
namespace recces { IntegerOptionKey const dump_freq( "recces:dump_freq" );  }
namespace recces { IntegerOptionKey const n_intermediate_dump( "recces:n_intermediate_dump" );  }
namespace recces { BooleanOptionKey const output_min_pose( "recces:output_min_pose" );  }
namespace recces { BooleanOptionKey const accept_no_op_moves( "recces:accept_no_op_moves" );  }
namespace recces { RealOptionKey const histogram_min( "recces:histogram_min" );  }
namespace recces { RealOptionKey const histogram_max( "recces:histogram_max" );  }
namespace recces { RealOptionKey const histogram_spacing( "recces:histogram_spacing" );  }
namespace recces { namespace base_pair { BooleanOptionKey const base_pair( "recces:base_pair" );  } }
namespace recces { namespace base_pair { RealOptionKey const rmsd_cutoff( "recces:base_pair:rmsd_cutoff" );  } }
namespace recces { namespace base_pair { RealOptionKey const translation_mag( "recces:base_pair:translation_mag" );  } }
namespace recces { namespace base_pair { RealOptionKey const rotation_mag( "recces:base_pair:rotation_mag" );  } }
namespace recces { namespace base_pair { BooleanOptionKey const recces( "recces:base_pair:recces" );  } }
namespace recces { namespace base_pair { BooleanOptionKey const block_stack( "recces:base_pair:block_stack" );  } }
namespace recces { namespace base_pair { BooleanOptionKey const sample_jump( "recces:base_pair:sample_jump" );  } }
namespace recces { namespace thermal_sampling { BooleanOptionKey const thermal_sampling( "recces:thermal_sampling" );  } }
namespace recces { namespace thermal_sampling { IntegerVectorOptionKey const sample_residues( "recces:thermal_sampling:sample_residues" );  } }
namespace recces { namespace thermal_sampling { IntegerVectorOptionKey const free_residues( "recces:thermal_sampling:free_residues" );  } }
namespace recces { namespace thermal_sampling { IntegerVectorOptionKey const loop_residues( "recces:thermal_sampling:loop_residues" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const angle_range_bb( "recces:thermal_sampling:angle_range_bb" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const angle_range_free_bb( "recces:thermal_sampling:angle_range_free_bb" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const angle_range_chi( "recces:thermal_sampling:angle_range_chi" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const angle_range_free_chi( "recces:thermal_sampling:angle_range_free_chi" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const chi_stdev( "recces:thermal_sampling:chi_stdev" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const bb_stdev( "recces:thermal_sampling:bb_stdev" );  } }
namespace recces { namespace thermal_sampling { RealOptionKey const standard_bb_stdev( "recces:thermal_sampling:standard_bb_stdev" );  } }
namespace recces { namespace thermal_sampling { BooleanOptionKey const setup_base_pair_constraints( "recces:thermal_sampling:setup_base_pair_constraints" );  } }
namespace strand_assembly { BooleanOptionKey const strand_assembly( "strand_assembly" );  }
namespace strand_assembly { IntegerOptionKey const min_num_strands_to_deal( "strand_assembly:min_num_strands_to_deal" );  }
namespace strand_assembly { IntegerOptionKey const max_num_strands_to_deal( "strand_assembly:max_num_strands_to_deal" );  }
namespace strand_assembly { BooleanOptionKey const extract_native_only( "strand_assembly:extract_native_only" );  }
namespace strand_assembly { IntegerOptionKey const min_res_in_strand( "strand_assembly:min_res_in_strand" );  }
namespace strand_assembly { IntegerOptionKey const max_res_in_strand( "strand_assembly:max_res_in_strand" );  }
namespace strand_assembly { RealOptionKey const min_O_N_dis( "strand_assembly:min_O_N_dis" );  }
namespace strand_assembly { RealOptionKey const max_O_N_dis( "strand_assembly:max_O_N_dis" );  }
namespace strand_assembly { RealOptionKey const min_sheet_dis( "strand_assembly:min_sheet_dis" );  }
namespace strand_assembly { RealOptionKey const max_sheet_dis( "strand_assembly:max_sheet_dis" );  }
namespace strand_assembly { RealOptionKey const min_sheet_torsion( "strand_assembly:min_sheet_torsion" );  }
namespace strand_assembly { RealOptionKey const max_sheet_torsion( "strand_assembly:max_sheet_torsion" );  }
namespace strand_assembly { RealOptionKey const min_sheet_angle( "strand_assembly:min_sheet_angle" );  }
namespace strand_assembly { RealOptionKey const max_sheet_angle( "strand_assembly:max_sheet_angle" );  }
namespace strand_assembly { RealOptionKey const min_shortest_dis_sidechain_inter_sheet( "strand_assembly:min_shortest_dis_sidechain_inter_sheet" );  }
namespace TailSegment { BooleanOptionKey const TailSegment( "TailSegment" );  }
namespace TailSegment { IntegerOptionKey const refine_cycles( "TailSegment:refine_cycles" );  }
namespace TailSegment { IntegerOptionKey const refine_repack_cycles( "TailSegment:refine_repack_cycles" );  }
namespace templates { BooleanOptionKey const templates( "templates" );  }
namespace templates { FileOptionKey const config( "templates:config" );  }
namespace templates { BooleanOptionKey const fix_aligned_residues( "templates:fix_aligned_residues" );  }
namespace templates { FileOptionKey const fix_frag_file( "templates:fix_frag_file" );  }
namespace templates { IntegerOptionKey const fix_margin( "templates:fix_margin" );  }
namespace templates { IntegerOptionKey const min_nr_large_frags( "templates:min_nr_large_frags" );  }
namespace templates { IntegerOptionKey const min_nr_small_frags( "templates:min_nr_small_frags" );  }
namespace templates { BooleanOptionKey const no_pick_fragments( "templates:no_pick_fragments" );  }
namespace templates { IntegerOptionKey const nr_large_copies( "templates:nr_large_copies" );  }
namespace templates { IntegerOptionKey const nr_small_copies( "templates:nr_small_copies" );  }
namespace templates { BooleanOptionKey const pairings( "templates:pairings" );  }
namespace templates { BooleanOptionKey const pick_multiple_sizes( "templates:pick_multiple_sizes" );  }
namespace templates { BooleanOptionKey const strand_constraint( "templates:strand_constraint" );  }
namespace templates { BooleanOptionKey const vary_frag_size( "templates:vary_frag_size" );  }
namespace templates { BooleanOptionKey const no_culling( "templates:no_culling" );  }
namespace templates { FileOptionKey const helix_pairings( "templates:helix_pairings" );  }
namespace templates { FileOptionKey const prefix( "templates:prefix" );  }
namespace templates { IntegerOptionKey const change_movemap( "templates:change_movemap" );  }
