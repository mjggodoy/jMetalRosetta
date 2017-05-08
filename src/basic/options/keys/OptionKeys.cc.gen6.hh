namespace boinc { FileOptionKey const score_cut_fl( "boinc:score_cut_fl" );  }
namespace boinc { BooleanOptionKey const score_cut_smart_throttle( "boinc:score_cut_smart_throttle" );  }
namespace broker { BooleanOptionKey const broker( "broker" );  }
namespace broker { FileVectorOptionKey const setup( "broker:setup" );  }
namespace broker { RealOptionKey const rb_mover_stage1_weight( "broker:rb_mover_stage1_weight" );  }
namespace broker { RealOptionKey const large_frag_mover_stage1_weight( "broker:large_frag_mover_stage1_weight" );  }
namespace broker { RealOptionKey const small_frag_mover_stage1_weight( "broker:small_frag_mover_stage1_weight" );  }
namespace bunsat_calc2 { BooleanOptionKey const bunsat_calc2( "bunsat_calc2" );  }
namespace bunsat_calc2 { RealOptionKey const sasa_burial_cutoff( "bunsat_calc2:sasa_burial_cutoff" );  }
namespace bunsat_calc2 { BooleanOptionKey const layered_sasa( "bunsat_calc2:layered_sasa" );  }
namespace bunsat_calc2 { BooleanOptionKey const generous_hbonds( "bunsat_calc2:generous_hbonds" );  }
namespace bunsat_calc2 { RealOptionKey const AHD_cutoff( "bunsat_calc2:AHD_cutoff" );  }
namespace bunsat_calc2 { RealOptionKey const dist_cutoff( "bunsat_calc2:dist_cutoff" );  }
namespace bunsat_calc2 { RealOptionKey const hxl_dist_cutoff( "bunsat_calc2:hxl_dist_cutoff" );  }
namespace bunsat_calc2 { RealOptionKey const sulph_dist_cutoff( "bunsat_calc2:sulph_dist_cutoff" );  }
namespace bunsat_calc2 { RealOptionKey const metal_dist_cutoff( "bunsat_calc2:metal_dist_cutoff" );  }
namespace canonical_sampling { BooleanOptionKey const canonical_sampling( "canonical_sampling" );  }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const probabilities( "canonical_sampling:probabilities" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc( "canonical_sampling:probabilities:sc" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const localbb( "canonical_sampling:probabilities:localbb" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_uniform( "canonical_sampling:probabilities:sc_prob_uniform" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_withinrot( "canonical_sampling:probabilities:sc_prob_withinrot" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_perturbcurrent( "canonical_sampling:probabilities:sc_prob_perturbcurrent" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const MPI_sync_pools( "canonical_sampling:probabilities:MPI_sync_pools" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const MPI_bcast( "canonical_sampling:probabilities:MPI_bcast" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const fast_sc_moves( "canonical_sampling:probabilities:fast_sc_moves" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const fast_sc_moves_ntrials( "canonical_sampling:probabilities:fast_sc_moves_ntrials" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const no_jd2_output( "canonical_sampling:probabilities:no_jd2_output" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const use_hierarchical_clustering( "canonical_sampling:probabilities:use_hierarchical_clustering" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const backrub( "canonical_sampling:probabilities:backrub" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const conrot( "canonical_sampling:probabilities:conrot" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const sampling( "canonical_sampling:sampling" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const no_detailed_balance( "canonical_sampling:sampling:no_detailed_balance" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const ntrials( "canonical_sampling:sampling:ntrials" );  } }
namespace canonical_sampling { namespace sampling { RealOptionKey const mc_kt( "canonical_sampling:sampling:mc_kt" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const interval_pose_dump( "canonical_sampling:sampling:interval_pose_dump" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const interval_data_dump( "canonical_sampling:sampling:interval_data_dump" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const output_only_cluster_transitions( "canonical_sampling:sampling:output_only_cluster_transitions" );  } }
namespace canonical_sampling { namespace sampling { RealOptionKey const transition_threshold( "canonical_sampling:sampling:transition_threshold" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const max_files_per_dir( "canonical_sampling:sampling:max_files_per_dir" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const save_loops_only( "canonical_sampling:sampling:save_loops_only" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const dump_loops_only( "canonical_sampling:sampling:dump_loops_only" );  } }
namespace canonical_sampling { namespace out { BooleanOptionKey const out( "canonical_sampling:out" );  } }
namespace canonical_sampling { namespace out { FileOptionKey const new_structures( "canonical_sampling:out:new_structures" );  } }
namespace casp { BooleanOptionKey const casp( "casp" );  }
namespace casp { RealOptionKey const opt_radius( "casp:opt_radius" );  }
namespace casp { BooleanOptionKey const repack( "casp:repack" );  }
namespace casp { BooleanOptionKey const sc_min( "casp:sc_min" );  }
namespace casp { BooleanOptionKey const sequential( "casp:sequential" );  }
namespace casp { RealOptionKey const num_iterations( "casp:num_iterations" );  }
namespace casp { StringOptionKey const refine_res( "casp:refine_res" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const chemically_conjugated_docking( "chemically_conjugated_docking" );  }
namespace chemically_conjugated_docking { FileOptionKey const UBQpdb( "chemically_conjugated_docking:UBQpdb" );  }
namespace chemically_conjugated_docking { FileOptionKey const E2pdb( "chemically_conjugated_docking:E2pdb" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const E2_residue( "chemically_conjugated_docking:E2_residue" );  }
namespace chemically_conjugated_docking { RealOptionKey const SASAfilter( "chemically_conjugated_docking:SASAfilter" );  }
namespace chemically_conjugated_docking { RealOptionKey const scorefilter( "chemically_conjugated_docking:scorefilter" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const publication( "chemically_conjugated_docking:publication" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const n_tail_res( "chemically_conjugated_docking:n_tail_res" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const two_ubiquitins( "chemically_conjugated_docking:two_ubiquitins" );  }
namespace chemically_conjugated_docking { FileVectorOptionKey const extra_bodies( "chemically_conjugated_docking:extra_bodies" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const UBQ2_lys( "chemically_conjugated_docking:UBQ2_lys" );  }
namespace chemically_conjugated_docking { FileOptionKey const UBQ2_pdb( "chemically_conjugated_docking:UBQ2_pdb" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const dont_minimize_omega( "chemically_conjugated_docking:dont_minimize_omega" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const pdz( "chemically_conjugated_docking:pdz" );  }
namespace chemically_conjugated_docking { FileOptionKey const GTPasepdb( "chemically_conjugated_docking:GTPasepdb" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const GTPase_residue( "chemically_conjugated_docking:GTPase_residue" );  }
namespace chunk { BooleanOptionKey const chunk( "chunk" );  }
namespace chunk { FileOptionKey const pdb2( "chunk:pdb2" );  }
namespace chunk { FileOptionKey const loop2( "chunk:loop2" );  }
namespace cluster { BooleanOptionKey const cluster( "cluster" );  }
namespace cluster { BooleanOptionKey const lite( "cluster:lite" );  }
namespace cluster { RealOptionKey const input_score_filter( "cluster:input_score_filter" );  }
namespace cluster { RealOptionKey const output_score_filter( "cluster:output_score_filter" );  }
namespace cluster { IntegerVectorOptionKey const exclude_res( "cluster:exclude_res" );  }
namespace cluster { RealOptionKey const thinout_factor( "cluster:thinout_factor" );  }
namespace cluster { RealOptionKey const radius( "cluster:radius" );  }
namespace cluster { IntegerOptionKey const limit_cluster_size( "cluster:limit_cluster_size" );  }
namespace cluster { RealOptionKey const limit_cluster_size_percent( "cluster:limit_cluster_size_percent" );  }
namespace cluster { RealOptionKey const random_limit_cluster_size_percent( "cluster:random_limit_cluster_size_percent" );  }
namespace cluster { IntegerOptionKey const limit_clusters( "cluster:limit_clusters" );  }
namespace cluster { IntegerOptionKey const limit_total_structures( "cluster:limit_total_structures" );  }
namespace cluster { IntegerOptionKey const max_total_cluster( "cluster:max_total_cluster" );  }
namespace cluster { BooleanOptionKey const gdtmm( "cluster:gdtmm" );  }
namespace cluster { BooleanOptionKey const skip_align( "cluster:skip_align" );  }
namespace cluster { IntegerOptionKey const max_rms_matrix( "cluster:max_rms_matrix" );  }
namespace cluster { BooleanOptionKey const rna_P( "cluster:rna_P" );  }
namespace cluster { BooleanOptionKey const sort_groups_by_energy( "cluster:sort_groups_by_energy" );  }
namespace cluster { BooleanOptionKey const sort_groups_by_size( "cluster:sort_groups_by_size" );  }
namespace cluster { BooleanOptionKey const remove_singletons( "cluster:remove_singletons" );  }
namespace cluster { BooleanOptionKey const export_only_low( "cluster:export_only_low" );  }
namespace cluster { BooleanOptionKey const remove_highest_energy_member( "cluster:remove_highest_energy_member" );  }
namespace cluster { BooleanOptionKey const idealize_final_structures( "cluster:idealize_final_structures" );  }
namespace cluster { IntegerOptionKey const limit_dist_matrix( "cluster:limit_dist_matrix" );  }
namespace cluster { BooleanOptionKey const make_ensemble_cst( "cluster:make_ensemble_cst" );  }
namespace cluster { BooleanOptionKey const hotspot_hash( "cluster:hotspot_hash" );  }
namespace cluster { BooleanOptionKey const loops( "cluster:loops" );  }
namespace cluster { RealOptionKey const population_weight( "cluster:population_weight" );  }
namespace cluster { StringOptionKey const template_scores( "cluster:template_scores" );  }
namespace cluster { IntegerOptionKey const K_level( "cluster:K_level" );  }
namespace cluster { RealVectorOptionKey const K_radius( "cluster:K_radius" );  }
namespace cluster { IntegerVectorOptionKey const K_n_cluster( "cluster:K_n_cluster" );  }
namespace cluster { StringVectorOptionKey const K_style( "cluster:K_style" );  }
namespace cluster { IntegerOptionKey const K_n_sub( "cluster:K_n_sub" );  }
namespace cluster { IntegerOptionKey const K_deque_size( "cluster:K_deque_size" );  }
namespace cluster { IntegerOptionKey const K_deque_level( "cluster:K_deque_level" );  }
namespace cluster { BooleanOptionKey const K_redundant( "cluster:K_redundant" );  }
namespace cluster { BooleanOptionKey const K_not_fit_xyz( "cluster:K_not_fit_xyz" );  }
namespace cluster { BooleanOptionKey const K_save_headers( "cluster:K_save_headers" );  }
namespace cluster { RealOptionKey const score_diff_cut( "cluster:score_diff_cut" );  }
namespace cluster { BooleanOptionKey const auto_tune( "cluster:auto_tune" );  }
namespace cluster { BooleanOptionKey const write_centers( "cluster:write_centers" );  }
namespace cm { BooleanOptionKey const cm( "cm" );  }
namespace cm { namespace sanitize { BooleanOptionKey const sanitize( "cm:sanitize" );  } }
namespace cm { namespace sanitize { IntegerOptionKey const num_fragments( "cm:sanitize:num_fragments" );  } }
namespace cm { namespace sanitize { RealOptionKey const cst_weight_pair( "cm:sanitize:cst_weight_pair" );  } }
namespace cm { namespace sanitize { RealOptionKey const cst_weight_coord( "cm:sanitize:cst_weight_coord" );  } }
namespace cm { BooleanOptionKey const start_models_only( "cm:start_models_only" );  }
namespace cm { StringOptionKey const aln_format( "cm:aln_format" );  }
namespace cm { BooleanOptionKey const recover_side_chains( "cm:recover_side_chains" );  }
namespace cm { FileVectorOptionKey const steal_extra_residues( "cm:steal_extra_residues" );  }
namespace cm { StringOptionKey const loop_mover( "cm:loop_mover" );  }
namespace cm { IntegerOptionKey const loop_close_level( "cm:loop_close_level" );  }
namespace cm { IntegerOptionKey const min_loop_size( "cm:min_loop_size" );  }
namespace cm { IntegerOptionKey const max_loop_rebuild( "cm:max_loop_rebuild" );  }
namespace cm { RealOptionKey const loop_rebuild_filter( "cm:loop_rebuild_filter" );  }
namespace cm { RealOptionKey const aln_length_filter_quantile( "cm:aln_length_filter_quantile" );  }
namespace cm { IntegerOptionKey const aln_length_filter( "cm:aln_length_filter" );  }
namespace cm { StringVectorOptionKey const seq_score( "cm:seq_score" );  }
namespace cm { StringOptionKey const aligner( "cm:aligner" );  }
namespace cm { RealOptionKey const min_gap_open( "cm:min_gap_open" );  }
namespace cm { RealOptionKey const max_gap_open( "cm:max_gap_open" );  }
namespace cm { RealOptionKey const min_gap_extend( "cm:min_gap_extend" );  }
namespace cm { RealOptionKey const max_gap_extend( "cm:max_gap_extend" );  }
namespace cm { IntegerOptionKey const nn( "cm:nn" );  }
namespace cm { FileVectorOptionKey const ev_map( "cm:ev_map" );  }
namespace cm { FileVectorOptionKey const hh_map( "cm:hh_map" );  }
namespace cm { namespace hybridize { BooleanOptionKey const hybridize( "cm:hybridize" );  } }
namespace cm { namespace hybridize { IntegerVectorOptionKey const starting_template( "cm:hybridize:starting_template" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const realign_domains( "cm:hybridize:realign_domains" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const realign_domains_stage2( "cm:hybridize:realign_domains_stage2" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const add_non_init_chunks( "cm:hybridize:add_non_init_chunks" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage1_increase_cycles( "cm:hybridize:stage1_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2_increase_cycles( "cm:hybridize:stage2_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2min_increase_cycles( "cm:hybridize:stage2min_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage1_probability( "cm:hybridize:stage1_probability" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const skip_stage2( "cm:hybridize:skip_stage2" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const no_global_frame( "cm:hybridize:no_global_frame" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const linmin_only( "cm:hybridize:linmin_only" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const relax( "cm:hybridize:relax" );  } }
namespace cm { namespace hybridize { RealOptionKey const frag_weight_aligned( "cm:hybridize:frag_weight_aligned" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const max_registry_shift( "cm:hybridize:max_registry_shift" );  } }
namespace cm { namespace hybridize { RealOptionKey const frag_1mer_insertion_weight( "cm:hybridize:frag_1mer_insertion_weight" );  } }
namespace cm { namespace hybridize { RealOptionKey const small_frag_insertion_weight( "cm:hybridize:small_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { RealOptionKey const big_frag_insertion_weight( "cm:hybridize:big_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const auto_frag_insertion_weight( "cm:hybridize:auto_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_1_cycles( "cm:hybridize:stage1_1_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_2_cycles( "cm:hybridize:stage1_2_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_3_cycles( "cm:hybridize:stage1_3_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_4_cycles( "cm:hybridize:stage1_4_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2_temperature( "cm:hybridize:stage2_temperature" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage1_4_cenrot_score( "cm:hybridize:stage1_4_cenrot_score" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const include_loop_ss_chunks( "cm:hybridize:include_loop_ss_chunks" );  } }
namespace contactMap { BooleanOptionKey const contactMap( "contactMap" );  }
namespace contactMap { StringOptionKey const prefix( "contactMap:prefix" );  }
namespace contactMap { RealOptionKey const distance_cutoff( "contactMap:distance_cutoff" );  }
namespace contactMap { StringOptionKey const region_def( "contactMap:region_def" );  }
namespace contactMap { BooleanOptionKey const row_format( "contactMap:row_format" );  }
namespace contactMap { BooleanOptionKey const distance_matrix( "contactMap:distance_matrix" );  }
namespace coupled_moves { BooleanOptionKey const coupled_moves( "coupled_moves" );  }
namespace coupled_moves { IntegerOptionKey const ntrials( "coupled_moves:ntrials" );  }
namespace coupled_moves { IntegerOptionKey const number_ligands( "coupled_moves:number_ligands" );  }
namespace coupled_moves { RealOptionKey const mc_kt( "coupled_moves:mc_kt" );  }
namespace coupled_moves { RealOptionKey const boltzmann_kt( "coupled_moves:boltzmann_kt" );  }
namespace coupled_moves { RealOptionKey const mm_bend_weight( "coupled_moves:mm_bend_weight" );  }
namespace coupled_moves { BooleanOptionKey const trajectory( "coupled_moves:trajectory" );  }
namespace coupled_moves { BooleanOptionKey const trajectory_gz( "coupled_moves:trajectory_gz" );  }
namespace coupled_moves { IntegerOptionKey const trajectory_stride( "coupled_moves:trajectory_stride" );  }
namespace coupled_moves { StringOptionKey const trajectory_file( "coupled_moves:trajectory_file" );  }
namespace coupled_moves { StringOptionKey const output_fasta( "coupled_moves:output_fasta" );  }
namespace coupled_moves { StringOptionKey const output_stats( "coupled_moves:output_stats" );  }
namespace coupled_moves { BooleanOptionKey const ligand_mode( "coupled_moves:ligand_mode" );  }
namespace coupled_moves { BooleanOptionKey const initial_repack( "coupled_moves:initial_repack" );  }
namespace coupled_moves { BooleanOptionKey const min_pack( "coupled_moves:min_pack" );  }
namespace coupled_moves { BooleanOptionKey const save_sequences( "coupled_moves:save_sequences" );  }
namespace coupled_moves { BooleanOptionKey const save_structures( "coupled_moves:save_structures" );  }
namespace coupled_moves { RealOptionKey const ligand_prob( "coupled_moves:ligand_prob" );  }
namespace coupled_moves { BooleanOptionKey const fix_backbone( "coupled_moves:fix_backbone" );  }
namespace coupled_moves { BooleanOptionKey const uniform_backrub( "coupled_moves:uniform_backrub" );  }
namespace coupled_moves { BooleanOptionKey const bias_sampling( "coupled_moves:bias_sampling" );  }
namespace coupled_moves { BooleanOptionKey const bump_check( "coupled_moves:bump_check" );  }
namespace coupled_moves { RealOptionKey const ligand_weight( "coupled_moves:ligand_weight" );  }
namespace coupled_moves { StringOptionKey const output_prefix( "coupled_moves:output_prefix" );  }
namespace cp { BooleanOptionKey const cp( "cp" );  }
namespace cp { RealOptionKey const cutoff( "cp:cutoff" );  }
namespace cp { StringOptionKey const relax_sfxn( "cp:relax_sfxn" );  }
namespace cp { StringOptionKey const pack_sfxn( "cp:pack_sfxn" );  }
namespace cp { StringOptionKey const minimizer_score_fxn( "cp:minimizer_score_fxn" );  }
namespace cp { StringOptionKey const output( "cp:output" );  }
namespace cp { IntegerOptionKey const ncycles( "cp:ncycles" );  }
namespace cp { IntegerOptionKey const max_failures( "cp:max_failures" );  }
namespace cp { BooleanOptionKey const print_reports( "cp:print_reports" );  }
namespace cp { StringOptionKey const vipReportFile( "cp:vipReportFile" );  }
namespace cp { StringOptionKey const exclude_file( "cp:exclude_file" );  }
namespace cp { StringOptionKey const relax_mover( "cp:relax_mover" );  }
namespace cp { BooleanOptionKey const skip_relax( "cp:skip_relax" );  }
namespace cp { BooleanOptionKey const local_relax( "cp:local_relax" );  }
namespace cp { BooleanOptionKey const print_intermediate_pdbs( "cp:print_intermediate_pdbs" );  }
namespace cp { BooleanOptionKey const use_unrelaxed_starting_points( "cp:use_unrelaxed_starting_points" );  }
namespace cp { BooleanOptionKey const easy_vip_acceptance( "cp:easy_vip_acceptance" );  }
namespace cryst { BooleanOptionKey const cryst( "cryst" );  }
namespace cryst { StringOptionKey const mtzfile( "cryst:mtzfile" );  }
namespace cryst { BooleanOptionKey const crystal_refine( "cryst:crystal_refine" );  }
namespace cryst { BooleanOptionKey const refinable_lattice( "cryst:refinable_lattice" );  }
namespace cryst { RealOptionKey const interaction_shell( "cryst:interaction_shell" );  }
namespace csa { BooleanOptionKey const csa( "csa" );  }
namespace csa { BooleanOptionKey const useZ( "csa:useZ" );  }
namespace cutoutdomain { BooleanOptionKey const cutoutdomain( "cutoutdomain" );  }
namespace cutoutdomain { IntegerOptionKey const start( "cutoutdomain:start" );  }
namespace cutoutdomain { IntegerOptionKey const end( "cutoutdomain:end" );  }
namespace cyclization { BooleanOptionKey const cyclization( "cyclization" );  }
namespace cyclization { IntegerVectorOptionKey const chains_to_cyclize( "cyclization:chains_to_cyclize" );  }
namespace cyclization { IntegerOptionKey const num_min_rebuild( "cyclization:num_min_rebuild" );  }
namespace cyclization { BooleanOptionKey const add_constraints( "cyclization:add_constraints" );  }
namespace cyclic_peptide { BooleanOptionKey const cyclic_peptide( "cyclic_peptide" );  }
namespace cyclic_peptide { StringOptionKey const rand_checkpoint_file( "cyclic_peptide:rand_checkpoint_file" );  }
namespace cyclic_peptide { StringOptionKey const checkpoint_file( "cyclic_peptide:checkpoint_file" );  }
namespace cyclic_peptide { StringOptionKey const checkpoint_job_identifier( "cyclic_peptide:checkpoint_job_identifier" );  }
namespace cyclic_peptide { StringOptionKey const default_rama_sampling_table( "cyclic_peptide:default_rama_sampling_table" );  }
namespace cyclic_peptide { StringVectorOptionKey const rama_sampling_table_by_res( "cyclic_peptide:rama_sampling_table_by_res" );  }
namespace cyclic_peptide { StringOptionKey const sequence_file( "cyclic_peptide:sequence_file" );  }
namespace cyclic_peptide { IntegerOptionKey const genkic_closure_attempts( "cyclic_peptide:genkic_closure_attempts" );  }
namespace cyclic_peptide { IntegerOptionKey const genkic_min_solution_count( "cyclic_peptide:genkic_min_solution_count" );  }
namespace cyclic_peptide { BooleanOptionKey const cyclic_permutations( "cyclic_peptide:cyclic_permutations" );  }
namespace cyclic_peptide { BooleanOptionKey const use_rama_filter( "cyclic_peptide:use_rama_filter" );  }
namespace cyclic_peptide { RealOptionKey const rama_cutoff( "cyclic_peptide:rama_cutoff" );  }
namespace cyclic_peptide { RealOptionKey const high_hbond_weight_multiplier( "cyclic_peptide:high_hbond_weight_multiplier" );  }
namespace cyclic_peptide { RealOptionKey const min_genkic_hbonds( "cyclic_peptide:min_genkic_hbonds" );  }
namespace cyclic_peptide { RealOptionKey const min_final_hbonds( "cyclic_peptide:min_final_hbonds" );  }
namespace cyclic_peptide { RealOptionKey const total_energy_cutoff( "cyclic_peptide:total_energy_cutoff" );  }
namespace cyclic_peptide { RealOptionKey const hbond_energy_cutoff( "cyclic_peptide:hbond_energy_cutoff" );  }
namespace cyclic_peptide { BooleanOptionKey const do_not_count_adjacent_res_hbonds( "cyclic_peptide:do_not_count_adjacent_res_hbonds" );  }
namespace cyclic_peptide { IntegerOptionKey const fast_relax_rounds( "cyclic_peptide:fast_relax_rounds" );  }
namespace cyclic_peptide { BooleanOptionKey const count_sc_hbonds( "cyclic_peptide:count_sc_hbonds" );  }
namespace cyclic_peptide { BooleanOptionKey const require_disulfides( "cyclic_peptide:require_disulfides" );  }
namespace cyclic_peptide { RealOptionKey const disulf_cutoff_prerelax( "cyclic_peptide:disulf_cutoff_prerelax" );  }
