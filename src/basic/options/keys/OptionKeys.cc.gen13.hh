namespace remodel { namespace design { BooleanOptionKey const skip_partial( "remodel:design:skip_partial" );  } }
namespace remodel { namespace design { BooleanOptionKey const design_neighbors( "remodel:design:design_neighbors" );  } }
namespace remodel { namespace design { BooleanOptionKey const find_neighbors( "remodel:design:find_neighbors" );  } }
namespace remodel { namespace design { BooleanOptionKey const include_current( "remodel:design:include_current" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const RemodelLoopMover( "remodel:RemodelLoopMover" );  } }
namespace remodel { namespace RemodelLoopMover { RealOptionKey const max_linear_chainbreak( "remodel:RemodelLoopMover:max_linear_chainbreak" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const randomize_loops( "remodel:RemodelLoopMover:randomize_loops" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const use_loop_hash( "remodel:RemodelLoopMover:use_loop_hash" );  } }
namespace remodel { namespace RemodelLoopMover { FileOptionKey const set_segment( "remodel:RemodelLoopMover:set_segment" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const allowed_closure_attempts( "remodel:RemodelLoopMover:allowed_closure_attempts" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const loophash_cycles( "remodel:RemodelLoopMover:loophash_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const simultaneous_cycles( "remodel:RemodelLoopMover:simultaneous_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const independent_cycles( "remodel:RemodelLoopMover:independent_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const boost_closure_cycles( "remodel:RemodelLoopMover:boost_closure_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { RealOptionKey const threshold_for_boost_closure( "remodel:RemodelLoopMover:threshold_for_boost_closure" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const force_cutting_index( "remodel:RemodelLoopMover:force_cutting_index" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const force_cutting_N( "remodel:RemodelLoopMover:force_cutting_N" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const bypass_closure( "remodel:RemodelLoopMover:bypass_closure" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const cyclic_peptide( "remodel:RemodelLoopMover:cyclic_peptide" );  } }
namespace remodel { namespace RemodelLoopMover { RealOptionKey const temperature( "remodel:RemodelLoopMover:temperature" );  } }
namespace resample { BooleanOptionKey const resample( "resample" );  }
namespace resample { FileOptionKey const silent( "resample:silent" );  }
namespace resample { StringOptionKey const tag( "resample:tag" );  }
namespace resample { BooleanOptionKey const stage1( "resample:stage1" );  }
namespace resample { BooleanOptionKey const stage2( "resample:stage2" );  }
namespace resample { BooleanOptionKey const jumps( "resample:jumps" );  }
namespace resample { RealVectorOptionKey const min_max_start_seq_sep( "resample:min_max_start_seq_sep" );  }
namespace rescore { BooleanOptionKey const rescore( "rescore" );  }
namespace rescore { BooleanOptionKey const assign_ss( "rescore:assign_ss" );  }
namespace rescore { BooleanOptionKey const skip( "rescore:skip" );  }
namespace rescore { BooleanOptionKey const verbose( "rescore:verbose" );  }
namespace rna { BooleanOptionKey const rna( "rna" );  }
namespace rna { BooleanOptionKey const corrected_geo( "rna:corrected_geo" );  }
namespace rna { BooleanOptionKey const rna_prot_erraser( "rna:rna_prot_erraser" );  }
namespace rna { BooleanOptionKey const vary_geometry( "rna:vary_geometry" );  }
namespace rna { StringOptionKey const data_file( "rna:data_file" );  }
namespace rna { BooleanOptionKey const cut_at_rna_chainbreak( "rna:cut_at_rna_chainbreak" );  }
namespace rna { BooleanOptionKey const evaluate_base_pairs( "rna:evaluate_base_pairs" );  }
namespace rna { namespace farna { BooleanOptionKey const farna( "rna:farna" );  } }
namespace rna { namespace farna { IntegerOptionKey const cycles( "rna:farna:cycles" );  } }
namespace rna { namespace farna { IntegerOptionKey const rna_protein_docking_freq( "rna:farna:rna_protein_docking_freq" );  } }
namespace rna { namespace farna { IntegerOptionKey const rounds( "rna:farna:rounds" );  } }
namespace rna { namespace farna { RealOptionKey const temperature( "rna:farna:temperature" );  } }
namespace rna { namespace farna { BooleanOptionKey const minimize_rna( "rna:farna:minimize_rna" );  } }
namespace rna { namespace farna { StringVectorOptionKey const sequence( "rna:farna:sequence" );  } }
namespace rna { namespace farna { StringOptionKey const secstruct( "rna:farna:secstruct" );  } }
namespace rna { namespace farna { StringOptionKey const secstruct_general( "rna:farna:secstruct_general" );  } }
namespace rna { namespace farna { StringOptionKey const secstruct_file( "rna:farna:secstruct_file" );  } }
namespace rna { namespace farna { StringOptionKey const secstruct_general_file( "rna:farna:secstruct_general_file" );  } }
namespace rna { namespace farna { StringOptionKey const secstruct_legacy( "rna:farna:secstruct_legacy" );  } }
namespace rna { namespace farna { StringOptionKey const lores_scorefxn( "rna:farna:lores_scorefxn" );  } }
namespace rna { namespace farna { StringOptionKey const params_file( "rna:farna:params_file" );  } }
namespace rna { namespace farna { BooleanOptionKey const filter_lores_base_pairs( "rna:farna:filter_lores_base_pairs" );  } }
namespace rna { namespace farna { BooleanOptionKey const filter_lores_base_pairs_early( "rna:farna:filter_lores_base_pairs_early" );  } }
namespace rna { namespace farna { BooleanOptionKey const filter_chain_closure( "rna:farna:filter_chain_closure" );  } }
namespace rna { namespace farna { BooleanOptionKey const filter_chain_closure_halfway( "rna:farna:filter_chain_closure_halfway" );  } }
namespace rna { namespace farna { RealOptionKey const filter_chain_closure_distance( "rna:farna:filter_chain_closure_distance" );  } }
namespace rna { namespace farna { BooleanOptionKey const relax_rna( "rna:farna:relax_rna" );  } }
namespace rna { namespace farna { BooleanOptionKey const simple_relax( "rna:farna:simple_relax" );  } }
namespace rna { namespace farna { BooleanOptionKey const ignore_secstruct( "rna:farna:ignore_secstruct" );  } }
namespace rna { namespace farna { RealOptionKey const jump_change_frequency( "rna:farna:jump_change_frequency" );  } }
namespace rna { namespace farna { BooleanOptionKey const close_loops( "rna:farna:close_loops" );  } }
namespace rna { namespace farna { BooleanOptionKey const close_loops_after_each_move( "rna:farna:close_loops_after_each_move" );  } }
namespace rna { namespace farna { BooleanOptionKey const heat( "rna:farna:heat" );  } }
namespace rna { namespace farna { BooleanOptionKey const staged_constraints( "rna:farna:staged_constraints" );  } }
namespace rna { namespace farna { StringOptionKey const jump_library_file( "rna:farna:jump_library_file" );  } }
namespace rna { namespace farna { StringOptionKey const vall_torsions( "rna:farna:vall_torsions" );  } }
namespace rna { namespace farna { BooleanOptionKey const use_1jj2_torsions( "rna:farna:use_1jj2_torsions" );  } }
namespace rna { namespace farna { RealOptionKey const rna_lores_chainbreak_weight( "rna:farna:rna_lores_chainbreak_weight" );  } }
namespace rna { namespace farna { RealOptionKey const rna_lores_linear_chainbreak_weight( "rna:farna:rna_lores_linear_chainbreak_weight" );  } }
namespace rna { namespace farna { BooleanOptionKey const fixed_stems( "rna:farna:fixed_stems" );  } }
namespace rna { namespace farna { BooleanOptionKey const allow_bulge( "rna:farna:allow_bulge" );  } }
namespace rna { namespace farna { IntegerVectorOptionKey const allowed_bulge_res( "rna:farna:allowed_bulge_res" );  } }
namespace rna { namespace farna { BooleanOptionKey const allow_consecutive_bulges( "rna:farna:allow_consecutive_bulges" );  } }
namespace rna { namespace farna { BooleanOptionKey const move_first_rigid_body( "rna:farna:move_first_rigid_body" );  } }
namespace rna { namespace farna { BooleanOptionKey const root_at_first_rigid_body( "rna:farna:root_at_first_rigid_body" );  } }
namespace rna { namespace farna { RealOptionKey const suppress_bp_constraint( "rna:farna:suppress_bp_constraint" );  } }
namespace rna { namespace farna { BooleanOptionKey const output_filters( "rna:farna:output_filters" );  } }
namespace rna { namespace farna { BooleanOptionKey const autofilter( "rna:farna:autofilter" );  } }
namespace rna { namespace farna { BooleanOptionKey const no_filters( "rna:farna:no_filters" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const output_res_num( "rna:farna:output_res_num" );  } }
namespace rna { namespace farna { IntegerOptionKey const offset( "rna:farna:offset" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const input_silent_res( "rna:farna:input_silent_res" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const virtual_anchor( "rna:farna:virtual_anchor" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const obligate_pair( "rna:farna:obligate_pair" );  } }
namespace rna { namespace farna { StringVectorOptionKey const obligate_pair_explicit( "rna:farna:obligate_pair_explicit" );  } }
namespace rna { namespace farna { StringVectorOptionKey const chain_connection( "rna:farna:chain_connection" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const remove_pair( "rna:farna:remove_pair" );  } }
namespace rna { namespace farna { ResidueChainVectorOptionKey const remove_obligate_pair( "rna:farna:remove_obligate_pair" );  } }
namespace rna { namespace farna { StringOptionKey const refine_silent_file( "rna:farna:refine_silent_file" );  } }
namespace rna { namespace farna { BooleanOptionKey const refine_native( "rna:farna:refine_native" );  } }
namespace rna { namespace farna { BooleanOptionKey const bps_moves( "rna:farna:bps_moves" );  } }
namespace rna { namespace farna { BooleanOptionKey const disallow_bps_at_extra_min_res( "rna:farna:disallow_bps_at_extra_min_res" );  } }
namespace rna { namespace farna { BooleanOptionKey const allow_fragment_moves_in_bps( "rna:farna:allow_fragment_moves_in_bps" );  } }
namespace rna { namespace farna { IntegerOptionKey const frag_size( "rna:farna:frag_size" );  } }
namespace rna { namespace farna { BooleanOptionKey const VDW_rep_screen_include_sidechains( "rna:farna:VDW_rep_screen_include_sidechains" );  } }
namespace rna { namespace farna { BooleanOptionKey const gradual_constraints( "rna:farna:gradual_constraints" );  } }
namespace rna { namespace farna { RealOptionKey const grid_vdw_weight( "rna:farna:grid_vdw_weight" );  } }
namespace rna { namespace farna { StringOptionKey const tag( "rna:farna:tag" );  } }
namespace rna { namespace farna { StringOptionKey const working_native( "rna:farna:working_native" );  } }
namespace rna { namespace farna { BooleanOptionKey const use_legacy_setup( "rna:farna:use_legacy_setup" );  } }
namespace rna { namespace farna { BooleanOptionKey const cst_gap( "rna:farna:cst_gap" );  } }
namespace rna { namespace farna { BooleanOptionKey const convert_protein_CEN( "rna:farna:convert_protein_CEN" );  } }
namespace rna { namespace farna { BooleanOptionKey const rna_protein_docking( "rna:farna:rna_protein_docking" );  } }
namespace rna { namespace farna { IntegerVectorOptionKey const exclude_fragments( "rna:farna:exclude_fragments" );  } }
namespace rna { namespace farna { StringOptionKey const exclusion_match_type( "rna:farna:exclusion_match_type" );  } }
namespace rna { namespace farna { RealOptionKey const fragment_homology_rmsd( "rna:farna:fragment_homology_rmsd" );  } }
namespace rna { namespace farna { BooleanOptionKey const exclude_native_fragments( "rna:farna:exclude_native_fragments" );  } }
namespace rna { namespace farna { StringVectorOptionKey const exclude_fragment_files( "rna:farna:exclude_fragment_files" );  } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const out( "rna:farna:out" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const output_lores_silent_file( "rna:farna:out:output_lores_silent_file" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const dump( "rna:farna:out:dump" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const binary_output( "rna:farna:out:binary_output" );  } } }
namespace rna { namespace farna { namespace out { StringOptionKey const output_score_file( "rna:farna:out:output_score_file" );  } } }
namespace rna { namespace farna { namespace out { IntegerOptionKey const output_score_frequency( "rna:farna:out:output_score_frequency" );  } } }
namespace rna { namespace farna { namespace out { IntegerVectorOptionKey const output_jump_res( "rna:farna:out:output_jump_res" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const output_jump_o3p_to_o5p( "rna:farna:out:output_jump_o3p_to_o5p" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const output_rotation_vector( "rna:farna:out:output_rotation_vector" );  } } }
namespace rna { namespace farna { namespace out { RealVectorOptionKey const target_xyz( "rna:farna:out:target_xyz" );  } } }
namespace rna { namespace farna { namespace out { BooleanOptionKey const save_jump_histogram( "rna:farna:out:save_jump_histogram" );  } } }
namespace rna { namespace farna { namespace out { StringOptionKey const output_histogram_file( "rna:farna:out:output_histogram_file" );  } } }
namespace rna { namespace farna { namespace out { RealOptionKey const jump_histogram_boxsize( "rna:farna:out:jump_histogram_boxsize" );  } } }
namespace rna { namespace farna { namespace out { RealOptionKey const jump_histogram_binwidth( "rna:farna:out:jump_histogram_binwidth" );  } } }
namespace rna { namespace farna { namespace out { RealOptionKey const jump_histogram_binwidth_rotvector( "rna:farna:out:jump_histogram_binwidth_rotvector" );  } } }
namespace rna { namespace farna { namespace db { BooleanOptionKey const db( "rna:farna:db" );  } } }
namespace rna { namespace farna { namespace db { BooleanOptionKey const jump_database( "rna:farna:db:jump_database" );  } } }
namespace rna { namespace farna { namespace db { BooleanOptionKey const bps_database( "rna:farna:db:bps_database" );  } } }
namespace rna { namespace farna { namespace erraser { BooleanOptionKey const erraser( "rna:farna:erraser" );  } } }
namespace rna { namespace farna { namespace erraser { BooleanOptionKey const constrain_P( "rna:farna:erraser:constrain_P" );  } } }
namespace rna { namespace farna { namespace erraser { IntegerVectorOptionKey const fixed_res( "rna:farna:erraser:fixed_res" );  } } }
namespace rna { namespace farna { namespace erraser { BooleanOptionKey const ready_set_only( "rna:farna:erraser:ready_set_only" );  } } }
namespace rna { namespace farna { namespace erraser { BooleanOptionKey const skip_minimize( "rna:farna:erraser:skip_minimize" );  } } }
namespace rna { namespace farna { namespace erraser { BooleanOptionKey const attempt_pyrimidine_flip( "rna:farna:erraser:attempt_pyrimidine_flip" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const minimize( "rna:farna:minimize" );  } } }
namespace rna { namespace farna { namespace minimize { IntegerOptionKey const minimize_rounds( "rna:farna:minimize:minimize_rounds" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const skip_coord_constraints( "rna:farna:minimize:skip_coord_constraints" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const skip_o2prime_trials( "rna:farna:minimize:skip_o2prime_trials" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const deriv_check( "rna:farna:minimize:deriv_check" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const minimizer_use_coordinate_constraints( "rna:farna:minimize:minimizer_use_coordinate_constraints" );  } } }
namespace rna { namespace farna { namespace minimize { StringOptionKey const min_type( "rna:farna:minimize:min_type" );  } } }
namespace rna { namespace farna { namespace minimize { BooleanOptionKey const minimize_bps( "rna:farna:minimize:minimize_bps" );  } } }
namespace rna { namespace farna { namespace minimize { ResidueChainVectorOptionKey const extra_minimize_res( "rna:farna:minimize:extra_minimize_res" );  } } }
namespace rna { namespace farna { namespace minimize { ResidueChainVectorOptionKey const extra_minimize_chi_res( "rna:farna:minimize:extra_minimize_chi_res" );  } } }
namespace rotamerdump { BooleanOptionKey const rotamerdump( "rotamerdump" );  }
namespace rotamerdump { BooleanOptionKey const xyz( "rotamerdump:xyz" );  }
namespace rotamerdump { BooleanOptionKey const one_body( "rotamerdump:one_body" );  }
namespace rotamerdump { BooleanOptionKey const two_body( "rotamerdump:two_body" );  }
namespace rotamerdump { BooleanOptionKey const annealer( "rotamerdump:annealer" );  }
namespace sample_around { BooleanOptionKey const sample_around( "sample_around" );  }
namespace sample_around { RealOptionKey const alpha_increment( "sample_around:alpha_increment" );  }
namespace sample_around { RealOptionKey const cosbeta_increment( "sample_around:cosbeta_increment" );  }
namespace sample_around { RealOptionKey const gamma_increment( "sample_around:gamma_increment" );  }
namespace sicdock { BooleanOptionKey const sicdock( "sicdock" );  }
namespace sicdock { RealOptionKey const clash_dis( "sicdock:clash_dis" );  }
namespace sicdock { RealOptionKey const contact_dis( "sicdock:contact_dis" );  }
namespace sicdock { RealOptionKey const hash_2D_vs_3D( "sicdock:hash_2D_vs_3D" );  }
namespace sewing { BooleanOptionKey const sewing( "sewing" );  }
namespace sewing { FileOptionKey const model_file_name( "sewing:model_file_name" );  }
namespace sewing { FileOptionKey const score_file_name( "sewing:score_file_name" );  }
namespace sewing { FileOptionKey const new_model_file_name( "sewing:new_model_file_name" );  }
namespace sewing { StringOptionKey const remove_any_dssp( "sewing:remove_any_dssp" );  }
namespace sewing { StringOptionKey const remove_all_dssp( "sewing:remove_all_dssp" );  }
namespace sewing { IntegerOptionKey const min_helix_length( "sewing:min_helix_length" );  }
namespace sewing { IntegerOptionKey const max_helix_length( "sewing:max_helix_length" );  }
namespace sewing { IntegerOptionKey const min_loop_length( "sewing:min_loop_length" );  }
namespace sewing { IntegerOptionKey const max_loop_length( "sewing:max_loop_length" );  }
namespace sewing { IntegerOptionKey const min_strand_length( "sewing:min_strand_length" );  }
namespace sewing { IntegerOptionKey const max_strand_length( "sewing:max_strand_length" );  }
namespace sewing { BooleanOptionKey const leave_models_by_ss_num( "sewing:leave_models_by_ss_num" );  }
namespace sewing { IntegerOptionKey const model_should_have_this_num_of_ss( "sewing:model_should_have_this_num_of_ss" );  }
namespace sewing { BooleanOptionKey const model_should_have_at_least_1_E_at_terminal_segment( "sewing:model_should_have_at_least_1_E_at_terminal_segment" );  }
namespace sewing { BooleanOptionKey const model_should_have_at_least_1_E( "sewing:model_should_have_at_least_1_E" );  }
namespace sewing { BooleanOptionKey const model_should_have_at_least_1_H( "sewing:model_should_have_at_least_1_H" );  }
namespace sewing { BooleanOptionKey const leave_models_with_E_terminal_ss( "sewing:leave_models_with_E_terminal_ss" );  }
namespace sewing { BooleanOptionKey const leave_models_with_H_terminal_ss( "sewing:leave_models_with_H_terminal_ss" );  }
namespace sewing { BooleanOptionKey const leave_antiparallel_way_H_bonded_models_by_terminal_strands( "sewing:leave_antiparallel_way_H_bonded_models_by_terminal_strands" );  }
namespace sewing { BooleanOptionKey const leave_parallel_way_H_bonded_models_by_terminal_strands( "sewing:leave_parallel_way_H_bonded_models_by_terminal_strands" );  }
namespace sewing { BooleanOptionKey const leave_certain_model_ids( "sewing:leave_certain_model_ids" );  }
namespace sewing { StringOptionKey const leave_these_model_ids( "sewing:leave_these_model_ids" );  }
namespace sewing { IntegerOptionKey const box_length( "sewing:box_length" );  }
namespace sewing { StringOptionKey const mode( "sewing:mode" );  }
namespace sewing { BooleanOptionKey const disregard_num_segment_matches( "sewing:disregard_num_segment_matches" );  }
namespace sewing { BooleanOptionKey const do_not_remove_connection_inconsistencies( "sewing:do_not_remove_connection_inconsistencies" );  }
namespace sewing { BooleanOptionKey const score_between_opposite_terminal_segments( "sewing:score_between_opposite_terminal_segments" );  }
namespace sewing { IntegerOptionKey const num_models_to_dump( "sewing:num_models_to_dump" );  }
namespace sewing { IntegerVectorOptionKey const models_to_dump( "sewing:models_to_dump" );  }
namespace sewing { IntegerOptionKey const min_hash_score( "sewing:min_hash_score" );  }
namespace sewing { IntegerOptionKey const max_clash_score( "sewing:max_clash_score" );  }
namespace sewing { IntegerOptionKey const num_segments_to_match( "sewing:num_segments_to_match" );  }
namespace sewing { IntegerVectorOptionKey const match_segments( "sewing:match_segments" );  }
namespace sewing { IntegerOptionKey const max_models( "sewing:max_models" );  }
namespace sewing { IntegerOptionKey const starting_model( "sewing:starting_model" );  }
namespace sewing { IntegerOptionKey const num_procs( "sewing:num_procs" );  }
namespace sewing { IntegerOptionKey const rank( "sewing:rank" );  }
namespace sewing { BooleanOptionKey const hash_tag_only_terminal_Es( "sewing:hash_tag_only_terminal_Es" );  }
namespace sewing { StringOptionKey const assembly_type( "sewing:assembly_type" );  }
namespace sewing { IntegerOptionKey const num_edges_to_follow( "sewing:num_edges_to_follow" );  }
namespace sewing { RealOptionKey const base_native_bonus( "sewing:base_native_bonus" );  }
namespace sewing { IntegerOptionKey const neighbor_cutoff( "sewing:neighbor_cutoff" );  }
namespace sewing { BooleanOptionKey const dump_pdbs( "sewing:dump_pdbs" );  }
namespace sewing { BooleanOptionKey const skip_refinement( "sewing:skip_refinement" );  }
namespace sewing { BooleanOptionKey const skip_filters( "sewing:skip_filters" );  }
namespace sewing { RealOptionKey const min_motif_score( "sewing:min_motif_score" );  }
namespace sewing { BooleanOptionKey const may_add_already_added_model( "sewing:may_add_already_added_model" );  }
namespace sewing { RealOptionKey const offset_bump_dsq( "sewing:offset_bump_dsq" );  }
namespace sewing { IntegerOptionKey const num_repeats( "sewing:num_repeats" );  }
namespace sewing { BooleanOptionKey const repeat( "sewing:repeat" );  }
namespace sewing { IntegerVectorOptionKey const pose_segment_starts( "sewing:pose_segment_starts" );  }
namespace sewing { IntegerVectorOptionKey const pose_segment_ends( "sewing:pose_segment_ends" );  }
namespace sewing { BooleanOptionKey const keep_source_segments( "sewing:keep_source_segments" );  }
namespace sewing { FileOptionKey const partner_pdb( "sewing:partner_pdb" );  }
namespace sewing { IntegerVectorOptionKey const keep_model_residues( "sewing:keep_model_residues" );  }
namespace sewing { IntegerOptionKey const min_lh_fragments( "sewing:min_lh_fragments" );  }
namespace sewing { BooleanOptionKey const skip_loop_generation( "sewing:skip_loop_generation" );  }
namespace sewing { IntegerOptionKey const max_ss_num( "sewing:max_ss_num" );  }
namespace sewing { BooleanOptionKey const dump_every_model( "sewing:dump_every_model" );  }
namespace SSrbrelax { BooleanOptionKey const SSrbrelax( "SSrbrelax" );  }
namespace SSrbrelax { FileOptionKey const rb_file( "SSrbrelax:rb_file" );  }
namespace SSrbrelax { FileOptionKey const rb_param_file( "SSrbrelax:rb_param_file" );  }
namespace SSrbrelax { IntegerVectorOptionKey const frag_sizes( "SSrbrelax:frag_sizes" );  }
namespace SSrbrelax { FileVectorOptionKey const frag_files( "SSrbrelax:frag_files" );  }
namespace stepwise { BooleanOptionKey const stepwise( "stepwise" );  }
namespace stepwise { StringVectorOptionKey const s1( "stepwise:s1" );  }
namespace stepwise { StringVectorOptionKey const s2( "stepwise:s2" );  }
namespace stepwise { StringVectorOptionKey const silent1( "stepwise:silent1" );  }
namespace stepwise { StringVectorOptionKey const silent2( "stepwise:silent2" );  }
namespace stepwise { StringVectorOptionKey const tags1( "stepwise:tags1" );  }
namespace stepwise { StringVectorOptionKey const tags2( "stepwise:tags2" );  }
namespace stepwise { IntegerVectorOptionKey const slice_res1( "stepwise:slice_res1" );  }
namespace stepwise { IntegerVectorOptionKey const slice_res2( "stepwise:slice_res2" );  }
namespace stepwise { IntegerVectorOptionKey const input_res1( "stepwise:input_res1" );  }
namespace stepwise { IntegerVectorOptionKey const input_res2( "stepwise:input_res2" );  }
namespace stepwise { BooleanOptionKey const backbone_only1( "stepwise:backbone_only1" );  }
namespace stepwise { BooleanOptionKey const backbone_only2( "stepwise:backbone_only2" );  }
namespace stepwise { IntegerVectorOptionKey const fixed_res( "stepwise:fixed_res" );  }
namespace stepwise { BooleanOptionKey const test_encapsulation( "stepwise:test_encapsulation" );  }
namespace stepwise { BooleanOptionKey const choose_random( "stepwise:choose_random" );  }
namespace stepwise { IntegerOptionKey const num_random_samples( "stepwise:num_random_samples" );  }
namespace stepwise { IntegerOptionKey const max_tries_multiplier_for_ccd( "stepwise:max_tries_multiplier_for_ccd" );  }
namespace stepwise { IntegerOptionKey const num_pose_minimize( "stepwise:num_pose_minimize" );  }
namespace stepwise { BooleanOptionKey const atr_rep_screen( "stepwise:atr_rep_screen" );  }
namespace stepwise { BooleanOptionKey const atr_rep_screen_for_docking( "stepwise:atr_rep_screen_for_docking" );  }
namespace stepwise { StringOptionKey const align_pdb( "stepwise:align_pdb" );  }
namespace stepwise { BooleanOptionKey const enumerate( "stepwise:enumerate" );  }
namespace stepwise { BooleanOptionKey const preminimize( "stepwise:preminimize" );  }
namespace stepwise { BooleanOptionKey const skip_preminimize( "stepwise:skip_preminimize" );  }
