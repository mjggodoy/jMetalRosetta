namespace pocket_grid { BooleanOptionKey const dump_connollySurface( "pocket_grid:dump_connollySurface" );  }
namespace pocket_grid { RealOptionKey const esp_buffer_dist( "pocket_grid:esp_buffer_dist" );  }
namespace pocket_grid { BooleanOptionKey const round_pocketGrid_center( "pocket_grid:round_pocketGrid_center" );  }
namespace gen_pharmacophore { BooleanOptionKey const gen_pharmacophore( "gen_pharmacophore" );  }
namespace gen_pharmacophore { RealOptionKey const clash_distance_cutoff( "gen_pharmacophore:clash_distance_cutoff" );  }
namespace gen_pharmacophore { IntegerOptionKey const ring_sasa_cutoff( "gen_pharmacophore:ring_sasa_cutoff" );  }
namespace gen_pharmacophore { IntegerOptionKey const min_num_ring( "gen_pharmacophore:min_num_ring" );  }
namespace gen_pharmacophore { RealOptionKey const ring_ring_dist_cutoff( "gen_pharmacophore:ring_ring_dist_cutoff" );  }
namespace gen_pharmacophore { RealOptionKey const ring_atm_dist_cutoff( "gen_pharmacophore:ring_atm_dist_cutoff" );  }
namespace fingerprint { BooleanOptionKey const fingerprint( "fingerprint" );  }
namespace fingerprint { BooleanOptionKey const print_eggshell( "fingerprint:print_eggshell" );  }
namespace fingerprint { RealOptionKey const atom_radius_scale( "fingerprint:atom_radius_scale" );  }
namespace fingerprint { RealOptionKey const atom_radius_buffer( "fingerprint:atom_radius_buffer" );  }
namespace fingerprint { RealOptionKey const packing_weight( "fingerprint:packing_weight" );  }
namespace fingerprint { RealOptionKey const dist_cut_off( "fingerprint:dist_cut_off" );  }
namespace fingerprint { BooleanOptionKey const include_extrashell_to_set_origin( "fingerprint:include_extrashell_to_set_origin" );  }
namespace fingerprint { BooleanOptionKey const include_hydrogens( "fingerprint:include_hydrogens" );  }
namespace fingerprint { BooleanOptionKey const multiple_origin( "fingerprint:multiple_origin" );  }
namespace fingerprint { BooleanOptionKey const use_DARC_gpu( "fingerprint:use_DARC_gpu" );  }
namespace fingerprint { BooleanOptionKey const square_score( "fingerprint:square_score" );  }
namespace fingerprint { BooleanOptionKey const darc_components( "fingerprint:darc_components" );  }
namespace fingerprint { IntegerOptionKey const set_origin( "fingerprint:set_origin" );  }
namespace fingerprint { StringOptionKey const origin_res_num( "fingerprint:origin_res_num" );  }
namespace fingerprint { BooleanOptionKey const add_esp( "fingerprint:add_esp" );  }
namespace fingerprint { BooleanOptionKey const darc_shape_only( "fingerprint:darc_shape_only" );  }
namespace fingerprint { BooleanOptionKey const darc_elsts_only( "fingerprint:darc_elsts_only" );  }
namespace fingerprint { RealOptionKey const esp_weight( "fingerprint:esp_weight" );  }
namespace fingerprint { RealOptionKey const esp_protein_wt( "fingerprint:esp_protein_wt" );  }
namespace fingerprint { RealOptionKey const esp_surface_wt( "fingerprint:esp_surface_wt" );  }
namespace fingerprint { BooleanOptionKey const delphi_grid( "fingerprint:delphi_grid" );  }
namespace fingerprint { RealOptionKey const cap_e_potential( "fingerprint:cap_e_potential" );  }
namespace fingerprint { BooleanOptionKey const return_zero_darc_score( "fingerprint:return_zero_darc_score" );  }
namespace fingerprint { BooleanOptionKey const set_surface_esp_to_zero( "fingerprint:set_surface_esp_to_zero" );  }
namespace fingerprint { BooleanOptionKey const set_protein_esp_to_zero( "fingerprint:set_protein_esp_to_zero" );  }
namespace fingerprint { StringOptionKey const inp_lig( "fingerprint:inp_lig" );  }
namespace fingerprint { StringOptionKey const ref_lig( "fingerprint:ref_lig" );  }
namespace ProQ { BooleanOptionKey const ProQ( "ProQ" );  }
namespace ProQ { IntegerOptionKey const svmmodel( "ProQ:svmmodel" );  }
namespace ProQ { StringOptionKey const basename( "ProQ:basename" );  }
namespace ProQ { BooleanOptionKey const membrane( "ProQ:membrane" );  }
namespace ProQ { BooleanOptionKey const prof_bug( "ProQ:prof_bug" );  }
namespace ProQ { BooleanOptionKey const output_feature_vector( "ProQ:output_feature_vector" );  }
namespace ProQ { BooleanOptionKey const output_local_prediction( "ProQ:output_local_prediction" );  }
namespace ProQ { StringOptionKey const prefix( "ProQ:prefix" );  }
namespace ProQ { BooleanOptionKey const use_gzip( "ProQ:use_gzip" );  }
namespace ProQ { RealOptionKey const normalize( "ProQ:normalize" );  }
namespace qsar { BooleanOptionKey const qsar( "qsar" );  }
namespace qsar { StringOptionKey const weights( "qsar:weights" );  }
namespace qsar { StringOptionKey const grid_dir( "qsar:grid_dir" );  }
namespace qsar { IntegerOptionKey const max_grid_cache_size( "qsar:max_grid_cache_size" );  }
namespace rbe { BooleanOptionKey const rbe( "rbe" );  }
namespace rbe { StringOptionKey const server_url( "rbe:server_url" );  }
namespace rbe { StringOptionKey const server_port( "rbe:server_port" );  }
namespace rbe { RealOptionKey const poll_frequency( "rbe:poll_frequency" );  }
namespace RBSegmentRelax { BooleanOptionKey const RBSegmentRelax( "RBSegmentRelax" );  }
namespace RBSegmentRelax { FileOptionKey const input_pdb( "RBSegmentRelax:input_pdb" );  }
namespace RBSegmentRelax { FileOptionKey const rb_file( "RBSegmentRelax:rb_file" );  }
namespace RBSegmentRelax { IntegerOptionKey const nrbmoves( "RBSegmentRelax:nrbmoves" );  }
namespace RBSegmentRelax { IntegerOptionKey const nrboutercycles( "RBSegmentRelax:nrboutercycles" );  }
namespace RBSegmentRelax { StringOptionKey const rb_scorefxn( "RBSegmentRelax:rb_scorefxn" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_fragment_moves( "RBSegmentRelax:skip_fragment_moves" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_seqshift_moves( "RBSegmentRelax:skip_seqshift_moves" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_rb_moves( "RBSegmentRelax:skip_rb_moves" );  }
namespace RBSegmentRelax { RealVectorOptionKey const helical_movement_params( "RBSegmentRelax:helical_movement_params" );  }
namespace RBSegmentRelax { RealVectorOptionKey const strand_movement_params( "RBSegmentRelax:strand_movement_params" );  }
namespace RBSegmentRelax { RealVectorOptionKey const default_movement_params( "RBSegmentRelax:default_movement_params" );  }
namespace RBSegmentRelax { IntegerOptionKey const cst_seqwidth( "RBSegmentRelax:cst_seqwidth" );  }
namespace rdc { BooleanOptionKey const rdc( "rdc" );  }
namespace rdc { BooleanOptionKey const correct_NH_length( "rdc:correct_NH_length" );  }
namespace rdc { BooleanOptionKey const reduced_couplings( "rdc:reduced_couplings" );  }
namespace rdc { FileOptionKey const weights( "rdc:weights" );  }
namespace rdc { RealOptionKey const iterate_weights( "rdc:iterate_weights" );  }
namespace rdc { FileOptionKey const segment_file( "rdc:segment_file" );  }
namespace rdc { StringOptionKey const segment_scoring_mode( "rdc:segment_scoring_mode" );  }
namespace rdc { RealOptionKey const total_weight( "rdc:total_weight" );  }
namespace rdc { RealOptionKey const tensor_weight( "rdc:tensor_weight" );  }
namespace rdc { FileOptionKey const print_rdc_values( "rdc:print_rdc_values" );  }
namespace rdc { RealOptionKey const iterate_tol( "rdc:iterate_tol" );  }
namespace rdc { BooleanOptionKey const iterate_reset( "rdc:iterate_reset" );  }
namespace rdc { FileOptionKey const dump_weight_trajectory( "rdc:dump_weight_trajectory" );  }
namespace rdc { RealVectorOptionKey const fix_normAzz( "rdc:fix_normAzz" );  }
namespace rdc { StringOptionKey const fit_method( "rdc:fit_method" );  }
namespace rdc { RealVectorOptionKey const fixDa( "rdc:fixDa" );  }
namespace rdc { RealVectorOptionKey const fixR( "rdc:fixR" );  }
namespace rdc { IntegerOptionKey const nlsrepeat( "rdc:nlsrepeat" );  }
namespace relax { BooleanOptionKey const relax( "relax" );  }
namespace relax { BooleanOptionKey const fast( "relax:fast" );  }
namespace relax { BooleanOptionKey const thorough( "relax:thorough" );  }
namespace relax { BooleanOptionKey const centroid_mode( "relax:centroid_mode" );  }
namespace relax { IntegerOptionKey const default_repeats( "relax:default_repeats" );  }
namespace relax { BooleanOptionKey const dualspace( "relax:dualspace" );  }
namespace relax { BooleanOptionKey const cyclic_peptide( "relax:cyclic_peptide" );  }
namespace relax { namespace range { BooleanOptionKey const range( "relax:range" );  } }
namespace relax { namespace range { BooleanOptionKey const set_tm_helical( "relax:range:set_tm_helical" );  } }
namespace relax { namespace range { RealOptionKey const kT( "relax:range:kT" );  } }
namespace relax { namespace range { RealOptionKey const angle_max( "relax:range:angle_max" );  } }
namespace relax { namespace range { StringOptionKey const nmoves( "relax:range:nmoves" );  } }
namespace relax { namespace range { BooleanOptionKey const spherical_wave( "relax:range:spherical_wave" );  } }
namespace relax { namespace range { BooleanOptionKey const repack_again( "relax:range:repack_again" );  } }
namespace relax { namespace range { IntegerOptionKey const cycles( "relax:range:cycles" );  } }
namespace relax { namespace range { IntegerOptionKey const min_cycles( "relax:range:min_cycles" );  } }
namespace relax { namespace range { BooleanOptionKey const idealize( "relax:range:idealize" );  } }
namespace relax { namespace range { BooleanOptionKey const skip_relax( "relax:range:skip_relax" );  } }
namespace relax { BooleanOptionKey const ramady( "relax:ramady" );  }
namespace relax { RealOptionKey const ramady_rms_limit( "relax:ramady_rms_limit" );  }
namespace relax { RealOptionKey const ramady_cutoff( "relax:ramady_cutoff" );  }
namespace relax { IntegerOptionKey const ramady_max_rebuild( "relax:ramady_max_rebuild" );  }
namespace relax { BooleanOptionKey const ramady_force( "relax:ramady_force" );  }
namespace relax { FileOptionKey const script( "relax:script" );  }
namespace relax { IntegerOptionKey const script_max_accept( "relax:script_max_accept" );  }
namespace relax { BooleanOptionKey const superimpose_to_native( "relax:superimpose_to_native" );  }
namespace relax { FileOptionKey const superimpose_to_file( "relax:superimpose_to_file" );  }
namespace relax { BooleanOptionKey const constrain_relax_to_native_coords( "relax:constrain_relax_to_native_coords" );  }
namespace relax { BooleanOptionKey const constrain_relax_to_start_coords( "relax:constrain_relax_to_start_coords" );  }
namespace relax { BooleanOptionKey const coord_constrain_sidechains( "relax:coord_constrain_sidechains" );  }
namespace relax { RealOptionKey const sc_cst_maxdist( "relax:sc_cst_maxdist" );  }
namespace relax { BooleanOptionKey const limit_aroma_chi2( "relax:limit_aroma_chi2" );  }
namespace relax { BooleanOptionKey const respect_resfile( "relax:respect_resfile" );  }
namespace relax { BooleanOptionKey const bb_move( "relax:bb_move" );  }
namespace relax { BooleanOptionKey const chi_move( "relax:chi_move" );  }
namespace relax { BooleanOptionKey const jump_move( "relax:jump_move" );  }
namespace relax { BooleanOptionKey const dna_move( "relax:dna_move" );  }
namespace relax { BooleanOptionKey const fix_omega( "relax:fix_omega" );  }
namespace relax { BooleanOptionKey const minimize_bond_lengths( "relax:minimize_bond_lengths" );  }
namespace relax { BooleanOptionKey const minimize_bond_angles( "relax:minimize_bond_angles" );  }
namespace relax { IntegerOptionKey const minimize_bondlength_subset( "relax:minimize_bondlength_subset" );  }
namespace relax { IntegerOptionKey const minimize_bondangle_subset( "relax:minimize_bondangle_subset" );  }
namespace relax { StringOptionKey const min_type( "relax:min_type" );  }
namespace relax { BooleanOptionKey const cartesian( "relax:cartesian" );  }
namespace relax { RealOptionKey const chainbreak_weight( "relax:chainbreak_weight" );  }
namespace relax { RealOptionKey const linear_chainbreak_weight( "relax:linear_chainbreak_weight" );  }
namespace relax { RealOptionKey const overlap_chainbreak_weight( "relax:overlap_chainbreak_weight" );  }
namespace relax { BooleanOptionKey const classic( "relax:classic" );  }
namespace relax { FileOptionKey const sequence_file( "relax:sequence_file" );  }
namespace relax { BooleanOptionKey const quick( "relax:quick" );  }
namespace relax { BooleanOptionKey const sequence( "relax:sequence" );  }
namespace relax { IntegerOptionKey const minirelax_repeats( "relax:minirelax_repeats" );  }
namespace relax { RealOptionKey const minirelax_sdev( "relax:minirelax_sdev" );  }
namespace relax { BooleanOptionKey const wobblemoves( "relax:wobblemoves" );  }
namespace relax { FileOptionKey const constrain_relax_segments( "relax:constrain_relax_segments" );  }
namespace relax { RealOptionKey const coord_cst_width( "relax:coord_cst_width" );  }
namespace relax { RealOptionKey const coord_cst_stdev( "relax:coord_cst_stdev" );  }
namespace relax { BooleanOptionKey const ramp_constraints( "relax:ramp_constraints" );  }
namespace relax { RealOptionKey const energycut( "relax:energycut" );  }
namespace relax { BooleanOptionKey const mini( "relax:mini" );  }
namespace relax { IntegerOptionKey const stage1_ramp_cycles( "relax:stage1_ramp_cycles" );  }
namespace relax { IntegerOptionKey const stage1_ramp_inner_cycles( "relax:stage1_ramp_inner_cycles" );  }
namespace relax { IntegerOptionKey const stage2_repack_period( "relax:stage2_repack_period" );  }
namespace relax { IntegerOptionKey const stage2_cycles( "relax:stage2_cycles" );  }
namespace relax { RealOptionKey const min_tolerance( "relax:min_tolerance" );  }
namespace relax { IntegerOptionKey const stage3_cycles( "relax:stage3_cycles" );  }
namespace relax { RealOptionKey const cycle_ratio( "relax:cycle_ratio" );  }
namespace relax { RealOptionKey const filter_stage2_beginning( "relax:filter_stage2_beginning" );  }
namespace relax { RealOptionKey const filter_stage2_quarter( "relax:filter_stage2_quarter" );  }
namespace relax { RealOptionKey const filter_stage2_half( "relax:filter_stage2_half" );  }
namespace relax { RealOptionKey const filter_stage2_end( "relax:filter_stage2_end" );  }
namespace relax { namespace centroid { BooleanOptionKey const centroid( "relax:centroid" );  } }
namespace relax { namespace centroid { StringOptionKey const weights( "relax:centroid:weights" );  } }
namespace relax { namespace centroid { BooleanOptionKey const ramp_vdw( "relax:centroid:ramp_vdw" );  } }
namespace relax { namespace centroid { BooleanOptionKey const ramp_rama( "relax:centroid:ramp_rama" );  } }
namespace relax { namespace centroid { StringOptionKey const parameters( "relax:centroid:parameters" );  } }
namespace relax { namespace centroid { BooleanOptionKey const do_final_repack( "relax:centroid:do_final_repack" );  } }
namespace relax { namespace centroid { BooleanOptionKey const increase_vdw_radii( "relax:centroid:increase_vdw_radii" );  } }
namespace remodel { BooleanOptionKey const remodel( "remodel" );  }
namespace remodel { FileOptionKey const blueprint( "remodel:blueprint" );  }
namespace remodel { FileOptionKey const cstfile( "remodel:cstfile" );  }
namespace remodel { BooleanOptionKey const cst_fa_only( "remodel:cst_fa_only" );  }
namespace remodel { IntegerOptionKey const cstfilter( "remodel:cstfilter" );  }
namespace remodel { StringOptionKey const cen_sfxn( "remodel:cen_sfxn" );  }
namespace remodel { BooleanOptionKey const check_scored_centroid( "remodel:check_scored_centroid" );  }
namespace remodel { IntegerOptionKey const num_trajectory( "remodel:num_trajectory" );  }
namespace remodel { IntegerOptionKey const save_top( "remodel:save_top" );  }
namespace remodel { BooleanOptionKey const swap_refine_confirm_protocols( "remodel:swap_refine_confirm_protocols" );  }
namespace remodel { BooleanOptionKey const bypass_fragments( "remodel:bypass_fragments" );  }
namespace remodel { BooleanOptionKey const use_same_length_fragments( "remodel:use_same_length_fragments" );  }
namespace remodel { BooleanOptionKey const no_jumps( "remodel:no_jumps" );  }
namespace remodel { IntegerOptionKey const reroot_tree( "remodel:reroot_tree" );  }
namespace remodel { BooleanOptionKey const use_blueprint_sequence ( "remodel:use_blueprint_sequence " );  }
namespace remodel { BooleanOptionKey const quick_and_dirty ( "remodel:quick_and_dirty " );  }
namespace remodel { BooleanOptionKey const checkpoint ( "remodel:checkpoint " );  }
namespace remodel { BooleanOptionKey const use_pose_relax ( "remodel:use_pose_relax " );  }
namespace remodel { BooleanOptionKey const use_cart_relax ( "remodel:use_cart_relax " );  }
namespace remodel { BooleanOptionKey const free_relax ( "remodel:free_relax " );  }
namespace remodel { StringOptionKey const generic_aa( "remodel:generic_aa" );  }
namespace remodel { RealOptionKey const cluster_radius( "remodel:cluster_radius" );  }
namespace remodel { BooleanOptionKey const use_clusters( "remodel:use_clusters" );  }
namespace remodel { BooleanOptionKey const run_confirmation( "remodel:run_confirmation" );  }
namespace remodel { BooleanOptionKey const cluster_on_entire_pose( "remodel:cluster_on_entire_pose" );  }
namespace remodel { IntegerOptionKey const dr_cycles( "remodel:dr_cycles" );  }
namespace remodel { IntegerOptionKey const two_chain_tree( "remodel:two_chain_tree" );  }
namespace remodel { IntegerOptionKey const repeat_structure( "remodel:repeat_structure" );  }
namespace remodel { IntegerOptionKey const lh_ex_limit( "remodel:lh_ex_limit" );  }
namespace remodel { StringVectorOptionKey const lh_filter_string( "remodel:lh_filter_string" );  }
namespace remodel { IntegerOptionKey const lh_cbreak_selection( "remodel:lh_cbreak_selection" );  }
namespace remodel { BooleanOptionKey const lh_closure_filter( "remodel:lh_closure_filter" );  }
namespace remodel { BooleanOptionKey const cen_minimize( "remodel:cen_minimize" );  }
namespace remodel { IntegerOptionKey const core_cutoff( "remodel:core_cutoff" );  }
namespace remodel { IntegerOptionKey const boundary_cutoff( "remodel:boundary_cutoff" );  }
namespace remodel { StringOptionKey const coreAA( "remodel:coreAA" );  }
namespace remodel { StringOptionKey const boundaryAA( "remodel:boundaryAA" );  }
namespace remodel { StringOptionKey const surfaceAA( "remodel:surfaceAA" );  }
namespace remodel { BooleanOptionKey const design_around_ligand( "remodel:design_around_ligand" );  }
namespace remodel { BooleanOptionKey const move_ligand( "remodel:move_ligand" );  }
namespace remodel { BooleanOptionKey const resclass_by_sasa( "remodel:resclass_by_sasa" );  }
namespace remodel { RealOptionKey const helical_rise( "remodel:helical_rise" );  }
namespace remodel { RealOptionKey const helical_radius( "remodel:helical_radius" );  }
namespace remodel { RealOptionKey const helical_omega( "remodel:helical_omega" );  }
namespace remodel { RealVectorOptionKey const filter_rise( "remodel:filter_rise" );  }
namespace remodel { RealVectorOptionKey const filter_radius( "remodel:filter_radius" );  }
namespace remodel { RealVectorOptionKey const filter_omega( "remodel:filter_omega" );  }
namespace remodel { RealOptionKey const COM_sd( "remodel:COM_sd" );  }
namespace remodel { RealOptionKey const COM_tolerance( "remodel:COM_tolerance" );  }
namespace remodel { RealOptionKey const vdw( "remodel:vdw" );  }
namespace remodel { RealOptionKey const rama( "remodel:rama" );  }
namespace remodel { RealOptionKey const cbeta( "remodel:cbeta" );  }
namespace remodel { RealOptionKey const cenpack( "remodel:cenpack" );  }
namespace remodel { RealOptionKey const rg_local( "remodel:rg_local" );  }
namespace remodel { RealOptionKey const hb_lrbb( "remodel:hb_lrbb" );  }
namespace remodel { RealOptionKey const hb_srbb( "remodel:hb_srbb" );  }
namespace remodel { RealOptionKey const rg( "remodel:rg" );  }
namespace remodel { RealOptionKey const rsigma( "remodel:rsigma" );  }
namespace remodel { RealOptionKey const ss_pair( "remodel:ss_pair" );  }
namespace remodel { BooleanOptionKey const build_disulf( "remodel:build_disulf" );  }
namespace remodel { RealOptionKey const match_rt_limit( "remodel:match_rt_limit" );  }
namespace remodel { IntegerVectorOptionKey const disulf_landing_range( "remodel:disulf_landing_range" );  }
namespace remodel { BooleanOptionKey const rank_by_bsasa( "remodel:rank_by_bsasa" );  }
namespace remodel { namespace staged_sampling { BooleanOptionKey const staged_sampling( "remodel:staged_sampling" );  } }
namespace remodel { namespace staged_sampling { FileOptionKey const residues_to_sample( "remodel:staged_sampling:residues_to_sample" );  } }
namespace remodel { namespace staged_sampling { StringOptionKey const starting_sequence( "remodel:staged_sampling:starting_sequence" );  } }
namespace remodel { namespace staged_sampling { FileOptionKey const starting_pdb( "remodel:staged_sampling:starting_pdb" );  } }
namespace remodel { namespace staged_sampling { StringOptionKey const starting_non_canonical( "remodel:staged_sampling:starting_non_canonical" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const require_frags_match_blueprint( "remodel:staged_sampling:require_frags_match_blueprint" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const start_w_ideal_helices( "remodel:staged_sampling:start_w_ideal_helices" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const sample_over_loops( "remodel:staged_sampling:sample_over_loops" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const small_moves( "remodel:staged_sampling:small_moves" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const fa_mode( "remodel:staged_sampling:fa_mode" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const fa_relax_moves( "remodel:staged_sampling:fa_relax_moves" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const sym_move( "remodel:staged_sampling:sym_move" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const loop_btw_parametric_components( "remodel:staged_sampling:loop_btw_parametric_components" );  } }
namespace remodel { namespace staged_sampling { BooleanOptionKey const pre_centroid( "remodel:staged_sampling:pre_centroid" );  } }
namespace remodel { namespace domainFusion { BooleanOptionKey const domainFusion( "remodel:domainFusion" );  } }
namespace remodel { namespace domainFusion { FileOptionKey const insert_segment_from_pdb( "remodel:domainFusion:insert_segment_from_pdb" );  } }
namespace remodel { namespace domainFusion { FileOptionKey const insert_segment2_from_pdb( "remodel:domainFusion:insert_segment2_from_pdb" );  } }
namespace remodel { namespace design { BooleanOptionKey const design( "remodel:design" );  } }
namespace remodel { namespace design { BooleanOptionKey const no_design ( "remodel:design:no_design " );  } }
namespace remodel { namespace design { BooleanOptionKey const design_all( "remodel:design:design_all" );  } }
