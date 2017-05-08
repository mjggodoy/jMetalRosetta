namespace loops { BooleanOptionKey const derive_torsion_string_from_native_pose( "loops:derive_torsion_string_from_native_pose" );  }
namespace loops { BooleanOptionKey const always_remodel_full_loop( "loops:always_remodel_full_loop" );  }
namespace loops { BooleanOptionKey const taboo_sampling( "loops:taboo_sampling" );  }
namespace loops { BooleanOptionKey const taboo_in_fa( "loops:taboo_in_fa" );  }
namespace loops { BooleanOptionKey const ramp_fa_rep( "loops:ramp_fa_rep" );  }
namespace loops { BooleanOptionKey const ramp_rama( "loops:ramp_rama" );  }
namespace loops { BooleanOptionKey const kic_rama2b( "loops:kic_rama2b" );  }
namespace loops { BooleanOptionKey const kic_leave_centroid_after_initial_closure( "loops:kic_leave_centroid_after_initial_closure" );  }
namespace loops { BooleanOptionKey const legacy_kic( "loops:legacy_kic" );  }
namespace loops { BooleanOptionKey const alternative_closure_protocol( "loops:alternative_closure_protocol" );  }
namespace loops { RealOptionKey const chainbreak_max_accept( "loops:chainbreak_max_accept" );  }
namespace loops { BooleanOptionKey const debug_loop_closure( "loops:debug_loop_closure" );  }
namespace loops { BooleanOptionKey const non_ideal_loop_closing( "loops:non_ideal_loop_closing" );  }
namespace loops { RealOptionKey const scored_frag_cycles( "loops:scored_frag_cycles" );  }
namespace loops { RealOptionKey const short_frag_cycles( "loops:short_frag_cycles" );  }
namespace loops { RealOptionKey const rmsd_tol( "loops:rmsd_tol" );  }
namespace loops { RealOptionKey const chain_break_tol( "loops:chain_break_tol" );  }
namespace loops { BooleanOptionKey const random_loop( "loops:random_loop" );  }
namespace loops { FileVectorOptionKey const stealfrags( "loops:stealfrags" );  }
namespace loops { IntegerOptionKey const stealfrags_times( "loops:stealfrags_times" );  }
namespace loops { RealOptionKey const coord_cst( "loops:coord_cst" );  }
namespace loops { RealOptionKey const skip_1mers( "loops:skip_1mers" );  }
namespace loops { RealOptionKey const skip_3mers( "loops:skip_3mers" );  }
namespace loops { RealOptionKey const skip_9mers( "loops:skip_9mers" );  }
namespace loops { BooleanOptionKey const loop_model( "loops:loop_model" );  }
namespace loops { RealOptionKey const score_filter_cutoff( "loops:score_filter_cutoff" );  }
namespace loops { BooleanOptionKey const ccd_closure( "loops:ccd_closure" );  }
namespace loops { BooleanOptionKey const skip_ccd_moves( "loops:skip_ccd_moves" );  }
namespace loops { BooleanOptionKey const no_randomize_loop( "loops:no_randomize_loop" );  }
namespace loops { BooleanOptionKey const loops_subset( "loops:loops_subset" );  }
namespace loops { IntegerOptionKey const num_desired_loops( "loops:num_desired_loops" );  }
namespace loops { RealOptionKey const loop_combine_rate( "loops:loop_combine_rate" );  }
namespace loops { RealOptionKey const final_score_filter( "loops:final_score_filter" );  }
namespace loops { BooleanOptionKey const no_combine_if_fail( "loops:no_combine_if_fail" );  }
namespace loops { BooleanOptionKey const shorten_long_terminal_loop( "loops:shorten_long_terminal_loop" );  }
namespace loops { IntegerOptionKey const backrub_trials( "loops:backrub_trials" );  }
namespace loops { RealOptionKey const looprlx_cycle_ratio( "loops:looprlx_cycle_ratio" );  }
namespace loops { RealOptionKey const extended_beta( "loops:extended_beta" );  }
namespace loops { BooleanOptionKey const no_looprebuild( "loops:no_looprebuild" );  }
namespace loops { BooleanOptionKey const allow_lig_move( "loops:allow_lig_move" );  }
namespace loops { FileOptionKey const keep_natro( "loops:keep_natro" );  }
namespace loops { IntegerOptionKey const refine_design_iterations( "loops:refine_design_iterations" );  }
namespace loops { namespace ccd { BooleanOptionKey const ccd( "loops:ccd" );  } }
namespace loops { namespace ccd { RealOptionKey const max_rama_score_increase( "loops:ccd:max_rama_score_increase" );  } }
namespace loops { namespace ccd { RealVectorOptionKey const max_torsion_delta_per_move( "loops:ccd:max_torsion_delta_per_move" );  } }
namespace loops { namespace ccd { RealVectorOptionKey const max_torsion_delta( "loops:ccd:max_torsion_delta" );  } }
namespace loops { namespace ccd { RealOptionKey const tolerance( "loops:ccd:tolerance" );  } }
namespace loops { namespace ccd { IntegerOptionKey const max_cycles( "loops:ccd:max_cycles" );  } }
namespace loops { namespace ccd { BooleanOptionKey const check_rama_scores( "loops:ccd:check_rama_scores" );  } }
namespace loops { namespace ccd { BooleanOptionKey const rama_2b( "loops:ccd:rama_2b" );  } }
namespace loops { namespace ccd { RealOptionKey const temperature( "loops:ccd:temperature" );  } }
namespace intf { BooleanOptionKey const intf( "intf" );  }
namespace intf { StringOptionKey const chains( "intf:chains" );  }
namespace mp { BooleanOptionKey const mp( "mp" );  }
namespace mp { RealOptionKey const thickness( "mp:thickness" );  }
namespace mp { RealOptionKey const steepness( "mp:steepness" );  }
namespace mp { RealVectorOptionKey const center_start( "mp:center_start" );  }
namespace mp { RealOptionKey const center_delta( "mp:center_delta" );  }
namespace mp { RealOptionKey const center_search_cycles( "mp:center_search_cycles" );  }
namespace mp { RealVectorOptionKey const normal_start( "mp:normal_start" );  }
namespace mp { RealOptionKey const normal_angle_start( "mp:normal_angle_start" );  }
namespace mp { RealOptionKey const normal_angle_delta( "mp:normal_angle_delta" );  }
namespace mp { RealOptionKey const normal_search_cycles( "mp:normal_search_cycles" );  }
namespace mp { RealOptionKey const chain_normal_angle_max( "mp:chain_normal_angle_max" );  }
namespace mp { RealOptionKey const pose_normal_angle_max( "mp:pose_normal_angle_max" );  }
namespace mp { BooleanOptionKey const no_interpolate_Mpair( "mp:no_interpolate_Mpair" );  }
namespace mp { BooleanOptionKey const Hbond_depth_correction( "mp:Hbond_depth_correction" );  }
namespace mp { BooleanOptionKey const TMprojection( "mp:TMprojection" );  }
namespace mp { RealOptionKey const wt_TMprojection( "mp:wt_TMprojection" );  }
namespace mp { BooleanOptionKey const non_helix( "mp:non_helix" );  }
namespace mp { RealOptionKey const wt_non_helix( "mp:wt_non_helix" );  }
namespace mp { BooleanOptionKey const termini( "mp:termini" );  }
namespace mp { RealOptionKey const wt_termini( "mp:wt_termini" );  }
namespace mp { BooleanOptionKey const secstruct( "mp:secstruct" );  }
namespace mp { RealOptionKey const wt_secstruct( "mp:wt_secstruct" );  }
namespace mp { BooleanOptionKey const spanning( "mp:spanning" );  }
namespace mp { RealOptionKey const wt_spanning( "mp:wt_spanning" );  }
namespace mp { namespace viewer { BooleanOptionKey const viewer( "mp:viewer" );  } }
namespace mp { namespace viewer { RealOptionKey const thickness( "mp:viewer:thickness" );  } }
namespace mp { namespace viewer { IntegerOptionKey const num_points( "mp:viewer:num_points" );  } }
namespace mp { namespace visualize { BooleanOptionKey const visualize( "mp:visualize" );  } }
namespace mp { namespace visualize { BooleanOptionKey const embedding( "mp:visualize:embedding" );  } }
namespace mp { namespace visualize { RealOptionKey const spacing( "mp:visualize:spacing" );  } }
namespace mp { namespace visualize { RealOptionKey const width( "mp:visualize:width" );  } }
namespace mp { namespace visualize { RealOptionKey const thickness( "mp:visualize:thickness" );  } }
namespace mp { namespace visualize { RealOptionKey const plane_radius( "mp:visualize:plane_radius" );  } }
namespace mp { namespace scoring { BooleanOptionKey const scoring( "mp:scoring" );  } }
namespace mp { namespace scoring { BooleanOptionKey const hbond( "mp:scoring:hbond" );  } }
namespace mp { namespace setup { BooleanOptionKey const setup( "mp:setup" );  } }
namespace mp { namespace setup { StringVectorOptionKey const spanfiles( "mp:setup:spanfiles" );  } }
namespace mp { namespace setup { BooleanOptionKey const spans_from_structure( "mp:setup:spans_from_structure" );  } }
namespace mp { namespace setup { StringOptionKey const lipsfile( "mp:setup:lipsfile" );  } }
namespace mp { namespace setup { RealVectorOptionKey const center( "mp:setup:center" );  } }
namespace mp { namespace setup { RealVectorOptionKey const normal( "mp:setup:normal" );  } }
namespace mp { namespace setup { RealOptionKey const membrane_rsd( "mp:setup:membrane_rsd" );  } }
namespace mp { namespace setup { BooleanOptionKey const transform_into_membrane( "mp:setup:transform_into_membrane" );  } }
namespace mp { namespace setup { BooleanOptionKey const position_from_topo( "mp:setup:position_from_topo" );  } }
namespace mp { namespace setup { BooleanOptionKey const optimize1( "mp:setup:optimize1" );  } }
namespace mp { namespace setup { BooleanOptionKey const optimize2( "mp:setup:optimize2" );  } }
namespace mp { namespace lipid_acc { BooleanOptionKey const lipid_acc( "mp:lipid_acc" );  } }
namespace mp { namespace lipid_acc { RealOptionKey const angle_cutoff( "mp:lipid_acc:angle_cutoff" );  } }
namespace mp { namespace lipid_acc { RealOptionKey const slice_width( "mp:lipid_acc:slice_width" );  } }
namespace mp { namespace lipid_acc { RealOptionKey const shell_radius( "mp:lipid_acc:shell_radius" );  } }
namespace mp { namespace lipid_acc { RealOptionKey const dist_cutoff( "mp:lipid_acc:dist_cutoff" );  } }
namespace mp { namespace lipid_acc { BooleanOptionKey const tm_alpha( "mp:lipid_acc:tm_alpha" );  } }
namespace mp { namespace transform { BooleanOptionKey const transform( "mp:transform" );  } }
namespace mp { namespace transform { BooleanOptionKey const optimize_embedding( "mp:transform:optimize_embedding" );  } }
namespace mp { namespace dock { BooleanOptionKey const dock( "mp:dock" );  } }
namespace mp { namespace dock { StringOptionKey const weights_cen( "mp:dock:weights_cen" );  } }
namespace mp { namespace dock { StringOptionKey const weights_fa( "mp:dock:weights_fa" );  } }
namespace mp { namespace dock { BooleanOptionKey const lowres( "mp:dock:lowres" );  } }
namespace mp { namespace dock { BooleanOptionKey const allow_flips( "mp:dock:allow_flips" );  } }
namespace mp { namespace dock { BooleanOptionKey const flexible_bb( "mp:dock:flexible_bb" );  } }
namespace mp { namespace dock { BooleanOptionKey const flexible_sc( "mp:dock:flexible_sc" );  } }
namespace mp { namespace dock { RealOptionKey const slide_threshold( "mp:dock:slide_threshold" );  } }
namespace mp { namespace quickrelax { BooleanOptionKey const quickrelax( "mp:quickrelax" );  } }
namespace mp { namespace quickrelax { RealOptionKey const angle_max( "mp:quickrelax:angle_max" );  } }
namespace mp { namespace quickrelax { StringOptionKey const nmoves( "mp:quickrelax:nmoves" );  } }
namespace mp { namespace quickrelax { BooleanOptionKey const repack_again( "mp:quickrelax:repack_again" );  } }
namespace mp { namespace mutate_relax { BooleanOptionKey const mutate_relax( "mp:mutate_relax" );  } }
namespace mp { namespace mutate_relax { StringOptionKey const mutation( "mp:mutate_relax:mutation" );  } }
namespace mp { namespace mutate_relax { StringOptionKey const mutant_file( "mp:mutate_relax:mutant_file" );  } }
namespace mp { namespace mutate_relax { IntegerOptionKey const nmodels( "mp:mutate_relax:nmodels" );  } }
namespace mp { namespace mutate_relax { BooleanOptionKey const repack_mutation_only( "mp:mutate_relax:repack_mutation_only" );  } }
namespace mp { namespace mutate_relax { RealOptionKey const repack_radius( "mp:mutate_relax:repack_radius" );  } }
namespace mp { namespace mutate_relax { BooleanOptionKey const relax( "mp:mutate_relax:relax" );  } }
namespace mp { namespace benchmark { BooleanOptionKey const benchmark( "mp:benchmark" );  } }
namespace mp { namespace benchmark { namespace ideal_helix { BooleanOptionKey const ideal_helix( "mp:benchmark:ideal_helix" );  } } }
namespace mp { namespace benchmark { namespace ideal_helix { RealOptionKey const helix_start( "mp:benchmark:ideal_helix:helix_start" );  } } }
namespace mp { namespace benchmark { namespace ideal_helix { RealOptionKey const helix_end( "mp:benchmark:ideal_helix:helix_end" );  } } }
namespace mp { namespace benchmark { namespace tilt_angle { BooleanOptionKey const tilt_angle( "mp:benchmark:tilt_angle" );  } } }
namespace mp { namespace benchmark { namespace tilt_angle { StringOptionKey const output( "mp:benchmark:tilt_angle:output" );  } } }
namespace mp { namespace output { BooleanOptionKey const output( "mp:output" );  } }
namespace mp { namespace output { BooleanOptionKey const normalize_to_thk( "mp:output:normalize_to_thk" );  } }
namespace membrane { BooleanOptionKey const membrane( "membrane" );  }
namespace membrane { IntegerOptionKey const normal_cycles( "membrane:normal_cycles" );  }
namespace membrane { RealOptionKey const normal_mag( "membrane:normal_mag" );  }
namespace membrane { RealOptionKey const center_mag( "membrane:center_mag" );  }
namespace membrane { RealOptionKey const smooth_move_frac( "membrane:smooth_move_frac" );  }
namespace membrane { BooleanOptionKey const no_interpolate_Mpair( "membrane:no_interpolate_Mpair" );  }
namespace membrane { BooleanOptionKey const Menv_penalties( "membrane:Menv_penalties" );  }
namespace membrane { BooleanOptionKey const Membed_init( "membrane:Membed_init" );  }
namespace membrane { BooleanOptionKey const Fa_Membed_update( "membrane:Fa_Membed_update" );  }
namespace membrane { BooleanOptionKey const center_search( "membrane:center_search" );  }
namespace membrane { BooleanOptionKey const normal_search( "membrane:normal_search" );  }
namespace membrane { IntegerOptionKey const center_max_delta( "membrane:center_max_delta" );  }
namespace membrane { IntegerOptionKey const normal_start_angle( "membrane:normal_start_angle" );  }
namespace membrane { IntegerOptionKey const normal_delta_angle( "membrane:normal_delta_angle" );  }
namespace membrane { IntegerOptionKey const normal_max_angle( "membrane:normal_max_angle" );  }
namespace membrane { BooleanOptionKey const debug( "membrane:debug" );  }
namespace membrane { BooleanOptionKey const fixed_membrane( "membrane:fixed_membrane" );  }
namespace membrane { RealVectorOptionKey const membrane_center( "membrane:membrane_center" );  }
namespace membrane { RealVectorOptionKey const membrane_normal( "membrane:membrane_normal" );  }
namespace membrane { BooleanOptionKey const view( "membrane:view" );  }
namespace membrane { BooleanOptionKey const Mhbond_depth( "membrane:Mhbond_depth" );  }
namespace membrane { RealOptionKey const thickness( "membrane:thickness" );  }
namespace mistakes { BooleanOptionKey const mistakes( "mistakes" );  }
namespace mistakes { BooleanOptionKey const restore_pre_talaris_2013_behavior( "mistakes:restore_pre_talaris_2013_behavior" );  }
namespace mistakes { namespace chemical { BooleanOptionKey const chemical( "mistakes:chemical" );  } }
namespace mistakes { namespace chemical { BooleanOptionKey const pre_talaris2013_geometries( "mistakes:chemical:pre_talaris2013_geometries" );  } }
namespace MonteCarlo { BooleanOptionKey const MonteCarlo( "MonteCarlo" );  }
namespace MonteCarlo { RealOptionKey const temp_initial( "MonteCarlo:temp_initial" );  }
namespace MonteCarlo { RealOptionKey const temp_final( "MonteCarlo:temp_final" );  }
namespace optimization { BooleanOptionKey const optimization( "optimization" );  }
namespace optimization { IntegerOptionKey const default_max_cycles( "optimization:default_max_cycles" );  }
namespace optimization { RealOptionKey const armijo_min_stepsize( "optimization:armijo_min_stepsize" );  }
namespace optimization { RealOptionKey const scale_normalmode_dampen( "optimization:scale_normalmode_dampen" );  }
namespace optimization { IntegerOptionKey const lbfgs_M( "optimization:lbfgs_M" );  }
namespace optimization { RealOptionKey const scale_d( "optimization:scale_d" );  }
namespace optimization { RealOptionKey const scale_theta( "optimization:scale_theta" );  }
namespace optimization { RealOptionKey const scale_rb( "optimization:scale_rb" );  }
namespace optimization { RealOptionKey const scale_rbangle( "optimization:scale_rbangle" );  }
namespace optimization { BooleanOptionKey const scmin_nonideal( "optimization:scmin_nonideal" );  }
namespace optimization { BooleanOptionKey const scmin_cartesian( "optimization:scmin_cartesian" );  }
namespace optimization { BooleanOptionKey const nonideal( "optimization:nonideal" );  }
namespace optimization { BooleanOptionKey const old_sym_min( "optimization:old_sym_min" );  }
namespace optimization { BooleanOptionKey const debug_inaccurate_G( "optimization:debug_inaccurate_G" );  }
namespace orbitals { BooleanOptionKey const orbitals( "orbitals" );  }
namespace orbitals { BooleanOptionKey const Hpol( "orbitals:Hpol" );  }
namespace orbitals { BooleanOptionKey const Haro( "orbitals:Haro" );  }
namespace orbitals { BooleanOptionKey const bb_stats( "orbitals:bb_stats" );  }
namespace orbitals { BooleanOptionKey const sc_stats( "orbitals:sc_stats" );  }
namespace pose_metrics { BooleanOptionKey const pose_metrics( "pose_metrics" );  }
namespace pose_metrics { RealOptionKey const atomic_burial_cutoff( "pose_metrics:atomic_burial_cutoff" );  }
namespace pose_metrics { RealOptionKey const sasa_calculator_probe_radius( "pose_metrics:sasa_calculator_probe_radius" );  }
namespace pose_metrics { RealOptionKey const interface_cutoff( "pose_metrics:interface_cutoff" );  }
namespace pose_metrics { IntegerOptionKey const min_sequence_separation( "pose_metrics:min_sequence_separation" );  }
namespace pose_metrics { RealOptionKey const contact_cutoffE( "pose_metrics:contact_cutoffE" );  }
namespace pose_metrics { RealOptionKey const neighbor_by_distance_cutoff( "pose_metrics:neighbor_by_distance_cutoff" );  }
namespace pose_metrics { RealOptionKey const inter_group_neighbors_cutoff( "pose_metrics:inter_group_neighbors_cutoff" );  }
namespace pose_metrics { RealOptionKey const semiex_water_burial_cutoff( "pose_metrics:semiex_water_burial_cutoff" );  }
namespace pose_metrics { namespace shobuns { BooleanOptionKey const shobuns( "pose_metrics:shobuns" );  } }
namespace pose_metrics { namespace shobuns { StringOptionKey const tgt_amino( "pose_metrics:shobuns:tgt_amino" );  } }
namespace pose_metrics { namespace shobuns { StringOptionKey const tgt_atom( "pose_metrics:shobuns:tgt_atom" );  } }
namespace pose_metrics { namespace shobuns { FileOptionKey const tgt_res( "pose_metrics:shobuns:tgt_res" );  } }
namespace pose_metrics { namespace shobuns { RealOptionKey const sho_cutoff( "pose_metrics:shobuns:sho_cutoff" );  } }
namespace rigid { BooleanOptionKey const rigid( "rigid" );  }
namespace rigid { RealOptionKey const chainbreak_bias( "rigid:chainbreak_bias" );  }
namespace rigid { BooleanOptionKey const close_loops( "rigid:close_loops" );  }
namespace rigid { IntegerOptionKey const fragment_cycles( "rigid:fragment_cycles" );  }
namespace rigid { BooleanOptionKey const log_accepted_moves( "rigid:log_accepted_moves" );  }
namespace rigid { RealOptionKey const max_ca_ca_dist( "rigid:max_ca_ca_dist" );  }
namespace rigid { FileOptionKey const patch( "rigid:patch" );  }
namespace rigid { IntegerOptionKey const residues_backbone_move( "rigid:residues_backbone_move" );  }
namespace rigid { RealOptionKey const rotation( "rigid:rotation" );  }
namespace rigid { FileOptionKey const sampling_prob( "rigid:sampling_prob" );  }
namespace rigid { StringOptionKey const score( "rigid:score" );  }
namespace rigid { IntegerOptionKey const sequence_separation( "rigid:sequence_separation" );  }
namespace rigid { IntegerOptionKey const small_cycles( "rigid:small_cycles" );  }
namespace rigid { IntegerOptionKey const stages( "rigid:stages" );  }
namespace rigid { RealOptionKey const temperature( "rigid:temperature" );  }
namespace rigid { RealOptionKey const translation( "rigid:translation" );  }
namespace sasa { BooleanOptionKey const sasa( "sasa" );  }
namespace sasa { StringOptionKey const method( "sasa:method" );  }
namespace sasa { BooleanOptionKey const include_hydrogens_explicitly( "sasa:include_hydrogens_explicitly" );  }
namespace sasa { RealOptionKey const probe_radius( "sasa:probe_radius" );  }
namespace sasa { BooleanOptionKey const include_probe_radius_in_atom_radii( "sasa:include_probe_radius_in_atom_radii" );  }
namespace sasa { BooleanOptionKey const include_only_C_S_in_hsasa( "sasa:include_only_C_S_in_hsasa" );  }
namespace sasa { BooleanOptionKey const exclude_polar_atoms_by_charge_in_hsasa( "sasa:exclude_polar_atoms_by_charge_in_hsasa" );  }
namespace sasa { RealOptionKey const polar_charge_cutoff( "sasa:polar_charge_cutoff" );  }
namespace sasa { StringOptionKey const implicit_hydrogen_radii_set( "sasa:implicit_hydrogen_radii_set" );  }
namespace sasa { StringOptionKey const explicit_hydrogen_radii_set( "sasa:explicit_hydrogen_radii_set" );  }
namespace symmetry { BooleanOptionKey const symmetry( "symmetry" );  }
namespace symmetry { StringOptionKey const symmetry_definition( "symmetry:symmetry_definition" );  }
namespace symmetry { RealOptionKey const reweight_symm_interactions( "symmetry:reweight_symm_interactions" );  }
namespace symmetry { BooleanOptionKey const initialize_rigid_body_dofs( "symmetry:initialize_rigid_body_dofs" );  }
namespace symmetry { BooleanOptionKey const detect_bonds( "symmetry:detect_bonds" );  }
namespace symmetry { RealVectorOptionKey const perturb_rigid_body_dofs( "symmetry:perturb_rigid_body_dofs" );  }
namespace symmetry { BooleanOptionKey const symmetric_rmsd( "symmetry:symmetric_rmsd" );  }
namespace abinitio { BooleanOptionKey const abinitio( "abinitio" );  }
namespace abinitio { BooleanOptionKey const membrane( "abinitio:membrane" );  }
namespace abinitio { BooleanOptionKey const use_loophash_filter( "abinitio:use_loophash_filter" );  }
namespace abinitio { RealOptionKey const loophash_filter_acceptance_rate( "abinitio:loophash_filter_acceptance_rate" );  }
namespace abinitio { FileOptionKey const kill_hairpins( "abinitio:kill_hairpins" );  }
namespace abinitio { RealOptionKey const kill_hairpins_frequency( "abinitio:kill_hairpins_frequency" );  }
namespace abinitio { BooleanOptionKey const smooth_cycles_only( "abinitio:smooth_cycles_only" );  }
namespace abinitio { BooleanOptionKey const relax( "abinitio:relax" );  }
namespace abinitio { BooleanOptionKey const final_clean_relax( "abinitio:final_clean_relax" );  }
namespace abinitio { BooleanOptionKey const fastrelax( "abinitio:fastrelax" );  }
namespace abinitio { BooleanOptionKey const multifastrelax( "abinitio:multifastrelax" );  }
namespace abinitio { BooleanOptionKey const debug( "abinitio:debug" );  }
namespace abinitio { BooleanOptionKey const clear_pose_cache( "abinitio:clear_pose_cache" );  }
namespace abinitio { BooleanOptionKey const explicit_pdb_debug( "abinitio:explicit_pdb_debug" );  }
namespace abinitio { BooleanOptionKey const use_filters( "abinitio:use_filters" );  }
namespace abinitio { RealOptionKey const increase_cycles( "abinitio:increase_cycles" );  }
namespace abinitio { StringOptionKey const jMetal_strategy( "abinitio:jMetal_strategy" );  }
