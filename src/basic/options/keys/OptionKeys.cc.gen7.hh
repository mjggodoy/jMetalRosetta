namespace cyclic_peptide { BooleanOptionKey const count_sc_hbonds( "cyclic_peptide:count_sc_hbonds" );  }
namespace cyclic_peptide { BooleanOptionKey const require_disulfides( "cyclic_peptide:require_disulfides" );  }
namespace cyclic_peptide { RealOptionKey const disulf_cutoff_prerelax( "cyclic_peptide:disulf_cutoff_prerelax" );  }
namespace cyclic_peptide { RealOptionKey const disulf_cutoff_postrelax( "cyclic_peptide:disulf_cutoff_postrelax" );  }
namespace cyclic_peptide { RealVectorOptionKey const user_set_alpha_dihedrals( "cyclic_peptide:user_set_alpha_dihedrals" );  }
namespace cyclic_peptide { RealOptionKey const user_set_alpha_dihedral_perturbation( "cyclic_peptide:user_set_alpha_dihedral_perturbation" );  }
namespace cyclic_peptide { BooleanOptionKey const filter_oversaturated_hbond_acceptors( "cyclic_peptide:filter_oversaturated_hbond_acceptors" );  }
namespace cyclic_peptide { RealOptionKey const hbond_acceptor_energy_cutoff( "cyclic_peptide:hbond_acceptor_energy_cutoff" );  }
namespace cyclic_peptide { RealOptionKey const sample_cis_pro_frequency( "cyclic_peptide:sample_cis_pro_frequency" );  }
namespace cyclic_peptide { BooleanOptionKey const design_peptide( "cyclic_peptide:design_peptide" );  }
namespace cyclic_peptide { FileOptionKey const allowed_residues_by_position( "cyclic_peptide:allowed_residues_by_position" );  }
namespace cyclic_peptide { BooleanOptionKey const prohibit_D_at_negative_phi( "cyclic_peptide:prohibit_D_at_negative_phi" );  }
namespace cyclic_peptide { BooleanOptionKey const prohibit_L_at_positive_phi( "cyclic_peptide:prohibit_L_at_positive_phi" );  }
namespace cyclic_peptide { FileOptionKey const L_alpha_comp_file( "cyclic_peptide:L_alpha_comp_file" );  }
namespace cyclic_peptide { FileOptionKey const D_alpha_comp_file( "cyclic_peptide:D_alpha_comp_file" );  }
namespace cyclic_peptide { FileOptionKey const L_beta_comp_file( "cyclic_peptide:L_beta_comp_file" );  }
namespace cyclic_peptide { FileOptionKey const D_beta_comp_file( "cyclic_peptide:D_beta_comp_file" );  }
namespace cyclic_peptide { IntegerOptionKey const angle_relax_rounds( "cyclic_peptide:angle_relax_rounds" );  }
namespace cyclic_peptide { IntegerOptionKey const angle_length_relax_rounds( "cyclic_peptide:angle_length_relax_rounds" );  }
namespace cyclic_peptide { IntegerOptionKey const cartesian_relax_rounds( "cyclic_peptide:cartesian_relax_rounds" );  }
namespace cyclic_peptide { BooleanOptionKey const use_classic_rama_for_sampling( "cyclic_peptide:use_classic_rama_for_sampling" );  }
namespace cyclic_peptide { IntegerVectorOptionKey const n_methyl_positions( "cyclic_peptide:n_methyl_positions" );  }
namespace cyclic_peptide { IntegerVectorOptionKey const TBMB_positions( "cyclic_peptide:TBMB_positions" );  }
namespace cyclic_peptide { BooleanOptionKey const use_TBMB_filters( "cyclic_peptide:use_TBMB_filters" );  }
namespace cyclic_peptide { RealOptionKey const TBMB_sidechain_distance_filter_multiplier( "cyclic_peptide:TBMB_sidechain_distance_filter_multiplier" );  }
namespace cyclic_peptide { RealOptionKey const TBMB_constraints_energy_filter_multiplier( "cyclic_peptide:TBMB_constraints_energy_filter_multiplier" );  }
namespace cyclic_peptide { BooleanOptionKey const link_all_cys_with_TBMB( "cyclic_peptide:link_all_cys_with_TBMB" );  }
namespace cyclic_peptide { IntegerOptionKey const require_symmetry_repeats( "cyclic_peptide:require_symmetry_repeats" );  }
namespace cyclic_peptide { BooleanOptionKey const require_symmetry_mirroring( "cyclic_peptide:require_symmetry_mirroring" );  }
namespace cyclic_peptide { RealOptionKey const require_symmetry_angle_threshold( "cyclic_peptide:require_symmetry_angle_threshold" );  }
namespace cyclic_peptide { RealOptionKey const require_symmetry_perturbation( "cyclic_peptide:require_symmetry_perturbation" );  }
namespace cyclic_peptide { IntegerVectorOptionKey const MPI_processes_by_level( "cyclic_peptide:MPI_processes_by_level" );  }
namespace cyclic_peptide { IntegerVectorOptionKey const MPI_batchsize_by_level( "cyclic_peptide:MPI_batchsize_by_level" );  }
namespace cyclic_peptide { StringOptionKey const MPI_sort_by( "cyclic_peptide:MPI_sort_by" );  }
namespace cyclic_peptide { BooleanOptionKey const MPI_choose_highest( "cyclic_peptide:MPI_choose_highest" );  }
namespace cyclic_peptide { RealOptionKey const MPI_output_fraction( "cyclic_peptide:MPI_output_fraction" );  }
namespace cyclic_peptide { IntegerOptionKey const MPI_stop_after_time( "cyclic_peptide:MPI_stop_after_time" );  }
namespace cyclic_peptide { RealOptionKey const MPI_pnear_lambda( "cyclic_peptide:MPI_pnear_lambda" );  }
namespace cyclic_peptide { RealOptionKey const MPI_pnear_kbt( "cyclic_peptide:MPI_pnear_kbt" );  }
namespace dc { BooleanOptionKey const dc( "dc" );  }
namespace dc { BooleanOptionKey const useZ( "dc:useZ" );  }
namespace ddg { BooleanOptionKey const ddg( "ddg" );  }
namespace ddg { BooleanOptionKey const opt_input_structure( "ddg:opt_input_structure" );  }
namespace ddg { BooleanOptionKey const pack_until_converge( "ddg:pack_until_converge" );  }
namespace ddg { BooleanOptionKey const no_constraints( "ddg:no_constraints" );  }
namespace ddg { BooleanOptionKey const use_rotamer_constraints_to_native( "ddg:use_rotamer_constraints_to_native" );  }
namespace ddg { BooleanOptionKey const suppress_checkpointing( "ddg:suppress_checkpointing" );  }
namespace ddg { BooleanOptionKey const wt_only( "ddg:wt_only" );  }
namespace ddg { BooleanOptionKey const mut_only( "ddg:mut_only" );  }
namespace ddg { BooleanOptionKey const output_silent( "ddg:output_silent" );  }
namespace ddg { StringOptionKey const minimization_scorefunction( "ddg:minimization_scorefunction" );  }
namespace ddg { StringOptionKey const minimization_patch( "ddg:minimization_patch" );  }
namespace ddg { BooleanOptionKey const min_cst( "ddg:min_cst" );  }
namespace ddg { IntegerOptionKey const lowest_x_decoys( "ddg:lowest_x_decoys" );  }
namespace ddg { BooleanOptionKey const local_opt_only( "ddg:local_opt_only" );  }
namespace ddg { BooleanOptionKey const print_per_res_diff( "ddg:print_per_res_diff" );  }
namespace ddg { BooleanOptionKey const mean( "ddg:mean" );  }
namespace ddg { BooleanOptionKey const min( "ddg:min" );  }
namespace ddg { BooleanOptionKey const no_bb_movement( "ddg:no_bb_movement" );  }
namespace ddg { BooleanOptionKey const initial_repack( "ddg:initial_repack" );  }
namespace ddg { IntegerOptionKey const interface_ddg( "ddg:interface_ddg" );  }
namespace ddg { RealOptionKey const ens_variation( "ddg:ens_variation" );  }
namespace ddg { BooleanOptionKey const sc_min_only( "ddg:sc_min_only" );  }
namespace ddg { StringOptionKey const min_cst_weights( "ddg:min_cst_weights" );  }
namespace ddg { RealOptionKey const opt_radius( "ddg:opt_radius" );  }
namespace ddg { StringOptionKey const output_dir( "ddg:output_dir" );  }
namespace ddg { StringOptionKey const last_accepted_pose_dir( "ddg:last_accepted_pose_dir" );  }
namespace ddg { BooleanOptionKey const min_with_cst( "ddg:min_with_cst" );  }
namespace ddg { RealOptionKey const temperature( "ddg:temperature" );  }
namespace ddg { BooleanOptionKey const ramp_repulsive( "ddg:ramp_repulsive" );  }
namespace ddg { StringOptionKey const mut_file( "ddg:mut_file" );  }
namespace ddg { StringOptionKey const out_pdb_prefix( "ddg:out_pdb_prefix" );  }
namespace ddg { RealOptionKey const constraint_weight( "ddg:constraint_weight" );  }
namespace ddg { RealOptionKey const harmonic_ca_tether( "ddg:harmonic_ca_tether" );  }
namespace ddg { IntegerOptionKey const iterations( "ddg:iterations" );  }
namespace ddg { StringOptionKey const out( "ddg:out" );  }
namespace ddg { BooleanOptionKey const debug_output( "ddg:debug_output" );  }
namespace ddg { BooleanOptionKey const dump_pdbs( "ddg:dump_pdbs" );  }
namespace ddg { StringOptionKey const weight_file( "ddg:weight_file" );  }
namespace DenovoProteinDesign { BooleanOptionKey const DenovoProteinDesign( "DenovoProteinDesign" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_core( "DenovoProteinDesign:redesign_core" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_loops( "DenovoProteinDesign:redesign_loops" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_surface( "DenovoProteinDesign:redesign_surface" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_complete( "DenovoProteinDesign:redesign_complete" );  }
namespace DenovoProteinDesign { BooleanOptionKey const optimize_loops( "DenovoProteinDesign:optimize_loops" );  }
namespace DenovoProteinDesign { FileOptionKey const secondary_structure_file( "DenovoProteinDesign:secondary_structure_file" );  }
namespace DenovoProteinDesign { FileOptionKey const hydrophobic_polar_pattern( "DenovoProteinDesign:hydrophobic_polar_pattern" );  }
namespace DenovoProteinDesign { BooleanOptionKey const use_template_sequence( "DenovoProteinDesign:use_template_sequence" );  }
namespace DenovoProteinDesign { BooleanOptionKey const use_template_topology( "DenovoProteinDesign:use_template_topology" );  }
namespace DenovoProteinDesign { FileOptionKey const create_from_template_pdb( "DenovoProteinDesign:create_from_template_pdb" );  }
namespace DenovoProteinDesign { BooleanOptionKey const create_from_secondary_structure( "DenovoProteinDesign:create_from_secondary_structure" );  }
namespace dna { BooleanOptionKey const dna( "dna" );  }
namespace dna { namespace specificity { BooleanOptionKey const specificity( "dna:specificity" );  } }
namespace dna { namespace specificity { BooleanOptionKey const exclude_dna_dna( "dna:specificity:exclude_dna_dna" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const params( "dna:specificity:params" );  } }
namespace dna { namespace specificity { FileVectorOptionKey const frag_files( "dna:specificity:frag_files" );  } }
namespace dna { namespace specificity { BooleanOptionKey const only_repack( "dna:specificity:only_repack" );  } }
namespace dna { namespace specificity { BooleanOptionKey const design_DNA( "dna:specificity:design_DNA" );  } }
namespace dna { namespace specificity { BooleanOptionKey const soft_rep( "dna:specificity:soft_rep" );  } }
namespace dna { namespace specificity { BooleanOptionKey const dump_pdbs( "dna:specificity:dump_pdbs" );  } }
namespace dna { namespace specificity { BooleanOptionKey const fast( "dna:specificity:fast" );  } }
namespace dna { namespace specificity { BooleanOptionKey const randomize_motif( "dna:specificity:randomize_motif" );  } }
namespace dna { namespace specificity { RealOptionKey const Wfa_elec( "dna:specificity:Wfa_elec" );  } }
namespace dna { namespace specificity { RealOptionKey const Wdna_bs( "dna:specificity:Wdna_bs" );  } }
namespace dna { namespace specificity { RealOptionKey const Wdna_bp( "dna:specificity:Wdna_bp" );  } }
namespace dna { namespace specificity { RealOptionKey const minimize_tolerance( "dna:specificity:minimize_tolerance" );  } }
namespace dna { namespace specificity { StringOptionKey const weights_tag( "dna:specificity:weights_tag" );  } }
namespace dna { namespace specificity { StringOptionKey const weights_tag_list( "dna:specificity:weights_tag_list" );  } }
namespace dna { namespace specificity { StringOptionKey const min_type( "dna:specificity:min_type" );  } }
namespace dna { namespace specificity { StringOptionKey const mode( "dna:specificity:mode" );  } }
namespace dna { namespace specificity { StringOptionKey const score_function( "dna:specificity:score_function" );  } }
namespace dna { namespace specificity { BooleanOptionKey const pre_minimize( "dna:specificity:pre_minimize" );  } }
namespace dna { namespace specificity { BooleanOptionKey const post_minimize( "dna:specificity:post_minimize" );  } }
namespace dna { namespace specificity { BooleanOptionKey const pre_pack( "dna:specificity:pre_pack" );  } }
namespace dna { namespace specificity { IntegerOptionKey const nloop( "dna:specificity:nloop" );  } }
namespace dna { namespace specificity { IntegerOptionKey const n_inner( "dna:specificity:n_inner" );  } }
namespace dna { namespace specificity { IntegerOptionKey const n_outer( "dna:specificity:n_outer" );  } }
namespace dna { namespace specificity { IntegerOptionKey const nstep_water( "dna:specificity:nstep_water" );  } }
namespace dna { namespace specificity { IntegerOptionKey const moving_jump( "dna:specificity:moving_jump" );  } }
namespace dna { namespace specificity { IntegerOptionKey const motif_begin( "dna:specificity:motif_begin" );  } }
namespace dna { namespace specificity { IntegerOptionKey const motif_size( "dna:specificity:motif_size" );  } }
namespace dna { namespace specificity { StringVectorOptionKey const pdb_pos( "dna:specificity:pdb_pos" );  } }
namespace dna { namespace specificity { StringVectorOptionKey const methylate( "dna:specificity:methylate" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const dna_backbone_torsion_sdevs( "dna:specificity:dna_backbone_torsion_sdevs" );  } }
namespace dna { namespace specificity { RealOptionKey const dna_sugar_torsion_sdev( "dna:specificity:dna_sugar_torsion_sdev" );  } }
namespace dna { namespace specificity { RealOptionKey const dna_chi_torsion_sdev( "dna:specificity:dna_chi_torsion_sdev" );  } }
namespace dna { namespace specificity { StringOptionKey const lk_ball_wtd_tag( "dna:specificity:lk_ball_wtd_tag" );  } }
namespace dna { namespace specificity { BooleanOptionKey const lk_ball_for_bb( "dna:specificity:lk_ball_for_bb" );  } }
namespace dna { namespace specificity { RealOptionKey const lk_ball_ramp_width_A2( "dna:specificity:lk_ball_ramp_width_A2" );  } }
namespace dna { namespace specificity { RealOptionKey const lk_ball_overlap_width_A2( "dna:specificity:lk_ball_overlap_width_A2" );  } }
namespace dna { namespace specificity { RealOptionKey const lk_ball_water_fade( "dna:specificity:lk_ball_water_fade" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const lk_ball_wtd_prefactors( "dna:specificity:lk_ball_wtd_prefactors" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const lk_ball_waters_sp2( "dna:specificity:lk_ball_waters_sp2" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const lk_ball_waters_sp3( "dna:specificity:lk_ball_waters_sp3" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const lk_ball_waters_ring( "dna:specificity:lk_ball_waters_ring" );  } }
namespace dna { namespace specificity { RealOptionKey const lk_ball_waters_donor( "dna:specificity:lk_ball_waters_donor" );  } }
namespace dna { namespace specificity { RealOptionKey const lk_ball_bridge_angle_widthscale( "dna:specificity:lk_ball_bridge_angle_widthscale" );  } }
namespace dna { namespace design { BooleanOptionKey const design( "dna:design" );  } }
namespace dna { namespace design { BooleanOptionKey const output_unbound_pdb( "dna:design:output_unbound_pdb" );  } }
namespace dna { namespace design { RealOptionKey const z_cutoff( "dna:design:z_cutoff" );  } }
namespace dna { namespace design { StringOptionKey const protein_scan( "dna:design:protein_scan" );  } }
namespace dna { namespace design { StringOptionKey const checkpoint( "dna:design:checkpoint" );  } }
namespace dna { namespace design { BooleanOptionKey const minimize( "dna:design:minimize" );  } }
namespace dna { namespace design { StringVectorOptionKey const dna_defs( "dna:design:dna_defs" );  } }
namespace dna { namespace design { StringOptionKey const dna_defs_file( "dna:design:dna_defs_file" );  } }
namespace dna { namespace design { BooleanOptionKey const nopdb( "dna:design:nopdb" );  } }
namespace dna { namespace design { BooleanOptionKey const designable_second_shell( "dna:design:designable_second_shell" );  } }
namespace dna { namespace design { BooleanOptionKey const base_contacts_only( "dna:design:base_contacts_only" );  } }
namespace dna { namespace design { IntegerOptionKey const probe_specificity( "dna:design:probe_specificity" );  } }
namespace dna { namespace design { BooleanOptionKey const reversion_scan( "dna:design:reversion_scan" );  } }
namespace dna { namespace design { BooleanOptionKey const binding( "dna:design:binding" );  } }
namespace dna { namespace design { RealOptionKey const Boltz_temp( "dna:design:Boltz_temp" );  } }
namespace dna { namespace design { BooleanOptionKey const repack_only( "dna:design:repack_only" );  } }
namespace dna { namespace design { BooleanOptionKey const sparse_pdb_output( "dna:design:sparse_pdb_output" );  } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const specificity( "dna:design:specificity" );  } } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const output_structures( "dna:design:specificity:output_structures" );  } } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const include_dna_potentials( "dna:design:specificity:include_dna_potentials" );  } } }
namespace dna { namespace design { namespace reversion { BooleanOptionKey const reversion( "dna:design:reversion" );  } } }
namespace dna { namespace design { namespace reversion { RealOptionKey const dscore_cutoff( "dna:design:reversion:dscore_cutoff" );  } } }
namespace dna { namespace design { namespace reversion { RealOptionKey const dspec_cutoff( "dna:design:reversion:dspec_cutoff" );  } } }
namespace docking { BooleanOptionKey const kick_relax( "docking:kick_relax" );  }
namespace docking { BooleanOptionKey const docking( "docking" );  }
namespace docking { BooleanOptionKey const view( "docking:view" );  }
namespace docking { BooleanOptionKey const no_filters( "docking:no_filters" );  }
namespace docking { StringVectorOptionKey const design_chains( "docking:design_chains" );  }
namespace docking { FileOptionKey const recover_sidechains( "docking:recover_sidechains" );  }
namespace docking { StringOptionKey const partners( "docking:partners" );  }
namespace docking { BooleanOptionKey const docking_local_refine( "docking:docking_local_refine" );  }
namespace docking { BooleanOptionKey const low_res_protocol_only( "docking:low_res_protocol_only" );  }
namespace docking { BooleanOptionKey const randomize1( "docking:randomize1" );  }
namespace docking { BooleanOptionKey const randomize2( "docking:randomize2" );  }
namespace docking { BooleanOptionKey const use_ellipsoidal_randomization( "docking:use_ellipsoidal_randomization" );  }
namespace docking { BooleanOptionKey const spin( "docking:spin" );  }
namespace docking { RealVectorOptionKey const tilt( "docking:tilt" );  }
namespace docking { StringOptionKey const tilt1_center( "docking:tilt1_center" );  }
namespace docking { StringOptionKey const tilt2_center( "docking:tilt2_center" );  }
namespace docking { RealVectorOptionKey const dock_pert( "docking:dock_pert" );  }
namespace docking { RealOptionKey const uniform_trans( "docking:uniform_trans" );  }
namespace docking { BooleanOptionKey const center_at_interface( "docking:center_at_interface" );  }
namespace docking { IntegerOptionKey const dock_mcm_first_cycles( "docking:dock_mcm_first_cycles" );  }
namespace docking { IntegerOptionKey const dock_mcm_second_cycles( "docking:dock_mcm_second_cycles" );  }
namespace docking { IntegerOptionKey const docking_centroid_outer_cycles( "docking:docking_centroid_outer_cycles" );  }
namespace docking { IntegerOptionKey const docking_centroid_inner_cycles( "docking:docking_centroid_inner_cycles" );  }
namespace docking { BooleanOptionKey const dock_min( "docking:dock_min" );  }
namespace docking { StringOptionKey const flexible_bb_docking( "docking:flexible_bb_docking" );  }
namespace docking { RealOptionKey const flexible_bb_docking_interface_dist( "docking:flexible_bb_docking_interface_dist" );  }
namespace docking { StringOptionKey const ensemble1( "docking:ensemble1" );  }
namespace docking { StringOptionKey const ensemble2( "docking:ensemble2" );  }
namespace docking { RealOptionKey const dock_mcm_trans_magnitude( "docking:dock_mcm_trans_magnitude" );  }
namespace docking { RealOptionKey const dock_mcm_rot_magnitude( "docking:dock_mcm_rot_magnitude" );  }
namespace docking { RealOptionKey const minimization_threshold( "docking:minimization_threshold" );  }
namespace docking { RealOptionKey const temperature( "docking:temperature" );  }
namespace docking { IntegerOptionKey const repack_period( "docking:repack_period" );  }
namespace docking { BooleanOptionKey const extra_rottrial( "docking:extra_rottrial" );  }
namespace docking { BooleanOptionKey const dock_rtmin( "docking:dock_rtmin" );  }
namespace docking { BooleanOptionKey const sc_min( "docking:sc_min" );  }
namespace docking { BooleanOptionKey const norepack1( "docking:norepack1" );  }
namespace docking { BooleanOptionKey const norepack2( "docking:norepack2" );  }
namespace docking { IntegerVectorOptionKey const bb_min_res( "docking:bb_min_res" );  }
namespace docking { IntegerVectorOptionKey const sc_min_res( "docking:sc_min_res" );  }
namespace docking { BooleanOptionKey const dock_ppk( "docking:dock_ppk" );  }
namespace docking { IntegerOptionKey const max_repeats( "docking:max_repeats" );  }
namespace docking { RealVectorOptionKey const dock_lowres_filter( "docking:dock_lowres_filter" );  }
namespace docking { IntegerVectorOptionKey const multibody( "docking:multibody" );  }
namespace docking { BooleanOptionKey const ignore_default_docking_task( "docking:ignore_default_docking_task" );  }
namespace docking { StringOptionKey const low_patch( "docking:low_patch" );  }
namespace docking { StringOptionKey const high_patch( "docking:high_patch" );  }
namespace docking { StringOptionKey const high_min_patch( "docking:high_min_patch" );  }
namespace docking { StringOptionKey const pack_patch( "docking:pack_patch" );  }
namespace docking { BooleanOptionKey const use_legacy_protocol( "docking:use_legacy_protocol" );  }
namespace docking { RealOptionKey const docklowres_trans_magnitude( "docking:docklowres_trans_magnitude" );  }
namespace docking { RealOptionKey const docklowres_rot_magnitude( "docking:docklowres_rot_magnitude" );  }
namespace docking { namespace ligand { BooleanOptionKey const ligand( "docking:ligand" );  } }
namespace docking { namespace ligand { StringOptionKey const protocol( "docking:ligand:protocol" );  } }
namespace docking { namespace ligand { BooleanOptionKey const soft_rep( "docking:ligand:soft_rep" );  } }
namespace docking { namespace ligand { BooleanOptionKey const tweak_sxfn( "docking:ligand:tweak_sxfn" );  } }
namespace docking { namespace ligand { BooleanOptionKey const old_estat( "docking:ligand:old_estat" );  } }
namespace docking { namespace ligand { BooleanOptionKey const random_conformer( "docking:ligand:random_conformer" );  } }
namespace docking { namespace ligand { IntegerOptionKey const improve_orientation( "docking:ligand:improve_orientation" );  } }
namespace docking { namespace ligand { BooleanOptionKey const mutate_same_name3( "docking:ligand:mutate_same_name3" );  } }
namespace docking { namespace ligand { RealOptionKey const subset_to_keep( "docking:ligand:subset_to_keep" );  } }
namespace docking { namespace ligand { RealOptionKey const min_rms( "docking:ligand:min_rms" );  } }
namespace docking { namespace ligand { IntegerOptionKey const max_poses( "docking:ligand:max_poses" );  } }
namespace docking { namespace ligand { BooleanOptionKey const minimize_ligand( "docking:ligand:minimize_ligand" );  } }
namespace docking { namespace ligand { RealOptionKey const harmonic_torsions( "docking:ligand:harmonic_torsions" );  } }
namespace docking { namespace ligand { BooleanOptionKey const use_ambig_constraints( "docking:ligand:use_ambig_constraints" );  } }
namespace docking { namespace ligand { IntegerOptionKey const shear_moves( "docking:ligand:shear_moves" );  } }
namespace docking { namespace ligand { BooleanOptionKey const minimize_backbone( "docking:ligand:minimize_backbone" );  } }
namespace docking { namespace ligand { RealOptionKey const harmonic_Calphas( "docking:ligand:harmonic_Calphas" );  } }
namespace docking { namespace ligand { RealOptionKey const tether_ligand( "docking:ligand:tether_ligand" );  } }
namespace docking { namespace ligand { RealVectorOptionKey const start_from( "docking:ligand:start_from" );  } }
namespace docking { namespace ligand { StringOptionKey const option_file( "docking:ligand:option_file" );  } }
namespace docking { namespace ligand { RealOptionKey const ligand_ensemble( "docking:ligand:ligand_ensemble" );  } }
namespace docking { namespace ligand { namespace grid { BooleanOptionKey const grid( "docking:ligand:grid" );  } } }
namespace docking { namespace ligand { namespace grid { FileOptionKey const grid_kin( "docking:ligand:grid:grid_kin" );  } } }
namespace docking { namespace ligand { namespace grid { FileOptionKey const grid_map( "docking:ligand:grid:grid_map" );  } } }
namespace DomainAssembly { BooleanOptionKey const DomainAssembly( "DomainAssembly" );  }
namespace DomainAssembly { BooleanOptionKey const da_setup( "DomainAssembly:da_setup" );  }
namespace DomainAssembly { FileOptionKey const da_setup_option_file( "DomainAssembly:da_setup_option_file" );  }
namespace DomainAssembly { FileOptionKey const da_setup_output_pdb( "DomainAssembly:da_setup_output_pdb" );  }
namespace DomainAssembly { FileOptionKey const da_linker_file( "DomainAssembly:da_linker_file" );  }
namespace DomainAssembly { FileOptionKey const da_require_buried( "DomainAssembly:da_require_buried" );  }
namespace DomainAssembly { FileOptionKey const da_start_pdb( "DomainAssembly:da_start_pdb" );  }
namespace DomainAssembly { BooleanOptionKey const run_fullatom( "DomainAssembly:run_fullatom" );  }
namespace DomainAssembly { BooleanOptionKey const run_centroid( "DomainAssembly:run_centroid" );  }
namespace DomainAssembly { BooleanOptionKey const run_centroid_abinitio( "DomainAssembly:run_centroid_abinitio" );  }
