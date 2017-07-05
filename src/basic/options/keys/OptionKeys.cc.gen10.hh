namespace matdes { namespace design { BooleanOptionKey const design( "matdes:design" );  } }
namespace matdes { namespace design { RealOptionKey const contact_dist( "matdes:design:contact_dist" );  } }
namespace matdes { namespace design { RealOptionKey const grid_size_angle( "matdes:design:grid_size_angle" );  } }
namespace matdes { namespace design { RealOptionKey const grid_size_radius( "matdes:design:grid_size_radius" );  } }
namespace matdes { namespace design { IntegerOptionKey const grid_nsamp_angle( "matdes:design:grid_nsamp_angle" );  } }
namespace matdes { namespace design { IntegerOptionKey const grid_nsamp_radius( "matdes:design:grid_nsamp_radius" );  } }
namespace matdes { namespace design { RealOptionKey const fav_nat_bonus( "matdes:design:fav_nat_bonus" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const mutalyze( "matdes:mutalyze" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const calc_rot_boltz( "matdes:mutalyze:calc_rot_boltz" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const ala_scan( "matdes:mutalyze:ala_scan" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const revert_scan( "matdes:mutalyze:revert_scan" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const min_rb( "matdes:mutalyze:min_rb" );  } }
namespace mc { BooleanOptionKey const mc( "mc" );  }
namespace mc { BooleanOptionKey const log_scores_in_MC( "mc:log_scores_in_MC" );  }
namespace mc { StringOptionKey const hierarchical_pool( "mc:hierarchical_pool" );  }
namespace mc { FileOptionKey const read_structures_into_pool( "mc:read_structures_into_pool" );  }
namespace mc { IntegerOptionKey const convergence_check_frequency( "mc:convergence_check_frequency" );  }
namespace mc { FileOptionKey const known_structures( "mc:known_structures" );  }
namespace mc { RealOptionKey const max_rmsd_against_known_structures( "mc:max_rmsd_against_known_structures" );  }
namespace mc { IntegerVectorOptionKey const excluded_residues_from_rmsd( "mc:excluded_residues_from_rmsd" );  }
namespace mc { IntegerOptionKey const heat_convergence_check( "mc:heat_convergence_check" );  }
namespace mh { BooleanOptionKey const mh( "mh" );  }
namespace mh { StringOptionKey const motif_out_file( "mh:motif_out_file" );  }
namespace mh { FileVectorOptionKey const harvest_motifs( "mh:harvest_motifs" );  }
namespace mh { FileVectorOptionKey const print_motifs( "mh:print_motifs" );  }
namespace mh { FileVectorOptionKey const remove_duplicates( "mh:remove_duplicates" );  }
namespace mh { FileVectorOptionKey const dump_motif_pdbs( "mh:dump_motif_pdbs" );  }
namespace mh { FileVectorOptionKey const merge_motifs( "mh:merge_motifs" );  }
namespace mh { FileVectorOptionKey const merge_scores( "mh:merge_scores" );  }
namespace mh { BooleanOptionKey const merge_motifs_one_per_bin( "mh:merge_motifs_one_per_bin" );  }
namespace mh { BooleanOptionKey const gen_reverse_motifs_on_load( "mh:gen_reverse_motifs_on_load" );  }
namespace mh { FileVectorOptionKey const dump_input_pdb( "mh:dump_input_pdb" );  }
namespace mh { FileVectorOptionKey const score_pdbs( "mh:score_pdbs" );  }
namespace mh { FileVectorOptionKey const sequence_recovery( "mh:sequence_recovery" );  }
namespace mh { FileVectorOptionKey const explicit_motif_score( "mh:explicit_motif_score" );  }
namespace mh { FileVectorOptionKey const harvest_scores( "mh:harvest_scores" );  }
namespace mh { FileOptionKey const print_scores( "mh:print_scores" );  }
namespace mh { FileVectorOptionKey const dump_matching_motifs( "mh:dump_matching_motifs" );  }
namespace mh { BooleanOptionKey const score_across_chains_only( "mh:score_across_chains_only" );  }
namespace mh { BooleanOptionKey const normalize_score_ncontact( "mh:normalize_score_ncontact" );  }
namespace mh { IntegerOptionKey const harvest_motifs_min_hh_ends( "mh:harvest_motifs_min_hh_ends" );  }
namespace mh { BooleanOptionKey const ignore_io_errors( "mh:ignore_io_errors" );  }
namespace mh { RealOptionKey const motif_match_radius( "mh:motif_match_radius" );  }
namespace mh { RealVectorOptionKey const merge_similar_motifs( "mh:merge_similar_motifs" );  }
namespace mh { namespace score { BooleanOptionKey const score( "mh:score" );  } }
namespace mh { namespace score { RealOptionKey const background_weight( "mh:score:background_weight" );  } }
namespace mh { namespace score { RealOptionKey const ca_cb_clash_weight( "mh:score:ca_cb_clash_weight" );  } }
namespace mh { namespace score { BooleanOptionKey const noloops( "mh:score:noloops" );  } }
namespace mh { namespace score { BooleanOptionKey const nosheets( "mh:score:nosheets" );  } }
namespace mh { namespace score { BooleanOptionKey const nohelix( "mh:score:nohelix" );  } }
namespace mh { namespace score { BooleanOptionKey const spread_ss_element( "mh:score:spread_ss_element" );  } }
namespace mh { namespace score { RealOptionKey const min_cover_fraction( "mh:score:min_cover_fraction" );  } }
namespace mh { namespace score { RealOptionKey const strand_pair_weight( "mh:score:strand_pair_weight" );  } }
namespace mh { namespace score { RealOptionKey const anti_polar_weight( "mh:score:anti_polar_weight" );  } }
namespace mh { namespace score { RealOptionKey const min_contact_pairs( "mh:score:min_contact_pairs" );  } }
namespace mh { namespace score { RealOptionKey const max_contact_pairs( "mh:score:max_contact_pairs" );  } }
namespace mh { namespace score { RealOptionKey const max_cb_dis( "mh:score:max_cb_dis" );  } }
namespace mh { namespace score { RealOptionKey const coverage_pow( "mh:score:coverage_pow" );  } }
namespace mh { namespace score { BooleanOptionKey const use_ss1( "mh:score:use_ss1" );  } }
namespace mh { namespace score { BooleanOptionKey const use_ss2( "mh:score:use_ss2" );  } }
namespace mh { namespace score { BooleanOptionKey const use_aa1( "mh:score:use_aa1" );  } }
namespace mh { namespace score { BooleanOptionKey const use_aa2( "mh:score:use_aa2" );  } }
namespace mh { namespace score { BooleanOptionKey const use_log( "mh:score:use_log" );  } }
namespace mh { namespace path { BooleanOptionKey const path( "mh:path" );  } }
namespace mh { namespace path { StringVectorOptionKey const biounit( "mh:path:biounit" );  } }
namespace mh { namespace path { StringVectorOptionKey const biounit_ideal( "mh:path:biounit_ideal" );  } }
namespace mh { namespace path { StringVectorOptionKey const pdb( "mh:path:pdb" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs( "mh:path:motifs" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_SC_SC( "mh:path:motifs_SC_SC" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_SC_BB( "mh:path:motifs_SC_BB" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_BB_BB( "mh:path:motifs_BB_BB" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_BB_PH( "mh:path:motifs_BB_PH" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_BB_PO( "mh:path:motifs_BB_PO" );  } }
namespace mh { namespace path { StringVectorOptionKey const motifs_PH_PO( "mh:path:motifs_PH_PO" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores( "mh:path:scores" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_SC_SC( "mh:path:scores_SC_SC" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_SC_BB( "mh:path:scores_SC_BB" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_BB_BB( "mh:path:scores_BB_BB" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_BB_PH( "mh:path:scores_BB_PH" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_BB_PO( "mh:path:scores_BB_PO" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_PH_PO( "mh:path:scores_PH_PO" );  } }
namespace mh { namespace path { StringVectorOptionKey const scores_frags( "mh:path:scores_frags" );  } }
namespace mh { namespace harvest { BooleanOptionKey const harvest( "mh:harvest" );  } }
namespace mh { namespace harvest { RealOptionKey const hash_cart_resl( "mh:harvest:hash_cart_resl" );  } }
namespace mh { namespace harvest { RealOptionKey const hash_angle_resl( "mh:harvest:hash_angle_resl" );  } }
namespace mh { namespace harvest { RealOptionKey const smoothing_factor( "mh:harvest:smoothing_factor" );  } }
namespace mh { namespace harvest { BooleanOptionKey const idealize( "mh:harvest:idealize" );  } }
namespace mh { namespace harvest { BooleanOptionKey const dump( "mh:harvest:dump" );  } }
namespace mh { namespace harvest { RealOptionKey const min_bin_val( "mh:harvest:min_bin_val" );  } }
namespace mh { namespace harvest { BooleanOptionKey const sep_aa( "mh:harvest:sep_aa" );  } }
namespace mh { namespace harvest { BooleanOptionKey const sep_aa1( "mh:harvest:sep_aa1" );  } }
namespace mh { namespace harvest { BooleanOptionKey const sep_aa2( "mh:harvest:sep_aa2" );  } }
namespace mh { namespace harvest { BooleanOptionKey const sep_ss( "mh:harvest:sep_ss" );  } }
namespace mh { namespace harvest { BooleanOptionKey const sep_dssp( "mh:harvest:sep_dssp" );  } }
namespace mh { namespace harvest { RealVectorOptionKey const sep_lj( "mh:harvest:sep_lj" );  } }
namespace mh { namespace harvest { RealVectorOptionKey const sep_hb( "mh:harvest:sep_hb" );  } }
namespace mh { namespace harvest { RealVectorOptionKey const sep_nbrs( "mh:harvest:sep_nbrs" );  } }
namespace mh { namespace harvest { RealVectorOptionKey const sep_bfac( "mh:harvest:sep_bfac" );  } }
namespace mh { namespace harvest { RealVectorOptionKey const sep_dist( "mh:harvest:sep_dist" );  } }
namespace mh { namespace harvest { BooleanOptionKey const weight_by_energy( "mh:harvest:weight_by_energy" );  } }
namespace mh { namespace harvest { RealOptionKey const max_rmsd( "mh:harvest:max_rmsd" );  } }
namespace mh { namespace harvest { IntegerOptionKey const max_res( "mh:harvest:max_res" );  } }
namespace mh { namespace harvest { BooleanOptionKey const agg_with_max( "mh:harvest:agg_with_max" );  } }
namespace mh { namespace harvest { RealOptionKey const multiplier( "mh:harvest:multiplier" );  } }
namespace mh { namespace match { BooleanOptionKey const match( "mh:match" );  } }
namespace mh { namespace match { BooleanOptionKey const interface_only( "mh:match:interface_only" );  } }
namespace mh { namespace match { BooleanOptionKey const ss( "mh:match:ss" );  } }
namespace mh { namespace match { BooleanOptionKey const ss1( "mh:match:ss1" );  } }
namespace mh { namespace match { BooleanOptionKey const ss2( "mh:match:ss2" );  } }
namespace mh { namespace match { BooleanOptionKey const aa( "mh:match:aa" );  } }
namespace mh { namespace match { BooleanOptionKey const aa1( "mh:match:aa1" );  } }
namespace mh { namespace match { BooleanOptionKey const aa2( "mh:match:aa2" );  } }
namespace mh { namespace dump { BooleanOptionKey const dump( "mh:dump" );  } }
namespace mh { namespace dump { IntegerOptionKey const limit_per_pair( "mh:dump:limit_per_pair" );  } }
namespace mh { namespace dump { IntegerOptionKey const max_per_res( "mh:dump:max_per_res" );  } }
namespace mh { namespace dump { RealOptionKey const max_ca_dis( "mh:dump:max_ca_dis" );  } }
namespace mh { namespace dump { RealOptionKey const max_rms( "mh:dump:max_rms" );  } }
namespace mh { namespace dump { RealOptionKey const resfile_min_pair_score( "mh:dump:resfile_min_pair_score" );  } }
namespace mh { namespace dump { RealOptionKey const resfile_min_tot_score( "mh:dump:resfile_min_tot_score" );  } }
namespace mh { namespace dump { BooleanOptionKey const resfile_dump( "mh:dump:resfile_dump" );  } }
namespace mh { namespace dump { BooleanOptionKey const symmetric_motifs( "mh:dump:symmetric_motifs" );  } }
namespace mh { namespace filter { BooleanOptionKey const filter( "mh:filter" );  } }
namespace mh { namespace filter { BooleanOptionKey const filter_harvest( "mh:filter:filter_harvest" );  } }
namespace mh { namespace filter { BooleanOptionKey const filter_io( "mh:filter:filter_io" );  } }
namespace mh { namespace filter { StringOptionKey const pdb( "mh:filter:pdb" );  } }
namespace mh { namespace filter { StringOptionKey const lig( "mh:filter:lig" );  } }
namespace mh { namespace filter { StringOptionKey const motif_type( "mh:filter:motif_type" );  } }
namespace mh { namespace filter { StringOptionKey const restype1( "mh:filter:restype1" );  } }
namespace mh { namespace filter { StringOptionKey const restype2( "mh:filter:restype2" );  } }
namespace mh { namespace filter { StringOptionKey const restype( "mh:filter:restype" );  } }
namespace mh { namespace filter { StringOptionKey const restype_one( "mh:filter:restype_one" );  } }
namespace mh { namespace filter { StringOptionKey const not_restype( "mh:filter:not_restype" );  } }
namespace mh { namespace filter { StringOptionKey const not_restype_one( "mh:filter:not_restype_one" );  } }
namespace mh { namespace filter { IntegerOptionKey const seqsep( "mh:filter:seqsep" );  } }
namespace mh { namespace filter { IntegerOptionKey const max_seqsep( "mh:filter:max_seqsep" );  } }
namespace mh { namespace filter { BooleanOptionKey const no_hb_bb( "mh:filter:no_hb_bb" );  } }
namespace mh { namespace filter { RealOptionKey const mindist2( "mh:filter:mindist2" );  } }
namespace mh { namespace filter { RealOptionKey const maxdist2( "mh:filter:maxdist2" );  } }
namespace mh { namespace filter { StringOptionKey const ss1( "mh:filter:ss1" );  } }
namespace mh { namespace filter { StringOptionKey const ss2( "mh:filter:ss2" );  } }
namespace mh { namespace filter { StringOptionKey const dssp1( "mh:filter:dssp1" );  } }
namespace mh { namespace filter { StringOptionKey const dssp2( "mh:filter:dssp2" );  } }
namespace mh { namespace filter { StringOptionKey const aa1( "mh:filter:aa1" );  } }
namespace mh { namespace filter { StringOptionKey const aa2( "mh:filter:aa2" );  } }
namespace mh { namespace filter { RealOptionKey const sasa( "mh:filter:sasa" );  } }
namespace mh { namespace filter { RealOptionKey const faatr( "mh:filter:faatr" );  } }
namespace mh { namespace filter { RealOptionKey const hb_sc( "mh:filter:hb_sc" );  } }
namespace mh { namespace filter { RealOptionKey const hb_bb_sc( "mh:filter:hb_bb_sc" );  } }
namespace mh { namespace filter { RealOptionKey const hb_bb( "mh:filter:hb_bb" );  } }
namespace mh { namespace filter { RealOptionKey const occupancy( "mh:filter:occupancy" );  } }
namespace mh { namespace filter { RealOptionKey const coorderr( "mh:filter:coorderr" );  } }
namespace mh { namespace filter { BooleanOptionKey const uniformfrag( "mh:filter:uniformfrag" );  } }
namespace mh { namespace filter { RealOptionKey const faatr_or_hbbb( "mh:filter:faatr_or_hbbb" );  } }
namespace mh { namespace filter { RealOptionKey const faatr_or_hb( "mh:filter:faatr_or_hb" );  } }
namespace mh { namespace filter { BooleanOptionKey const noloops( "mh:filter:noloops" );  } }
namespace mh { namespace filter { BooleanOptionKey const oneloop( "mh:filter:oneloop" );  } }
namespace mh { namespace filter { BooleanOptionKey const nodisulf( "mh:filter:nodisulf" );  } }
namespace mh { namespace filter { RealOptionKey const score( "mh:filter:score" );  } }
namespace magnesium { BooleanOptionKey const magnesium( "magnesium" );  }
namespace magnesium { BooleanOptionKey const scan( "magnesium:scan" );  }
namespace magnesium { IntegerVectorOptionKey const mg_res( "magnesium:mg_res" );  }
namespace magnesium { BooleanOptionKey const minimize_during_scoring( "magnesium:minimize_during_scoring" );  }
namespace magnesium { ResidueChainVectorOptionKey const ligand_res( "magnesium:ligand_res" );  }
namespace magnesium { IntegerVectorOptionKey const pose_ligand_res( "magnesium:pose_ligand_res" );  }
namespace magnesium { BooleanOptionKey const lores_scan( "magnesium:lores_scan" );  }
namespace magnesium { RealOptionKey const xyz_step( "magnesium:xyz_step" );  }
namespace magnesium { RealOptionKey const score_cut( "magnesium:score_cut" );  }
namespace magnesium { RealOptionKey const score_cut_PDB( "magnesium:score_cut_PDB" );  }
namespace magnesium { BooleanOptionKey const integration_test( "magnesium:integration_test" );  }
namespace magnesium { BooleanOptionKey const tether_to_closest_res( "magnesium:tether_to_closest_res" );  }
namespace magnesium { BooleanOptionKey const fixup( "magnesium:fixup" );  }
namespace magnesium { BooleanOptionKey const pack_water_hydrogens( "magnesium:pack_water_hydrogens" );  }
namespace magnesium { BooleanOptionKey const hydrate( "magnesium:hydrate" );  }
namespace magnesium { BooleanOptionKey const monte_carlo( "magnesium:monte_carlo" );  }
namespace magnesium { BooleanOptionKey const scored_hydrogen_sampling( "magnesium:scored_hydrogen_sampling" );  }
namespace magnesium { BooleanOptionKey const all_hydration_frames( "magnesium:all_hydration_frames" );  }
namespace magnesium { BooleanOptionKey const leave_other_waters( "magnesium:leave_other_waters" );  }
namespace magnesium { BooleanOptionKey const minimize( "magnesium:minimize" );  }
namespace magnesium { RealOptionKey const minimize_mg_coord_constraint_distance( "magnesium:minimize_mg_coord_constraint_distance" );  }
namespace magnesium { namespace montecarlo { BooleanOptionKey const montecarlo( "magnesium:montecarlo" );  } }
namespace magnesium { namespace montecarlo { RealOptionKey const temperature( "magnesium:montecarlo:temperature" );  } }
namespace magnesium { namespace montecarlo { IntegerOptionKey const cycles( "magnesium:montecarlo:cycles" );  } }
namespace magnesium { namespace montecarlo { BooleanOptionKey const dump( "magnesium:montecarlo:dump" );  } }
namespace magnesium { namespace montecarlo { RealOptionKey const add_delete_frequency( "magnesium:montecarlo:add_delete_frequency" );  } }
namespace motifs { BooleanOptionKey const motifs( "motifs" );  }
namespace motifs { RealOptionKey const close_enough( "motifs:close_enough" );  }
namespace motifs { IntegerOptionKey const max_depth( "motifs:max_depth" );  }
namespace motifs { BooleanOptionKey const keep_motif_xtal_location( "motifs:keep_motif_xtal_location" );  }
namespace motifs { RealOptionKey const pack_score_cutoff( "motifs:pack_score_cutoff" );  }
namespace motifs { RealOptionKey const hb_score_cutoff( "motifs:hb_score_cutoff" );  }
namespace motifs { RealOptionKey const water_score_cutoff( "motifs:water_score_cutoff" );  }
namespace motifs { RealOptionKey const pack_min_threshold( "motifs:pack_min_threshold" );  }
namespace motifs { RealOptionKey const pack_max_threshold( "motifs:pack_max_threshold" );  }
namespace motifs { RealOptionKey const hbond_min_threshold( "motifs:hbond_min_threshold" );  }
namespace motifs { RealOptionKey const hbond_max_threshold( "motifs:hbond_max_threshold" );  }
namespace motifs { RealOptionKey const elec_min_threshold( "motifs:elec_min_threshold" );  }
namespace motifs { RealOptionKey const elec_max_threshold( "motifs:elec_max_threshold" );  }
namespace motifs { RealOptionKey const duplicate_dist_cutoff( "motifs:duplicate_dist_cutoff" );  }
namespace motifs { RealOptionKey const duplicate_angle_cutoff( "motifs:duplicate_angle_cutoff" );  }
namespace motifs { StringOptionKey const motif_output_directory( "motifs:motif_output_directory" );  }
namespace motifs { BooleanOptionKey const eliminate_weak_motifs( "motifs:eliminate_weak_motifs" );  }
namespace motifs { RealOptionKey const duplicate_motif_cutoff( "motifs:duplicate_motif_cutoff" );  }
namespace motifs { BooleanOptionKey const preminimize_motif_pdbs( "motifs:preminimize_motif_pdbs" );  }
namespace motifs { BooleanOptionKey const preminimize_motif_pdbs_sconly( "motifs:preminimize_motif_pdbs_sconly" );  }
namespace motifs { BooleanOptionKey const place_adduct_waters( "motifs:place_adduct_waters" );  }
namespace motifs { FileVectorOptionKey const list_motifs( "motifs:list_motifs" );  }
namespace motifs { StringOptionKey const motif_filename( "motifs:motif_filename" );  }
namespace motifs { StringOptionKey const file_prefix( "motifs:file_prefix" );  }
namespace motifs { StringOptionKey const build_residue_file( "motifs:build_residue_file" );  }
namespace motifs { StringOptionKey const motif_flexible_loop_file( "motifs:motif_flexible_loop_file" );  }
namespace motifs { StringOptionKey const residue_trim_file( "motifs:residue_trim_file" );  }
namespace motifs { StringOptionKey const BPData( "motifs:BPData" );  }
namespace motifs { FileVectorOptionKey const list_dnaconformers( "motifs:list_dnaconformers" );  }
namespace motifs { StringVectorOptionKey const target_dna_defs( "motifs:target_dna_defs" );  }
namespace motifs { StringVectorOptionKey const motif_build_defs( "motifs:motif_build_defs" );  }
namespace motifs { IntegerVectorOptionKey const motif_build_positions( "motifs:motif_build_positions" );  }
namespace motifs { RealOptionKey const r1( "motifs:r1" );  }
namespace motifs { RealOptionKey const r2( "motifs:r2" );  }
namespace motifs { RealOptionKey const z1( "motifs:z1" );  }
namespace motifs { RealOptionKey const z2( "motifs:z2" );  }
namespace motifs { RealOptionKey const dtest( "motifs:dtest" );  }
namespace motifs { IntegerOptionKey const rotlevel( "motifs:rotlevel" );  }
namespace motifs { IntegerOptionKey const num_repacks( "motifs:num_repacks" );  }
namespace motifs { BooleanOptionKey const minimize( "motifs:minimize" );  }
namespace motifs { BooleanOptionKey const minimize_dna( "motifs:minimize_dna" );  }
namespace motifs { BooleanOptionKey const run_motifs( "motifs:run_motifs" );  }
namespace motifs { BooleanOptionKey const expand_motifs( "motifs:expand_motifs" );  }
namespace motifs { BooleanOptionKey const aromatic_motifs( "motifs:aromatic_motifs" );  }
namespace motifs { BooleanOptionKey const dump_motifs( "motifs:dump_motifs" );  }
namespace motifs { BooleanOptionKey const quick_and_dirty( "motifs:quick_and_dirty" );  }
namespace motifs { RealOptionKey const special_rotweight( "motifs:special_rotweight" );  }
namespace motifs { StringOptionKey const output_file( "motifs:output_file" );  }
namespace motifs { StringOptionKey const data_file( "motifs:data_file" );  }
namespace motifs { StringOptionKey const target_aa( "motifs:target_aa" );  }
namespace motifs { BooleanOptionKey const flex_sugar( "motifs:flex_sugar" );  }
namespace motifs { BooleanOptionKey const clear_bprots( "motifs:clear_bprots" );  }
namespace motifs { IntegerOptionKey const rots2add( "motifs:rots2add" );  }
namespace motifs { BooleanOptionKey const restrict_to_wt( "motifs:restrict_to_wt" );  }
namespace motifs { BooleanOptionKey const rerun_motifsearch( "motifs:rerun_motifsearch" );  }
namespace motifs { BooleanOptionKey const no_rotamer_bump( "motifs:no_rotamer_bump" );  }
namespace motifs { RealOptionKey const ligand_motif_sphere( "motifs:ligand_motif_sphere" );  }
namespace ms { BooleanOptionKey const ms( "ms" );  }
namespace ms { IntegerOptionKey const pop_from_ss( "ms:pop_from_ss" );  }
namespace ms { IntegerOptionKey const pop_size( "ms:pop_size" );  }
namespace ms { IntegerOptionKey const generations( "ms:generations" );  }
namespace ms { IntegerOptionKey const num_packs( "ms:num_packs" );  }
