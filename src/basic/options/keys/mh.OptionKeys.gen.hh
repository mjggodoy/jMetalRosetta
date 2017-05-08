// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/mh.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mh_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_mh_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mh { extern BooleanOptionKey const mh; }
namespace mh { extern StringOptionKey const motif_out_file; }
namespace mh { extern FileVectorOptionKey const harvest_motifs; }
namespace mh { extern FileVectorOptionKey const print_motifs; }
namespace mh { extern FileVectorOptionKey const remove_duplicates; }
namespace mh { extern FileVectorOptionKey const dump_motif_pdbs; }
namespace mh { extern FileVectorOptionKey const merge_motifs; }
namespace mh { extern FileVectorOptionKey const merge_scores; }
namespace mh { extern BooleanOptionKey const merge_motifs_one_per_bin; }
namespace mh { extern BooleanOptionKey const gen_reverse_motifs_on_load; }
namespace mh { extern FileVectorOptionKey const dump_input_pdb; }
namespace mh { extern FileVectorOptionKey const score_pdbs; }
namespace mh { extern FileVectorOptionKey const sequence_recovery; }
namespace mh { extern FileVectorOptionKey const explicit_motif_score; }
namespace mh { extern FileVectorOptionKey const harvest_scores; }
namespace mh { extern FileOptionKey const print_scores; }
namespace mh { extern FileVectorOptionKey const dump_matching_motifs; }
namespace mh { extern BooleanOptionKey const score_across_chains_only; }
namespace mh { extern BooleanOptionKey const normalize_score_ncontact; }
namespace mh { extern IntegerOptionKey const harvest_motifs_min_hh_ends; }
namespace mh { extern BooleanOptionKey const ignore_io_errors; }
namespace mh { extern RealOptionKey const motif_match_radius; }
namespace mh { extern RealVectorOptionKey const merge_similar_motifs; }
namespace mh { namespace score { extern BooleanOptionKey const score; } }
namespace mh { namespace score { extern RealOptionKey const background_weight; } }
namespace mh { namespace score { extern RealOptionKey const ca_cb_clash_weight; } }
namespace mh { namespace score { extern BooleanOptionKey const noloops; } }
namespace mh { namespace score { extern BooleanOptionKey const nosheets; } }
namespace mh { namespace score { extern BooleanOptionKey const nohelix; } }
namespace mh { namespace score { extern BooleanOptionKey const spread_ss_element; } }
namespace mh { namespace score { extern RealOptionKey const min_cover_fraction; } }
namespace mh { namespace score { extern RealOptionKey const strand_pair_weight; } }
namespace mh { namespace score { extern RealOptionKey const anti_polar_weight; } }
namespace mh { namespace score { extern RealOptionKey const min_contact_pairs; } }
namespace mh { namespace score { extern RealOptionKey const max_contact_pairs; } }
namespace mh { namespace score { extern RealOptionKey const max_cb_dis; } }
namespace mh { namespace score { extern RealOptionKey const coverage_pow; } }
namespace mh { namespace score { extern BooleanOptionKey const use_ss1; } }
namespace mh { namespace score { extern BooleanOptionKey const use_ss2; } }
namespace mh { namespace score { extern BooleanOptionKey const use_aa1; } }
namespace mh { namespace score { extern BooleanOptionKey const use_aa2; } }
namespace mh { namespace score { extern BooleanOptionKey const use_log; } }
namespace mh { namespace path { extern BooleanOptionKey const path; } }
namespace mh { namespace path { extern StringVectorOptionKey const biounit; } }
namespace mh { namespace path { extern StringVectorOptionKey const biounit_ideal; } }
namespace mh { namespace path { extern StringVectorOptionKey const pdb; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_SC_SC; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_SC_BB; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_BB_BB; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_BB_PH; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_BB_PO; } }
namespace mh { namespace path { extern StringVectorOptionKey const motifs_PH_PO; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_SC_SC; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_SC_BB; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_BB_BB; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_BB_PH; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_BB_PO; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_PH_PO; } }
namespace mh { namespace path { extern StringVectorOptionKey const scores_frags; } }
namespace mh { namespace harvest { extern BooleanOptionKey const harvest; } }
namespace mh { namespace harvest { extern RealOptionKey const hash_cart_resl; } }
namespace mh { namespace harvest { extern RealOptionKey const hash_angle_resl; } }
namespace mh { namespace harvest { extern RealOptionKey const smoothing_factor; } }
namespace mh { namespace harvest { extern BooleanOptionKey const idealize; } }
namespace mh { namespace harvest { extern BooleanOptionKey const dump; } }
namespace mh { namespace harvest { extern RealOptionKey const min_bin_val; } }
namespace mh { namespace harvest { extern BooleanOptionKey const sep_aa; } }
namespace mh { namespace harvest { extern BooleanOptionKey const sep_aa1; } }
namespace mh { namespace harvest { extern BooleanOptionKey const sep_aa2; } }
namespace mh { namespace harvest { extern BooleanOptionKey const sep_ss; } }
namespace mh { namespace harvest { extern BooleanOptionKey const sep_dssp; } }
namespace mh { namespace harvest { extern RealVectorOptionKey const sep_lj; } }
namespace mh { namespace harvest { extern RealVectorOptionKey const sep_hb; } }
namespace mh { namespace harvest { extern RealVectorOptionKey const sep_nbrs; } }
namespace mh { namespace harvest { extern RealVectorOptionKey const sep_bfac; } }
namespace mh { namespace harvest { extern RealVectorOptionKey const sep_dist; } }
namespace mh { namespace harvest { extern BooleanOptionKey const weight_by_energy; } }
namespace mh { namespace harvest { extern RealOptionKey const max_rmsd; } }
namespace mh { namespace harvest { extern IntegerOptionKey const max_res; } }
namespace mh { namespace harvest { extern BooleanOptionKey const agg_with_max; } }
namespace mh { namespace harvest { extern RealOptionKey const multiplier; } }
namespace mh { namespace match { extern BooleanOptionKey const match; } }
namespace mh { namespace match { extern BooleanOptionKey const interface_only; } }
namespace mh { namespace match { extern BooleanOptionKey const ss; } }
namespace mh { namespace match { extern BooleanOptionKey const ss1; } }
namespace mh { namespace match { extern BooleanOptionKey const ss2; } }
namespace mh { namespace match { extern BooleanOptionKey const aa; } }
namespace mh { namespace match { extern BooleanOptionKey const aa1; } }
namespace mh { namespace match { extern BooleanOptionKey const aa2; } }
namespace mh { namespace dump { extern BooleanOptionKey const dump; } }
namespace mh { namespace dump { extern IntegerOptionKey const limit_per_pair; } }
namespace mh { namespace dump { extern IntegerOptionKey const max_per_res; } }
namespace mh { namespace dump { extern RealOptionKey const max_ca_dis; } }
namespace mh { namespace dump { extern RealOptionKey const max_rms; } }
namespace mh { namespace dump { extern RealOptionKey const resfile_min_pair_score; } }
namespace mh { namespace dump { extern RealOptionKey const resfile_min_tot_score; } }
namespace mh { namespace dump { extern BooleanOptionKey const resfile_dump; } }
namespace mh { namespace dump { extern BooleanOptionKey const symmetric_motifs; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter_harvest; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter_io; } }
namespace mh { namespace filter { extern StringOptionKey const pdb; } }
namespace mh { namespace filter { extern StringOptionKey const lig; } }
namespace mh { namespace filter { extern StringOptionKey const motif_type; } }
namespace mh { namespace filter { extern StringOptionKey const restype1; } }
namespace mh { namespace filter { extern StringOptionKey const restype2; } }
namespace mh { namespace filter { extern StringOptionKey const restype; } }
namespace mh { namespace filter { extern StringOptionKey const restype_one; } }
namespace mh { namespace filter { extern StringOptionKey const not_restype; } }
namespace mh { namespace filter { extern StringOptionKey const not_restype_one; } }
namespace mh { namespace filter { extern IntegerOptionKey const seqsep; } }
namespace mh { namespace filter { extern IntegerOptionKey const max_seqsep; } }
namespace mh { namespace filter { extern BooleanOptionKey const no_hb_bb; } }
namespace mh { namespace filter { extern RealOptionKey const mindist2; } }
namespace mh { namespace filter { extern RealOptionKey const maxdist2; } }
namespace mh { namespace filter { extern StringOptionKey const ss1; } }
namespace mh { namespace filter { extern StringOptionKey const ss2; } }
namespace mh { namespace filter { extern StringOptionKey const dssp1; } }
namespace mh { namespace filter { extern StringOptionKey const dssp2; } }
namespace mh { namespace filter { extern StringOptionKey const aa1; } }
namespace mh { namespace filter { extern StringOptionKey const aa2; } }
namespace mh { namespace filter { extern RealOptionKey const sasa; } }
namespace mh { namespace filter { extern RealOptionKey const faatr; } }
namespace mh { namespace filter { extern RealOptionKey const hb_sc; } }
namespace mh { namespace filter { extern RealOptionKey const hb_bb_sc; } }
namespace mh { namespace filter { extern RealOptionKey const hb_bb; } }
namespace mh { namespace filter { extern RealOptionKey const occupancy; } }
namespace mh { namespace filter { extern RealOptionKey const coorderr; } }
namespace mh { namespace filter { extern BooleanOptionKey const uniformfrag; } }
namespace mh { namespace filter { extern RealOptionKey const faatr_or_hbbb; } }
namespace mh { namespace filter { extern RealOptionKey const faatr_or_hb; } }
namespace mh { namespace filter { extern BooleanOptionKey const noloops; } }
namespace mh { namespace filter { extern BooleanOptionKey const oneloop; } }
namespace mh { namespace filter { extern BooleanOptionKey const nodisulf; } }
namespace mh { namespace filter { extern RealOptionKey const score; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
