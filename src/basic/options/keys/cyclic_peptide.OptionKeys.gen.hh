// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/cyclic_peptide.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_cyclic_peptide_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_cyclic_peptide_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace cyclic_peptide { extern BooleanOptionKey const cyclic_peptide; }
namespace cyclic_peptide { extern StringOptionKey const rand_checkpoint_file; }
namespace cyclic_peptide { extern StringOptionKey const checkpoint_file; }
namespace cyclic_peptide { extern StringOptionKey const checkpoint_job_identifier; }
namespace cyclic_peptide { extern StringOptionKey const default_rama_sampling_table; }
namespace cyclic_peptide { extern StringVectorOptionKey const rama_sampling_table_by_res; }
namespace cyclic_peptide { extern StringOptionKey const sequence_file; }
namespace cyclic_peptide { extern IntegerOptionKey const genkic_closure_attempts; }
namespace cyclic_peptide { extern IntegerOptionKey const genkic_min_solution_count; }
namespace cyclic_peptide { extern BooleanOptionKey const cyclic_permutations; }
namespace cyclic_peptide { extern BooleanOptionKey const use_rama_filter; }
namespace cyclic_peptide { extern RealOptionKey const rama_cutoff; }
namespace cyclic_peptide { extern RealOptionKey const high_hbond_weight_multiplier; }
namespace cyclic_peptide { extern RealOptionKey const min_genkic_hbonds; }
namespace cyclic_peptide { extern RealOptionKey const min_final_hbonds; }
namespace cyclic_peptide { extern RealOptionKey const total_energy_cutoff; }
namespace cyclic_peptide { extern RealOptionKey const hbond_energy_cutoff; }
namespace cyclic_peptide { extern BooleanOptionKey const do_not_count_adjacent_res_hbonds; }
namespace cyclic_peptide { extern IntegerOptionKey const fast_relax_rounds; }
namespace cyclic_peptide { extern BooleanOptionKey const count_sc_hbonds; }
namespace cyclic_peptide { extern BooleanOptionKey const require_disulfides; }
namespace cyclic_peptide { extern RealOptionKey const disulf_cutoff_prerelax; }
namespace cyclic_peptide { extern RealOptionKey const disulf_cutoff_postrelax; }
namespace cyclic_peptide { extern RealVectorOptionKey const user_set_alpha_dihedrals; }
namespace cyclic_peptide { extern RealOptionKey const user_set_alpha_dihedral_perturbation; }
namespace cyclic_peptide { extern BooleanOptionKey const filter_oversaturated_hbond_acceptors; }
namespace cyclic_peptide { extern RealOptionKey const hbond_acceptor_energy_cutoff; }
namespace cyclic_peptide { extern RealOptionKey const sample_cis_pro_frequency; }
namespace cyclic_peptide { extern BooleanOptionKey const design_peptide; }
namespace cyclic_peptide { extern FileOptionKey const allowed_residues_by_position; }
namespace cyclic_peptide { extern BooleanOptionKey const prohibit_D_at_negative_phi; }
namespace cyclic_peptide { extern BooleanOptionKey const prohibit_L_at_positive_phi; }
namespace cyclic_peptide { extern FileOptionKey const L_alpha_comp_file; }
namespace cyclic_peptide { extern FileOptionKey const D_alpha_comp_file; }
namespace cyclic_peptide { extern FileOptionKey const L_beta_comp_file; }
namespace cyclic_peptide { extern FileOptionKey const D_beta_comp_file; }
namespace cyclic_peptide { extern IntegerOptionKey const angle_relax_rounds; }
namespace cyclic_peptide { extern IntegerOptionKey const angle_length_relax_rounds; }
namespace cyclic_peptide { extern IntegerOptionKey const cartesian_relax_rounds; }
namespace cyclic_peptide { extern BooleanOptionKey const use_classic_rama_for_sampling; }
namespace cyclic_peptide { extern IntegerVectorOptionKey const n_methyl_positions; }
namespace cyclic_peptide { extern IntegerVectorOptionKey const TBMB_positions; }
namespace cyclic_peptide { extern BooleanOptionKey const use_TBMB_filters; }
namespace cyclic_peptide { extern RealOptionKey const TBMB_sidechain_distance_filter_multiplier; }
namespace cyclic_peptide { extern RealOptionKey const TBMB_constraints_energy_filter_multiplier; }
namespace cyclic_peptide { extern BooleanOptionKey const link_all_cys_with_TBMB; }
namespace cyclic_peptide { extern IntegerOptionKey const require_symmetry_repeats; }
namespace cyclic_peptide { extern BooleanOptionKey const require_symmetry_mirroring; }
namespace cyclic_peptide { extern RealOptionKey const require_symmetry_angle_threshold; }
namespace cyclic_peptide { extern RealOptionKey const require_symmetry_perturbation; }
namespace cyclic_peptide { extern IntegerVectorOptionKey const MPI_processes_by_level; }
namespace cyclic_peptide { extern IntegerVectorOptionKey const MPI_batchsize_by_level; }
namespace cyclic_peptide { extern StringOptionKey const MPI_sort_by; }
namespace cyclic_peptide { extern BooleanOptionKey const MPI_choose_highest; }
namespace cyclic_peptide { extern RealOptionKey const MPI_output_fraction; }
namespace cyclic_peptide { extern IntegerOptionKey const MPI_stop_after_time; }
namespace cyclic_peptide { extern RealOptionKey const MPI_pnear_lambda; }
namespace cyclic_peptide { extern RealOptionKey const MPI_pnear_kbt; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
