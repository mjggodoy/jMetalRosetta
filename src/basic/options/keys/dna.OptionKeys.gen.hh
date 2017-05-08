// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/dna.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_dna_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_dna_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace dna { extern BooleanOptionKey const dna; }
namespace dna { namespace specificity { extern BooleanOptionKey const specificity; } }
namespace dna { namespace specificity { extern BooleanOptionKey const exclude_dna_dna; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const params; } }
namespace dna { namespace specificity { extern FileVectorOptionKey const frag_files; } }
namespace dna { namespace specificity { extern BooleanOptionKey const only_repack; } }
namespace dna { namespace specificity { extern BooleanOptionKey const design_DNA; } }
namespace dna { namespace specificity { extern BooleanOptionKey const soft_rep; } }
namespace dna { namespace specificity { extern BooleanOptionKey const dump_pdbs; } }
namespace dna { namespace specificity { extern BooleanOptionKey const fast; } }
namespace dna { namespace specificity { extern BooleanOptionKey const randomize_motif; } }
namespace dna { namespace specificity { extern RealOptionKey const Wfa_elec; } }
namespace dna { namespace specificity { extern RealOptionKey const Wdna_bs; } }
namespace dna { namespace specificity { extern RealOptionKey const Wdna_bp; } }
namespace dna { namespace specificity { extern RealOptionKey const minimize_tolerance; } }
namespace dna { namespace specificity { extern StringOptionKey const weights_tag; } }
namespace dna { namespace specificity { extern StringOptionKey const weights_tag_list; } }
namespace dna { namespace specificity { extern StringOptionKey const min_type; } }
namespace dna { namespace specificity { extern StringOptionKey const mode; } }
namespace dna { namespace specificity { extern StringOptionKey const score_function; } }
namespace dna { namespace specificity { extern BooleanOptionKey const pre_minimize; } }
namespace dna { namespace specificity { extern BooleanOptionKey const post_minimize; } }
namespace dna { namespace specificity { extern BooleanOptionKey const pre_pack; } }
namespace dna { namespace specificity { extern IntegerOptionKey const nloop; } }
namespace dna { namespace specificity { extern IntegerOptionKey const n_inner; } }
namespace dna { namespace specificity { extern IntegerOptionKey const n_outer; } }
namespace dna { namespace specificity { extern IntegerOptionKey const nstep_water; } }
namespace dna { namespace specificity { extern IntegerOptionKey const moving_jump; } }
namespace dna { namespace specificity { extern IntegerOptionKey const motif_begin; } }
namespace dna { namespace specificity { extern IntegerOptionKey const motif_size; } }
namespace dna { namespace specificity { extern StringVectorOptionKey const pdb_pos; } }
namespace dna { namespace specificity { extern StringVectorOptionKey const methylate; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const dna_backbone_torsion_sdevs; } }
namespace dna { namespace specificity { extern RealOptionKey const dna_sugar_torsion_sdev; } }
namespace dna { namespace specificity { extern RealOptionKey const dna_chi_torsion_sdev; } }
namespace dna { namespace specificity { extern StringOptionKey const lk_ball_wtd_tag; } }
namespace dna { namespace specificity { extern BooleanOptionKey const lk_ball_for_bb; } }
namespace dna { namespace specificity { extern RealOptionKey const lk_ball_ramp_width_A2; } }
namespace dna { namespace specificity { extern RealOptionKey const lk_ball_overlap_width_A2; } }
namespace dna { namespace specificity { extern RealOptionKey const lk_ball_water_fade; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const lk_ball_wtd_prefactors; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const lk_ball_waters_sp2; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const lk_ball_waters_sp3; } }
namespace dna { namespace specificity { extern RealVectorOptionKey const lk_ball_waters_ring; } }
namespace dna { namespace specificity { extern RealOptionKey const lk_ball_waters_donor; } }
namespace dna { namespace specificity { extern RealOptionKey const lk_ball_bridge_angle_widthscale; } }
namespace dna { namespace design { extern BooleanOptionKey const design; } }
namespace dna { namespace design { extern BooleanOptionKey const output_unbound_pdb; } }
namespace dna { namespace design { extern RealOptionKey const z_cutoff; } }
namespace dna { namespace design { extern StringOptionKey const protein_scan; } }
namespace dna { namespace design { extern StringOptionKey const checkpoint; } }
namespace dna { namespace design { extern BooleanOptionKey const minimize; } }
namespace dna { namespace design { extern StringVectorOptionKey const dna_defs; } }
namespace dna { namespace design { extern StringOptionKey const dna_defs_file; } }
namespace dna { namespace design { extern BooleanOptionKey const nopdb; } }
namespace dna { namespace design { extern BooleanOptionKey const designable_second_shell; } }
namespace dna { namespace design { extern BooleanOptionKey const base_contacts_only; } }
namespace dna { namespace design { extern IntegerOptionKey const probe_specificity; } }
namespace dna { namespace design { extern BooleanOptionKey const reversion_scan; } }
namespace dna { namespace design { extern BooleanOptionKey const binding; } }
namespace dna { namespace design { extern RealOptionKey const Boltz_temp; } }
namespace dna { namespace design { extern BooleanOptionKey const repack_only; } }
namespace dna { namespace design { extern BooleanOptionKey const sparse_pdb_output; } }
namespace dna { namespace design { namespace specificity { extern BooleanOptionKey const specificity; } } }
namespace dna { namespace design { namespace specificity { extern BooleanOptionKey const output_structures; } } }
namespace dna { namespace design { namespace specificity { extern BooleanOptionKey const include_dna_potentials; } } }
namespace dna { namespace design { namespace reversion { extern BooleanOptionKey const reversion; } } }
namespace dna { namespace design { namespace reversion { extern RealOptionKey const dscore_cutoff; } } }
namespace dna { namespace design { namespace reversion { extern RealOptionKey const dspec_cutoff; } } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
