// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/corrections.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_corrections_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_corrections_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace corrections { extern BooleanOptionKey const corrections; }
namespace corrections { extern BooleanOptionKey const beta; }
namespace corrections { extern BooleanOptionKey const beta_cart; }
namespace corrections { extern BooleanOptionKey const beta_nov16; }
namespace corrections { extern BooleanOptionKey const beta_nov16_cart; }
namespace corrections { extern BooleanOptionKey const beta_nov15; }
namespace corrections { extern BooleanOptionKey const beta_nov15_cart; }
namespace corrections { extern BooleanOptionKey const beta_july15; }
namespace corrections { extern BooleanOptionKey const beta_july15_cart; }
namespace corrections { extern BooleanOptionKey const beta_patch; }
namespace corrections { extern BooleanOptionKey const beta_nov15_patch; }
namespace corrections { extern BooleanOptionKey const newdna; }
namespace corrections { extern BooleanOptionKey const correct; }
namespace corrections { extern BooleanOptionKey const hbond_sp2_correction; }
namespace corrections { extern BooleanOptionKey const facts_default; }
namespace corrections { namespace score { extern BooleanOptionKey const score; } }
namespace corrections { namespace score { extern BooleanOptionKey const bbdep_omega; } }
namespace corrections { namespace score { extern BooleanOptionKey const bbdep_bond_params; } }
namespace corrections { namespace score { extern BooleanOptionKey const bbdep_bond_devs; } }
namespace corrections { namespace score { extern BooleanOptionKey const no_his_his_pairE; } }
namespace corrections { namespace score { extern BooleanOptionKey const no_his_DE_pairE; } }
namespace corrections { namespace score { extern StringOptionKey const p_aa_pp; } }
namespace corrections { namespace score { extern BooleanOptionKey const p_aa_pp_nogridshift; } }
namespace corrections { namespace score { extern BooleanOptionKey const rama_not_squared; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map; } }
namespace corrections { namespace score { extern FileOptionKey const rama_pp_map; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_average_L_flat; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_average_L_flat; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_gly_flat; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_pro_flat; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_average_L_flat_stringent; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_average_L_flat_stringent; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_gly_flat_stringent; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map_sym_pro_flat_stringent; } }
namespace corrections { namespace score { extern BooleanOptionKey const rama_prepro_steep; } }
namespace corrections { namespace score { extern BooleanOptionKey const rama_prepro_nobidentate; } }
namespace corrections { namespace score { extern BooleanOptionKey const cenrot; } }
namespace corrections { namespace score { extern BooleanOptionKey const dun10; } }
namespace corrections { namespace score { extern StringOptionKey const dun10_dir; } }
namespace corrections { namespace score { extern StringOptionKey const dun02_file; } }
namespace corrections { namespace score { extern StringOptionKey const ch_o_bond_potential; } }
namespace corrections { namespace score { extern RealOptionKey const lj_hbond_hdis; } }
namespace corrections { namespace score { extern RealOptionKey const lj_hbond_OH_donor_dis; } }
namespace corrections { namespace score { extern BooleanOptionKey const score12prime; } }
namespace corrections { namespace score { extern RealOptionKey const hbond_energy_shift; } }
namespace corrections { namespace score { extern RealOptionKey const hb_sp2_BAH180_rise; } }
namespace corrections { namespace score { extern RealOptionKey const hb_sp2_outer_width; } }
namespace corrections { namespace score { extern BooleanOptionKey const hb_sp2_chipen; } }
namespace corrections { namespace score { extern BooleanOptionKey const hbond_measure_sp3acc_BAH_from_hvy; } }
namespace corrections { namespace score { extern BooleanOptionKey const hb_fade_energy; } }
namespace corrections { namespace score { extern BooleanOptionKey const hb_cen_soft; } }
namespace corrections { namespace score { extern BooleanOptionKey const use_bicubic_interpolation; } }
namespace corrections { namespace score { extern BooleanOptionKey const dun_normsd; } }
namespace corrections { namespace score { extern BooleanOptionKey const dun_entropy_correction; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const chemical; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const icoor_05_2009; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const parse_charge; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const expand_st_chi2sampling; } }
namespace corrections { extern BooleanOptionKey const shapovalov_lib_fixes_enable; }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shapovalov_lib; } }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shap_dun10_enable; } }
namespace corrections { namespace shapovalov_lib { extern StringOptionKey const shap_dun10_smooth_level; } }
namespace corrections { namespace shapovalov_lib { extern StringOptionKey const shap_dun10_dir; } }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shap_dun10_use_minus_log_P_ignore_P; } }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shap_rama_enable; } }
namespace corrections { namespace shapovalov_lib { extern StringOptionKey const shap_rama_smooth_level; } }
namespace corrections { namespace shapovalov_lib { extern FileOptionKey const shap_rama_map; } }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shap_rama_nogridshift; } }
namespace corrections { namespace shapovalov_lib { extern BooleanOptionKey const shap_p_aa_pp_enable; } }
namespace corrections { namespace shapovalov_lib { extern StringOptionKey const shap_p_aa_pp_smooth_level; } }
namespace corrections { namespace shapovalov_lib { extern StringOptionKey const shap_p_aa_pp; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
