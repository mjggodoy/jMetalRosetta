// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/fingerprint.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_fingerprint_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_fingerprint_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace fingerprint { extern BooleanOptionKey const fingerprint; }
namespace fingerprint { extern BooleanOptionKey const print_eggshell; }
namespace fingerprint { extern RealOptionKey const atom_radius_scale; }
namespace fingerprint { extern RealOptionKey const atom_radius_buffer; }
namespace fingerprint { extern RealOptionKey const packing_weight; }
namespace fingerprint { extern RealOptionKey const dist_cut_off; }
namespace fingerprint { extern BooleanOptionKey const include_extrashell_to_set_origin; }
namespace fingerprint { extern BooleanOptionKey const include_hydrogens; }
namespace fingerprint { extern BooleanOptionKey const multiple_origin; }
namespace fingerprint { extern BooleanOptionKey const use_DARC_gpu; }
namespace fingerprint { extern BooleanOptionKey const square_score; }
namespace fingerprint { extern BooleanOptionKey const darc_components; }
namespace fingerprint { extern IntegerOptionKey const set_origin; }
namespace fingerprint { extern StringOptionKey const origin_res_num; }
namespace fingerprint { extern BooleanOptionKey const add_esp; }
namespace fingerprint { extern BooleanOptionKey const darc_shape_only; }
namespace fingerprint { extern BooleanOptionKey const darc_elsts_only; }
namespace fingerprint { extern RealOptionKey const esp_weight; }
namespace fingerprint { extern RealOptionKey const esp_protein_wt; }
namespace fingerprint { extern RealOptionKey const esp_surface_wt; }
namespace fingerprint { extern BooleanOptionKey const delphi_grid; }
namespace fingerprint { extern RealOptionKey const cap_e_potential; }
namespace fingerprint { extern BooleanOptionKey const return_zero_darc_score; }
namespace fingerprint { extern BooleanOptionKey const set_surface_esp_to_zero; }
namespace fingerprint { extern BooleanOptionKey const set_protein_esp_to_zero; }
namespace fingerprint { extern StringOptionKey const inp_lig; }
namespace fingerprint { extern StringOptionKey const ref_lig; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
