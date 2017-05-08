// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/peptide_deriver.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_peptide_deriver_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_peptide_deriver_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace peptide_deriver { extern BooleanOptionKey const peptide_deriver; }
namespace peptide_deriver { extern IntegerVectorOptionKey const pep_lengths; }
namespace peptide_deriver { extern BooleanOptionKey const skip_zero_isc; }
namespace peptide_deriver { extern BooleanOptionKey const dump_peptide_pose; }
namespace peptide_deriver { extern BooleanOptionKey const dump_cyclic_poses; }
namespace peptide_deriver { extern BooleanOptionKey const dump_prepared_pose; }
namespace peptide_deriver { extern BooleanOptionKey const dump_report_file; }
namespace peptide_deriver { extern StringVectorOptionKey const restrict_receptors_to_chains; }
namespace peptide_deriver { extern StringVectorOptionKey const restrict_partners_to_chains; }
namespace peptide_deriver { extern BooleanOptionKey const do_minimize; }
namespace peptide_deriver { extern RealOptionKey const optimize_cyclic_threshold; }
namespace peptide_deriver { extern StringOptionKey const report_format; }
namespace peptide_deriver { extern BooleanOptionKey const report_gzip; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
