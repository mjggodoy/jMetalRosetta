// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/jd3.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_jd3_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_jd3_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace jd3 { extern BooleanOptionKey const jd3; }
namespace jd3 { extern BooleanOptionKey const mpi_work_partition_job_distributor; }
namespace jd3 { extern FileOptionKey const job_definition_schema; }
namespace jd3 { extern BooleanOptionKey const mpi_fast_nonblocking_output; }
namespace jd3 { extern IntegerOptionKey const n_archive_nodes; }
namespace jd3 { extern BooleanOptionKey const do_not_archive_on_node0; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
