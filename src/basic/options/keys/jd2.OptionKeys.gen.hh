// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/jd2.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_jd2_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_jd2_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace jd2 { extern BooleanOptionKey const jd2; }
namespace jd2 { extern BooleanOptionKey const pose_input_stream; }
namespace jd2 { extern BooleanOptionKey const lazy_silent_file_reader; }
namespace jd2 { extern BooleanOptionKey const mpi_nowait_for_remaining_jobs; }
namespace jd2 { extern RealOptionKey const mpi_timeout_factor; }
namespace jd2 { extern BooleanOptionKey const mpi_work_partition_job_distributor; }
namespace jd2 { extern BooleanOptionKey const mpi_file_buf_job_distributor; }
namespace jd2 { extern BooleanOptionKey const mpi_filebuf_jobdistributor; }
namespace jd2 { extern BooleanOptionKey const mpi_fast_nonblocking_output; }
namespace jd2 { extern BooleanOptionKey const dd_parser; }
namespace jd2 { extern IntegerOptionKey const ntrials; }
namespace jd2 { extern StringOptionKey const generic_job_name; }
namespace jd2 { extern BooleanOptionKey const no_output; }
namespace jd2 { extern BooleanOptionKey const enzdes_out; }
namespace jd2 { extern IntegerOptionKey const buffer_silent_output; }
namespace jd2 { extern RealOptionKey const buffer_flush_frequency; }
namespace jd2 { extern BooleanOptionKey const delete_old_poses; }
namespace jd2 { extern FileVectorOptionKey const resource_definition_files; }
namespace jd2 { extern FileOptionKey const checkpoint_file; }
namespace jd2 { extern BooleanOptionKey const failed_job_exception; }
namespace jd2 { extern IntegerOptionKey const max_nstruct_in_memory; }
namespace jd2 { extern BooleanOptionKey const sequential_mpi_job_distribution; }
namespace jd2 { extern BooleanOptionKey const grid_ensemble; }
namespace jd2 { extern BooleanOptionKey const seed_ensemble; }
namespace jd2 { extern RealVectorOptionKey const seed_ensemble_weights; }
namespace jd2 { extern FileOptionKey const seed_ensemble_weights_file; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
