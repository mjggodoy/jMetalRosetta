// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJob.fwd.hh
/// @brief  Forward header file for MakeRotLibJob class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibJob_fwd_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibJob_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace make_rot_lib {

class MakeRotLibJob;
typedef utility::pointer::shared_ptr< MakeRotLibJob > MakeRotLibJobOP;
typedef utility::pointer::shared_ptr< MakeRotLibJob const > MakeRotLibJobCOP;

}//make_rot_lib
}//protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibJob_fwd_hh
