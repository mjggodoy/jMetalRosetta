// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_outputters/SecondaryPoseOutputterCreator.fwd.hh
/// @brief  Forward declaration of the %SecondaryPoseOutputterCreator class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_FWD_HH
#define INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

class SecondaryPoseOutputterCreator;

typedef utility::pointer::shared_ptr< SecondaryPoseOutputterCreator > SecondaryPoseOutputterCreatorOP;
typedef utility::pointer::shared_ptr< SecondaryPoseOutputterCreator const > SecondaryPoseOutputterCreatorCOP;

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_HH
