// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseFromPDBLoader.fwd.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_import_pose_PoseFromPDBLoader_FWD_HH
#define INCLUDED_core_import_pose_PoseFromPDBLoader_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace import_pose {

class PoseFromPDBLoader;
typedef utility::pointer::shared_ptr< PoseFromPDBLoader > PoseFromPDBLoaderOP;
typedef utility::pointer::shared_ptr< PoseFromPDBLoader const > PoseFromPDBLoaderCOP;

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_PoseFromPDBLoader_FWD_HH

