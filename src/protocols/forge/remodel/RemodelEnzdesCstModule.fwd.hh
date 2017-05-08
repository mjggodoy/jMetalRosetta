// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelEnzdesCstModule.hh
///
/// @brief this file handles merging constraint defined by enzdes type cstfile
/// @brief and blueprint definition of positions and add them to the pose
/// @author Possu Huang, possu@u.washington.edu, Jan 2010

#ifndef INCLUDED_protocols_forge_remodel_RemodelEnzdesCstModule_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelEnzdesCstModule_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace forge {
namespace remodel {

class RemodelEnzdesCstModule;
typedef utility::pointer::shared_ptr< RemodelEnzdesCstModule > RemodelEnzdesCstModuleOP;

} // namespace remodel
} // namespace forge
} // namespace protocols

#endif
