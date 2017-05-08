// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints during forge loop remodelling
/// @author Florian Richter, floric@u.washington.edu, april 2009


#ifndef INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_fwd_hh


// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace forge {
namespace remodel {

class RemodelConstraintGenerator;
typedef utility::pointer::shared_ptr< RemodelConstraintGenerator > RemodelConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< RemodelConstraintGenerator const > RemodelConstraintGeneratorCOP;
typedef utility::pointer::weak_ptr< RemodelConstraintGenerator > RemodelConstraintGeneratorAP;
typedef utility::pointer::weak_ptr< RemodelConstraintGenerator const > RemodelConstraintGeneratorCAP;

class GenericRemodelConstraintGenerator;
typedef utility::pointer::shared_ptr< GenericRemodelConstraintGenerator > GenericRemodelConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< GenericRemodelConstraintGenerator const > GenericRemodelConstraintGeneratorCOP;

}
}
}


#endif
