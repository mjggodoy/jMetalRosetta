// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_protocols_constraint_generator_ConstraintGenerator_fwd_hh
#define INCLUDED_protocols_constraint_generator_ConstraintGenerator_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace constraint_generator {

class ConstraintGenerator;
typedef utility::pointer::shared_ptr< ConstraintGenerator > ConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< ConstraintGenerator const > ConstraintGeneratorCOP;
typedef utility::pointer::weak_ptr< ConstraintGenerator > ConstraintGeneratorAP;
typedef utility::pointer::weak_ptr< ConstraintGenerator const > ConstraintGeneratorCAP;

typedef utility::vector1< ConstraintGeneratorOP > ConstraintGeneratorOPs;
typedef utility::vector1< ConstraintGeneratorCOP > ConstraintGeneratorCOPs;

}
}

#endif
