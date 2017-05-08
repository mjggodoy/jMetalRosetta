// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/AtomPairConstraintGenerator.fwd.hh
/// @brief Generates atom pair constraints for a set of residues from the current or reference pose
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_fwd_hh
#define INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace constraint_generator {

class AtomPairConstraintGenerator;

typedef utility::pointer::shared_ptr< AtomPairConstraintGenerator > AtomPairConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< AtomPairConstraintGenerator const > AtomPairConstraintGeneratorCOP;

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_fwd_hh
