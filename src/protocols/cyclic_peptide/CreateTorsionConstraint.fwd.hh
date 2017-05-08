// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cyclic_peptide/CreateTorsionConstraint.fwd.hh
/// @brief  Defines owning pointers for CreateTorsionConstraint mover class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide {

class CreateTorsionConstraint; // fwd declaration
typedef utility::pointer::shared_ptr< CreateTorsionConstraint > CreateTorsionConstraintOP;
typedef utility::pointer::shared_ptr< CreateTorsionConstraint const > CreateTorsionConstraintCOP;
typedef utility::vector1<CreateTorsionConstraintOP> CreateTorsionConstraintOPs;
typedef utility::vector1<CreateTorsionConstraintCOP> CreateTorsionConstraintCOPs;

} // namespace cyclic_peptide
} // namespace protocols

#endif // INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_fwd_hh
