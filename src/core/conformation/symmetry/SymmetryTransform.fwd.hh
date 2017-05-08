// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmetryTransform.fwd.hh
/// @brief  Forward declarations and owning pointer declarations for the SymmetryTransform
/// class, which stores information about a transform between two symmetry subunits.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_conformation_symmetry_SymmetryTransform_fwd_hh
#define INCLUDED_core_conformation_symmetry_SymmetryTransform_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace symmetry {

class SymmetryTransform;
typedef utility::pointer::shared_ptr< SymmetryTransform > SymmetryTransformOP;
typedef utility::pointer::shared_ptr< SymmetryTransform const > SymmetryTransformCOP;

} // symmetry
} // conformation
} // core

#endif

