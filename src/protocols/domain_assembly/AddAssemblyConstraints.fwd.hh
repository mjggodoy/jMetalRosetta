// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/domain_assembly/AddAssemblyConstraints.fwd.hh
/// @brief  Adds constraints for the assembling process
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_AddAssemblyConstraints_fwd_hh
#define INCLUDED_protocols_domain_assembly_AddAssemblyConstraints_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace domain_assembly {

class AddAssemblyConstraints;
typedef utility::pointer::shared_ptr< AddAssemblyConstraints > AddAssemblyConstraintsOP;
typedef utility::pointer::shared_ptr< AddAssemblyConstraints const > AddAssemblyConstraintsCOP;


}//namespace domain_assembly
}//namespace protocols

#endif

