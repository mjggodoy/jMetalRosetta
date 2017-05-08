// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Interface.fwd.hh
/// @author Steven Lewis


#ifndef INCLUDED_protocols_scoring_Interface_fwd_hh
#define INCLUDED_protocols_scoring_Interface_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
namespace scoring {
class Interface;
typedef  utility::pointer::shared_ptr< Interface >  InterfaceOP;
typedef  utility::pointer::shared_ptr< Interface const >  InterfaceCOP;

} // namespace scoring
} // namespace protocols

#endif // INCLUDED_core_conformation_Interface_FWD_HH
