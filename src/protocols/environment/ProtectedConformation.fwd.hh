// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ProtectedConformation.fwd.hh
/// @brief definition of the ProtectedConformation class
/// @author

#ifndef INCLUDED_protocols_environment_ProtectedConformation_fwd_hh
#define INCLUDED_protocols_environment_ProtectedConformation_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace environment {

class ProtectedConformation;
typedef utility::pointer::shared_ptr< ProtectedConformation > ProtectedConformationOP;
typedef utility::pointer::shared_ptr< ProtectedConformation const > ProtectedConformationCOP;

typedef utility::pointer::weak_ptr< ProtectedConformation > ProtectedConformationAP;
typedef utility::pointer::weak_ptr< ProtectedConformation const > ProtectedConformationCAP;

} // environment
} // protocols

#endif //INCLUDED_protocols_moves_ProtectedConformation_fwd_HH
