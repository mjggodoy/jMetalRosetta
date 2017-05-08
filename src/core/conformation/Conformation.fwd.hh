// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  forward declaration for Conformation.hh
/// @file   core/conformation/Conformation.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Conformation_fwd_hh
#define INCLUDED_core_conformation_Conformation_fwd_hh


// Unit headers

// Project headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// C++ headers

namespace core {
namespace conformation {

class Conformation;
typedef utility::pointer::shared_ptr< Conformation       > ConformationOP;
typedef utility::pointer::shared_ptr< Conformation const > ConformationCOP;
typedef utility::pointer::weak_ptr< Conformation       > ConformationAP;
typedef utility::pointer::weak_ptr< Conformation const > ConformationCAP;

} // conformation
} // core


#endif
