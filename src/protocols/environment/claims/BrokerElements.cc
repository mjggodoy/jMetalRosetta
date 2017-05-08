// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BrokerElements.cc
/// @brief some implementation details of the BrokerElements (such as they are).
/// @details responsibilities:
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/BrokerElements.hh>

// Package Headers

// Project Headers

// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace environment {
namespace claims {

std::string const ResidueElement::type = "ResidueElement";
std::string const JumpElement::type    = "JumpElement";
std::string const CutElement::type     = "CutElement";
std::string const CutBiasElement::type = "CutBiasElement";
std::string const DOFElement::type     = "DOFElement";

std::ostream& operator<<( std::ostream& os, ControlStrength const& cstr ) {
	assert( cstr >= DOES_NOT_CONTROL && cstr <= EXCLUSIVE );
	if ( cstr == DOES_NOT_CONTROL ) {
		os << "DOES_NOT_CONTROL";
	} else if ( cstr == CAN_CONTROL ) {
		os << "CAN_CONTROL";
	} else if ( cstr == MUST_CONTROL ) {
		os << "MUST_CONTROL";
	} else {
		os << "EXCLUSIVE";
	}
	return os;
}


} //claims
} //environment
} //protocols
