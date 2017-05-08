// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/SinXOverX.cc
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// Unit headers
#include <core/scoring/saxs/SinXOverX.hh>

// C++
#include <cmath>

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace scoring {
namespace saxs {

SinXOverX::SinXOverX() {
	fill_sin_x_over_x_table();
}

void SinXOverX::fill_sin_x_over_x_table() {
	sin_x_over_x_.push_back(2.0);
	for ( Real x=0.001; x<=1000; x+=0.01 ) {
		sin_x_over_x_.push_back( 2.0*sin(x)/x );
	}
}


} // saxs
} // scoring
} // core
