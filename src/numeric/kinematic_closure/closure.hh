// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// This file is contains all the same functions as `bridgeObjects.cc', but
// radians are assumed for all angles.  Although this is easier to use in most
// cases, it is not backwards compatible.  That's why a new file was created.


/// @file   bridgeObjects.hh
/// @brief  Header file for bridgeObjects code for ALC
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

#ifndef INCLUDED_numeric_kinematic_closure_hh
#define INCLUDED_numeric_kinematic_closure_hh

// Rosetta Headers
#include <numeric/types.hh>

// Utility headers

#include <utility/vector1.fwd.hh>


namespace numeric {
namespace kinematic_closure {
namespace radians {

void bridge_objects (
	const utility::vector1<utility::vector1<numeric::Real> >& atoms,
	const utility::vector1<numeric::Real> & dt,
	const utility::vector1<numeric::Real> & da,
	const utility::vector1<numeric::Real> & db,
	const utility::vector1<int>& pivots,
	const utility::vector1<int>& order,
	utility::vector1<utility::vector1<numeric::Real> >& t_ang,
	utility::vector1<utility::vector1<numeric::Real> >& b_ang,
	utility::vector1<utility::vector1<numeric::Real> >& b_len,
	int& nsol);

void chainTORS (
	const int& n,
	const utility::vector1<utility::vector1<numeric::Real> >& atoms,
	utility::vector1<numeric::Real>& t_ang,
	utility::vector1<numeric::Real>& b_ang,
	utility::vector1<numeric::Real>& b_len,
	utility::vector1<numeric::Real>& R0,
	utility::vector1<utility::vector1<numeric::Real> >& Q);

numeric::Real torsion(
	const utility::vector1<numeric::Real>& a,
	const utility::vector1<numeric::Real>& b,
	const utility::vector1<numeric::Real>& c,
	const utility::vector1<numeric::Real>& d);

void chainXYZ  (
	const int& n,
	const utility::vector1<numeric::Real>& b_len,
	const utility::vector1<numeric::Real>& b_ang,
	const utility::vector1<numeric::Real>& t_ang,
	const bool space,
	const utility::vector1<numeric::Real>& R0,
	const utility::vector1<utility::vector1<numeric::Real> >& Q,
	utility::vector1<utility::vector1<numeric::Real> >& atoms);

void chainXYZ  (
	const int& n,
	const utility::vector1<numeric::Real>& b_len,
	const utility::vector1<numeric::Real>& b_ang,
	const utility::vector1<numeric::Real>& t_ang,
	utility::vector1<utility::vector1<numeric::Real> >& atoms);

numeric::Real bondangle(
	const utility::vector1<numeric::Real>& a,
	const utility::vector1<numeric::Real>& b,
	const utility::vector1<numeric::Real>& c);

void to_radians(
	utility::vector1<Real> & degrees);

void to_degrees(
	utility::vector1<Real> & radians);

} // end namespace radians
} // end namespace kinematic_closure
} // end namespace numeric

#endif
