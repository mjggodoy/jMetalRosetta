// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file HelixParams.cc.
/// @brief Functions implementing the Crick equations for a straight helix (NOT a helical bundle).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//C++ headers
#include <cmath>

// Unit headers
#include <numeric/types.hh>

//Project headers
#include <numeric/crick_equations/HelixParams.hh>

namespace numeric {
namespace crick_equations {

/// @brief Returns the x-coordinate of a point on a helix given r1, omega1, and t.
/// Note that this is not translated or rotated.
Real x (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return r1*cos( omega1*t + delta_omega1 );
}

/// @brief Returns the y-coordinate of a point on a helix given r1, omega1, and t.
/// Note that this is not translated or rotated.
Real y (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return r1*sin( omega1*t + delta_omega1);
}

/// @brief Returns the z-coordinate of a point on a helix given r1, omega1, and t.
/// Note that this is not translated or rotated.
Real z (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &t,
	Real const &dz1,
	Real const &/*delta_omega1*/,
	Real const &delta_z1
) {
	return dz1*t + delta_z1;
}

/// @brief Returns the x-, y-, and z-coordinates of a point on a helix given r1, omega1, and t.
/// Note that this is not translated or rotated.
xyzVector <Real> xyz(
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &dz1,
	Real const &delta_omega1,
	Real const &delta_z1
) {
	return xyzVector<Real> ( x(r1,omega1,t,dz1,delta_omega1,delta_z1), y(r1,omega1,t,dz1,delta_omega1,delta_z1), z(r1,omega1,t,dz1,delta_omega1,delta_z1) );
}

/// @brief Returns the derivative of x with respect to r1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dx_dr1 (
	Real const &/*r1*/,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return cos(omega1*t + delta_omega1);
}

/// @brief Returns the derivative of y with respect to r1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dy_dr1 (
	Real const &/*r1*/,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return sin(omega1*t + delta_omega1);
}

/// @brief Returns the derivative of z with respect to r1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dz_dr1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/ ,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) {
	return 0;
}

/// @brief Returns the derivative of x with respect to omega1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dx_domega1 (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return -r1*t*sin(omega1*t + delta_omega1);
}

/// @brief Returns the derivative of y with respect to omega1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dy_domega1 (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/ ,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) {
	return r1*t*cos(omega1*t + delta_omega1);
}

/// @brief Returns the derivative of z with respect to omega1 for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dz_domega1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) {
	return 0;
}

/// @brief Returns the derivative of x with respect to dz1 (rise per residue) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dx_ddz1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) {return 0;}

/// @brief Returns the derivative of y with respect to dz1 (rise per residue) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dy_ddz1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) {return 0;}

/// @brief Returns the derivative of z with respect to dz1 (rise per residue) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dz_ddz1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &t,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) {return t;}

/// @brief Returns the derivative of x with respect to delta_omega1 (omega-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dx_ddelta_omega1 (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) { return -r1*sin(omega1*t+delta_omega1); }

/// @brief Returns the derivative of y with respect to delta_omega1 (omega-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dy_ddelta_omega1 (
	Real const &r1,
	Real const &omega1,
	Real const &t,
	Real const &/*dz1*/,
	Real const &delta_omega1,
	Real const &/*delta_z1*/
) { return r1*cos(omega1*t+delta_omega1); }

/// @brief Returns the derivative of z with respect to delta_omega1 (omega-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dz_ddelta_omega1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) { return 0; }

/// @brief Returns the derivative of x with respect to delta_z1 (z-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dx_ddelta_z1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) { return 0; }

/// @brief Returns the derivative of y with respect to delta_z1 (z-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dy_ddelta_z1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) { return 0; }

/// @brief Returns the derivative of z with respect to delta_z1 (z-offset) for a given value of r1, omega1, and t.
/// Note that this is not translated or rotated.
Real dz_ddelta_z1 (
	Real const &/*r1*/,
	Real const &/*omega1*/,
	Real const &/*t*/,
	Real const &/*dz1*/,
	Real const &/*delta_omega1*/,
	Real const &/*delta_z1*/
) { return 1; }

} //namespace crick_equations
} //namespace numeric
