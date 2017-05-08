// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/Func.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Andrew Leaver-Fay
/// @author James Thompson
/// @author Oliver Lange

#ifndef INCLUDED_core_scoring_func_Func_hh
#define INCLUDED_core_scoring_func_Func_hh

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <core/scoring/func/Func.fwd.hh>
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Func is an abstract base class representing a function used to define
/// constraints, in which func(r) gives the constraint score for the given value
/// r.
class Func : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Func();

	/// @brief This method must return a deep copy of this %Func, meaning if this %Func holds pointers
	/// to other %Func objects, that it must clone those %Func objects as well.
	virtual
	FuncOP clone() const = 0;

	/// @brief Equality operator.  Looks for strict equality.  Floating-point comparison is the
	/// rule rather than the exception.
	virtual
	bool
	operator == ( Func const & other ) const = 0;

	/// @brief inequality operator -- simply the negation of the (virtual) equality operator
	bool
	operator != ( Func const & other ) const;

	/// @brief Does the input Func, "other", have the same type as me?  Necessary for the
	/// equality operator to function correctly.  All derived Func classes must implement
	/// this function.
	virtual
	bool
	same_type_as_me( Func const & other ) const = 0;

	/// @brief initialize this Func from the given std::istream.
	virtual
	void read_data( std::istream & );

	/// @brief Returns a value representing this function evaluated at a given
	/// point.
	virtual
	Real
	func( Real const ) const = 0;

	/// @brief Returns a value representing the derivative of this function
	/// evaluated at a given point.
	virtual
	Real
	dfunc( Real const ) const = 0;

	/// @brief Estimates the derivative of this function at a given radius by
	/// calculating the slope of the secant line from func(r) and func(r+1e-05).
	virtual
	Real estimate_dfunc( Real const r ) const;

	/// @brief Estimates the derivative of this function at a given radius by
	/// calculating the slope of the secant line from func(r) and func(r+h).
	virtual
	Real estimate_dfunc( Real const r, Real const h ) const;

	/// @brief Prints out space-delimited columns for r, func, dfunc and
	/// dfunc_est. The values for func, dfunc and dfunc_est are plotted as a
	/// function of r, which is varied from 2-20 in steps of 0.5. The value for
	/// dfunc_est is the estimated by the method estimate_dfunc( r ).
	virtual void show( std::ostream & out ) const;

	/// @brief shows the definition of this function, usually the string type of
	/// function and the parameters passed in to the constructor.
	virtual void show_definition( std::ostream & out ) const;

	/// @brief show some sort of stringified representation of the violations for
	/// this constraint.
	virtual Size show_violations(
		std::ostream& out, Real r, Size verbose_level, Real threshold = 1
	) const;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class Func

std::ostream & operator << (std::ostream & out, Func const & f );

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_Func )
#endif // SERIALIZATION


#endif
