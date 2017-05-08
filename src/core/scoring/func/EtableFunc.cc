// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/EtableFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <core/scoring/func/EtableFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <core/types.hh>

#include <utility/exit.hh>

// C++ Headers

#include <iostream>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

using namespace core::scoring::constraints;

bool EtableFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	EtableFunc const & other_downcast( static_cast< EtableFunc const & > (other) );
	if ( func_     != other_downcast.func_     ) return false;
	if ( min_      != other_downcast.min_      ) return false;
	if ( max_      != other_downcast.max_      ) return false;
	if ( stepsize_ != other_downcast.stepsize_ ) return false;

	return true;
}

bool EtableFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< EtableFunc const * > ( &other );
}

void
EtableFunc::read_data( std::istream& in ) {
	in  >> min_ >> max_;
	stepsize_ = 0.1;
	for ( Real r = min_; r <= max_; r += stepsize_ ) {
		core::Real func_temp;
		in >> func_temp;
		func_.push_back( func_temp  );
	}
}

Real
EtableFunc::func( Real const x ) const {
	// find the appopriate index into func_
	Real index = ( x - min_ ) / stepsize_;
	Size x_lower_idx = (Size) (index);
	Size x_upper_idx = x_lower_idx + 1;

	Real x_lower = min_ + (x_lower_idx * stepsize_);
	Real x_upper = min_ + (x_upper_idx * stepsize_);
	return linear_interpolate( x, x_lower, x_upper, func_[x_lower_idx], func_[x_upper_idx] );
} // func

Real
EtableFunc::dfunc( Real const ) const {
	utility_exit_with_message( "dfunc not implemented!\n" );
	return -1;
} // dfunc

void EtableFunc::show_definition( std::ostream& out ) const {
	out << "ETABLEFUNC " << min_ << ' ' << max_;
	for ( utility::vector1< core::Real >::const_iterator f_it = func_.begin(), f_end = func_.end();
			f_it != f_end;
			++f_it
			) {
		out << ' ' << *f_it;
	} // for func_ and dfunc_

	out << "\n";
} // show_definition

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::EtableFunc::EtableFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::EtableFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( func_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( min_ ) ); // core::Real
	arc( CEREAL_NVP( max_ ) ); // core::Real
	arc( CEREAL_NVP( stepsize_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::EtableFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( func_ ); // utility::vector1<core::Real>
	arc( min_ ); // core::Real
	arc( max_ ); // core::Real
	arc( stepsize_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::EtableFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::EtableFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_EtableFunc )
#endif // SERIALIZATION
