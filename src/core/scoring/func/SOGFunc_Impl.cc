// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SOGFunc_Impl.cc
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <numeric/util.hh>

#include <core/scoring/func/SOGFunc_Impl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// C++ Headers

#include <iostream>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

using namespace core::scoring::constraints;

SOGFunc_Impl::~SOGFunc_Impl() {}

void
SOGFunc_Impl::read_data( std::istream & in ) {
	Size n_funcs;

	clear_(); // don't forget to remove old data!

	Real lowest_sdev(100);
	in >> n_funcs;
	for ( Size i = 1; i <= n_funcs; ++i ) {
		Real mean(0.0), sdev(0.0), weight(0.0);
		in >> mean >> sdev >> weight;

		if ( sdev == 0 ) {
			sdev = 4.0;
		}

		means_.  push_back( mean );
		sdevs_.  push_back( sdev );
		weights_.push_back( weight );

		if ( sdev < lowest_sdev ) lowest_sdev = sdev;
	}

	bool renormalize(
		weights_.size() != sdevs_.size()
	);

	if ( !renormalize ) {
		for ( Size ii = 1; ii <= weights_.size(); ++ii ) {
			//bool renormalize(false);
			if ( numeric::isinf( weights_[ii] ) || numeric::isnan( weights_[ii] ) ) {
				renormalize = true;
			}
		}
	}

	if ( renormalize ) renormalize_weights();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ OptionKeys::constraints::sog_upper_bound ].user() ) {
		upper_bound(option[ OptionKeys::constraints::sog_upper_bound ]());
	} else {
		upper_bound(10);
	}
}

void SOGFunc_Impl::renormalize_weights() {
	using core::Size;
	using core::Real;

	Real total(0.0);
	for ( Size ii = 1; ii <= sdevs_.size(); ++ii ) {
		total += std::pow( sdevs_[ii], -10 );
	}

	weights_.resize( sdevs_.size(), 0.0 );
	for ( Size ii = 1; ii <= sdevs_.size(); ++ii ) {
		weights_[ii] = std::pow( sdevs_[ii], -10 ) / total;
	}

	runtime_assert( sdevs_.size() == means_.size() );
	runtime_assert( sdevs_.size() == weights_.size() );
}

bool SOGFunc_Impl::smooth_to_zero() const {
	return smooth_;
}

void SOGFunc_Impl::smooth_to_zero( bool const setting ) {
	smooth_ = setting;
}

void SOGFunc_Impl::upper_bound( Real const r ) {
	debug_assert( r > 0 );
	upper_bound_ = r;
	score_upper_ = -1 * std::log( prob_sum_of_gaussians(r) );
}

Real SOGFunc_Impl::upper_bound() const {
	return upper_bound_;
}

Real SOGFunc_Impl::upper_bound_score() const {
	return score_upper_;
}

void
SOGFunc_Impl::clear_() {
	means_.  clear();
	sdevs_.  clear();
	weights_.clear();
}

core::Real
SOGFunc_Impl::get_alt_score_( Real const x ) const {
	Real const alt_mean( 18.991 );
	Real const alt_sdev(  7.353 );
	return -1 * std::log( dgaussian( x, alt_mean, alt_sdev, 1.0 ) );
}

Real SOGFunc_Impl::prob_sum_of_gaussians( Real const x ) const {
	Real score( 0.0 );
	for ( utility::vector1< Real >::const_iterator
			w = weights_.begin(), w_end = weights_.end(),
			m = means_.begin(), m_end = means_.end(),
			s = sdevs_.begin(), s_end = sdevs_.end(); // iterators
			w != w_end && m != m_end && s != s_end; // condition
			++w, ++m, ++s  // iteration
			) {
		Real temp_sc = dgaussian( x, *m, *s, *w );
		score += temp_sc;
	}

	return score;
}

Real
SOGFunc_Impl::func( Real const x ) const {
	if ( x > upper_bound() ) {
		return 0.0;
	}
	Real score = prob_sum_of_gaussians(x);
	check_bounds( x, score );
	if ( score <= 1e-50 ) {  // avoid floating point comparison to 0
		return 0.0;
	} else {
		Real sc = -std::log(score);
		sc = sc - upper_bound_score();
		return sc;
	}
} // func

Real
SOGFunc_Impl::dfunc( Real const x ) const {
	Real const h( 1e-6 );
	Real const df( (func(x+h) - func(x-h)) / (2*h) );
	check_bounds( x, df );
	return df;
} // dfunc

void SOGFunc_Impl::check_bounds( Real const x, Real const val ) const {
	if ( numeric::isinf( val ) || numeric::isnan( val ) ) {
		std::cerr << "bounds error (radius = " << x << ", val = " << val << "), def = ";
		show_definition( std::cerr );
#ifdef BOINC
		utility_exit_with_message( "Fatal SOGFunc_Impl error." );
#endif
	}
}

void SOGFunc_Impl::show_definition( std::ostream & out ) const {
	debug_assert( weights_.size() == means_.size() );
	debug_assert( weights_.size() == sdevs_.size() );

	out << "SOGFUNC " << weights_.size();
	for ( Size i = 1; i <= weights_.size(); ++i ) {
		out << " " << means_[i] << " " << sdevs_[i] << " " << weights_[i];
	}
	out << '\n';
} // show_definition

core::Real SOGFunc_Impl::sog_cst_param() const {
	return sog_cst_param_;
}

void SOGFunc_Impl::sog_cst_param( core::Real const param ) {
	sog_cst_param_ = param;
}

void SOGFunc_Impl::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ basic::options::OptionKeys::constraints::sog_cst_param ].user() ) {
		sog_cst_param( option[ basic::options::OptionKeys::constraints::sog_cst_param ]() );
	}
}

bool SOGFunc_Impl::operator == ( SOGFunc_Impl const & other ) const
{
	if ( means_ != other.means_ ) return false;
	if ( sdevs_ != other.sdevs_ ) return false;
	if ( weights_ != other.weights_ ) return false;
	if ( upper_bound_ != other.upper_bound_ ) return false;
	if ( score_upper_ != other.score_upper_ ) return false;
	if ( sog_cst_param_ != other.sog_cst_param_ ) return false;
	if ( smooth_ != other.smooth_ ) return false;
	return true;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SOGFunc_Impl::save( Archive & arc ) const {
	arc( CEREAL_NVP( means_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( sdevs_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( weights_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( upper_bound_ ) ); // core::Real
	arc( CEREAL_NVP( score_upper_ ) ); // core::Real
	arc( CEREAL_NVP( sog_cst_param_ ) ); // core::Real
	arc( CEREAL_NVP( smooth_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SOGFunc_Impl::load( Archive & arc ) {
	arc( means_ ); // utility::vector1<core::Real>
	arc( sdevs_ ); // utility::vector1<core::Real>
	arc( weights_ ); // utility::vector1<core::Real>
	arc( upper_bound_ ); // core::Real
	arc( score_upper_ ); // core::Real
	arc( sog_cst_param_ ); // core::Real
	arc( smooth_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SOGFunc_Impl );
CEREAL_REGISTER_TYPE( core::scoring::func::SOGFunc_Impl )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SOGFunc_Impl )
#endif // SERIALIZATION
