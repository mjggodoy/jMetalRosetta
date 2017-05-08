// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainDoubleFunc.hh
/// @brief Definition for functions used in definition of loop closure energies.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_func_GaussianChainDoubleFunc_HH
#define INCLUDED_core_scoring_func_GaussianChainDoubleFunc_HH

#include <core/scoring/func/GaussianChainDoubleFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// See GaussianChainFunc.cc for more information, including link to mathematical derivation.

// C++ Headers
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class GaussianChainDoubleFunc : public Func {
public:

	GaussianChainDoubleFunc( Real const gaussian_variance,
		Real const loop_fixed_cost,
		Real const D2 );

	FuncOP
	clone() const;

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

private:

	void initialize_parameters();
	void recompute_parameters();

private:

	Real gaussian_variance_;
	Real loop_fixed_cost_;
	Real D2_;

	// following have nice default values
	Real kB_T_;

	// derived from above.
	Real loop_fixed_cost_total_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	GaussianChainDoubleFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_GaussianChainDoubleFunc )
#endif // SERIALIZATION


#endif
