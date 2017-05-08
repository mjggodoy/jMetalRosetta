// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SquareWellFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_func_SquareWellFunc_hh
#define INCLUDED_core_scoring_func_SquareWellFunc_hh

#include <core/scoring/func/SquareWellFunc.fwd.hh>

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


// C++ Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class SquareWellFunc : public Func {
public:
	SquareWellFunc( Real const x0_in, Real const well_depth_in ): x0_( x0_in ), well_depth_( well_depth_in ){}

	FuncOP
	clone() const { return FuncOP( new SquareWellFunc( *this ) ); }
	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

	Real x0() const { return x0_; }
	Real well_depth() const { return well_depth_; }
	void x0( Real x ) { x0_ = x; }
	void well_depth( Real well_depth ){ well_depth_ = well_depth;}

	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	Real x0_;
	Real well_depth_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SquareWellFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_SquareWellFunc )
#endif // SERIALIZATION


#endif
