// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/hbtrie/HBAtom.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBAtom.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {


HBAtom::HBAtom() :
	xyz_(0.0,0.0,0.0),
	base_xyz_(0.0,0.0,0.0),
	base2_xyz_(0.0,0.0,0.0),
	is_hydrogen_( false ),
	is_backbone_( false ),
	is_protein_( false ),
	is_dna_( false ),
	hb_chem_type_( 0 )
	//,seqpos_( 0 )
{}

HBAtom::~HBAtom() {}

/// @brief send a description of the atom to standard out
void
HBAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
HBAtom::print( std::ostream & os ) const
{
	os << "HBAtom" <<  " ";
	os << "(" << xyz().x();
	os << ", " << xyz().y();
	os << ", " << xyz().z() << "), abase";
	os << "(" << base_xyz().x();
	os << ", " << base_xyz().y();
	os << ", " << base_xyz().z() << "), abase2";
	os << "(" << base2_xyz().x();
	os << ", " << base2_xyz().y();
	os << ", " << base2_xyz().z() << "), ";
	os << hb_chem_type_ << ";";
}

std::ostream & operator << ( std::ostream & os, HBAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::hbtrie::HBAtom::save( Archive & arc ) const {
	arc( CEREAL_NVP( xyz_ ) ); // Vector
	arc( CEREAL_NVP( base_xyz_ ) ); // Vector
	arc( CEREAL_NVP( base2_xyz_ ) ); // Vector
	arc( CEREAL_NVP( is_hydrogen_ ) ); // _Bool
	arc( CEREAL_NVP( is_backbone_ ) ); // _Bool
	arc( CEREAL_NVP( is_protein_ ) ); // _Bool
	arc( CEREAL_NVP( is_dna_ ) ); // _Bool
	arc( CEREAL_NVP( hb_chem_type_ ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::hbtrie::HBAtom::load( Archive & arc ) {
	arc( xyz_ ); // Vector
	arc( base_xyz_ ); // Vector
	arc( base2_xyz_ ); // Vector
	arc( is_hydrogen_ ); // _Bool
	arc( is_backbone_ ); // _Bool
	arc( is_protein_ ); // _Bool
	arc( is_dna_ ); // _Bool
	arc( hb_chem_type_ ); // int
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::hbtrie::HBAtom );
#endif // SERIALIZATION
