// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Stub.cc
/// @brief  Stub class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <numeric/conversions.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

#endif // SERIALIZATION


namespace core {
namespace kinematics {

Stub::Stub( RT const & rt ):
	M( rt.get_rotation() ),
	v( rt.get_translation() )
{}

/// @brief output operator, 3x3 matrix followed by an xyzVector
std::ostream &
operator<<( std::ostream & os, Stub const & a )
{
	os << "STUB";
	for ( int i=1; i<= 3; ++i ) {
		for ( int j=1; j<= 3; ++j ) {
			os << ' ' << a.M(i,j);
		}
	}
	for ( int i=1; i<= 3; ++i ) {
		os << ' ' << a.v(i);
	}
	return os;
}

/////////////////////////////////////////////////////////////////////////////
/// @brief build a stub from a center points and other three points a, b, c
///
/// @details orthogonal coord frame M contains three unit vectors by column,
///the first one is the unit vector from b pointing to a, the second one is the
///unit vector which is in the plane defined by vector b->a and b->c and
///perpendicular to b->a, the third one is the cross product of the first two.
void
Stub::from_four_points(
	Vector const & center,
	Vector const & a,
	Vector const & b,
	Vector const & c
)
{
	Vector e1( a - b);
	e1.normalize();

	Vector e3( cross( e1, c - b ) );
	e3.normalize();

	Vector e2( cross( e3,e1) );
	M.col_x( e1 ).col_y( e2 ).col_z( e3 );
	v = center;
}


Vector
Stub::build_fake_xyz( Size const index ) const
{
	Real const length( 1.4 );
	Real const angle( numeric::conversions::radians( 120.0 ) );
	Vector const xyz1( v );
	Vector const xyz2( xyz1 + length * M * Vector( -1, 0, 0 ) );
	Vector const xyz3( xyz2 + length * M * Vector( std::cos( angle ), std::sin( angle ), 0 ) );

	switch( index ) {
	case 1 :
		return xyz1;
	case 2 :
		return xyz2;
	case 3 :
		return xyz3;
	default :
		utility_exit_with_message( "Stub::build_fake_xyz must be called with 1<= index <= 3" );
	}

	return Vector( 0.0 ); // wont get here

}

/// @details Floating point comparison is notoriously fickle.  Perform the same
/// series of rotations and translations with ever so slightly different code, and
/// you will produce two different stubs.  The comparison operator is looking for
/// exact matches.  The only reason to expect that stub A and stub B are the same
/// is because A was copied from B; if neither have changed since the copy, then
/// you will have the same stub.
bool
Stub::operator == ( Stub const & rhs ) const {
	return M == rhs.M && v == rhs.v;
}

/// @details a stub center at 0.0 with lab frame(identity matrix)
Stub /* const */ default_stub( Stub::Matrix::identity(), Stub::Vector( 0.0 ) );


} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::Stub::save( Archive & arc ) const {
	arc( CEREAL_NVP( M ) ); // Matrix
	arc( CEREAL_NVP( v ) ); // Vector
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::Stub::load( Archive & arc ) {
	arc( M ); // Matrix
	arc( v ); // Vector
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::Stub );
#endif // SERIALIZATION
