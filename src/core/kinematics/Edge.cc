// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Edge.cc
/// @brief  Fold tree edge class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/Edge.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <string>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {

static THREAD_LOCAL basic::Tracer tr( "core.kinematics" );

// PHIL -- need to update these routines to include atom info?
// ANDREW -- made a stab at adding atom info to comparison operators and ostream operators; no guarantee

/////////////////////////////////////////////////////////////////////////////
// these two should be inverses:

std::ostream &
operator <<( std::ostream & os, const Edge & e )
{
	std::string tag ( "EDGE" );
	if ( e.is_jump() && e.has_atom_info() ) tag = "JEDGE";
	os << " " << tag << " " << e.start() << ' ' << e.stop() << ' ' << e.label() << ' ';
	if ( e.label() == Edge::CHEMICAL ) os << e.start_atom() << ' ' << e.stop_atom() << ' ';
	if ( e.is_jump() ) {
		if ( e.start_atom().size() ) {
			os << e.start_atom() << ' ' << e.stop_atom() << ' ';
		} else {
			//  os << " X X "; not-necessary with JEDGE tag //otherwise reading becomes difficult
		}
	}
	if ( e.is_jump() && e.has_atom_info() ) {
		if ( e.keep_stub_in_residue() ) {
			os << " INTRA_RES_STUB ";
			debug_assert( e.start_atom().size() );
		} else {
			os << " END ";
		}
	}
	return os;
}

/////////////////////////////////////////////////////////////////////////////
std::istream &
operator >>( std::istream & is, Edge & e )
{
	std::string tag;


	//Removed temp foldtree hack from 08 - JAB
	is >> tag;
	if ( ! (tag == "EDGE" || tag == "JEDGE") ) {
		tr.Trace << "failed reading EDGE tag --- found instead: " << tag << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}
	is >> e.start_ >> e.stop_ >> e.label_;

	e.start_atom_ = ""; e.stop_atom_ = "";
	if ( e.label() == Edge::CHEMICAL ) is >> e.start_atom_ >> e.stop_atom_;
	if ( e.is_jump() && tag == "JEDGE" ) {
		is >> e.start_atom_;
		is >> e.stop_atom_;
	} else return is;
	if ( e.start_atom_ == "X" ) e.start_atom_ = "";
	if ( e.stop_atom_ == "X" ) e.stop_atom_ = "";

	// allow either both atoms set or both unset
	debug_assert( ( (e.start_atom_.size() && e.stop_atom_.size()) )
		|| (( e.start_atom_.size() == 0) && (e.stop_atom_.size() ==0 )));

	is >> tag;
	e.bKeepStubInResidue_ = false;
	if ( tag == "END" ) return is;
	debug_assert( tag == "INTRA_RES_STUB" ); //only allowed if also atoms are specified;
	e.bKeepStubInResidue_ = true;
	return is;
}


/////////////////////////////////////////////////////////////////////////////
/// @details compare start residue number first, then stop residue number, then label
/// index number, then start_atom, then stop_atom
bool
operator <( Edge const & a, Edge const & b )
{
	//return ( a.start() <  b.start() ||
	// a.start() == b.start() && a.stop() <  b.stop() ||
	// a.start() == b.start() && a.stop() == b.stop() && a.label() < b.label() ||
	// a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label() && a.start_atom() < b.start_atom() ||
	// a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label() && a.start_atom() == b.start_atom() && a.stop_atom() < b.stop_atom() );
	//);
	return ( a.start() == b.start() ? ( a.stop() == b.stop() ? ( a.label() == b.label() ?
		( a.start_atom() == b.start_atom() ? a.stop_atom() < b.stop_atom() : a.start_atom() < b.start_atom() ) :
		a.label() < b.label() ) : a.stop() < b.stop() ) : a.start() < b.start() );
}


/////////////////////////////////////////////////////////////////////////////
/// @details when start residue number, stop residue number and label index number are all equal
bool
operator ==( Edge const & a, Edge const & b )
{
	return ( a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label()
		&& a.start_atom() == b.start_atom() && a.stop_atom() == b.stop_atom() );
}


/////////////////////////////////////////////////////////////////////////////
/// @details when any of start residue number, stop residue number and label index number is not equal
bool
operator !=( Edge const & a, Edge const & b )
{
	return ( a.start() != b.start() || a.stop() != b.stop() || a.label() != b.label() );
}

const int Edge::PEPTIDE;
const int Edge::CHEMICAL;

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::Edge::save( Archive & arc ) const {
	arc( CEREAL_NVP( start_ ) ); // int
	arc( CEREAL_NVP( stop_ ) ); // int
	arc( CEREAL_NVP( label_ ) ); // int
	arc( CEREAL_NVP( start_atom_ ) ); // std::string
	arc( CEREAL_NVP( stop_atom_ ) ); // std::string
	arc( CEREAL_NVP( bKeepStubInResidue_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::Edge::load( Archive & arc ) {
	arc( start_ ); // int
	arc( stop_ ); // int
	arc( label_ ); // int
	arc( start_atom_ ); // std::string
	arc( stop_atom_ ); // std::string
	arc( bKeepStubInResidue_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::Edge );
#endif // SERIALIZATION
