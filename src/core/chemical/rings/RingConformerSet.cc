// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/RingConformerSet.cc
/// @brief   Method definitions for RingConformerSet.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/rings/RingConformerManager.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>

// C++ headers
#include <cmath>
#include <iostream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.chemical.rings.RingConformerSet" );


namespace core {
namespace chemical {
namespace rings {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////////
// Standard constructor
/// @param    <ring_size>: an unsigned integer expressing the size of the ring
/// @param    <lowest_conformer>: IUPAC name for the lowest-energy ring conformer, if known
/// @param    <low_conformers>: IUPAC name for other low-energy ring conformers
RingConformerSet::RingConformerSet(
	core::Size const ring_size,
	std::string const & lowest_conformer,
	utility::vector1< std::string > const & low_conformers ) : utility::pointer::ReferenceCount()
{
	init( ring_size, lowest_conformer, low_conformers );
}

// Copy constructor
RingConformerSet::RingConformerSet( RingConformerSet const & object_to_copy ) :
	utility::pointer::ReferenceCount( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Assignment operator
RingConformerSet &
RingConformerSet::operator=( RingConformerSet const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}

// Destructor
RingConformerSet::~RingConformerSet() {}


// Standard Rosetta methods ///////////////////////////////////////////////////////////////////////////////////////////
// General methods
void
RingConformerSet::show( std::ostream & output ) const
{
	using namespace std;

	output << "Possible ring conformers:" << endl;

	Size const n_conformers( nondegenerate_conformers_.size() );
	for ( uint i( 1 ); i <= n_conformers; ++i ) {
		output << "   " << nondegenerate_conformers_[ i ] << endl;
	}

	output << endl;
}


// Accessors/Mutators
// Return the size of the conformer set.
core::Size
RingConformerSet::size() const
{
	return nondegenerate_conformers_.size();
}

// Are the low-energy conformers known for this set?
bool
RingConformerSet::low_energy_conformers_are_known() const
{
	// energy_minima_conformers_ should always include at least energy_minimum_conformer_, but then there is only one
	// option to select from, so return false.
	return ( energy_minima_conformers_.size() > 1 );
}


// Return a list of all nondegenerate conformers in the set.
utility::vector1< RingConformer > const &
RingConformerSet::get_all_nondegenerate_conformers() const
{
	return nondegenerate_conformers_;
}


// Return the conformer corresponding to the requested name.
/// @param    <name>: the IUPAC name for a specific ring conformation, e.g., "1C4"
/// @details  For a saccharide residue, the provided name assumes a ring with the anomeric carbon labeled 1.  That is,
/// for a 2-ketopyranose in the 2C5 chair form, provide 1C4.
/// @return   matching RingConformer or exits, if no match is found
/// @note     This is slow, but it should not be called by most protocols, which will pull randomly from various
/// subsets.
// AMW: cppcheck wants you to change to pass by reference; DO NOT
RingConformer const &
RingConformerSet::get_ideal_conformer_by_name( std::string const & name ) const
{
	using namespace std;

	Size const n_conformers( size() );
	for ( uint i( 1 ); i <= n_conformers; ++i ) {
		RingConformer const & conformer( nondegenerate_conformers_[ i ] );
		if ( conformer.specific_name == name ) {
			return conformer;
		}
	}

	utility_exit_with_message( "No conformer with given name found in this set; exiting." );
	{
		static RingConformer DUMMY_CONFORMER;
		return DUMMY_CONFORMER;  // will never be reached
	}
}

// Return the conformer that is the best fit for the provided Cremer-Pople parameters.
/// @param    <parameter>: an appropriate, ordered list of C-P parameters, with angles in degrees:\n
/// @details  For 4-membered rings, provide q.\n
/// For 5-membered rings, provide q, phi.\n
/// For 6-membered rings, provide q, phi, theta.
/// For ideal conformers, q is ignored, if non-zero, except for 4-membered rings, where only the sign matters.
/// @return   matching RingConformer or exits, if no match is found
/// @note     This is slow, but it should not be called by most protocols, which will pull randomly from various
/// subsets.
// AMW: cppcheck wants you to change to pass by reference; DO NOT
RingConformer const &
RingConformerSet::get_ideal_conformer_by_CP_parameters( utility::vector1< core::Real > const & parameters ) const
{
	using namespace std;
	using namespace numeric;

	// Drop out if this is not a 4- to 6-membered ring set.
	if ( ring_size_ == 3 ) {
		utility_exit_with_message( "A 3-membered ring can only be planar; exiting." );
	}
	if ( ring_size_ > 6 ) {
		utility_exit_with_message(
			"Rosetta does not currently handle C-P parameters for rings larger than size 6; exiting." );
	}

	Size const n_parameters( parameters.size() );

	// Check for reasonable values.
	if ( n_parameters != ring_size_ - 3 ) {
		utility_exit_with_message( "An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
			"yet a different number was provided; exiting." );
	}
	if ( parameters[ q ] == 0.0 ) {
		utility_exit_with_message( "Planar ring conformations are not handled by Rosetta; "
			"please specify a non-zero q value; exiting." );
	}

	Size const n_conformers( size() );
	uint adjusted_phi, adjusted_theta;

	// Interpret parameters as appropriate for this ring size.
	switch ( ring_size_ ) {
	case 4 :
		utility_exit_with_message( "4-membered rings not yet handled in Rosetta; exiting." );
		//if ( parameters[ q ] > 0 ) {
		// TODO
		//} else /* q < 0 */ {
		// TODO
		//}
		break;
	case 5 :
		// Adjust input to the nearest 18th degree.
		adjusted_phi = uint( ( nonnegative_principal_angle_degrees( parameters[ PHI ] ) + 18/2 ) / 18 ) * 18;

		if ( TR.Debug.visible() ) {
			TR.Debug << "Searching for phi = " << adjusted_phi << "..." << endl;
		}

		for ( uint i( 1 ); i <= n_conformers; ++i ) {
			RingConformer const & conformer( nondegenerate_conformers_[ i ] );
			if ( uint( conformer.CP_parameters[ PHI ] ) == adjusted_phi ) {
				return conformer;
			}
		}
		break;
	case 6 :
		// Adjust input to the nearest 30th and 45th degrees, respectively.
		// TODO: Use numeric::nearest_angle_degrees.
		adjusted_phi = uint( ( nonnegative_principal_angle_degrees( parameters[ PHI ] ) + 30/2 ) / 30 ) * 30;
		adjusted_theta = uint( ( nonnegative_principal_angle_degrees( parameters[ THETA ] ) + 45/2 ) / 45 ) * 45;

		// Theta must be between 0 and 180.
		if ( adjusted_theta > 180 ) {
			adjusted_theta = 360 - adjusted_theta;
		}

		// If at the poles, phi is meaningless, so set to arbitrary value 180, as listed in database.
		if ( adjusted_theta == 180 || adjusted_theta == 0 ) {
			adjusted_phi = 180;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Searching for phi = " << adjusted_phi;
			TR.Debug << " and theta = " << adjusted_theta << "..." << endl;
		}

		for ( uint i( 1 ); i <= n_conformers; ++i ) {
			RingConformer const & conformer( nondegenerate_conformers_[ i ] );
			if ( uint( conformer.CP_parameters[ PHI ] ) == adjusted_phi ) {
				if ( uint( conformer.CP_parameters[ THETA ] ) == adjusted_theta ) {
					return conformer;
				}
			}
		}
		break;
	}

	utility_exit_with_message( "No conformer with given parameters found in this set; exiting." );
	{
		static RingConformer DUMMY_CONFORMER;
		return DUMMY_CONFORMER;  // will never be reached
	}
}

// Return the conformer that is the best fit for the provided list of nu angles.
/// @param    <angles>: an appropriate, ordered list of nu angles in degrees, one less than the ring size.
/// @return   matching RingConformer or exits, if no match is found
/// @note     This is slow, but it should not be called by most protocols, which will pull randomly from various
/// subsets.
RingConformer const &
RingConformerSet::get_ideal_conformer_from_nus( utility::vector1< core::Angle > const & angles, core::Real limit /*=90.0*/ ) const
{
	using namespace std;
	using namespace numeric;

	// Drop out if this is not a 4- to 6-membered ring set.
	if ( ring_size_ == 3 ) {
		utility_exit_with_message( "A 3-membered ring can only be planar; exiting." );
	}
	if ( ring_size_ > 6 ) {
		utility_exit_with_message(
			"Rosetta does not currently handle rings larger than size 6; exiting." );
	}

	Size const n_angles( angles.size() );

	if ( n_angles != ring_size_ - 1 ) {
		TR.Error << "A " << ring_size_ << "-membered ring requires ";
		TR.Error << ring_size_ - 1 << " nu angles; exiting." << endl;
		utility_exit();
	}

	// Now search for best match.
	if ( TR.Debug.visible() ) {
		TR.Debug << "Searching for";
		for ( uint i( 1 ); i <= n_angles; ++i ) {
			TR.Debug << " nu" << i << " = " << angles[ i ];
			if ( i != n_angles ) {
				TR.Debug << ',';
			}
		}
		TR.Debug << "..." << endl;
	}
	Size const n_conformers( size() );
	Angle best_match_score( limit * n_angles );  // the highest possible sum of angle deviations halved
	uint best_conformer( 0 );
	for ( uint i( 1 ); i <= n_conformers; ++i ) {
		RingConformer const & conformer( nondegenerate_conformers_[ i ] );
		Angle match_score( 0.0 );
		for ( uint j( 1 ); j <= n_angles; ++j ) {
			// The match score is simply the sum of the angle deviations. 0 is a perfect match.
			match_score += abs( conformer.nu_angles[ j ] - principal_angle_degrees( angles[ j ] ) );
		}
		if ( match_score < best_match_score ) {
			best_match_score = match_score;
			best_conformer = i;
		}
	}

	if ( ! best_conformer ) {
		utility_exit_with_message( "No conformer with given nu angles found in this set; exiting." );
	}

	RingConformer const & conformer( nondegenerate_conformers_[ best_conformer ] );
	return conformer;
}


// Return the conformer that is known from studies (if available) to be the lowest energy ring conformer.
RingConformer const &
RingConformerSet::get_lowest_energy_conformer() const
{
	return energy_minimum_conformer_;
}

// Return a random conformer from the set.
RingConformer const &
RingConformerSet::get_random_conformer() const
{
	uint const i( uint( numeric::random::rg().uniform() * degenerate_conformers_.size() + 1 ) );

	return degenerate_conformers_[ i ];
}

// Return a random conformer from the subset of conformers that are local minima.
RingConformer const &
RingConformerSet::get_random_local_min_conformer() const
{
	uint const i( uint( numeric::random::rg().uniform() * energy_minima_conformers_.size() + 1 ) );

	return energy_minima_conformers_[ i ];
}

// Return the subset of conformers that are local minima.
utility::vector1< RingConformer >
RingConformerSet::get_local_min_conformers() const
{
	return energy_minima_conformers_;
}

// Private methods ////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize data members for the given ring size.
void
RingConformerSet::init(
	core::uint const ring_size,
	std::string const & lowest_conformer_in,
	utility::vector1< std::string > const & low_conformers )
{
	using namespace utility;

	ring_size_ = ring_size;
	std::string lowest_conformer( lowest_conformer_in );

	// For 6-membered rings, assume 4C1, if lowest conformer is not provided or known.
	if ( lowest_conformer == "" && ring_size_ == 6 ) {
		lowest_conformer = "4C1";
	}

	Size const n_conformers( RingConformerManager::conformers_for_ring_size( ring_size_ ).size() );
	for ( uint i( 1 ); i <= n_conformers; ++i ) {
		RingConformer const & conformer( RingConformerManager::conformers_for_ring_size( ring_size_ ).at( i ) );
		nondegenerate_conformers_.push_back( conformer );
		for ( uint copy_num( 1 ); copy_num <= conformer.degeneracy; ++copy_num ) {
			degenerate_conformers_.push_back( conformer );
		}

		// Other subsets
		if ( conformer.specific_name == lowest_conformer ) {
			energy_minimum_conformer_ = conformer;
		}
		if ( conformer.specific_name == lowest_conformer || low_conformers.contains( conformer.specific_name ) ) {
			energy_minima_conformers_.push_back( conformer );
		}
	}
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformerSet::copy_data(
	RingConformerSet & object_to_copy_to,
	RingConformerSet const & object_to_copy_from )
{
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.nondegenerate_conformers_ = object_to_copy_from.nondegenerate_conformers_;
	object_to_copy_to.degenerate_conformers_ = object_to_copy_from.degenerate_conformers_;
	object_to_copy_to.energy_minimum_conformer_ = object_to_copy_from.energy_minimum_conformer_;
	object_to_copy_to.energy_minima_conformers_ = object_to_copy_from.energy_minima_conformers_;
}


// Helper methods /////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Two conformers are equivalent if their CP parameters are the same.
bool
operator==( RingConformer const & ring1, RingConformer const & ring2 )
{
	return ( ring1.CP_parameters == ring2.CP_parameters );
}

/// @details  Two conformers are not equivalent if their CP parameters are not the same.
bool
operator!=( RingConformer const & ring1, RingConformer const & ring2 )
{
	return ( ring1.CP_parameters != ring2.CP_parameters );
}

// Insertion operators (overloaded so that RingConformer and RingConformerSet can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, RingConformer const & object_to_output )
{
	output << object_to_output.specific_name << " (" << object_to_output.general_name << "): ";
	Size const n_CP_parameters( object_to_output.CP_parameters.size() );
	if ( n_CP_parameters == 2 /* (for 5-membered rings) */ ||
			n_CP_parameters == 3 /* (for 6-membered rings) */ ) {
		output << "C-P parameters (q, phi, theta): ";
		for ( uint parameter( 1 ); parameter <= n_CP_parameters; ++parameter ) {
			if ( parameter > 1 ) {
				output << ", ";
			}
			output << object_to_output.CP_parameters[ parameter ];
		}
		output << "; ";
	}
	output << "nu angles (degrees): ";
	Size const n_angles( object_to_output.nu_angles.size() );
	for ( uint i( 1 ); i <= n_angles; ++i ) {
		if ( i > 1 ) {
			output << ", ";
		}
		output << object_to_output.nu_angles[ i ];
	}
	return output;
}

std::ostream &
operator<<( std::ostream & output, RingConformerSet const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace rings
}  // namespace chemical
}  // namespace core
