// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfo.cc
/// @brief   Method definitions for CarbohydrateInfo.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

// Utility headers
#include <utility/py/PyAssert.hh>
#include <utility/exit.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <algorithm>
#include <iostream>

#include <boost/algorithm/string/case_conv.hpp> // we don't need all string

// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.chemical.carbohydrates.CarbohydrateInfo" );


namespace core {
namespace chemical {
namespace carbohydrates {

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Standard constructor
/// @param    <residue_type>: the ResidueType object containing this CarbohydrateInfo
CarbohydrateInfo::CarbohydrateInfo( core::chemical::ResidueTypeCAP residue_type ) : utility::pointer::ReferenceCount()
{
	init( residue_type );
}

// "Copy constructor"
CarbohydrateInfo::CarbohydrateInfo( CarbohydrateInfo const & object_to_copy,
	core::chemical::ResidueTypeCAP new_owner ) :
	utility::pointer::ReferenceCount( object_to_copy )
{
	residue_type_ = new_owner;
	copy_data( *this, object_to_copy );
}

// Destructor
CarbohydrateInfo::~CarbohydrateInfo() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods
void
CarbohydrateInfo::show( std::ostream & output ) const
{
	using namespace std;

	string prefix, suffix, ring_form, modifications;

	// Parse properties.
	if ( is_aldose() ) {
		prefix = "aldo";
	} else /*is ketose*/ {
		char num = '0' + anomeric_carbon_;
		prefix = string( 1, num ) + string( "-keto" );
	}
	switch ( n_carbons_ ) {
	case 3 :
		suffix = "triose";
		break;
	case 4 :
		suffix = "tetrose";
		break;
	case 5 :
		suffix = "pentose";
		break;
	case 6 :
		suffix = "hexose";
		break;
	case 7 :
		suffix = "heptose";
		break;
	case 8 :
		suffix = "octose";
		break;
	case 9 :
		suffix = "nonose";
		break;
	}

	ring_form = CarbohydrateInfoManager::morpheme_from_ring_size( ring_size_ ) + "ose";

	for ( uint position = 1; position <= n_carbons_; ++position ) {
		if ( modifications_[ position ] != "" ) {
			modifications += string( "  " );
			modifications += modifications_[ position ];
			modifications += string( "\n" );
		}
	}
	if ( modifications == "" ) {
		modifications = "  none\n";
	}

	// Produce output.
	output << "Carbohydrate Properties for this Residue:" << endl;
	output << " Basic Name: " << base_name() << endl;
	output << " IUPAC Name: " << full_name_ << endl;
	output << " Abbreviation: " << short_name_ << endl;
	output << " Classification: " << prefix << suffix << endl;
	output << " Stereochemistry: " << stereochem_ << endl;
	if ( ring_size_ != 0 ) {
		output << " Ring Form: " << ring_form << endl;
		output << " Anomeric Form: " << anomer_ << endl;
	}
	output << " Modifications: " << endl << modifications;
	output << " Polymeric Information:" << endl;
	if ( mainchain_glycosidic_bond_acceptor_ ) {
		output << "  Main chain connection: (_->" << mainchain_glycosidic_bond_acceptor_ << ')' << endl;
	} else {
		output << "  Main chain connection: N/A" << endl;
	}
	output << "  Branch connections: ";
	if ( n_branches() == 0 ) {
		output << "none" << endl;
	} else {
		for ( uint i = 1; i <= n_branches(); ++i ) {
			output << "(_->" << branch_points_[i] << ')';
			if ( i != n_branches() ) { output << "; "; }
		}
		output << endl;
	}
}


// Accessors/Mutators
// Return the standard/common, non-residue, short name of the monosaccharide.
std::string
CarbohydrateInfo::base_name() const
{
	using namespace std;

	core::chemical::ResidueTypeCOP residue_type( residue_type_ );
	string const & code = residue_type->name3();
	string const & root = CarbohydrateInfoManager::root_from_code( code );

	// First cover accepted trivial names.
	// TODO: Remove hardcoding!
	if ( code == "Gly" ) {
		return root + "aldehyde";
	} else if ( code == "DHA" ) {
		return root + "one";
	} else if ( code == "Neu" || code == "Kdn" ) {
		return root + "ic acid";
	}

	// Then cover other names systematically.
	// The order here matters. Follow IUPAC priority rules.
	if ( is_uronic_acid() ) {
		return root + "uronic acid";
	} else if ( is_amino_sugar() ) {
		return root + "osamine";
	} else {
		return root + "ose";
	}
}


// Return the attachment point of the downstream saccharide residue attached to ith branch off of this residue.
/// @param    <i>: the branch point index
/// @return   an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point of the
/// downstream monosaccharide residue; e.g., 4 specifies O4
/// @details  A monosaccharide with a group linked to it at one position is a distinct residue type from the same
/// monosaccharide with the same group linked to it at another position.  For example, Rosetta treats (1->4)-beta-
/// D-glucopyranose as an entirely distinct residue type from (1->3)-beta-D-glucopyranose, with separate .params
/// files for each.\n
/// \n
/// See also:\n
///  CarbohydrateInfo.mainchain_glycosidic_bond_acceptor()\n
///  CarbohydrateInfo.n_branches()
core::uint
CarbohydrateInfo::branch_point( core::uint const i ) const
{
	debug_assert( ( i > 0 ) && ( i <= n_branches() ) );
	PyAssert( (i > 0) && ( i <= n_branches() ),
		"CarbohydrateInfo::branch_point( core::uint i ): "
		"There is no ith branch point on this carbohydrate residue.");

	return branch_points_[ i ];
}


// Side-chain modifications
// Return true if any hydroxyl group has been modified to an acetylated amino group.
bool
CarbohydrateInfo::is_N_acetylated() const {
	return residue_type_.lock()->properties().has_property( ACETYLAMINO_SUGAR );
}

// Return true if any hydroxyl group has been modified by acetylation.
bool
CarbohydrateInfo::is_O_acetylated() const {
	return residue_type_.lock()->properties().has_property( ACETYL_SUGAR );
}

// Return true if the sugar has been acetylated at any position.
bool
CarbohydrateInfo::is_acetylated() const {
	return is_N_acetylated() || is_O_acetylated();
}

// Return true if any hydroxyl group has been modified to an amino group or an acetylated amino group.
bool
CarbohydrateInfo::is_amino_sugar() const {
	return residue_type_.lock()->properties().has_property( AMINO_SUGAR ) || is_N_acetylated();
}

// Return true if the primary hydroxyl group is oxidized to the acid.
bool
CarbohydrateInfo::is_uronic_acid() const {
	return residue_type_.lock()->properties().has_property( URONIC_ACID );
}


// Private methods /////////////////////////////////////////////////////////////
// Empty constructor
CarbohydrateInfo::CarbohydrateInfo() : utility::pointer::ReferenceCount()
{
	init( core::chemical::ResidueTypeCAP() );
}

// Initialize data members from properties.
void
CarbohydrateInfo::init( core::chemical::ResidueTypeCAP residue_type_in )
{
	// Set default values.
	residue_type_ = residue_type_in;
	cyclic_oxygen_ = 0;  // assumes linear
	cyclic_oxygen_name_ = "";
	cyclic_oxygen_index_ = 0;
	virtual_cyclic_oxygen_index_ = 0;

	ResidueTypeCOP residue_type( residue_type_ );
	if ( residue_type->is_lower_terminus() ) {
		is_glycoside_ = false;  // can be overridden later
	} else {
		is_glycoside_ = true;
	}

	read_and_set_properties();

	determine_polymer_connections();

	determine_IUPAC_names();
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
CarbohydrateInfo::copy_data(
	CarbohydrateInfo & object_to_copy_to,
	CarbohydrateInfo const & object_to_copy_from)
{
	object_to_copy_to.full_name_ = object_to_copy_from.full_name_;
	object_to_copy_to.short_name_ = object_to_copy_from.short_name_;
	object_to_copy_to.anomeric_carbon_ = object_to_copy_from.anomeric_carbon_;
	object_to_copy_to.anomeric_carbon_name_ = object_to_copy_from.anomeric_carbon_;
	object_to_copy_to.anomeric_carbon_index_ = object_to_copy_from.anomeric_carbon_index_;
	object_to_copy_to.cyclic_oxygen_ = object_to_copy_from.cyclic_oxygen_;
	object_to_copy_to.cyclic_oxygen_name_ = object_to_copy_from.cyclic_oxygen_name_;
	object_to_copy_to.cyclic_oxygen_index_ = object_to_copy_from.cyclic_oxygen_index_;
	object_to_copy_to.virtual_cyclic_oxygen_index_ = object_to_copy_from.virtual_cyclic_oxygen_index_;
	object_to_copy_to.n_carbons_ = object_to_copy_from.n_carbons_;
	object_to_copy_to.stereochem_ = object_to_copy_from.stereochem_;
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.anomer_ = object_to_copy_from.anomer_;
	object_to_copy_to.is_glycoside_ = object_to_copy_from.is_glycoside_;
	object_to_copy_to.modifications_ = object_to_copy_from.modifications_;
	object_to_copy_to.mainchain_glycosidic_bond_acceptor_ = object_to_copy_from.mainchain_glycosidic_bond_acceptor_;
	object_to_copy_to.branch_points_ = object_to_copy_from.branch_points_;
	object_to_copy_to.has_exocyclic_linkage_to_child_mainchain_ = object_to_copy_from.has_exocyclic_linkage_to_child_mainchain_;
}


// Read through all the properties.  Check for impossible cases.  If any property type is not set, the default
// value will be maintained.
void
CarbohydrateInfo::read_and_set_properties()
{
	using namespace std;
	using namespace utility;

	ResidueTypeCOP residue_type( residue_type_ );
	ResidueProperties const & properties( residue_type->properties() );


	// Start with general sugar properties:
	// Oxidation type
	if ( properties.has_property( ALDOSE ) && ! properties.has_property( KETOSE ) ) {
		anomeric_carbon_ = 1;
		anomeric_carbon_name_ = "C1";
	} else if ( properties.has_property( KETOSE ) && ! properties.has_property( ALDOSE ) ) {
		anomeric_carbon_ = 2;  // TODO: Provide method for dealing with non-ulose ketoses.
		anomeric_carbon_name_ = "C2";
	} else {
		utility_exit_with_message( "A sugar must be EITHER an aldose OR a ketose; check the .params file." );
	}

	// Number of carbons
	Size const n_chain_length_properties( properties.has_property( TRIOSE ) + properties.has_property( TETROSE ) +
		properties.has_property( PENTOSE ) + properties.has_property( HEXOSE ) +
		properties.has_property( HEPTOSE ) + properties.has_property( OCTOSE ) +
		properties.has_property( NONOSE ) );
	if ( n_chain_length_properties == 1 ) {
		if ( properties.has_property( TRIOSE ) ) {
			n_carbons_ = 3;
		} else if ( properties.has_property( TETROSE ) ) {
			n_carbons_ = 4;
		} else if ( properties.has_property( PENTOSE ) ) {
			n_carbons_ = 5;
		} else if ( properties.has_property( HEXOSE ) ) {
			n_carbons_ = 6;
		} else if ( properties.has_property( HEPTOSE ) ) {
			n_carbons_ = 7;
		} else if ( properties.has_property( OCTOSE ) ) {
			n_carbons_ = 8;
		} else /* NONOSE */ {
			n_carbons_ = 9;
		}
	} else {
		utility_exit_with_message( "A sugar must have a single property designating length; check the .params file." );
	}

	// Stereochemistry
	if ( properties.has_property( L_SUGAR ) && ! properties.has_property( D_SUGAR ) ) {
		stereochem_ = 'L';
	} else if ( properties.has_property( D_SUGAR ) && ! properties.has_property( L_SUGAR ) ) {
		stereochem_ = 'D';
	} else {
		utility_exit_with_message( "A sugar must have EITHER L OR D stereochemistry; check the .params file." );
	}

	// Ring Size
	Size const n_ring_size_properties( properties.has_property( OXIROSE ) + properties.has_property( OXETOSE ) +
		properties.has_property( FURANOSE ) + properties.has_property( PYRANOSE ) +
		properties.has_property( SEPTANOSE ) );
	if ( n_ring_size_properties == 0 ) {
		ring_size_ = 0;  // assumes linear
	} else if ( n_ring_size_properties == 1 ) {
		if ( properties.has_property( OXIROSE ) ) {
			ring_size_ = 3;
		} else if ( properties.has_property( OXETOSE ) ) {
			ring_size_ = 4;
		} else if ( properties.has_property( FURANOSE ) ) {
			ring_size_ = 5;
		} else if ( properties.has_property( PYRANOSE ) ) {
			ring_size_ = 6;
		} else /* SEPTANOSE */ {
			ring_size_ = 7;
		}
	} else {
		utility_exit_with_message( "A sugar cannot have multiple ring sizes; check the .params file." );
	}

	// Anomer
	if ( properties.has_property( ALPHA_SUGAR ) && ! properties.has_property( BETA_SUGAR ) ) {
		anomer_ = "alpha";
	} else if ( properties.has_property( BETA_SUGAR ) && ! properties.has_property( ALPHA_SUGAR ) ) {
		anomer_ = "beta";
	} else if ( ! properties.has_property( ALPHA_SUGAR ) && ! properties.has_property( BETA_SUGAR ) ) {
		anomer_ = "";  // assumes linear
	} else {
		utility_exit_with_message( "A sugar cannot be both alpha and beta; check the .params file." );
	}

	// Next, look for modifications:
	modifications_.resize( n_carbons_ );

	// Modifications for Which the Position Is Inherent
	if ( properties.has_property( GLYCOSIDE ) ) {
		is_glycoside_ = true;
	}
	if ( properties.has_property( URONIC_ACID ) ) {
		modifications_[ n_carbons_ ] = "uronic acid";
	}

	// Modifications with Positions
	// TODO: This is crap; I am going to refactor this with a VariantType to position, string map. ~Labonte
	vector1< string > const & variants( residue_type->properties().get_list_of_variants() );
	Size const n_variants( variants.size() );
	string variant;
	uint position;  // location of modification; 0 for a property that does not have an associated position

	for ( uint i( 1 ); i <= n_variants; ++i ) {
		// If the 2nd character (index 1) of ith property is a number, it is a modification.
		// Otherwise, it is a regular property or a modification for which the position is inherent, such as uronic
		// acid.
		variant = variants[ i ];
		position = atoi( &variant[ 1 ] );
		if ( position ) {
			if ( modifications_[ position ] != "" ) {
				utility_exit_with_message(
					"A sugar cannot have multiple modifications at the same position; check the .params file.");
			} else {
				variant = variant.substr( 3 );  // assumes 1st character is a "C" and the 3rd an underscore
				if ( variant != "BRANCH_POINT" ) {
					boost::algorithm::to_lower( variant );
					replace( variant.begin(), variant.end(), '_', ' ' );
					modifications_[ position ] = variant;
				}
			}
		}
	}

	// Double-check for inconsistencies.
	if ( ( ring_size_ != 0 ) && ( anomer_ == "" ) ) {
		utility_exit_with_message(
			"A cyclic sugar must have its anomeric property declared; check the .params file." );
	}
	if ( ( ring_size_ == 0 ) && ( anomer_ != "" ) ) {
		utility_exit_with_message( "An acyclic sugar cannot be alpha or beta; check the .params file." );
	}

	anomeric_carbon_index_ = residue_type->atom_index( anomeric_carbon_name_ );

	// Determine cyclic oxygen from "cut bond" neighbor to the anomeric carbon, if applicable.
	if ( ring_size_ != 0 ) {
		cyclic_oxygen_index_ = residue_type->cut_bond_neighbor( anomeric_carbon_index_ )[ 1 ];
		cyclic_oxygen_name_ = residue_type->atom_name( cyclic_oxygen_index_ );
		string const cyclic_oxygen_position = cyclic_oxygen_name_.substr( 2, 1 );  // 3rd col. (index 2) is the atom #
		cyclic_oxygen_ = atoi( &cyclic_oxygen_position[ 0 ] );
		// TODO: There will likely be a better way to do this once Andrew finishes the virtual shadowing code.
		virtual_cyclic_oxygen_index_ = residue_type->atom_index( "VO" + cyclic_oxygen_position );
	}
}

// Get connection data from the residue type.
void
CarbohydrateInfo::determine_polymer_connections()
{
	using namespace std;
	using namespace id;

	ResidueTypeCOP residue_type( residue_type_ );

	// Main chain connections
	if ( ! residue_type->is_upper_terminus() ) {
		uint upper_atom_index = residue_type->upper_connect_atom();
		string atom_name = residue_type->atom_name( upper_atom_index );
		mainchain_glycosidic_bond_acceptor_ = atoi( &atom_name[ 2 ] );  // 3rd column (index 2) is the atom number
	} else {
		mainchain_glycosidic_bond_acceptor_ = 0;
	}

	// Branch points
	Size const n_connections( residue_type->n_possible_residue_connections() );
	for ( uint i = 1; i <= n_connections; ++i ) {
		if ( i == residue_type->lower_connect_id() || i == residue_type->upper_connect_id() ) { continue; }
		uint const branch_atom_index( residue_type->residue_connect_atom_index( i ) );
		string const & branch_atom_name( residue_type->atom_name( branch_atom_index ) );
		if ( branch_atom_name[ 1 ] == 'O' ) {  // 2nd column (index 1) is the element; must be oxygen
			branch_points_.push_back( atoi( &branch_atom_name[ 2 ] ) );  // 3rd column (index 2) is the atom number
		}
	}

	// Exocyclic linkage?
	if ( mainchain_glycosidic_bond_acceptor_ > last_carbon_in_ring() ) {
		has_exocyclic_linkage_to_child_mainchain_ = true;
	} else {
		has_exocyclic_linkage_to_child_mainchain_ = false;
	}
}

// Determine and set the full and abbreviated IUPAC names.
// The NAME property in the .params file is actually the standard IUPAC abbreviation (of an internal/unpatched
// residue), not the full name.  It, combined with any patches, is the Rosetta name for the ResidueType.  The IUPAC
// names will change depending on the residue's place in the sequence and/or any patches.
void
CarbohydrateInfo::determine_IUPAC_names()
{
	using namespace std;

	core::chemical::ResidueTypeCOP residue_type( residue_type_ );

	// Determine root.
	string const & code = residue_type->name3();
	string const & root = CarbohydrateInfoManager::root_from_code( code );

	// Flag for special cases.
	// TODO: Find a way to remove this hardcoding.
	bool const is_Neu( code == "Neu" );

	// Determine prefixes.
	stringstream long_prefixes( stringstream::out );
	stringstream short_prefixes( stringstream::out );

	// Connectivity
	if ( ! residue_type->is_upper_terminus() ) {
		long_prefixes << "->" << mainchain_glycosidic_bond_acceptor_ << ")-";
	}
	if ( ! residue_type->is_lower_terminus() ) {
		long_prefixes << anomer_ << '-';
	}
	short_prefixes << long_prefixes.str();

	// Substitutions
	// TODO: How do I alphabetize and group substitutions?  I'll need a vector to sort.  For now, order by position.
	for ( uint position = 1; position <= n_carbons_; ++position ) {
		if ( modifications_[ position ] == "amino sugar" ) {
			if ( ! is_Neu ) {  // (Neu is by definition an amino sugar.)
				long_prefixes << position << "-amino-" << position << "-deoxy-";
			}
		}
		if ( modifications_[ position ] == "acetylamino sugar" ) {
			if ( is_Neu ) {  // (Neu is by definition an amino sugar.)
				long_prefixes << position << "-acetyl-";
			} else {
				long_prefixes << position << "-(N-acetylamino)-" << position << "-deoxy-";
			}
		}
		if ( modifications_[ position ] == "acetyl sugar" ) {
			long_prefixes << position << "-acetyl-";
		}
	}

	// Stereochemistry
	if ( CarbohydrateInfoManager::default_stereochem_from_code( code ) != '*' ) {  // (* implies inherent stereochem.)
		long_prefixes << stereochem_ << '-';
		short_prefixes << stereochem_ << '-';
	}

	// Determine suffix.
	stringstream long_suffix( stringstream::out );
	stringstream short_suffix( stringstream::out );

	// Ring form
	if ( ! is_Neu ) {  // (For some odd reason, "opyran" is not used with Neu, though "p" is....)
		long_suffix << 'o' << CarbohydrateInfoManager::morpheme_from_ring_size( ring_size_ );
	}
	short_suffix << CarbohydrateInfoManager::ring_affix_from_ring_size( ring_size_ );

	if ( residue_type->is_lower_terminus() ) {
		if ( is_glycoside_ ) {  // TODO: Extract name of R-group.
			if ( is_uronic_acid() ) {
				long_suffix << "uronoside";
				short_suffix << "A";
			} else {
				long_suffix << "oside";
				if ( is_amino_sugar() ) {
					short_suffix << "N";
				}
			}
		} else /* is not a glycoside */ {
			if ( is_uronic_acid() ) {
				long_suffix << "uronate";
				short_suffix << "A";
			} else {
				if ( is_Neu ) {
					long_suffix << "ate";  // TODO: Rework such that all acids get "ate" or "ic acid".
				} else {
					long_suffix << "ose";
				}
				if ( is_amino_sugar() ) {
					if ( ! is_Neu ) {  // (Neu is by definition an amino sugar.)
						short_suffix << "N";
					}
				}
			}
		}
	} else {
		if ( is_uronic_acid() ) {
			long_suffix << "uronoyl";
			short_suffix << "A";
		} else {
			if ( is_Neu ) {
				long_suffix << "yl";
			} else {
				long_suffix << "osyl";
			}
			if ( is_amino_sugar() ) {
				if ( ! is_Neu ) {  // (Neu is by definition an amino sugar.)
					short_suffix << "N";
				}
			}
			if ( is_acetylated() ) {
				short_suffix << "Ac";
			}
		}
		short_suffix << '-';
	}

	full_name_ = long_prefixes.str() + root + long_suffix.str();
	short_name_ = short_prefixes.str() + code + short_suffix.str();
}



// Helper methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that CarbohydrateInfo can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, CarbohydrateInfo const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
