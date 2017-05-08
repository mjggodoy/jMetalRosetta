// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
/// @brief
/// A class for holding bond information
///
/// @author
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
///
/////////////////////////////////////////////////////////////////////////

// Rosetta headers
#include <core/chemical/Bond.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.Bond" );

BondName convert_to_BondName( std::string const & id ) {
	int id_int( utility::string2int( id ) );
	// Will return -1 if it's not an integer
	if ( id_int >= 0 && id_int <= 4 ) {
		return BondName( id_int );
	}
	if ( id_int >= 5 && id_int <= 8 ) { // SDF "Any" and degenerate query bonds.
		return UnknownBond;
	}

	// String conversions - allow for upper, lower, mixed and abreviated forms
	std::string const converted( ObjexxFCL::lowercased( id ) );
	if ( utility::startswith(converted, "ar") ) {
		return AromaticBond;
	}
	if ( utility::startswith(converted, "unk") ) {
		return UnknownBond;
	}
	if ( utility::startswith(converted, "sing") ) {
		return SingleBond;
	}
	if ( utility::startswith(converted, "doub") ) {
		return DoubleBond;
	}
	if ( utility::startswith(converted, "tri") ) {
		return TripleBond;
	}

	// Aliases and others
	if ( utility::startswith(converted, "am") ) {
		return AromaticBond; // Aromatic bond doubles for amide bond
	}
	if ( utility::startswith(converted, "carboxy") ) {
		return AromaticBond; // Aromatic bond doubles for delocalized carboxy bond
	}
	if ( utility::startswith(converted, "delocal") ) {
		return AromaticBond; // Aromatic bond doubles for general delocalized bond
	}
	if ( utility::startswith(converted, "pseu") ) {
		return UnknownBond; // PseudoBond
	}
	if ( utility::startswith(converted, "orb") ) {
		return OrbitalBond;
	}

	utility_exit_with_message("Unable to convert '" + id + "' to a bond type description.");
	return UnknownBond;
}

Bond::Bond():
	distance_(0),
	cut_bond_(false),
	bond_name_(UnknownBond),
	order_(UnknownBondOrder),
	conjug_(UnknownConjugability),
	aroma_(UnknownAromaticity),
	ring_(UnknownRingness),
	isometry_(UnknownIsometry)
{}

Bond::Bond(Real d, BondName name, bool cut_bond /*= false*/):
	distance_(d),
	cut_bond_(cut_bond),
	bond_name_(name),
	order_(UnknownBondOrder),
	conjug_(UnknownConjugability),
	aroma_(UnknownAromaticity),
	ring_(UnknownRingness),
	isometry_(UnknownIsometry)
{
	bond_name(name);
}

Bond::Bond(Real d, BondOrder order, BondConjugability conj,  BondAromaticity aroma,
	BondRingness ring /*= UnknownRingness*/, BondIsometry isom /*= UnknownIsometry*/, bool cut_bond /*= false*/):
	distance_(d),
	cut_bond_(cut_bond),
	order_(order),
	conjug_(conj),
	aroma_(aroma),
	ring_(ring),
	isometry_(isom)
{
	if ( order == OrbitalBondOrder ) { bond_name_ = OrbitalBond; }
	else if ( aroma == IsAromaticBond ) { bond_name_ = AromaticBond; }
	else if ( order == SingleBondOrder ) { bond_name_ = SingleBond; }
	else if ( order == DoubleBondOrder ) { bond_name_ = DoubleBond; }
	else if ( order == TripleBondOrder ) { bond_name_ = TripleBond; }
	else { bond_name_ = UnknownBond; }

	if ( aroma == IsAromaticBond && ring == UnknownRingness ) {
		ring_ = BondInRing; // As aromatic implies ring.
	}

	if ( isom == UnknownIsometry && ( order == SingleBondOrder || order == TripleBondOrder ) ) {
		isometry_ = NoBondIsometry;
	}

}

void
Bond::bond_name(BondName name){
	bond_name_ = name;
	// Attempt to keep the other parameters updated with the values, but only if we know it will change
	if ( name == SingleBond ) { order_ = SingleBondOrder; }
	else if ( name == DoubleBond ) { order_ = DoubleBondOrder; }
	else if ( name == TripleBond ) { order_ = TripleBondOrder; }
	else if ( name == OrbitalBond ) { order_ = OrbitalBondOrder; }

	if ( name == DoubleBond || name == TripleBond || name == AromaticBond ) { conjug_ = ConjugableBond; }

	if ( name == AromaticBond ) { aroma_ = IsAromaticBond; } // Probably not strictly speaking true - SDF bondtype 4 is often abused.

	if ( name != DoubleBond && name != AromaticBond ) { isometry_ = NoBondIsometry; }
}

void
Bond::order(BondOrder order){
	order_ = order;
	if ( order == OrbitalBondOrder ) { bond_name_ = OrbitalBond; }
	else if ( order == PseudoBondOrder ) { bond_name_ = UnknownBond; }
	else if ( aroma_ == NonaromaticBond ) {
		if ( order == SingleBondOrder ) { bond_name_ = SingleBond; }
		if ( order == DoubleBondOrder ) { bond_name_ = DoubleBond; }
		if ( order == TripleBondOrder ) { bond_name_ = TripleBond; }
	}
}

void
Bond::aromaticity(BondAromaticity aroma){
	aroma_ = aroma;

	if ( aroma == IsAromaticBond ) {
		bond_name_ = AromaticBond;
	} else if ( aroma == NonaromaticBond && bond_name_ == AromaticBond ) {
		bond_name_ = UnknownBond;
	} else if ( aroma == UnknownAromaticity && bond_name_ == AromaticBond ) {
		bond_name_ = UnknownBond;
	}
}


void
Bond::print( std::ostream & out ) const {
	out << distance_ << std::endl;
}

std::ostream &
operator<< (std::ostream & out, Bond const & bond){
	bond.print( out );
	return out;
}

void Bond::SetSDFType( const core::Size SDF_ID) {
	if ( SDF_ID == 0 || SDF_ID > 5 )  {
		TR.Warning << "Inappropriate SDF bond type: " << SDF_ID << std::endl;
	}
	bond_name_ = UnknownBond;
	order_ = UnknownBondOrder;
	conjug_ = UnknownConjugability;
	aroma_ = UnknownAromaticity;
	ring_ = UnknownRingness;
	isometry_ = UnknownIsometry;
	switch( SDF_ID ) {
	case 0 :
		break;
	case 1 :
		bond_name_ = SingleBond;
		order_ = SingleBondOrder;
		conjug_ = NotConjugableBond;
		isometry_ = NoBondIsometry;
		break;
	case 2 :
		bond_name_ = DoubleBond;
		order_ = DoubleBondOrder;
		conjug_ = ConjugableBond;
		break;
	case 3 :
		bond_name_ = TripleBond;
		order_ = TripleBondOrder;
		conjug_ = ConjugableBond;
		isometry_ = NoBondIsometry;
		break;
	case 4 :
		bond_name_ = AromaticBond;
		conjug_ = ConjugableBond;
		aroma_ = IsAromaticBond;
		break;
	case 5 :
		conjug_ = ConjugableBond;
		break;
	default :
		break;
	}
}

core::Size
Bond::GetNumberOfElectrons() const
{
	switch( order_ ) {
	case UnknownBondOrder :
		// Technically, a triple bond could be in the mix, but in practice
		// when we reach here we're talking about single/double unspecified.
		// (Besides, the BCL doesn't consider triple bond degeneracy.)
		if ( conjug_ == ConjugableBond || aroma_ == IsAromaticBond ) {
			return 3;
		} else {
			TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
			return 0;
		}
	case SingleBondOrder :
		return 2;
	case DoubleBondOrder :
		return 4;
	case TripleBondOrder :
		return 6;
	case OrbitalBondOrder:
	case PseudoBondOrder :
		TR.Warning << "Number of electrons for pseudobond requested." << std::endl;
		return 0;
	default :
		TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
		return 0;
	}
}

core::Size Bond::GetMinimumElectrons() const {
	switch( order_ ) {
	case UnknownBondOrder :
		if ( conjug_ == ConjugableBond || aroma_ == IsAromaticBond ) {
			return 2;
		} else {
			TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
			return 0;
		}
	case SingleBondOrder :
		return 2;
	case DoubleBondOrder :
		return 4;
	case TripleBondOrder :
		return 6;
	case OrbitalBondOrder:
	case PseudoBondOrder :
		TR.Warning << "Number of electrons for pseudobond requested." << std::endl;
		return 0;
	default :
		TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
		return 0;
	}
}

core::Size Bond::GetMaximumElectrons() const {
	switch( order_ ) {
	case UnknownBondOrder :
		// Technically, a triple bond could be in the mix, but in practice
		// when we reach here we're talking about single/double unspecified.
		// (Besides, the BCL doesn't consider triple bond degeneracy.)
		if ( conjug_ == ConjugableBond || aroma_ == IsAromaticBond ) {
			return 4;
		} else {
			TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
			return 0;
		}
	case SingleBondOrder :
		return 2;
	case DoubleBondOrder :
		return 4;
	case TripleBondOrder :
		return 6;
	case OrbitalBondOrder:
	case PseudoBondOrder :
		TR.Warning << "Number of electrons for pseudobond requested." << std::endl;
		return 0;
	default :
		TR.Warning << "Number of electrons for unspecified bond requested." << std::endl;
		return 0;
	}
}

core::Size Bond::GetSDFileID() const {
	switch( order_ ) {
	case UnknownBondOrder :
		if ( aroma_ == IsAromaticBond ) { return 4; }
		else if ( conjug_ == ConjugableBond ) { return 5; }
		else {
			TR.Warning << "SDF designation for unspecified bond requested." << std::endl;
			return 0;
		}
	case SingleBondOrder :
		return 1;
	case DoubleBondOrder :
		return 2;
	case TripleBondOrder :
		return 3;
	case OrbitalBondOrder :
		TR.Warning << "SDF designation for an Orbital bond requested." << std::endl;
		return 5; // TODO: Check if scombs actually used a literal 5 for orbitals.
	case PseudoBondOrder :
		TR.Warning << "SDF designation for a pseudo bond requested." << std::endl;
		return 0;
	default :
		TR.Warning << "SDF designation for unspecified bond requested." << std::endl;
		return 0;
	}
}


core::Size Bond::GetSDAltFileID() const {
	if ( order_ == TripleBondOrder ) { return 3; }
	else if ( aroma_ == IsAromaticBond && order_ == SingleBondOrder ) { return 6; }
	else if ( aroma_ == IsAromaticBond && order_ == DoubleBondOrder ) { return 7; }
	else if ( conjug_ == ConjugableBond ) { return 5; }
	else if ( order_ == SingleBondOrder ) { return 1; }
	else if ( order_ == DoubleBondOrder ) { return 2; }
	else {
		TR.Warning << "SDF Alt designation for unspecified bond requested." << std::endl;
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::Bond::save( Archive & arc ) const {
	arc( CEREAL_NVP( distance_ ) ); // Real
	arc( CEREAL_NVP( cut_bond_ ) ); // _Bool
	arc( CEREAL_NVP( bond_name_ ) ); // enum core::chemical::BondName
	arc( CEREAL_NVP( order_ ) ); // enum core::chemical::BondOrder
	arc( CEREAL_NVP( conjug_ ) ); // enum core::chemical::BondConjugability
	arc( CEREAL_NVP( aroma_ ) ); // enum core::chemical::BondAromaticity
	arc( CEREAL_NVP( ring_ ) ); // enum core::chemical::BondRingness
	arc( CEREAL_NVP( isometry_ ) ); // enum core::chemical::BondIsometry
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::Bond::load( Archive & arc ) {
	arc( distance_ ); // Real
	arc( cut_bond_ ); // _Bool
	arc( bond_name_ ); // enum core::chemical::BondName
	arc( order_ ); // enum core::chemical::BondOrder
	arc( conjug_ ); // enum core::chemical::BondConjugability
	arc( aroma_ ); // enum core::chemical::BondAromaticity
	arc( ring_ ); // enum core::chemical::BondRingness
	arc( isometry_ ); // enum core::chemical::BondIsometry
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::Bond );
#endif // SERIALIZATION
