// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for reading in the orbital type properties
///
/// @details
/// This class reads in the orbital_properties.txt file which contains the "chemical" information for orbitals.
/// This does not contain the actual properties, but sets the properties through the OrbitalType class.
/// This class is called by the ChemicalManager. Modeled off of atomtypeset.
///
///
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeSet_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeSet_hh


// Project headers
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <core/types.hh>
//#include <utility/vector1_bool.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace orbitals {


class OrbitalTypeSet : public utility::pointer::ReferenceCount {


public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~OrbitalTypeSet();
	OrbitalTypeSet(std::string const & directory, std::string const & name="" );

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string const &
	name() const { return name_; }

	void read_file(std::string const & filename);

	/// @brief [ ] operator, simulating vector index behavior
	///
	/// @details look up an OrbitalTypeSet by 1-based indexing
	///
	OrbitalType const &
	operator[](core::Size const index) const
	{
		return *(orbitals_[index]);
	}


	/// @brief lookup the orbital type by the orbital type name string
	int
	orbital_type_index( std::string const & orbital_type_name ) const;

	/// @brief lookup the orbital type by the orbital type name string
	int
	orbital_type_index( std::string & orbital_type_name ) const;

private:

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string name_;

	/// lookup map: get orbital_type_index by orbital_type_name
	std::map< std::string, int > orbital_type_index_;


	/// @brief  Save the directory name for future use
	std::string directory_;

	/// @brief a collection of OrbitalTypes,
	///
	/// @details OrbitalType has data of atom properties, and it can be
	/// looked up by orbital_type_index.
	utility::vector1< OrbitalTypeOP > orbitals_;


};


#ifdef    SERIALIZATION
/// @brief Serialize an OrbitalTypeSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void serialize_orbital_type_set( Archive & arc, OrbitalTypeSetCOP ptr );

/// @brief Deserialize an OrbitalTypeSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void deserialize_orbital_type_set( Archive & arc, OrbitalTypeSetCOP & ptr );
#endif // SERIALIZATION

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_orbitals_OrbitalTypeSet )

SPECIAL_COP_SERIALIZATION_HANDLING( core::chemical::orbitals::OrbitalTypeSet, core::chemical::orbitals::serialize_orbital_type_set, core::chemical::orbitals::deserialize_orbital_type_set )
#endif // SERIALIZATION

#endif /* ORBITALTYPESET_HH_ */
