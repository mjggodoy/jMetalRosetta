// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/bcl/ElementSet.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_ElementSet_hh
#define INCLUDED_core_chemical_ElementSet_hh


// Unit headers
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/Element.fwd.hh>
#include <core/chemical/Elements.hh>

// Project headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <map>

#include <utility/vector1.hh>

#include <core/types.hh>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @brief A set of Bcl Elements
///
/// @details This class contains a vector of pointers each of which points to an
/// Element and the vector index is looked up by an element_name string
/// in a map.
///
class ElementSet : public utility::pointer::ReferenceCount {

public:
	ElementSet( std::string const & name = "");
	~ElementSet() override;

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string const &
	name() const { return name_; }

	/// @brief Number of elements in the set
	Size
	n_elements() const
	{
		return elements_.size();
	}

	/// @brief Check if there is an element_type associated with an element_symbol string
	bool
	contains_element_type( std::string const & element_symbol ) const;

	/// @brief Lookup the element index by the element enum
	Size
	element_index( core::chemical::element::Elements ele ) const;

	/// @brief Lookup the element index by the element_symbol string
	Size
	element_index( std::string const & element_symbol ) const;

	/// @brief Lookup the element index by the element enum;
	ElementCOP
	element( core::chemical::element::Elements ele ) const;

	/// @brief Lookup the element object by the element_symbol string
	ElementCOP
	element( std::string const & element_symbol ) const;

	/// @brief Lookup an Element by 1-based indexing
	ElementCOP
	operator[] ( Size const index ) const;

	/// @brief Load the ElementSet from a file
	void
	read_file( std::string const & filename );

	// data
private:

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string name_;

	/// @brief element_index_ lookup map
	///
	/// @details element_index_ allows lookup of the element by its symbol
	std::map< std::string, core::Size > element_index_;

	/// @brief a collection of Elements,
	///
	/// @details Element has data of atom properties, and it can be
	/// looked up by element_index.
	utility::vector1< ElementOP > elements_;

};

#ifdef    SERIALIZATION
/// @brief Serialize an ElementSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void serialize_element_set( Archive & arc, ElementSetCOP ptr );

/// @brief Deserialize an AtomTypeSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void deserialize_element_set( Archive & arc, ElementSetCOP & ptr );
#endif // SERIALIZATION

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ElementSet )

SPECIAL_COP_SERIALIZATION_HANDLING( core::chemical::ElementSet, core::chemical::serialize_element_set, core::chemical::deserialize_element_set )
#endif // SERIALIZATION

#endif
