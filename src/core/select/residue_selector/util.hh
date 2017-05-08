// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/util.hh
/// @brief  Utility functions for the residue selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_util_HH
#define INCLUDED_core_select_residue_selector_util_HH

// Core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <list>
#include <algorithm>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Used to name the xs:complexType for a residue selector that is
/// created with the "rs_type" tag-name.  Does so by prepending "rs_" and
/// appending "Type" to the "rs_type".  E.g., "rs_AndType" would be the
/// name given to the complexType to describe the format of the
/// AndResidueSelector.
std::string
complex_type_name_for_residue_selector( std::string const & rs_type );

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains no other ResidueSelectors but may contain some number of
/// attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains a single ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	core::Size min_occurrence,
	core::Size max_occurrence,
	utility::tag::AttributeList const & attributes
);

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" (or whatever option_name is set to)
///  option in tag.
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
ResidueSelectorCOP
parse_residue_selector(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const &option_name="residue_selector"
);

/// @brief Companion function for parse_residue_selector
/// @brief This assumes the default residue selector option name ("residue_selector").
void
attributes_for_parse_residue_selector_default_option_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_residue_selector
void
attributes_for_parse_residue_selector(
	utility::tag::AttributeList & attlist,
	std::string const & option_name = "residue_selector",
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_residue_selector to be used when it is unacceptible
/// for the parse_residue_selector function to return a null pointer
void
attributes_for_parse_residue_selector_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & option_name /*= "residue_selector"*/,
	std::string const & documentation_string /*= ""*/
);

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data );

/// @brief returns a residue selector embeded in a tag.
/// for example, the Chain selector in the example below:
///
/// Example:
///<RESIDUE_SELECTORS>
///    <Neighborhood name="chAB_neighbors">
///        <Chain chains="A,B">
///    </Neighborhood>
///</RESIDUE_SELECTORS>
///
/// @thorws utility::excn::EXCN_Msg_Exception if selector is not embedded ( ! tag->size() > 1)
///
ResidueSelectorOP
get_embedded_residue_selector(  utility::tag::TagCOP tag, basic::datacache::DataMap & datamap );

/// @brief returns residue selectors embeded in a tag.
/// for example, the Chain selector in the example below:
///
/// Example:
///<RESIDUE_SELECTORS>
///    <Neighborhood name="chAB_neighbors">
///        <Chain chains="A,B">
///        <Glycan >
///    </Neighborhood>
///</RESIDUE_SELECTORS>
///
/// @thorws utility::excn::EXCN_Msg_Exception if no embedded selectors (tag->size() <= 1)
///
utility::vector1< ResidueSelectorOP >
get_embedded_residue_selectors( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap );


}
}
}

#endif
