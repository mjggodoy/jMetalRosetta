// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/select/residue_selector/ConstraintGeneratorFactory.hh
/// @brief  Class for instantiating arbitrary ResidueSelectors from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_protocols_constraint_generator_ConstraintGeneratorFactory_HH
#define INCLUDED_protocols_constraint_generator_ConstraintGeneratorFactory_HH

// Unit headers
#include <utility/SingletonBase.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// Package headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGeneratorCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace protocols {
namespace constraint_generator {

class ConstraintGeneratorFactory : public utility::SingletonBase< ConstraintGeneratorFactory > {
private:
	typedef std::map< std::string, ConstraintGeneratorCreatorOP > CreatorMap;

public:
	ConstraintGeneratorFactory();

	void factory_register( ConstraintGeneratorCreatorOP creator );

	bool has_type( std::string const & constraint_generator_name ) const;

	ConstraintGeneratorOP new_constraint_generator(
		std::string const & constraint_generator_name,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;


	void define_constraint_generator_xml_schema_group( utility::tag::XMLSchemaDefinition & xsd ) const;

	//All of these may move to a new file for schema generation

	static std::string constraint_generator_xml_schema_group_name();
	static std::string
	complex_type_name_for_constraint_generator( std::string const & constraint_name );

	static void
	xsd_constraint_generator_type_definition_w_attributes(
		utility::tag::XMLSchemaDefinition & xsd,
		std::string const & constraint_type,
		std::string const & description,
		utility::tag::AttributeList const & attributes);

private:
	CreatorMap creator_map_;
};


} //namespace constraint_generator
} //namespace protocols


#endif
