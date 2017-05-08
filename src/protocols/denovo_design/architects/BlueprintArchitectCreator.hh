// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/denovo_design/architects/BlueprintArchitectCreator.hh
/// @brief  Class for instantiating a particular DeNovoArchitect
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_denovo_design_architects_BlueprintArchitectCreator_HH
#define INCLUDED_protocols_denovo_design_architects_BlueprintArchitectCreator_HH

// Package headers
#include <protocols/denovo_design/architects/DeNovoArchitectCreator.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace protocols {
namespace denovo_design {
namespace architects {

class BlueprintArchitectCreator : public DeNovoArchitectCreator {
public:
	/// @brief Instantiate a particular DeNovoArchitect
	virtual DeNovoArchitectOP
	create_architect( std::string const & architect_id ) const;

	/// @brief Return a string that will be used to instantiate the particular DeNovoArchitect
	virtual std::string
	keyname() const;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace architects
} //namespace denovo_design
} //namespace protocols

#endif
