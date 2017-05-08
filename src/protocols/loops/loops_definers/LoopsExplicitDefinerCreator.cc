// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsExplicitDefinerCreator.cc
/// @brief  LoopsExplicitDefinerCreator for the LoosExplicitDefiner load-time factory registration scheme
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Unit Headers
#include <protocols/loops/loops_definers/LoopsExplicitDefinerCreator.hh>
#include <protocols/loops/loops_definers/LoopsExplicitDefiner.hh>

using std::string;

namespace protocols {
namespace loops {
namespace loops_definers {

LoopsExplicitDefinerCreator::LoopsExplicitDefinerCreator() {}
LoopsExplicitDefinerCreator::~LoopsExplicitDefinerCreator() {}

LoopsDefinerOP
LoopsExplicitDefinerCreator::create_loops_definer() const {
	return LoopsDefinerOP( new LoopsExplicitDefiner );
}

string
LoopsExplicitDefinerCreator::type_name() const {
	// return "Loops";
	return LoopsExplicitDefiner::class_name();
}

void LoopsExplicitDefinerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopsExplicitDefiner::provide_xml_schema( xsd );
}

} //namespace
} //namespace
} //namespace
