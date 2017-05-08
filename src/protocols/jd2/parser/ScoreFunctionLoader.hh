// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/parser/ScoreFunctionLoader.hh
/// @brief  Declartion of the XML parser's ScoreFunctionLoader class for adding named ScoreFunctions to the basic::datacache::DataMap
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_ScoreFunctionLoader_hh
#define INCLUDED_protocols_jd2_parser_ScoreFunctionLoader_hh

// Package Headers
#include <protocols/jd2/parser/DataLoader.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace parser {

/// @brief A class for loading arbitrary data into the XML parser's basic::datacache::DataMap.
class ScoreFunctionLoader : public DataLoader
{
public:
	ScoreFunctionLoader();
	virtual ~ScoreFunctionLoader();

	/// @brief The ScoreFunctionLoader will create named ScoreFunctions and load them into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

	static std::string loader_name();
	static std::string score_function_loader_ct_namer( std::string const & element_name );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} //namespace parser
} //namespace jd2
} //namespace protocols

#endif
