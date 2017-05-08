// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ConcurrencyTest.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_ConcurrencyTest_hh
#define INCLUDED_protocols_features_helixAssembly_ConcurrencyTest_hh

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <protocols/features/helixAssembly/ConcurrencyTest.fwd.hh>

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
//#include <basic/options/keys/sewing.OptionKeys.gen.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

class ConcurrencyTest : public protocols::features::FeaturesReporter
{

public:

	ConcurrencyTest(){}

	// XRW TEMP  virtual
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const  {
	// XRW TEMP   return "HelixBundleFeatures";
	// XRW TEMP  }

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

}
}
}
#endif
