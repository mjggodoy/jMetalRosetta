// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinResidueConformationFeatures.hh
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_ResidueConformationFeatures_hh
#define INCLUDED_protocols_features_ResidueConformationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueConformationFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

class ResidueConformationFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueConformationFeatures();

	ResidueConformationFeatures(
		ResidueConformationFeatures const & );

	~ResidueConformationFeatures() override= default;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief return string with class name
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	void
	delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP) override;

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose) override;

	void
	load_conformation(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void
	set_coords_for_residue(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::Size seqpos,
		core::pose::Pose & pose);

	void
	set_coords_for_residue_from_compact_schema(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::Size seqpos,
		core::pose::Pose & pose);

private:
	bool compact_residue_schema_;


};

} // features
} // protocols

#endif // include guard
