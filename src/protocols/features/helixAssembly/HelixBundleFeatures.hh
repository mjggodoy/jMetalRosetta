// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file HelixBundleFeatures.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh
#define INCLUDED_protocols_features_helixAssembly_HelixBundleFeatures_hh

//Unit
#include <protocols/features/helixAssembly/HelixBundleFeatures.fwd.hh>

//Protocols
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

//Devel
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

struct FragmentPair
{
	FragmentPair(HelicalFragmentOP frag1, HelicalFragmentOP frag2):
		fragment_1(frag1),
		fragment_2(frag2)
	{}

	HelicalFragmentOP fragment_1;
	HelicalFragmentOP fragment_2;

	core::Real end_1_distance;
	core::Real end_2_distance;
	core::Real fa_attr;
	core::Real fa_fraction;
	core::Real crossing_angle;

	bool
	operator == ( FragmentPair const & other ) const
	{
		return fragment_1 == other.fragment_1 && fragment_2 == other.fragment_2;
	}

	friend
	bool
	operator <(
		FragmentPair const & a,
		FragmentPair const & b
	)
	{
		return a.fragment_1 < b.fragment_1;
	}
};

typedef std::map< std::pair<core::Size, core::Size>, FragmentPair> PairMap;

class HelixBundleFeatures : public protocols::features::FeaturesReporter
{

public:

	HelixBundleFeatures();

	// XRW TEMP  virtual
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const  {
	// XRW TEMP   return "HelixBundleFeatures";
	// XRW TEMP  };

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session
	) const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	) override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	utility::vector1<HelicalFragmentOP>
	get_helices(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	bool
	validate_bundle(
		core::pose::Pose const & pose,
		utility::vector1<HelicalFragment> const & helices
	);

	void
	write_bundle_to_db(
		protocols::features::StructureID const & struct_id,
		utility::sql_database::sessionOP db_session,
		utility::vector1<HelicalFragment> const & bundle
	);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	//Maximum allowable distance, in angstroms, between any two helix ends
	core::Real helix_cap_dist_cutoff_;

	//number of helices in the bundle
	core::Size num_helices_per_bundle_;

	//number of residues in each helix
	core::Size min_helix_size_;

	//minimum residue_normalized fa_attr between two helix fragments in order to be considered
	//interacting
	core::Real min_per_residue_fa_attr_;

	//minimum fraction of
	core::Real min_interacting_set_fraction_;

	protocols::simple_filters::InterfaceSasaFilter sasa_filter_;

	core::scoring::ScoreFunctionOP scorefxn_;
};

} //namespace helixAssembly
} //namespace features
} //namespace protocols

#endif
