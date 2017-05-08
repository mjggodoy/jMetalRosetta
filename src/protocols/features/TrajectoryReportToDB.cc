// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/TrajectoryReportToDB.cc
///
/// @brief  Report features data to database multiple times per structure, creating a trajectory
/// @author Kyle Barlow (kb@kylebarlow.com)

// #ifdef USEMPI
// #include <mpi.h>
// #endif

#include <protocols/features/TrajectoryReportToDB.hh>
#include <protocols/features/TrajectoryMapFeatures.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
// Setup Mover
#include <protocols/features/TrajectoryReportToDBCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/database/sql_utils.hh>
#include <utility/tag/Tag.hh>

// Platform Headers
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace features {

// Constructors

TrajectoryReportToDB::TrajectoryReportToDB() :
	ReportToDB("TrajectoryReportToDB"),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	core::Size stride
) :
	ReportToDB("TrajectoryReportToDB"),
	stride_(stride),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	std::string const & name
) :
	ReportToDB(name),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	utility::sql_database::sessionOP db_session,
	std::string const & batch_name,
	std::string const & batch_description,
	bool use_transactions,
	core::Size cache_size
) :
	ReportToDB(
	"TrajectoryReportToDB",
	db_session, batch_name, batch_description,
	use_transactions, cache_size ),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
}

TrajectoryReportToDB::TrajectoryReportToDB(
	TrajectoryReportToDB const & src ) :
	ReportToDB(src),
	stride_(src.stride_),
	trajectory_map_features_reporter_(src.trajectory_map_features_reporter_)
{}

TrajectoryReportToDB::~TrajectoryReportToDB()= default;

// Functions

void
TrajectoryReportToDB::set_stride(
	core::Size stride
) {
	stride_ = stride;
}

core::Size
TrajectoryReportToDB::get_stride () const {
	return stride_;
}

std::map<std::string, core::Size>
TrajectoryReportToDB::get_cycle_counts() const {
	return cycle_counts_;
}

// XRW TEMP std::string
// XRW TEMP TrajectoryReportToDBCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return TrajectoryReportToDB::mover_name();
// XRW TEMP }

// XRW TEMP moves::MoverOP
// XRW TEMP TrajectoryReportToDBCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new TrajectoryReportToDB );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TrajectoryReportToDB::mover_name()
// XRW TEMP {
// XRW TEMP  return "TrajectoryReportToDB";
// XRW TEMP }

static THREAD_LOCAL basic::Tracer TR( "protocols.features.TrajectoryReportToDB" );

moves::MoverOP
TrajectoryReportToDB::fresh_instance() const { return moves::MoverOP( new TrajectoryReportToDB ); }

moves::MoverOP
TrajectoryReportToDB::clone() const
{
	return moves::MoverOP( new TrajectoryReportToDB( *this ) );
}

void
TrajectoryReportToDB::initialize_trajectory_reporter()
{
	trajectory_map_features_reporter_ = TrajectoryMapFeaturesOP( new TrajectoryMapFeatures() );
	ReportToDB::add_features_reporter( trajectory_map_features_reporter_ );
}

void
TrajectoryReportToDB::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose )
{
	ReportToDB::parse_my_tag(
		tag, data, filters, movers, pose
	);

	parse_stride_tag_item(tag);
}

void
TrajectoryReportToDB::parse_stride_tag_item(
	TagCOP const tag
) {
	if ( tag->hasOption("stride") ) {
		set_stride( tag->getOption<core::Size>("stride") );
	}
}

void
TrajectoryReportToDB::apply( Pose& pose )
{
	std::map<std::string, core::Size>::iterator tag_iter;
	std::string structure_tag;
	core::Size cycle_count;

	ReportToDB::ensure_structure_tags_are_ready();
	structure_tag = ReportToDB::get_structure_tag();

	tag_iter = cycle_counts_.find( structure_tag );
	if ( tag_iter == cycle_counts_.end() ) {
		// New job output tag - initialize cycle count at 0
		cycle_counts_[structure_tag] = 0;
		cycle_count = 0;
	} else {
		cycle_count = tag_iter->second;
	}

	trajectory_map_features_reporter_->set_current_cycle( cycle_count );

	if ( cycle_count % stride_ == 0 ) {
		ReportToDB::apply( pose );
	}

	// Increase cycle count
	cycle_counts_[structure_tag] += 1;
}

std::string TrajectoryReportToDB::get_name() const {
	return mover_name();
}

std::string TrajectoryReportToDB::mover_name() {
	return "TrajectoryReportToDB";
}

void TrajectoryReportToDB::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	FeaturesReporterFactory::get_instance()->define_features_reporter_xml_schema_group( xsd );
	AttributeList attlist;
	ReportToDB::attributes_for_report_to_db( attlist, xsd );
	attlist
		+ XMLSchemaAttribute( "stride", xsct_non_negative_integer, "Number of iterations between reports to database" );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & FeaturesReporterFactory::features_reporter_xml_schema_group_name );
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Reports to database every stride steps in a trajectory", attlist, subelements );
}

std::string TrajectoryReportToDBCreator::keyname() const {
	return TrajectoryReportToDB::mover_name();
}

protocols::moves::MoverOP
TrajectoryReportToDBCreator::create_mover() const {
	return protocols::moves::MoverOP( new TrajectoryReportToDB );
}

void TrajectoryReportToDBCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TrajectoryReportToDB::provide_xml_schema( xsd );
}


} // namespace
} // namespace
