// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/HolesFilter.cc
/// @brief filter structures by will's hole value
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/simple_filters/HolesFilter.hh>
#include <protocols/simple_filters/HolesFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/packing/compute_holes_score.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.filters.HolesFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
HolesFilter::HolesFilter():
	Filter( "Holes" ),
	filtered_value_( 2.0 ),
	cmd_( "" )
{}

// @brief copy constructor
HolesFilter::HolesFilter( HolesFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_value_( rval.filtered_value_ ),
	cmd_( rval.cmd_ )
{}

// @brief set filtered value
void HolesFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
HolesFilter::Real
HolesFilter::report_sm( Pose const & pose ) const
{
	return  compute( pose );
}

/// @brief
void
HolesFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "Holes: " <<  compute( pose ) << std::endl;
}

/// @brief
HolesFilter::Real
HolesFilter::compute( Pose const & pose ) const
{
	using core::scoring::packing::HolesResult;
	using core::scoring::packing::compute_holes_score;

	Size MAX_RES = 5000;
	runtime_assert( pose.size() <= MAX_RES );

	HolesResult result = compute_holes_score( pose, cmd_ );
	return result.dec15_score;
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool HolesFilter::apply( Pose const & pose ) const
{
	Real value = compute( pose );
	if ( value < filtered_value_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
HolesFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set filtered type
	cmd_ = tag->getOption<String>( "cmd", "" );
	if ( cmd_ == "" ) {
		tr << "cmd in xml file is emptry, so using -holes::dalphaball is expected. " << std::endl;
		tr << "if both are empty, this gonna be crash. " << std::endl;
	}

	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", 2.0 );
	tr << "Structures which have holes value less than " << filtered_value_ << " will be filtered." << std::endl;
}

// XRW TEMP filters::FilterOP
// XRW TEMP HolesFilterCreator::create_filter() const { return filters::FilterOP( new HolesFilter ); }

// XRW TEMP std::string
// XRW TEMP HolesFilterCreator::keyname() const { return "Holes"; }

std::string HolesFilter::name() const {
	return class_name();
}

std::string HolesFilter::class_name() {
	return "Holes";
}

void HolesFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("cmd", xs_string, "expects you to use -holes::dalphaball", "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "Threshold value for the filter", "2.0");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates holes value(?) and filters based upon topology", attlist );
}

std::string HolesFilterCreator::keyname() const {
	return HolesFilter::class_name();
}

protocols::filters::FilterOP
HolesFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new HolesFilter );
}

void HolesFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HolesFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
