// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/PoseCommentFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/PoseCommentFilter.hh>
#include <protocols/simple_filters/PoseCommentFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
//#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.PoseComment" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP PoseCommentFilterCreator::create_filter() const { return protocols::filters::FilterOP( new PoseComment ); }

// XRW TEMP std::string
// XRW TEMP PoseCommentFilterCreator::keyname() const { return "PoseComment"; }

//default ctor
PoseComment::PoseComment() :
	protocols::filters::Filter( "PoseComment" ),
	comment_name_( "" ),
	comment_value_( "" ),
	comment_exists_( false )
{
}

PoseComment::~PoseComment() = default;

void
PoseComment::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	comment_name( tag->getOption< std::string >( "comment_name","" ) );
	comment_value( tag->getOption< std::string >( "comment_value", "" ) );
	comment_exists( tag->getOption< bool > ( "comment_exists", false ) );

	TR<<"PoseComment with options: comment_name "<<comment_name()<<" comment_value "<<comment_value()<<" comment_exists: "<<comment_exists()<<std::endl;
}

bool
PoseComment::apply( core::pose::Pose const & pose ) const {
	core::Real const val ( compute( pose ) );
	TR<<"Pose comment "<<comment_name()<<":"<<comment_value();
	if ( val >= 0.9999 ) {
		TR<<" found"<<std::endl;
	} else {
		TR<<" not found"<<std::endl;
	}
	return( val >= 0.9999 );
}

void
PoseComment::report( std::ostream &o, core::pose::Pose const & pose ) const {
	bool const val = ( compute( pose ) >= 0.99999 );
	o << "PoseComment returns "<<val<<std::endl;
}

core::Real
PoseComment::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
PoseComment::compute(
	core::pose::Pose const & pose
) const {
	std::string val( "" );
	bool exists( false );

	using namespace core::pose;

	if ( comment_name() == "" ) {
		std::map< std::string, std::string > const comments = get_all_comments( pose );
		if ( comment_exists() && !comments.empty() ) { //size() > 0 )
			return 1.0;
		}
		if ( comments.empty() ) { //size() == 0 )
			return 0.0;
		}
		for ( auto const & comment : comments ) {// iterate over all comments and find the one with the value
			if ( comment.second == comment_value() ) {
				return 1.0;
			}
		}
		return 0.0;
	}

	exists = get_comment( pose, comment_name(), val );
	if ( comment_exists() ) {
		if ( exists ) {
			return 1.0;
		} else {
			return 0.0;
		}
	}

	if ( val == comment_value() ) {
		return 1.0;
	} else {
		return 0.0;
	}
}

std::string PoseComment::name() const {
	return class_name();
}

std::string PoseComment::class_name() {
	return "PoseComment";
}

void PoseComment::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"comment_name", xs_string,
		"the key value of the comment",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"comment_value", xs_string,
		"the comment's value",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"comment_exists", xsct_rosetta_bool,
		"check only whether the comment exists or not, regardless of its content",
		"false");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Test for the existence or the value of a comment in the pose. "
		"This is useful for controlling execution flow: if the pose comments "
		"have been modified you do one thing or another",
		attlist );
}

std::string PoseCommentFilterCreator::keyname() const {
	return PoseComment::class_name();
}

protocols::filters::FilterOP
PoseCommentFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new PoseComment );
}

void PoseCommentFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PoseComment::provide_xml_schema( xsd );
}


}
}
