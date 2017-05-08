// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ConstraintScoreFilter.cc
/// @brief Filter that computes scores of constraints generated by ConstraintGenerators
/// @author Tom Linsky (tlinsky@uw.edu)

#include <protocols/simple_filters/ConstraintScoreFilter.hh>
#include <protocols/simple_filters/ConstraintScoreFilterCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.ConstraintScoreFilter" );

namespace protocols {
namespace simple_filters {

ConstraintScoreFilter::ConstraintScoreFilter():
	protocols::filters::Filter( "ConstraintScoreFilter" ),
	cgs_(),
	threshold_( 0.0 )
{
}

ConstraintScoreFilter::~ConstraintScoreFilter()
= default;

void
ConstraintScoreFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );

	cgs_ = protocols::constraint_generator::parse_constraint_generators( tag, data );
	if ( cgs_.empty() ) {
		std::stringstream msg;
		msg << "ConstraintScoreFilter requires the constraint_generators' option. No constraint generators are current set."
			<< std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
}

protocols::filters::FilterOP
ConstraintScoreFilter::clone() const
{
	return protocols::filters::FilterOP( new ConstraintScoreFilter( *this ) );
}


protocols::filters::FilterOP
ConstraintScoreFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new ConstraintScoreFilter );
}

std::string
ConstraintScoreFilter::get_name() const
{
	return ConstraintScoreFilter::class_name();
}

bool
ConstraintScoreFilter::apply( core::pose::Pose const & pose ) const
{
	return ( report_sm( pose ) < threshold_ );
}

core::Real
ConstraintScoreFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::pose::Pose posecopy( pose );
	posecopy.remove_constraints();
	posecopy.energies().clear();

	core::scoring::ScoreFunctionOP sfx_op( new core::scoring::ScoreFunction );
	sfx_op->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfx_op->set_weight( core::scoring::angle_constraint, 1.0 );
	sfx_op->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
	sfx_op->set_weight( core::scoring::dihedral_constraint, 1.0 );
	sfx_op->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfx_op->set_weight( core::scoring::res_type_constraint, 1.0 );

	if ( core::pose::symmetry::is_symmetric( posecopy ) ) {
		// Why does this take an OP instead of a non-const &???
		core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( posecopy, sfx_op );
	}

	protocols::constraint_generator::AddConstraints( cgs_ ).apply( posecopy );
	sfx_op->show( TR, posecopy );
	TR.flush();
	return (*sfx_op)( posecopy );
}

void
ConstraintScoreFilter::report( std::ostream &, core::pose::Pose const & ) const
{

}

/////////////// Creator ///////////////

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ConstraintScoreFilterCreator::create_filter() const
// XRW TEMP {
// XRW TEMP  return protocols::filters::FilterOP( new ConstraintScoreFilter );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ConstraintScoreFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ConstraintScoreFilter::class_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ConstraintScoreFilter::class_name()
// XRW TEMP {
// XRW TEMP  return "ConstraintScore";
// XRW TEMP }

std::string ConstraintScoreFilter::name() const {
	return class_name();
}

std::string ConstraintScoreFilter::class_name() {
	return "ConstraintScore";
}

void ConstraintScoreFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "threshold" , xsct_real , "If the computed constraint score is less than or equal to this value, the filter returns true. Otherwise, it returns false." ) ;

	protocols::constraint_generator::attributes_for_parse_constraint_generators( attlist ) ; //Names of constraint generators used to generate constraints. Constraint Generators specified here must have been previously defined in an AddConstraints mover.

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the score of a specific set of constraints generated by Constraint Generators. This filter makes a copy of the pose and clears all constraints from the copy. It then uses the provided Constraint Generators to generate a set of constraints, and scores the pose using a basic scorefunction. This basic scorefunction consists only of constraint score terms." , attlist );
}

std::string ConstraintScoreFilterCreator::keyname() const {
	return ConstraintScoreFilter::class_name();
}

protocols::filters::FilterOP
ConstraintScoreFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ConstraintScoreFilter );
}

void ConstraintScoreFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConstraintScoreFilter::provide_xml_schema( xsd );
}


} //protocols
} //simple_filters

