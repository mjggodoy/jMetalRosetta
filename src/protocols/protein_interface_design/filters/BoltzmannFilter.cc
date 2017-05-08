// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/BoltzmannFilter.hh>
#include <protocols/protein_interface_design/filters/BoltzmannFilterCreator.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <math.h>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
using ObjexxFCL::format::F;


namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.BoltzmannFilter" );

/// @brief default ctor
BoltzmannFilter::BoltzmannFilter() :
	temperature_( 0.6 ),
	fitness_threshold_( 0.0 ),
	triage_threshold_( -9999 ),
	norm_neg_( false )
{
	positive_filters_.clear();
	negative_filters_.clear();
	anchors_.clear();
}


core::Real
BoltzmannFilter::fitness_threshold() const{
	return fitness_threshold_;
}

void
BoltzmannFilter::fitness_threshold( core::Real const f ){
	fitness_threshold_ = f;
}

core::Real
BoltzmannFilter::temperature() const{
	return temperature_;
}

void
BoltzmannFilter::temperature( core::Real const t ){
	temperature_ = t;
}

core::Real
BoltzmannFilter::triage_threshold() const{
	return triage_threshold_;
}

void
BoltzmannFilter::triage_threshold( core::Real const t ){
	triage_threshold_ = t;
}

void
BoltzmannFilter::norm_neg( bool const n ){
	norm_neg_ = n;
}

bool
BoltzmannFilter::norm_neg() const{
	return norm_neg_;
}


utility::vector1< protocols::filters::FilterOP >
BoltzmannFilter::get_positive_filters() const{
	return positive_filters_;
}

utility::vector1< protocols::filters::FilterOP >
BoltzmannFilter::get_negative_filters() const{
	return negative_filters_;
}

void
BoltzmannFilter::add_positive_filter( protocols::filters::FilterOP f ){
	positive_filters_.push_back( f );
}

void
BoltzmannFilter::add_negative_filter( protocols::filters::FilterOP f ){
	negative_filters_.push_back( f );
}

void
BoltzmannFilter::anchors( utility::vector1< core::Real > const & anchors ){
	anchors_ = anchors;
}

utility::vector1< core::Real >
BoltzmannFilter::anchors() const{
	return anchors_;
}

bool
BoltzmannFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const fitness( compute( pose ) );
	return( fitness <= fitness_threshold() );
}

/// NOTICE that this returns -Fitness [-1:0] for use in optimization
/// F = Sum_{+}( -E/T ) / [ Sum_{-}( -E/T ) + Sum_{+} ( -E/T ) + Sum_{+}(( E - anchor )/T) ]
/// This is the standard fitness function, except for anchor. Anchor can be (but doesn't have to be)
/// defined for each positive state and sets a threshold above which energy increases in the positive
/// state substantially reduce fitness, irrespective of what happened to all negative states.
/// Can be used to ensure that the stability of a target state is not compromised during design.
/// Set this to a very large number (99999) to eliminate the effects of the anchor, or specify no anchors at all
core::Real
BoltzmannFilter::compute( core::pose::Pose const & pose ) const{
	using protocols::filters::FilterCOP;

	core::Real positive_sum( 0.0 ), negative_sum( 0.0 );
	core::Size negative_counter( 0 );
	std::string s = "BOLTZ: ";
	for ( core::Size index = 1; index <= get_positive_filters().size(); ++index ) {
		core::Real const filter_val( get_positive_filters()[ index ]->report_sm( pose ));
		s += F(7,3,filter_val)+" ";
		positive_sum += exp( -filter_val / temperature() );
		if ( anchors_.size() >= index ) {
			negative_sum += exp( ( filter_val - anchors_[ index ] ) / temperature() );
		}
	}

	for ( FilterCOP filter : get_negative_filters() ) {
		core::Real filter_val = filter->report_sm( pose );
		s += F(7,3,filter_val)+" ";
		if ( filter_val >= triage_threshold() ) {
			negative_sum += exp( -filter_val / temperature() );
			negative_counter += 1;
			TR<<"Taken filter: "<<filter->get_user_defined_name()<<" filter val: "<< filter_val<<std::endl;
		}
	}
	TR << s << -positive_sum/(positive_sum+negative_sum) <<std::endl;
	if ( norm_neg() ) {
		TR<<"Negative counter: "<< negative_counter <<std::endl;
		// if there are no negative states then set a value that will definitely return fail in compute
		if ( !negative_counter ) {
			TR<<"Normalized fitness: 9999 "<<std::endl;
			return 9999;
		}
		TR<<"Normalized fitness: " << (( -(core::Real)positive_sum / ( positive_sum + negative_sum )) / ((core::Real) get_positive_filters().size()/(get_positive_filters().size()+negative_counter)))<<std::endl;
		TR<<"Number of positive states: " << get_positive_filters().size() << std::endl;
		return ( ( -(core::Real)positive_sum / ( positive_sum + negative_sum )) / ((core::Real) get_positive_filters().size()/(get_positive_filters().size()+negative_counter) ) );
	}
	TR<<"Fitness: " << ( -(core::Real)positive_sum / ( positive_sum + negative_sum )) << std::endl;
	return( -(core::Real)positive_sum / ( positive_sum + negative_sum ));
}

core::Real
BoltzmannFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
BoltzmannFilter::report( std::ostream &, core::pose::Pose const & ) const
{
}

void
BoltzmannFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	TR << "BoltzmannFilter"<<std::endl;
	runtime_assert( tag->hasOption( "anchors" ) || tag->hasOption( "negative_filters" ) );
	fitness_threshold( tag->getOption< core::Real >( "fitness_threshold", 0 ) );
	temperature( tag->getOption< core::Real >( "temperature", 0.6 ) );
	triage_threshold( tag->getOption< core::Real >( "triage_threshold", -9999 ) );
	norm_neg( tag->getOption< bool >( "norm_neg", false ) );
	utility::vector1< std::string > const positive_filter_names( utility::string_split( tag->getOption< std::string >( "positive_filters" ), ',' ) );
	utility::vector1< std::string > negative_filter_names, anchors_string;
	negative_filter_names.clear(); anchors_string.clear();
	if ( tag->hasOption( "negative_filters" ) ) {
		negative_filter_names = utility::string_split( tag->getOption< std::string >( "negative_filters" ), ',' );
	}
	if ( tag->hasOption( "anchors" ) ) {
		anchors_string = utility::string_split( tag->getOption< std::string >( "anchors"), ',' );
	}
	for ( std::string const & positive_filter_name : positive_filter_names ) {
		add_positive_filter( protocols::rosetta_scripts::parse_filter( positive_filter_name, filters ) );
	}
	for ( std::string const & negative_filter_name : negative_filter_names ) {
		add_negative_filter( protocols::rosetta_scripts::parse_filter( negative_filter_name, filters ) );
	}
	for ( std::string const & anchor_str : anchors_string ) {
		anchors_.push_back( (core::Real) utility::string2float( anchor_str ) );
	}

	TR<<"with options temperature: "<<temperature()<<"  triage_threshold "<<triage_threshold()<<" fitness_threshold "<<fitness_threshold()<<"  "<<get_positive_filters().size()<<" positive and "<<get_negative_filters().size()<<" negative filters."<<std::endl;
	if ( anchors().size() > 0 ) {
		TR<<"defined "<<anchors().size()<<" anchors"<<std::endl;
		runtime_assert( get_positive_filters().size() == anchors().size());
	}
}

protocols::filters::FilterOP
BoltzmannFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new BoltzmannFilter() );
}

BoltzmannFilter::~BoltzmannFilter(){}

protocols::filters::FilterOP
BoltzmannFilter::clone() const{
	return protocols::filters::FilterOP( new BoltzmannFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP BoltzmannFilterCreator::create_filter() const { return protocols::filters::FilterOP( new BoltzmannFilter ); }

// XRW TEMP std::string
// XRW TEMP BoltzmannFilterCreator::keyname() const { return "Boltzmann"; }

std::string BoltzmannFilter::name() const {
	return class_name();
}

std::string BoltzmannFilter::class_name() {
	return "Boltzmann";
}

void BoltzmannFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "fitness_threshold", xsct_real, "Threshold below which fitness must fall", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "Boltzmann temperature for simulation (kT)", "0.6" )
		+ XMLSchemaAttribute::attribute_w_default( "triage_threshold", xsct_real, "Above this value (e.g. delta score/delta ddg), a negative state will be counted in the Boltzmann fitness calculation. This prevents the dilution of negative states", "-9999" )
		+ XMLSchemaAttribute::attribute_w_default( "norm_neg", xsct_rosetta_bool, "normalize the fitness of the mutation state in relative to the original state. When triage_threshold is used the number of negative states is changed, therefore norm_neg is needed in order to compare mutations in the same position.", "false" )
		+ XMLSchemaAttribute( "positive_filter_names", xs_string, "Names of positive filters" )
		+ XMLSchemaAttribute( "negative_filter_names", xs_string, "Names of negative filters" )
		+ XMLSchemaAttribute( "anchors", xsct_real_cslist, "Threshold energy values above which positive states really hurt" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Compute the boltzmann-weighted sum of positive and negative filter values.", attlist );
}

std::string BoltzmannFilterCreator::keyname() const {
	return BoltzmannFilter::class_name();
}

protocols::filters::FilterOP
BoltzmannFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new BoltzmannFilter );
}

void BoltzmannFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BoltzmannFilter::provide_xml_schema( xsd );
}


} // filters
} // protein_interface_design
} // protocols
