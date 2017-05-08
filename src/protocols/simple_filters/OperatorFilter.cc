// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/OperatorFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/simple_filters/OperatorFilterCreator.hh>
#include <protocols/simple_filters/RelativePoseFilter.hh>
#include <protocols/simple_filters/SigmoidFilter.hh>
#include <protocols/simple_filters/MultipleSigmoidsFilter.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
//JD2 headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>
#include <protocols/filters/BasicFilters.hh>
#include <limits>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.Operator" );
using namespace protocols::filters;

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP OperatorFilterCreator::create_filter() const { return protocols::filters::FilterOP( new Operator ); }

// XRW TEMP std::string
// XRW TEMP OperatorFilterCreator::keyname() const { return "Operator"; }

//default ctor
Operator::Operator() :
	protocols::filters::Filter( "Operator" ),
	operation_( PRODUCT ),
	threshold_( 0.0 ),
	negate_( false ),
	logarithm_( false ),
	report_subvalues_( false )
{
	filters_.clear();
}

Operator::~Operator() = default;

void
Operator::reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ){
	for ( protocols::filters::FilterOP filter : filters() ) {
		if ( filter->get_type() == "Sigmoid" ) {
			SigmoidOP sigmoid_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Sigmoid > ( filter ) );
			runtime_assert( sigmoid_filter != nullptr );
			sigmoid_filter->reset_baseline( pose, attempt_read_from_checkpoint );
			TR<<"Resetting Sigmoid filter's baseline"<<std::endl;
		}
		if ( filter->get_type() == "MultipleSigmoids" ) {
			MultipleSigmoidsOP multisigmoid_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::MultipleSigmoids > ( filter ) );
			runtime_assert( multisigmoid_filter != nullptr );
			multisigmoid_filter->reset_baseline( pose, attempt_read_from_checkpoint );
			TR<<"Resetting MultipleSigmoids filter's baseline"<<std::endl;
		} else if ( filter->get_type() == "Operator" ) { //recursive call
			OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( filter ) );
			runtime_assert( operator_filter != nullptr );
			operator_filter->reset_baseline( pose, attempt_read_from_checkpoint );
			TR<<"Resetting Operator filter's baseline"<<std::endl;
		} else if ( filter->get_type() == "CompoundStatement" ) { ///all RosettaScripts user-defined filters with confidence!=1 are compoundstatements
			CompoundFilterOP comp_filt_op( utility::pointer::dynamic_pointer_cast< protocols::filters::CompoundFilter > ( filter ) );
			runtime_assert( comp_filt_op != nullptr );
			for ( auto & cs_it : *comp_filt_op ) {
				protocols::filters::FilterOP f( cs_it.first );
				if ( f->get_type() == "Sigmoid" ) {
					SigmoidOP sigmoid_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Sigmoid > ( f ) );
					runtime_assert( sigmoid_filter != nullptr );
					sigmoid_filter->reset_baseline( pose, attempt_read_from_checkpoint );
					TR<<"Resetting Sigmoid filter's baseline"<<std::endl;
				} else if ( f->get_type() == "Operator" ) { //recursive call
					OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( f ) );
					runtime_assert( operator_filter != nullptr );
					operator_filter->reset_baseline( pose, attempt_read_from_checkpoint );
					TR<<"Resetting Operator filter's baseline"<<std::endl;
				}
			}//for cs_it
		}//elseif CompoundStatement
	}//foreach
}

void
Operator::modify_relative_filters_pdb_names(){
	utility::vector1< FilterOP > erase_filters;
	erase_filters.clear();
	for ( FilterOP filter : filters_ ) {
		if ( filter->get_type() != "Sigmoid" ) continue;

		SigmoidOP sigmoid_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Sigmoid > ( filter ) );
		runtime_assert( sigmoid_filter != nullptr );
		if ( sigmoid_filter->filter()->get_type() != "RelativePose" ) continue;

		TR<<"Replicating and changing RelativePose's filter pdb fname. File names: ";
		for ( std::string const fname : relative_pose_names_ ) {
			SigmoidOP new_sigmoid( sigmoid_filter );
			RelativePoseFilterOP relative_pose( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::RelativePoseFilter > ( new_sigmoid->filter() ) );
			runtime_assert( relative_pose != nullptr );
			relative_pose->pdb_name( fname );
			add_filter( new_sigmoid );
			TR<<fname<<", ";
		}
		TR<<std::endl;
		erase_filters.push_back( filter );
	}
	for ( FilterOP erase_f : erase_filters ) {
		filters_.erase( std::find( filters_.begin(), filters_.end(), erase_f ) );
	}
}


void
Operator::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	std::string const op( tag->getOption< std::string >( "operation" ) );
	if ( op=="SUM" ) {
		operation( SUM );
	}
	if ( op=="PRODUCT" ) {
		operation( PRODUCT );
	}
	if ( op=="NORMALIZED_SUM" ) {
		operation( NORMALIZED_SUM );
	}
	if ( op=="MAX" ) {
		operation( MAX );
	}
	if ( op=="MIN" ) {
		operation( MIN );
	}
	if ( op=="SUBTRACT" ) {
		operation( SUBTRACT );
	}
	if ( op=="ABS" ) {
		operation( ABS );
	}
	if ( op=="XOR" ) {
		operation( XOR );
	}
	if ( op == "BOOLEAN_OR" ) {
		operation( BOOLEAN_OR );
	}
	if ( op != "SUM" && op != "PRODUCT" && op != "NORMALIZED_SUM" && op != "MAX" && op != "MIN" && op != "SUBTRACT" && op != "ABS" && op != "BOOLEAN_OR" && op != "XOR" ) {
		utility_exit_with_message( "Operation " + op + " not recognized" );
	}
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	negate( tag->getOption< bool >( "negate", false ) );

	utility::vector1< std::string > filter_names;
	filter_names.clear();

	report_subvalues( tag->getOption< bool >( "report_subvalues", false ) );

	if ( !tag->hasOption( "filters" ) ) {
		TR<<"filters parameter not set. I expect another mover/filter to set the filters, o/w you'll crash and burn in apply! See Shira"<<std::endl;
	} else {
		filter_names = utility::string_split( tag->getOption< std::string >( "filters" ), ',' );
	}
	for ( std::string const & fname : filter_names ) {
		add_filter( protocols::rosetta_scripts::parse_filter( fname, filters ) );
		TR<<"Adding filter "<<fname<<std::endl;
	}
	if ( ( operation() == SUBTRACT || operation() == BOOLEAN_OR ) && filters_.size() != 2 ) {
		utility_exit_with_message( "Operation "+op+" requested, but the number of filters provided is different than 2. I only know how to "+op+" one filter from another" );
	}
	if ( operation() == ABS && filters_.size() != 1 ) {
		utility_exit_with_message( "Operation ABS requested, but the number of filters provided is different than 1. I only know how to take the absolute value of one filter" );
	}
	multi_relative( tag->getOption< bool >( "multi_relative", false ) );
	if ( multi_relative() ) {
		TR<<"multi_relative is set. Duplicating filters to include relative_pose_names."<<std::endl;
		relative_pose_names_ = utility::string_split( tag->getOption< std::string >( "relative_pose_names" ), ',' );
		modify_relative_filters_pdb_names();
	}
	logarithm( tag->getOption< bool >( "logarithm", false ) );
	TR<<"setting logarithm to "<<logarithm()<<std::endl;
	TR<<" using operator "<<op<<" with "<< filters_.size()<<" filters "<<std::endl;
	TR<<"report_subvalues: "<<report_subvalues()<<std::endl;
}

bool
Operator::apply( core::pose::Pose const & pose ) const {
	core::Real const val ( compute( pose ) );
	TR<<"Filter returns "<<val<<std::endl;
	if ( negate() ) {
		return( val <= threshold() );
	} else {
		return( val >= threshold() );
	}
}

void
Operator::report( std::ostream &o, core::pose::Pose const & pose ) const {
	core::Real const val = compute( pose );
	o << "Operator returns "<<val<<std::endl;
}

core::Real
Operator::report_sm( core::pose::Pose const & pose ) const {
	core::Real const val( compute( pose ) );
	if ( report_subvalues() ) {
		//add sigmoid values to the scroring file
		using protocols::jd2::JobDistributor;
		protocols::jd2::JobOP job_me( protocols::jd2::JobDistributor::get_instance()->current_job() );
		TR<<"reporting operator subvalues for: ";
		for ( FilterOP filter : filters_ ) {
			core::Real const val_local( filter->report_sm( pose ) );
			TR<<filter->get_user_defined_name()<<" with value "<<val_local<<std::endl;
			job_me->add_string_real_pair(filter->get_user_defined_name(), val_local);
		}
	}//fi report_subvalues
	TR<<"Operator returning: "<<val<<std::endl;
	return( val );
}

core::Real
Operator::compute(
	core::pose::Pose const & pose
) const {
	core::Real val( 0.0 );
	runtime_assert( filters().size() > 0 ); /// you didn't set the filters parameter in the XML, or by another mover/filter; see Shira
	if ( operation() == SUBTRACT || operation() == BOOLEAN_OR ) {
		runtime_assert( filters().size() == 2 );
		core::Real const val1( filters()[ 1 ]->report_sm( pose ) );
		core::Real const val2( filters()[ 2 ]->report_sm( pose ) );
		core::Real const difference( negate() ? val2 - val1 : val1 - val2 );
		core::Real const boolean_or( negate() ? -( val1 + val2 - val1 * val2 ) : val1 + val2 - val1 * val2 );
		TR<<"Filters' values "<<val1<<" "<<val2;
		if ( operation() == SUBTRACT ) {
			TR<<" and the difference is "<<difference<<std::endl;
			return difference;
		} else if ( operation() == BOOLEAN_OR ) {
			TR<<" the boolean or value is "<<boolean_or<<std::endl;
			return boolean_or;
		}
	}
	if ( operation() == ABS ) {
		runtime_assert( filters().size() == 1 );
		core::Real const val( filters()[ 1 ]->report_sm( pose ) );
		core::Real const abs_val( negate() ? -std::abs( val ) : std::abs( val ) );
		TR<<"Filter returns "<<val<<" and its absolute value is "<<abs_val<<std::endl;
		return abs_val;
	}
	if ( operation() == PRODUCT ) {
		val = 1.0;
	}
	if ( operation() == MIN ) {
		val = 9999999.999;
	}
	if ( operation() == MAX ) {
		val = -99999999.999;
	}
	if ( operation() == XOR ) {
		runtime_assert( filters().size() == 2 );
		core::Real const val1( filters()[ 1 ]->report_sm( pose ) );
		core::Real const val2( filters()[ 2 ]->report_sm( pose ) );
		core::Real const xor_ret(  negate() ? -(val1*(1.0-val2)+(1.0-val1)*val2) : val1*(1.0-val2)+(1.0-val1)*val2 );
		//  TR<<"XOR. val1: "<<val1<<" vals: "<<val2<<" xor: "<<xor_ret<<std::endl;
		return xor_ret;
	}

	for ( protocols::filters::FilterOP f : filters() ) {
		core::Real const filter_val( f->report_sm( pose ) );
		TR<<"Filter "<<f->get_type()<<" return "<<filter_val<<std::endl;
		if ( operation() == SUM || operation() == NORMALIZED_SUM ) {
			val += filter_val;
		}
		if ( operation() == PRODUCT ) {
			val *= filter_val;
		}
		if ( operation() == MIN ) {
			if ( filter_val <= val ) {
				val = filter_val;
			}
		}
		if ( operation() == MAX ) {
			if ( filter_val >= val ) {
				val = filter_val;
			}
		}
	}
	if ( logarithm() ) {
		TR<<"value: "<<val<<" ";
		val = log10( val + std::numeric_limits< core::Real >::epsilon() );
		TR<<"log10(val) = "<<val<<std::endl;
	}
	if ( operation() == NORMALIZED_SUM ) {
		val /= (core::Real) filters().size();
	}
	if ( negate() ) {
		val = -1.0 * val;
	}

	TR<<"Operator returns "<<val<<std::endl;
	return( val );
}

utility::vector1< protocols::filters::FilterOP >
Operator::filters() const{ return filters_; }

void
Operator::add_filter( protocols::filters::FilterOP f ){ filters_.push_back( f ); }

std::string Operator::name() const {
	return class_name();
}

std::string Operator::class_name() {
	return "Operator";
}

void Operator::attributes( utility::tag::AttributeList & attlist ) {
	using namespace utility::tag;
	attlist + XMLSchemaAttribute::required_attribute("operation", xs_string, "operation to perform: sum, product, min, or max")
		+ XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "filter threshold", "0")
		+ XMLSchemaAttribute::attribute_w_default("negate", xsct_rosetta_bool, "multiply return value by -1. Useful in optimization", "false")
		+ XMLSchemaAttribute::attribute_w_default("report_subvalues", xsct_rosetta_bool, "report subvalues?", "false")
		+ XMLSchemaAttribute("filters", xs_string, "list of previously defined filters on which to carry out the operation")
		+ XMLSchemaAttribute::attribute_w_default("multi_relative", xsct_rosetta_bool, "If set, duplicates filters to include relative_pose_names.", "false")
		+ XMLSchemaAttribute("relative_pose_names", xs_string, "Comma seperated list of pose names")
		+ XMLSchemaAttribute::attribute_w_default("logarithm", xsct_rosetta_bool, " take the log(x) value of the resulting operator.", "false");
}

void Operator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	Operator::attributes( attlist );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Part of the fuzzy-logic design algorithm.", attlist );
}

std::string OperatorFilterCreator::keyname() const {
	return Operator::class_name();
}

protocols::filters::FilterOP
OperatorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new Operator );
}

void OperatorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Operator::provide_xml_schema( xsd );
}

}
}
