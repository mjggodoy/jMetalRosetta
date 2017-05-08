// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceSasaFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Jacob Bale (balej@.washington.edu)

#include <algorithm>
#include <boost/lexical_cast.hpp>

#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_filters/InterfaceSasaFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/MetricValue.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
// Jacob
#include <utility/string_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Project Headers
#include <utility/excn/Exceptions.hh>
#include <core/types.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.InterfaceSasaFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP InterfaceSasaFilterCreator::create_filter() const { return protocols::filters::FilterOP( new InterfaceSasaFilter ); }

// XRW TEMP std::string
// XRW TEMP InterfaceSasaFilterCreator::keyname() const { return "Sasa"; }


InterfaceSasaFilter::InterfaceSasaFilter() :
	Filter( "Sasa" ),
	lower_threshold_( 0.0 ),
	upper_threshold_(100000000000.0),
	hydrophobic_( false ),
	polar_( false ),
	jumps_(),
	sym_dof_names_()
{
	jump(1);
}

InterfaceSasaFilter::InterfaceSasaFilter( core::Real const lower_threshold, bool const hydrophobic/*=false*/, bool const polar/*=false*/, core::Real upper_threshold, std::string const & sym_dof_names ) :
	Filter( "Sasa" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold),
	hydrophobic_( hydrophobic ),
	polar_( polar ),
	jumps_(),
	sym_dof_names_()
{
	if ( sym_dof_names != "" ) {
		this->sym_dof_names(sym_dof_names);
	} else {
		jump(1);
	}
}

InterfaceSasaFilter::~InterfaceSasaFilter()= default;

filters::FilterOP
InterfaceSasaFilter::clone() const{
	return filters::FilterOP( new InterfaceSasaFilter( *this ) );
}

filters::FilterOP
InterfaceSasaFilter::fresh_instance() const{
	return filters::FilterOP( new InterfaceSasaFilter );
}

void
InterfaceSasaFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &,moves::Movers_map const &, core::pose::Pose const & pose )
{
	lower_threshold_ = tag->getOption<core::Real>( "threshold", 800 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 1000000);

	std::string specified_jumps = tag->getOption< std::string >( "jump", "" );
	std::string specified_sym_dof_names = tag->getOption< std::string >( "sym_dof_names", "" );

	if ( specified_jumps != "" && specified_sym_dof_names != "" ) {
		TR.Error << "Can not specify 'jump' and 'sym_dof_names' in InterfaceSasaFilter" << tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Can not specify 'jump' and 'sym_dof_names' in InterfaceSasaFilter" );
	} else if ( specified_jumps != "" ) {
		// Populate jumps_ with str->int converstions of the jump list.
		TR.Debug << "Reading jump list: " << specified_jumps << std::endl;

		jumps_.resize(0);

		utility::vector1<std::string> jump_strings = utility::string_split( specified_jumps, ',' );
		for ( core::Size i = 1; i <= jump_strings.size(); ++i ) {
			jumps_.push_back( boost::lexical_cast<core::Size>(jump_strings[i]));
		}

		sym_dof_names_.resize(0);
	} else if ( specified_sym_dof_names != "" ) {
		TR.Debug << "Reading sym_dof_name list: " << specified_sym_dof_names << std::endl;

		jumps_.resize(0);
		sym_dof_names_ = utility::string_split( specified_sym_dof_names, ',');
	} else if ( core::pose::symmetry::is_symmetric( pose ) && core::pose::symmetry::is_multicomponent( pose) ) {
		// get the sym_dof_names for all slidable dofs
		std::set<core::Size> sym_aware_jump_ids;
		Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
		for ( Size j = 1; j <= nslidedofs; j++ ) {
			sym_aware_jump_ids.insert(core::pose::symmetry::get_sym_aware_jump_num(pose, j ));
		}
		for ( core::Size sym_aware_jump_id : sym_aware_jump_ids ) {
			sym_dof_names_.push_back(core::pose::symmetry::jump_num_sym_dof(pose,sym_aware_jump_id));
		}
	} else {
		TR.Debug << "Defaulting to jump 1. " << std::endl;

		jump(1);
	}

	hydrophobic_ = tag->getOption<bool>( "hydrophobic", false );
	polar_ = tag->getOption<bool>( "polar", false );

	if ( polar_ && hydrophobic_ ) {
		TR.Error << "Polar and hydrophobic flags specified in Sasa filter: " << tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Polar and hydrophobic flags specified in Sasa filter." );
	}

	if ( ( polar_ || hydrophobic_ ) && (jumps_.size() != 1 || jumps_[1] != 1) ) {
		TR.Error << "Only total sasa is supported across a jump other than 1. Remove polar and hydrophobic flags and try again: " << tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Only total sasa is supported across a jump other than 1. Remove polar and hydrophobic flags and try again." );
	}

	TR.Debug << "Parsed Sasa Filter: <Sasa" <<
		" threshold=" << lower_threshold_ <<
		" upper_threshold=" << upper_threshold_ <<
		" jump=";
	std::copy(jumps_.begin(), jumps_.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
	TR.Debug <<
		" sym_dof_names=";
	std::copy(sym_dof_names_.begin(), sym_dof_names_.end(), std::ostream_iterator<std::string>(TR.Debug, ","));
	TR.Debug <<
		" hydrophobic=" << hydrophobic_ <<
		" polar=" << polar_ <<
		" />" << std::endl;
}

bool
InterfaceSasaFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ) );

	TR<<"sasa is "<<sasa<<". ";
	if ( sasa >= lower_threshold_ && sasa <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
InterfaceSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	out<<"Sasa= "<< sasa<<'\n';
}

core::Real
InterfaceSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	return( sasa );
}

void
InterfaceSasaFilter::jump( core::Size const jump )
{
	jumps_.resize(0);
	jumps_.push_back(jump);
	sym_dof_names_.resize(0);
}

void
InterfaceSasaFilter::add_jump( core::Size const jump )
{
	jumps_.push_back(jump);
}

void InterfaceSasaFilter::jumps( utility::vector1<core::Size> const & jumps )
{
	jumps_ = jumps;
}

void
InterfaceSasaFilter::sym_dof_names( std::string const & sym_dof_names )
{
	sym_dof_names_ = utility::string_split(sym_dof_names, ',');
	jumps_.resize(0);
}

void
InterfaceSasaFilter::add_sym_dof_name( std::string const & sym_dof_name )
{
	sym_dof_names_.push_back(sym_dof_name);
}

void
InterfaceSasaFilter::sym_dof_names( utility::vector1<std::string> const & sym_dof_names )
{
	sym_dof_names_ = sym_dof_names;
}

core::Real
InterfaceSasaFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose::metrics;
	using basic::MetricValue;
	using namespace core;
	using namespace protocols::moves;

	bool symm = core::pose::symmetry::is_symmetric( pose );
	if ( symm && core::pose::symmetry::is_multicomponent( pose) ) {
		Real buried_sasa = 0;
		for ( Size i = 1; i <= sym_dof_names_.size(); i++ ) {
			core::pose::Pose full_component_subpose = core::pose::symmetry::get_full_intracomponent_subpose(pose, sym_dof_names_[i]);
			core::pose::Pose full_component_neighbors_subpose = core::pose::symmetry::get_full_intracomponent_neighbor_subpose(pose, sym_dof_names_[i]);
			core::pose::Pose full_component_and_neighbors_subpose = core::pose::symmetry::get_full_intracomponent_and_neighbor_subpose(pose, sym_dof_names_[i]);
			//full_component_subpose.dump_pdb("full_component_subpose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");
			//full_component_neighbors_subpose.dump_pdb("full_component_neighbors_subpose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");
			//full_component_and_neighbors_subpose.dump_pdb("full_component_and_neighbors_subpose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");

			MetricValue< core::Real > mv_sasa;
			full_component_subpose.metric( "sasa", "total_sasa", mv_sasa);
			core::Real const full_component_sasa( mv_sasa.value() );
			TR << "sasa for component controlled by sym_dof " << sym_dof_names_[i] << " = " << full_component_sasa << std::endl;

			full_component_neighbors_subpose.metric( "sasa", "total_sasa", mv_sasa);
			core::Real const full_component_neighbors_sasa( mv_sasa.value() );
			TR << "sasa for the neighboring subunits of the component controlled by sym_dof " << sym_dof_names_[i] << " = " << full_component_neighbors_sasa << std::endl;

			full_component_and_neighbors_subpose.metric( "sasa", "total_sasa", mv_sasa);
			core::Real const full_component_and_neighbors_sasa( mv_sasa.value() );
			TR << "sasa for component controlled by sym_dof and its neighboring subunits" << sym_dof_names_[i] << " = " << full_component_and_neighbors_sasa << std::endl;
			utility::vector1<Size> sym_dof_subunits = core::pose::symmetry::get_jump_name_to_subunits(pose, sym_dof_names_[i]);
			TR << "number of subunits = " << sym_dof_subunits.size() << std::endl;
			core::Real individual_component_buried_sasa = ((full_component_sasa + full_component_neighbors_sasa) - full_component_and_neighbors_sasa) / (2*sym_dof_subunits.size());
			TR << "individual_component_buried_sasa = " << individual_component_buried_sasa << std::endl;
			buried_sasa = individual_component_buried_sasa + buried_sasa;
			TR << "buried_sasa = " << buried_sasa << std::endl;
		}
		return( buried_sasa );
	} else {
		core::pose::Pose split_pose( pose );
		std::set<core::Size> sym_aware_jump_ids;

		for ( Size i = 1; i <= sym_dof_names_.size(); i++ ) {
			TR.Debug << "getting sym_aware_jump_id from " << sym_dof_names_[i] << std::endl;
			sym_aware_jump_ids.insert(core::pose::symmetry::sym_dof_jump_num( split_pose, sym_dof_names_[i] ));
		}

		if ( !symm ) {
			for ( Size i = 1; i <= jumps_.size(); i++ ) {
				TR.Debug << "getting sym_aware_jump_id from " << jumps_[i] << std::endl;
				sym_aware_jump_ids.insert(jumps_[i]);
			}
		} else {
			// all slidable jumps
			Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
			for ( Size j = 1; j <= nslidedofs; j++ ) {
				sym_aware_jump_ids.insert( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
			}
		}

		runtime_assert( !sym_aware_jump_ids.empty() );

		for ( Size const sym_aware_jump_id : sym_aware_jump_ids ) {
			TR.Debug << "Moving jump id: " << sym_aware_jump_id << std::endl;
			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( split_pose, sym_aware_jump_id ) );
			translate->step_size( 1000.0 );
			translate->apply( split_pose );
		}
		//split_pose.dump_pdb("Interface_SASA_split.pdb");
		runtime_assert( !hydrophobic_ || !polar_ );
		if ( !hydrophobic_ && !polar_ ) {
			MetricValue< core::Real > mv_sasa;

			pose.metric( "sasa", "total_sasa", mv_sasa);
			core::Real const bound_sasa( mv_sasa.value() );
			TR.Debug << "Bound sasa: " << bound_sasa << std::endl;

			split_pose.metric( "sasa", "total_sasa", mv_sasa );
			core::Real const unbound_sasa( mv_sasa.value() );
			TR.Debug << "Unbound sasa: " << unbound_sasa << std::endl;

			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				//fpd  this logic is not correct for lattice symmetries
				core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
				core::Real const buried_sasa( (unbound_sasa - bound_sasa) /(sym_info->subunits()));

				TR.Debug << "Normalizing calculated sasa by subunits: " << sym_info->subunits() << std::endl;
				TR.Debug << "Buried sasa: " << buried_sasa << std::endl;

				return( buried_sasa );
			} else {
				core::Real const buried_sasa(unbound_sasa - bound_sasa);

				TR.Debug << "Buried sasa: " << buried_sasa << std::endl;
				return( buried_sasa );
			}
		} else {
			MetricValue< id::AtomID_Map< core::Real > > atom_sasa;
			pose.metric( "sasa_interface", "delta_atom_sasa", atom_sasa );
			core::Real polar_sasa( 0.0 ), hydrophobic_sasa( 0.0 );
			for ( core::Size pos(1); pos<=pose.size(); ++pos ) {
				for ( core::Size atomi( 1 ); atomi <= atom_sasa.value().n_atom( pos ); ++atomi ) {
					core::Real const atomi_delta_sasa( atom_sasa.value()( pos, atomi ) );
					core::conformation::Residue const pos_rsd( pose.residue( pos ) );
					core::chemical::AtomType const atom_type( pos_rsd.atom_type( atomi ) );
					bool const is_polar( atom_type.is_donor() || atom_type.is_acceptor() || atom_type.is_polar_hydrogen() );
					if ( is_polar ) polar_sasa += atomi_delta_sasa;
					else hydrophobic_sasa += atomi_delta_sasa;
				}
			}
			if ( hydrophobic_ ) return hydrophobic_sasa;
			if ( polar_ ) return polar_sasa;
		}
	}
	utility_exit_with_message( "Execution should never have reached this point." );
	return( 0 );
}

std::string InterfaceSasaFilter::name() const {
	return class_name();
}

std::string InterfaceSasaFilter::class_name() {
	return "Sasa";
}

void InterfaceSasaFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "lower level threshold for filter", "800")
		+ XMLSchemaAttribute::attribute_w_default("upper_threshold", xsct_real, "upper level threshold for filter", "1000000")
		+ XMLSchemaAttribute("jump", xs_string, "jump number")
		+ XMLSchemaAttribute("sym_dof_names", xs_string, "symmetric degrees of freedom")
		+ XMLSchemaAttribute::attribute_w_default("hydrophobic", xsct_rosetta_bool, "hydrophobic?", "false")
		+ XMLSchemaAttribute::attribute_w_default("polar", xsct_rosetta_bool, "polar?", "false");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filters based upon interface SASA using either jumps or symmetric degrees of freedom", attlist );
}

std::string InterfaceSasaFilterCreator::keyname() const {
	return InterfaceSasaFilter::class_name();
}

protocols::filters::FilterOP
InterfaceSasaFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new InterfaceSasaFilter );
}

void InterfaceSasaFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceSasaFilter::provide_xml_schema( xsd );
}


}
}
