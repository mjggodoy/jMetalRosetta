// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NetChargeFilter.cc
/// @brief
/// @author Dave La (davela@u.washington.edu)
// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/NetChargeFilter.hh>
#include <protocols/simple_filters/NetChargeFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <sstream>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.NetChargeFilter" );

NetChargeFilter::NetChargeFilter() :
	Filter( "NetCharge" ),
	chain_( 0 ),
	net_charge_max_( 100 ),
	net_charge_min_( -100 ),
	task_factory_( /* NULL */ ) {}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP NetChargeFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NetChargeFilter ); }

// XRW TEMP std::string
// XRW TEMP NetChargeFilterCreator::keyname() const { return "NetCharge"; }

NetChargeFilter::~NetChargeFilter()= default;

void
NetChargeFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	chain_ = tag->getOption<core::Size>( "chain", 0 );
	net_charge_max_ = tag->getOption<signed int>( "max", 100 );
	net_charge_min_ = tag->getOption<signed int>( "min", -100 );
	if ( tag->hasOption( "task_operations" ) ) {
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag->getOption< std::string >( "task_operations" ), data ) );
	}

	TR<<"Net charge will be caculated for chain " << chain_ << " with maximum cutoff " << net_charge_max_ << " and minimum cutoff " << net_charge_min_ << "." << std::endl;
}

bool
NetChargeFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	TR<<"Net Charge: "<<net_charge<<". " ;
	bool const status_max = (net_charge <= net_charge_max_) ? (true) : (false);
	bool const status_min = (net_charge >= net_charge_min_) ? (true) : (false);
	bool const status = (status_max && status_min) ? (true) : (false);
	if ( status_max && status_min ) {
		TR << "passing." << std::endl;
	} else {
		TR << "failing." << std::endl;
	}
	return status;
}

void
NetChargeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	out<<"Net Charge: "<< net_charge <<'\n';
}

core::Real
NetChargeFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	return( net_charge );
}

core::Real
NetChargeFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose copy_pose = pose;

	signed int net_charge = 0;


	utility::vector1< core::Size > target_res;
	target_res.clear();

	if ( !task_factory() ) {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			target_res.push_back( i );
		}
	} else {
		target_res = protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ );
	}

	for ( core::Size const i : target_res ) {
		core::Size const chain = copy_pose.residue( i ).chain();

		// Skip if current residue is not part of the chain specified.
		// Otherwise, the default chain=0 means consider all chains.
		if ( chain_ != 0 ) {
			if ( chain != chain_ ) continue;
		}

		std::string arg_res ("ARG");
		std::string lys_res ("LYS");

		std::string asp_res ("ASP");
		std::string glu_res ("GLU");

		std::stringstream out;
		std::string cur_res;

		out << pose.aa(i);
		cur_res = out.str();

		if ( arg_res.compare(cur_res) == 0 ) {
			TR << "AA:  +1  " << cur_res << " " << i << std::endl;
			net_charge++;
		} else if ( lys_res.compare(cur_res) == 0 ) {
			TR << "AA:  +1  " << cur_res << " " << i << std::endl;
			net_charge++;
		} else if ( asp_res.compare(cur_res) == 0 ) {
			TR << "AA:  -1  " << cur_res << " " << i << std::endl;
			net_charge--;
		} else if ( glu_res.compare(cur_res) == 0 ) {
			TR << "AA:  -1  " << cur_res << " " << i << std::endl;
			net_charge--;
		}

	}

	TR << "The net charge is: " << net_charge << std::endl;

	return( net_charge );
}


core::pack::task::TaskFactoryOP
NetChargeFilter::task_factory() const{ return task_factory_; }

void
NetChargeFilter::task_factory( core::pack::task::TaskFactoryOP tf ){ task_factory_ = tf; }

std::string NetChargeFilter::name() const {
	return class_name();
}

std::string NetChargeFilter::class_name() {
	return "NetCharge";
}

void NetChargeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("chain", xsct_non_negative_integer, "specify which chain you want to calculate the net charge (In the input PDB file, from top to bottom: 1 means first chain, 2 means the second chain, and so forth). Use the value 0 (default) if you want to consider all residues in the input PDB structure.", "0")
		+ XMLSchemaAttribute::attribute_w_default("max", xs_integer, "maximum net charge desired", "100")
		+ XMLSchemaAttribute::attribute_w_default("min", xs_integer, "minimum net charge desired", "-100");

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "This filter sums up all of the positively and negatively charged amino acids in your structure and reports a simplistic sequence-based net charge.", attlist );
}

std::string NetChargeFilterCreator::keyname() const {
	return NetChargeFilter::class_name();
}

protocols::filters::FilterOP
NetChargeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new NetChargeFilter );
}

void NetChargeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NetChargeFilter::provide_xml_schema( xsd );
}



}
}
