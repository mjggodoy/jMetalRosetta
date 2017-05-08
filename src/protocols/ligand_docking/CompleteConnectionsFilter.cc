// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/ligand_docking/CompleteConnectionsFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/CompleteConnectionsFilter.hh>
#include <protocols/ligand_docking/CompleteConnectionsFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer complete_connections_tracer( "protocols.ligand_docking.CompleteConnectionsFilter" );

bool
CompleteConnectionsFilter::apply( core::pose::Pose const & pose ) const {

	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	for ( ; start <= end; ++start ) {
		core::conformation::Residue const & res= pose.residue(start);
		complete_connections_tracer<< res.name();
		if ( res.has_incomplete_connection() ) {
			complete_connections_tracer<<" has incomplete connection"<< std::endl;
			return true;
		}
		complete_connections_tracer<<" completely connected"<< std::endl;
	}
	complete_connections_tracer<< "no more connections"<< std::endl;
	return false;
}

void
CompleteConnectionsFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "CompleteConnections" ) {
		complete_connections_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! tag->hasOption("chain") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("CompleteConnections filter needs a 'chain' option");
	}
	chain_ = tag->getOption<std::string>("chain");
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP CompleteConnectionsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new CompleteConnectionsFilter ); }

// XRW TEMP std::string
// XRW TEMP CompleteConnectionsFilterCreator::keyname() const { return "CompleteConnections"; }

std::string CompleteConnectionsFilter::name() const {
	return class_name();
}

std::string CompleteConnectionsFilter::class_name() {
	return "CompleteConnections";
}

void CompleteConnectionsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "chain", xsct_char, "XRW TO DO" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string CompleteConnectionsFilterCreator::keyname() const {
	return CompleteConnectionsFilter::class_name();
}

protocols::filters::FilterOP
CompleteConnectionsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new CompleteConnectionsFilter );
}

void CompleteConnectionsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CompleteConnectionsFilter::provide_xml_schema( xsd );
}



} // ligand_docking
} // protocols
