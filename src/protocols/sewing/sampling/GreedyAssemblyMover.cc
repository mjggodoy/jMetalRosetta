// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file GreedyAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/GreedyAssemblyMover.hh>
#include <protocols/sewing/sampling/GreedyAssemblyMoverCreator.hh>

// Package Headers
#include <protocols/sewing/conformation/AssemblyFactory.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.sampling.GreedyAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP GreedyAssemblyMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP(new GreedyAssemblyMover);
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP GreedyAssemblyMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return GreedyAssemblyMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP GreedyAssemblyMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "GreedyAssemblyMover";
// XRW TEMP }

protocols::moves::MoverOP
GreedyAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new GreedyAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
GreedyAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new GreedyAssemblyMover );
}

// XRW TEMP std::string
// XRW TEMP GreedyAssemblyMover::get_name() const {
// XRW TEMP  return "GreedyAssemblyMover";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  GreedyAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

GreedyAssemblyMover::GreedyAssemblyMover():
	best_complete_assembly_(0),
	best_score_(10000),
	cycles_(1000),
	max_edges_per_node_(300)
{}

AssemblyOP
GreedyAssemblyMover::generate_assembly(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Main loop
	core::Size starttime = time(NULL);
	core::Size cur_cycle=1;

	for ( cur_cycle=1; cur_cycle <= cycles_; ++cur_cycle ) {

		//Initialize starting Assembly for the current cycle
		AssemblyOP assembly;
		if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "discontinuous" ) {
			assembly = AssemblyFactory::create_assembly("discontinuous");
		} else if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "continuous" ) {
			assembly = AssemblyFactory::create_assembly("continuous");
		}
		add_starting_model(assembly);

		TR << "Cycle " << cur_cycle << std::endl;
		while ( requirement_set_->can_be_added_to(assembly) ) {
			TR << "Looking to add edge " << assembly->path().size()+1 << std::endl;
			core::Real best_edge_score_ = 10000;
			AssemblyOP best_edge_assembly_ = 0;

			ModelNode const * reference_node = graph_->get_model_node(assembly->get_next_reference_node(graph_));
			graph_->add_edges_from_binary(edge_file_, reference_node->get_node_index());
			core::Size num_edges = reference_node->num_edges();

			//If there are no edges from the selected node, abort this trajectory and go to the next cycle
			if ( num_edges == 0 ) {
				TR << "No edges for model " << reference_node->model().model_id_ << std::endl;
				break;
			}

			core::Size max_edge_attempts = std::min(num_edges, max_edges_per_node_);

			//Randomize the edge order
			utility::vector1<core::Size> edge_order(num_edges);
			for ( core::Size i = 0; i < num_edges; ++i ) {
				edge_order[i+1]=i;
			}
			numeric::random::random_permutation(edge_order, numeric::random::rg());

			for ( core::Size cur_edge_ind=1; cur_edge_ind<=max_edge_attempts; ++cur_edge_ind ) {

				AssemblyOP pre_edge_assembly = assembly->clone();

				utility::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
				for ( core::Size j=0; j<edge_order[cur_edge_ind]; ++j ) {
					++edge_it;
				}

				//Cast the edge to a proper HashEdge, check if the new model satisfies requirements, and follow it
				HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);
				assembly->follow_edge(graph_, cur_edge, reference_node->get_node_index());
				if ( !requirement_set_->violates(assembly) ) {
					core::Real edge_score = assembly_scorefxn_->score(assembly);
					if ( edge_score < best_edge_score_ ) {
						best_edge_score_ = edge_score;
						best_edge_assembly_ = assembly->clone();
					}
				}
				assembly = pre_edge_assembly;
			}

			//If we couldn't find an edge that satisfies node requirements, go to the next cycle
			if ( best_edge_assembly_ == 0 ) { break; }

			//Revert to the best scoring Assembly for the most recent edge addition, check to see if
			//this assembly is complete. If so, check to see if it's the best one and continue on
			assembly = best_edge_assembly_;
			if ( requirement_set_->satisfies(assembly) ) {
				core::Real complete_score = assembly_scorefxn_->score(assembly);
				if ( complete_score < best_score_ ) {
					best_score_ = complete_score;
					best_complete_assembly_ = assembly->clone();
					TR << "SAVING BEST " << best_score_ << std::endl;
				}
			}
		}// While can be added to
	}//Cycles

	core::Size endtime = time(NULL);
	TR << "Completed " << cur_cycle << " cycles in " << endtime - starttime << " seconds" << std::endl;
	return best_complete_assembly_;
}

void
GreedyAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){

	parent::parse_my_tag(tag, data, filters, movers, pose);

	if ( tag->hasOption("cycles") ) {
		cycles_ = tag->getOption<core::Size>("cycles");
	}

	if ( tag->hasOption("max_edges_per_node") ) {
		max_edges_per_node_ = tag->getOption<core::Size>("max_edges_per_node");
	}
}

std::string GreedyAssemblyMover::get_name() const {
	return mover_name();
}

std::string GreedyAssemblyMover::mover_name() {
	return "GreedyAssemblyMover";
}

void GreedyAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	// TO DO: perhaps this is not the right function to call? -- also, delete this comment
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string GreedyAssemblyMoverCreator::keyname() const {
	return GreedyAssemblyMover::mover_name();
}

protocols::moves::MoverOP
GreedyAssemblyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new GreedyAssemblyMover );
}

void GreedyAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GreedyAssemblyMover::provide_xml_schema( xsd );
}


} //sewing
} //protocols
