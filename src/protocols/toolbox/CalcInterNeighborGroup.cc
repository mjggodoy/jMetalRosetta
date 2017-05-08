// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/CalcInterNeighborGroup.cc
/// @brief This is meant for finding interfaces between protein domains - like protein-protein interfaces but within a protein.  It's more flexible than that, though.  You define groups of residues within a protein (say, the N and C terminal domains).  You then define which pairs of groups you are interested in.  This calculator returns the union of the sets of residues at the interfaces between these domains/groups.  This calculator contains a superset of the functionality of some of the other calculators, but is less efficient in simple cases.  The pose does NOT have to have been scored.
/// @author Steven Lewis
/// @author Jared Adolf-Bryfogle (split from IGNC pose calculator)


//Unit headers
#include <protocols/toolbox/CalcInterNeighborGroup.hh>


#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <utility/graph/Graph.hh>

//Utility headers
//#include <basic/options/option.hh>
//#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>


//C++ headers
//#include <set>
//#include <utility>

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.CalcInterNeighborGroup" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {


//typedef std::set< core::Size > one_group;
//typedef std::pair< onegroup > group_pair;
//typedef utility::vector1< group_pair > group_set;

CalcInterNeighborGroup::CalcInterNeighborGroup()
{
	num_neighbors_ = 0 ;
	dist_cutoff_ = 10.0;
}

CalcInterNeighborGroup::CalcInterNeighborGroup( group_set const & groups, core::Real dist_cutoff )
: groups_(groups), dist_cutoff_(dist_cutoff), num_neighbors_(0)

{}

CalcInterNeighborGroup::CalcInterNeighborGroup( CalcInterNeighborGroup const & calculator )
: groups_(calculator.groups_), dist_cutoff_(calculator.dist_cutoff_), num_neighbors_(calculator.num_neighbors_), neighbors_(calculator.neighbors_)
{}

CalcInterNeighborGroup::~CalcInterNeighborGroup() = default;


void
CalcInterNeighborGroup::compute( core::pose::Pose const & pose )
{
	//clear old data
	neighbors_.clear();
	num_neighbors_ = 0;

	//Might be a good idea to error-check that all group residues are within the pose? - can assert < nres later?
	core::Size const nres(pose.size());

	//this is the expensive part!
	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); //create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); //create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, dist_cutoff_ ); //create edges
	runtime_assert(nres == pg->num_vertices());

	//PointGraph is the one-way graph, but this is inefficient for group v group calculations - we do not want to iterate over the entire graph each time.  Instead we want to visit just the nodes in one group and see if its edges are in the second group, so we need a two-way graph to prevent reiterating the lower half every time.
	utility::graph::Graph neighborgraph(nres);
	for ( core::Size r(1); r <= nres; ++r ) {
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			neighborgraph.add_edge(r, edge_iter->upper_vertex());
		}
	}
	runtime_assert(nres == neighborgraph.num_nodes());
	runtime_assert(pg->num_edges() == neighborgraph.num_edges());

	//iterating through the graph is somewhat less expensive.  We will need to iterate once per group pair (domain pair)
	//for each group/domain pair
	for ( core::Size i(1), vecsize(groups_.size()); i <= vecsize; ++i ) {
		//for the first member of the group/domain, iterate through its residues
		for ( auto it(groups_[i].first.begin()), end(groups_[i].first.end()); it != end; ++it ) {
			//for all edges of that node
			for ( utility::graph::Graph::EdgeListConstIter edge_iter = neighborgraph.get_node(*it)->const_edge_list_begin(),
					edge_end_iter = neighborgraph.get_node(*it)->const_edge_list_end();
					edge_iter != edge_end_iter; ++edge_iter ) {
				core::Size const other = (*edge_iter)->get_other_ind(*it);
				//at this point, *it and other are neighbors.  *it is in the "first" group, we need to see if other is in the second.
				if ( groups_[i].second.find(other) != groups_[i].second.end() ) {
					// *it was in group 1 and other was in group 2 - store them!
					neighbors_.insert(*it);
					neighbors_.insert(other);
				} //if these are cross-group neighbors
			}//for all edges of a node
			//we also need to check if a residue is in both groups at once - it is its own neighbor, so this makes it part of the set
			if ( groups_[i].second.find(*it) != groups_[i].second.end() ) neighbors_.insert(*it);
		}// for all residues in a group
	}//for all group pairs

	num_neighbors_ = neighbors_.size();

	return;

} //compute



} //namespace toolbox
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::CalcInterNeighborGroup::save( Archive & arc ) const {
	arc( CEREAL_NVP( groups_ ) ); // group_set
	arc( CEREAL_NVP( dist_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( num_neighbors_ ) ); // core::Size
	arc( CEREAL_NVP( neighbors_ ) ); // std::set<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::CalcInterNeighborGroup::load( Archive & arc ) {
	arc( groups_ ); // group_set
	arc( dist_cutoff_ ); // core::Real
	arc( num_neighbors_ ); // core::Size
	arc( neighbors_ ); // std::set<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::CalcInterNeighborGroup );
CEREAL_REGISTER_TYPE( protocols::toolbox::CalcInterNeighborGroup )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_CalcInterNeighborGroup )
#endif // SERIALIZATION
