// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/InteractionGraphFactory.cc
/// @brief  Interation graph factory class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>

// Package headers
#include <utility/graph/Graph.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/MultiplexedAnnealableGraph.hh>

#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/DoubleLazyInteractionGraph.hh>
#include <core/pack/interaction_graph/LazyInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/SurfaceInteractionGraph.hh>
#include <core/pack/interaction_graph/HPatchInteractionGraph.hh>
#include <core/pack/interaction_graph/SymmLinMemInteractionGraph.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>

#include <core/pose/symmetry/util.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <utility/pointer/owning_ptr.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/BooleanVectorOption.hh>

namespace core {
namespace pack {
namespace interaction_graph {

static THREAD_LOCAL basic::Tracer T( "core.pack.interaction_graph.interaction_graph_factory", basic::t_info );

InteractionGraphBaseOP
InteractionGraphFactory::create_interaction_graph(
	task::PackerTask const & the_task,
	rotamer_set::RotamerSets const & rotsets,
	pose::Pose const & pose,
	scoring::ScoreFunction const & sfxn
)
{
	core::Real surface_weight( sfxn.get_weight( core::scoring::surface ) );
	core::Real hpatch_weight( sfxn.get_weight( core::scoring::hpatch ) );

	// don't use the surface or hpatch interaction graphs if we're not designing
	if ( ! the_task.design_any() ) { surface_weight = 0; hpatch_weight = 0; }

	if ( the_task.linmem_ig() ) {
		/// Symmetric OTFIGs are not currently capable of handling either the Surface or HPatch scores, so check
		/// for symmetry first and return a (pairwise-decomposable) SymmLinearMemoryInteractionGraph if requested.
		if ( pose::symmetry::is_symmetric( pose ) ) {
			T << "Instantiating SymmLinearMemoryInteractionGraph" << std::endl;
			SymmLinearMemoryInteractionGraphOP symlinmemig( new SymmLinearMemoryInteractionGraph( the_task.num_to_be_packed() ) );
			symlinmemig->set_pose( pose );
			symlinmemig->set_score_function( sfxn );
			symlinmemig->set_recent_history_size( the_task.linmem_ig_history_size() );
			return symlinmemig;
		}

		if ( surface_weight ) {
			T << "Instantiating LinearMemorySurfaceInteractionGraph" << std::endl;
			LinearMemorySurfaceInteractionGraphOP lmsolig( new LinearMemorySurfaceInteractionGraph( the_task.num_to_be_packed() ) );
			lmsolig->set_pose( pose );
			lmsolig->set_packer_task( the_task );
			lmsolig->set_score_function( sfxn );
			lmsolig->set_rotamer_sets( rotsets );
			lmsolig->set_surface_score_weight( surface_weight );
			lmsolig->set_recent_history_size( the_task.linmem_ig_history_size() );
			return lmsolig;
		}

		if ( hpatch_weight ) {
			T << "Instantiating LinearMemoryHPatchInteractionGraph" << std::endl;
			LinearMemoryHPatchInteractionGraphOP lmhig( new LinearMemoryHPatchInteractionGraph( the_task.num_to_be_packed() ) );
			lmhig->set_pose( pose );
			lmhig->set_packer_task( the_task );
			lmhig->set_score_function( sfxn );
			lmhig->set_rotamer_sets( rotsets );
			lmhig->set_score_weight( hpatch_weight );
			lmhig->set_recent_history_size( the_task.linmem_ig_history_size() );
			return lmhig;
		}

		T << "Instantiating LinearMemoryInteractionGraph" << std::endl;
		LinearMemoryInteractionGraphOP lmig( new LinearMemoryInteractionGraph( the_task.num_to_be_packed() ) );
		lmig->set_pose( pose );
		lmig->set_score_function( sfxn );
		lmig->set_recent_history_size( the_task.linmem_ig_history_size() );
		return lmig;

	} else if ( the_task.design_any() ) {

		if ( rotsets.nmoltenres() >= 1 ) { //we are altering at least one residue
			if ( rotsets.rotamer_set_for_moltenresidue(1)->num_rotamers() >= 1 ) { //and it has at least one rotamer
				if ( rotsets.rotamer_set_for_moltenresidue(1)->rotamer(1)->type().mode() != chemical::CENTROID_t ) { //and it's not centroid repacking
					if ( surface_weight ) { //Note that surface overrides lazy!
						T << "Instantiating PDSurfaceInteractionGraph" << std::endl;
						PDSurfaceInteractionGraphOP pdsig( new PDSurfaceInteractionGraph( the_task.num_to_be_packed() ) );
						pdsig->set_pose( pose );
						pdsig->set_packer_task( the_task );
						pdsig->set_rotamer_sets( rotsets );
						pdsig->set_surface_score_weight( surface_weight );
						return pdsig;

					} else if ( hpatch_weight ) {
						T << "Instantiating PDHPatchInteractionGraph" << std::endl;
						PDHPatchInteractionGraphOP hig( new PDHPatchInteractionGraph( the_task.num_to_be_packed() ) );
						hig->set_pose( pose );
						hig->set_packer_task( the_task );
						hig->set_rotamer_sets( rotsets );
						hig->set_score_weight( hpatch_weight );
						return hig;

					} else if ( the_task.lazy_ig() ) {
						T << "Instantiating LazyInteractionGraph" << std::endl;
						LazyInteractionGraphOP lazy_ig( new LazyInteractionGraph( the_task.num_to_be_packed() ) );
						lazy_ig->set_pose( pose );
						lazy_ig->set_score_function( sfxn );
						return lazy_ig;
					} else if ( the_task.double_lazy_ig() ) {
						T << "Instantiating DoubleLazyInteractionGraph" << std::endl;
						DoubleLazyInteractionGraphOP double_lazy_ig( new DoubleLazyInteractionGraph( the_task.num_to_be_packed() ) );
						double_lazy_ig->set_pose( pose );
						double_lazy_ig->set_score_function( sfxn );
						//T << "Setting DoubleLazyIngeractionGraph memory limit to " << the_task.double_lazy_ig_memlimit()  << std::endl;
						double_lazy_ig->set_memory_max_for_rpes( the_task.double_lazy_ig_memlimit() );
						return double_lazy_ig;
					} else {
						T << "Instantiating PDInteractionGraph" << std::endl;
						return InteractionGraphBaseOP( new PDInteractionGraph( the_task.num_to_be_packed() ) );
					}
				}
			}
		}
	}

	// either of the two below
	// 'linmem_ig flag is off and design is not being performed', or 'linmem_ig flag is off and centroid mode design is being performed'
	//This will also trigger if there are no rotamers
	T << "Instantiating DensePDInteractionGraph" << std::endl;
	return InteractionGraphBaseOP( new DensePDInteractionGraph( the_task.num_to_be_packed() ) );
}

InteractionGraphBaseOP
InteractionGraphFactory::create_and_initialize_two_body_interaction_graph(
	task::PackerTask const & packer_task,
	rotamer_set::RotamerSets & rotsets,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	utility::graph::GraphCOP packer_neighbor_graph)
{
	InteractionGraphBaseOP ig = create_interaction_graph( packer_task, rotsets, pose, scfxn);

	PROF_START( basic::GET_ENERGIES );
	rotsets.compute_energies( pose, scfxn, packer_neighbor_graph, ig );
	PROF_STOP( basic::GET_ENERGIES );
	T << "IG: " << ig->getTotalMemoryUsage() << " bytes" << std::endl;

	setup_IG_res_res_weights(pose, packer_task, rotsets, *ig);

	return ig;
}

AnnealableGraphBaseOP
InteractionGraphFactory::create_and_initialize_annealing_graph(
	task::PackerTask const & packer_task,
	rotamer_set::RotamerSets & rotsets,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	utility::graph::GraphCOP packer_neighbor_graph)
{
	//using namespace scoring::annealing;

	std::list< AnnealableGraphBaseOP > annealing_graphs;
	annealing_graphs.push_back(create_and_initialize_two_body_interaction_graph(packer_task, rotsets, pose, scfxn, packer_neighbor_graph));


	//TODO alexford This resolves target methods via dynamic cast of all activated WholeStructureMethods
	// in the score function. This should likely be replaced with some sort of registration mechanism.
	ResidueArrayAnnealingEvaluatorOP residue_array_evaluator( new ResidueArrayAnnealingEvaluator());
	residue_array_evaluator->initialize( scfxn, pose, rotsets, packer_neighbor_graph );
	if ( residue_array_evaluator->has_methods() ) {
		annealing_graphs.push_back( residue_array_evaluator );
	}

	if ( annealing_graphs.size() == 1 ) {
		return annealing_graphs.front();
	} else {
		return AnnealableGraphBaseOP(new MultiplexedAnnealableGraph(annealing_graphs));
	}
}

void
InteractionGraphFactory::setup_IG_res_res_weights(
	pose::Pose const & pose,
	task::PackerTask const & task,
	rotamer_set::RotamerSets const & rotsets,
	interaction_graph::InteractionGraphBase & ig
)
{
	task::IGEdgeReweightContainerCOP edge_reweights = task.IGEdgeReweights();

	if ( edge_reweights ) {

		for ( Size ii = 1; ii<= rotsets.nmoltenres(); ++ii ) {
			Size const res1id = rotsets.moltenres_2_resid( ii );

			for ( ig.reset_edge_list_iterator_for_node( ii ); !ig.edge_list_iterator_at_end(); ig.increment_edge_list_iterator() ) {

				interaction_graph::EdgeBase const & edge( ig.get_edge() );
				Size const other_node = edge.get_other_ind( ii );
				if ( other_node < ii ) {
					continue; // only deal with upper edges
				}
				Size const res2id = rotsets.moltenres_2_resid( other_node );
				ig.set_edge_weight( ii, other_node, edge_reweights->res_res_weight( pose, task, res1id, res2id ) );
			}
		}
	}

}//setup_IG_res_res_weights

}
}
}

