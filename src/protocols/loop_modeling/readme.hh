// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @namespace protocols::loop_modeling
///
/// @brief Framework for loop structure prediction.
///
/// @details The most important classes in this namespace are LoopModeler,
/// LoopProtocol and LoopMover.  LoopMover provides an interface for loop
/// sampling algorithms.  Methods are provided for specifying which loops need
/// to be sampled and for setting up a compatible fold tree.  Some important
/// LoopMover subclasses include RepackingRefiner, MinimizationRefiner, and
/// KicMover.  The LoopProtocol class provides a framework for combining any
/// number of loop movers into a single Monte Carlo simulation.  Aspects of the
/// simulation like the number of iterations, the score function, and the
/// temperature can be easily controlled.  LoopModeler combines initial loop
/// building, centroid refinement, and fullatom refinement into a complete loop
/// modeling simulation.  It also provides a lot of configurability, especially
/// from rosetta scripts.
///
/// Below is an excerpt from protocols::comparative_modeling::LoopRelaxMover
/// that shows some of the classes from this namespace in action.  This code
/// creates a protocol featuring KIC moves followed by a rotamer trials and
/// local minimization.  Full repacking is done once every 20 iterations.  The
/// number of iterations can be controlled from the command line, and a concise
/// progress update will be reported after every iteration.
///
/// @code
/// LoopProtocolOP protocol = new LoopProtocol;
///
/// Size repack_period = 20;
/// if (option[OptionKeys::loops::repack_period].user()) {
///   repack_period = option[OptionKeys::loops::repack_period]();
/// }
///
/// protocol->add_mover(new KicMover);
/// protocol->add_mover(new RepackingRefiner(repack_period));
/// protocol->add_mover(new RotamerTrialsRefiner);
/// protocol->add_mover(new FullatomMinimizationRefiner);
///
/// Size outer_cycles = 3;
/// Size inner_cycles = 10 * loops->loop_size();
///
/// if (option[OptionKeys::loops::outer_cycles].user()) {
///   outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
/// }
/// if (option[OptionKeys::loops::max_inner_cycles].user()) {
///   Size max_cycles = option[OptionKeys::loops::max_inner_cycles]();
///   inner_cycles = std::max(inner_cycles, max_cycles);
/// }
/// if (option[OptionKeys::loops::fast]) {
///   outer_cycles = 3;
///   inner_cycles = 12;
/// }
///
/// protocol->set_loops(*loops);
/// protocol->set_score_function(fa_scorefxn_);
/// protocol->set_iterations(outer_cycles, inner_cycles, 2);
/// protocol->apply(pose);
/// @endcode

