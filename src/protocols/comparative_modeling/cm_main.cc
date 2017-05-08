// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/cm_main.cc
/// @author James Thompson

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/comparative_modeling/LoopRelaxThreadingMover.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.comparative_modeling.cm_main" );

namespace protocols  {
namespace comparative_modeling {

void cm_main() {
	using protocols::comparative_modeling::LoopRelaxThreadingMover;
	using protocols::comparative_modeling::LoopRelaxThreadingMoverOP;

	// initialization
	LoopRelaxThreadingMoverOP lrt( new LoopRelaxThreadingMover() );
	lrt->setup();

	// run
	protocols::jd2::JobDistributor::get_instance()->go(lrt);
}

} // comparative_modeling
} // protocols
