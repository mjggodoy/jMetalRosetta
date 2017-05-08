// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/LoopRelaxMover.fwd.hh
/// @brief forward declaration for LoopRelaxMover class
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_LoopRelaxMover_fwd_hh
#define INCLUDED_protocols_comparative_modeling_LoopRelaxMover_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace comparative_modeling {

class LoopRelaxMover;

typedef  utility::pointer::shared_ptr< LoopRelaxMover >  LoopRelaxMoverOP;
typedef  utility::pointer::shared_ptr< LoopRelaxMover const >  LoopRelaxMoverCOP;


class LoopRelaxThreadingMover;

typedef  utility::pointer::shared_ptr< LoopRelaxThreadingMover >  LoopRelaxThreadingMoverOP;
typedef  utility::pointer::shared_ptr< LoopRelaxThreadingMover const >  LoopRelaxThreadingMoverCOP;


} // namespace comparative_modeling
} // namespace protocols

#endif
