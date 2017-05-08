// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StealLigandMover.fwd.hh
/// @brief forward declaration for StealLigandMover class.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_StealLigandMover_fwd_hh
#define INCLUDED_protocols_comparative_modeling_StealLigandMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace comparative_modeling {

class StealLigandMover;
typedef utility::pointer::shared_ptr< StealLigandMover > StealLigandMoverOP;
typedef utility::pointer::shared_ptr< StealLigandMover const > StealLigandMoverCOP;

} // comparative_modeling
} // protocols

#endif
