// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ShearMinCCDTrial.fwd.hh
///
/// @brief forward declaration
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_FWD_HH
#define INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

// Forward
class ShearMinCCDTrial;

typedef utility::pointer::shared_ptr< ShearMinCCDTrial > ShearMinCCDTrialOP;
typedef utility::pointer::shared_ptr< ShearMinCCDTrial const > ShearMinCCDTrialCOP;

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_FWD_HH

