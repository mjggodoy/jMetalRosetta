// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_mutate_relax.cc
/// @brief   Mutate a residue, then do quick relax for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPMutateRelaxMover_fwd_hh
#define INCLUDED_protocols_membrane_MPMutateRelaxMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class MPMutateRelaxMover;
typedef utility::pointer::shared_ptr< MPMutateRelaxMover > MPMutateRelaxMoverOP;
typedef utility::pointer::shared_ptr< MPMutateRelaxMover const > MPMutateRelaxMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPMutateRelaxMover_fwd_hh
