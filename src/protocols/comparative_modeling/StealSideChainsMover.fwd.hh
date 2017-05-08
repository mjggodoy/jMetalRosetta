// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StealSideChainsMover.fwd.hh
/// @brief forward declaration for StealSideChainsMover class.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_StealSideChainsMover_fwd_hh
#define INCLUDED_protocols_comparative_modeling_StealSideChainsMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace comparative_modeling {

class StealSideChainsMover;
typedef utility::pointer::shared_ptr< StealSideChainsMover > StealSideChainsMoverOP;
typedef utility::pointer::shared_ptr< StealSideChainsMover const > StealSideChainsMoverCOP;

} // comparative_modeling
} // protocols

#endif
