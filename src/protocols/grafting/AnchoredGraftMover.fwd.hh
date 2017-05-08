// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/grafting/AnchoredGraftMover.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_grafting_AnchoredGraftMover_fwd_hh
#define INCLUDED_protocols_grafting_AnchoredGraftMover_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace grafting {

// Forward
class AnchoredGraftMover;

// Types
typedef  utility::pointer::shared_ptr< AnchoredGraftMover >  AnchoredGraftMoverOP;
typedef  utility::pointer::shared_ptr< AnchoredGraftMover const >  AnchoredGraftMoverCOP;

typedef  utility::pointer::weak_ptr< AnchoredGraftMover >  AnchoredGraftMoverAP;
typedef  utility::pointer::weak_ptr< AnchoredGraftMover const >  AnchoredGraftMoverCAP;


} // namespace grafting
} // namespace protocols

#endif // #ifndef INCLUDED_protocols/grafting_AnchoredGraftMover_fwd_hh

