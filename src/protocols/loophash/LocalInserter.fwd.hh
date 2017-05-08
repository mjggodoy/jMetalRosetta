// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loophash/LoopHashLibrary.fwd.hh
/// @brief  LoopHashLibrary class forward declarations header
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_LocalInserter_fwd_hh
#define INCLUDED_protocols_loophash_LocalInserter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loophash {

class LocalInserter;
typedef utility::pointer::shared_ptr< LocalInserter > LocalInserterOP;
typedef utility::pointer::shared_ptr< LocalInserter const > LocalInserterCOP;

class LocalInserter_SimpleMin;
typedef utility::pointer::shared_ptr< LocalInserter_SimpleMin > LocalInserter_SimpleMinOP;
typedef utility::pointer::shared_ptr< LocalInserter_SimpleMin const > LocalInserter_SimpleMinCOP;

}
}

#endif

