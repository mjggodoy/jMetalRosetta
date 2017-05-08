// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/Interface_/DDPlookup.fwd.hh
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)

#ifndef INCLUDED_core_scoring_interface_DDPlookup_fwd_hh
#define INCLUDED_core_scoring_interface_DDPlookup_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace interface_ {

class DDPlookup;

typedef utility::pointer::shared_ptr< DDPlookup > DDPlookupOP;
typedef utility::pointer::shared_ptr< DDPlookup const > DDPlookupCOP;


} // Interface_
} // scoring
} // core

#endif /* INCLUDED_core_scoring_Interface_DDPLOOKUP_FWD_HH_ */
