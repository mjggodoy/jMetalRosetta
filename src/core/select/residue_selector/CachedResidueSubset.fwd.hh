// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   CachedResidueSubset.fwd.hh
/// @brief
/// @author Tom Linsky ( tlinsky at uw dot edu )


#ifndef INCLUDED_core_select_residue_selector_CachedResidueSubset_fwd_hh
#define INCLUDED_core_select_residue_selector_CachedResidueSubset_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace residue_selector {

class CachedResidueSubset;
typedef utility::pointer::shared_ptr< CachedResidueSubset > CachedResidueSubsetOP;
typedef utility::pointer::shared_ptr< CachedResidueSubset const > CachedResidueSubsetCOP;

} // residue_selector
} // select
} // core

#endif
