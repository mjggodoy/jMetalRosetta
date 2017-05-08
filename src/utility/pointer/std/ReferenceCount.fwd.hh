// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/ReferenceCount.fwd.hh
/// @brief  utility::pointer::ReferenceCount forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_pointer_std_ReferenceCount_fwd_hh
#define INCLUDED_utility_pointer_std_ReferenceCount_fwd_hh

#ifdef PTR_STD

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace pointer {


// Forward
class ReferenceCount;

typedef shared_ptr< ReferenceCount > ReferenceCountOP;
typedef shared_ptr< ReferenceCount const > ReferenceCountCOP;

} // namespace pointer
} // namespace utility

#endif // PTR_STD

#endif // INCLUDED_utility_pointer_std_ReferenceCount_FWD_HH
