// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundleHelix.fwd.hh
/// @brief  Defines owning pointers for MakeBundleHelix mover class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_MakeBundleHelix_fwd_hh
#define INCLUDED_protocols_helical_bundle_MakeBundleHelix_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace helical_bundle {

class MakeBundleHelix; // fwd declaration
typedef utility::pointer::shared_ptr< MakeBundleHelix > MakeBundleHelixOP;
typedef utility::pointer::shared_ptr< MakeBundleHelix const > MakeBundleHelixCOP;
typedef utility::vector1<MakeBundleHelixOP> MakeBundleHelixOPs;
typedef utility::vector1<MakeBundleHelixCOP> MakeBundleHelixCOPs;

} // namespace helical_bundle
} // namespace protocols

#endif // INCLUDED_protocols_helical_bundle_MakeBundleHelix_fwd_hh
