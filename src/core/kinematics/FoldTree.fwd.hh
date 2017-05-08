// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/FoldTree.fwd.hh
/// @brief  kinematics::FoldTree forward declarations header
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_FoldTree_fwd_hh
#define INCLUDED_core_kinematics_FoldTree_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace core {
namespace kinematics {


// Forward
class FoldTree;

// Types
typedef  utility::pointer::weak_ptr< FoldTree >  FoldTreeAP;
typedef  utility::pointer::weak_ptr< FoldTree const >  FoldTreeCAP;
typedef  utility::pointer::shared_ptr< FoldTree >  FoldTreeOP;
typedef  utility::pointer::shared_ptr< FoldTree const >  FoldTreeCOP;


} // namespace kinematics
} // namespace core

#endif // INCLUDED_core_kinematics_FoldTree_FWD_HH
