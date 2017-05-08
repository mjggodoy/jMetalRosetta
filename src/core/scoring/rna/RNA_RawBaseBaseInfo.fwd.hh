// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_RawBaseBasePotential.fwd.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_RawBaseBaseInfo_fwd_hh
#define INCLUDED_core_scoring_rna_RNA_RawBaseBaseInfo_fwd_hh

// C++

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.fwd.hh>

namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
class RNA_RawBaseBaseInfo;
typedef utility::pointer::shared_ptr< RNA_RawBaseBaseInfo > RNA_RawBaseBaseInfoOP;
typedef utility::pointer::weak_ptr< RNA_RawBaseBaseInfo > RNA_RawBaseBaseInfoAP;

} //rna
} //scoring
} //core

#endif
