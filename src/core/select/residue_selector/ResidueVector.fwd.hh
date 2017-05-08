// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueVector.fwd.hh
/// @brief  Forward declaration of a class that identifies continguous segments from a subset of residues from a Pose
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueVector_FWD_HH
#define INCLUDED_core_select_residue_selector_ResidueVector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace select {
namespace residue_selector {

class ResidueVector;

typedef utility::pointer::shared_ptr< ResidueVector > ResidueVectorOP;
typedef utility::pointer::shared_ptr< ResidueVector const > ResidueVectorCOP;
typedef utility::pointer::weak_ptr< ResidueVector > ResidueVectorAP;
typedef utility::pointer::weak_ptr< ResidueVector const > ResidueVectorCAP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
