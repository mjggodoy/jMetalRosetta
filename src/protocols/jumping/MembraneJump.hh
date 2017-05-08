// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MembraneJump
/// @brief read jump-definition file   setups fold tree an chainbreak variants
/// loop code didn't work because fold-tree to complicated ( overlapping loops )
/// @details
/// @author Bjorn Wallner


#ifndef INCLUDED_protocols_jumping_MembraneJump_hh
#define INCLUDED_protocols_jumping_MembraneJump_hh


// Unit Headers
#include <protocols/jumping/MembraneJump.fwd.hh>

// Package Headers
#include <protocols/jumping/PairingLibrary.hh>
#include <core/scoring/dssp/PairingsList.fwd.hh>


// Project Headers
#include <core/types.hh>
//#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
//#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.fwd.hh>
//#include <core/fragment/FrameList.fwd.hh>
//#include <core/fragment/FragSet.fwd.hh>
//#include <core/scoring/constraints/ConstraintForest.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


#ifdef PYROSETTA
#include <core/scoring/dssp/PairingsList.hh>
#endif


namespace protocols {
namespace jumping {


class MembraneJump : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~MembraneJump() override;
	MembraneJump();

	void
	init(std::string const& template_file,std::string const& pairings_file);

	bool
	defined() const {
		return(template_size_ > 0 && pairings_size_ > 0);
	}

	Size
	template_size() const {
		return template_size_;
	}

	Size
	pairings_size() const {
		return pairings_size_;
	}


	void
	setup_fold_tree(core::pose::Pose & pose, core::Size njumps);

	void
	rt_templates(core::pose::Pose & pose);

private:
	PairingLibrary templates_;
	core::scoring::dssp::PairingList pairings_;
	core::Size template_size_;
	core::Size pairings_size_;
	core::scoring::dssp::PairingList selected_pairings_;
};


} //jumping
} //protocols
#endif
