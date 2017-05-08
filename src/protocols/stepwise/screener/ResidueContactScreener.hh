// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/ResidueContactScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ResidueContactScreener_HH
#define INCLUDED_protocols_stepwise_screener_ResidueContactScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/ResidueContactScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class ResidueContactScreener: public StepWiseScreener {

public:

	//constructor
	// Undefined, commenting out to fix PyRosetta build  ResidueContactScreener();

	//constructor
	ResidueContactScreener(  core::pose::Pose & screening_pose,
		core::Size const last_append_res,
		core::Size const last_prepend_res,
		core::Distance const atom_atom_overlap_dist_cutoff );

	//destructor
	~ResidueContactScreener();

public:

	bool
	check_screen();

	std::string
	name() const { return "ResidueContactScreener"; }

	StepWiseScreenerType
	type() const { return RESIDUE_CONTACT; }

private:

	core::pose::Pose & screening_pose_;
	core::Size const last_append_res_;
	core::Size const last_prepend_res_;
	core::Distance const atom_atom_overlap_dist_cutoff_;

};

} //screener
} //stepwise
} //protocols

#endif
