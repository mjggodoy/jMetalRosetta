// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_pivot_pickers_StandardPivots_HH
#define INCLUDED_protocols_kinematic_closure_pivot_pickers_StandardPivots_HH

// Unit headers
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

/// @brief Randomly pick pivots in a way that may or may not span the whole
/// loop.
///
/// @details This is the default pivot picking algorithm used by
/// samplers::KicMover.  It works well for most applications, and should only
/// really need to be swapped out in in favor of custom pivot pickers that take
/// into account specific knowledge of the loop being sampled.  Earlier
/// versions of this algorithm were biased towards the fount half of the loop.
/// This bias has been mitigated, but not completely removed, in the current
/// version.

class StandardPivots : public PivotPicker {

public:
	/// @brief Default constructor.
	StandardPivots() : counter_(0) {}

public:
	/// @copydoc PivotPicker::get_name
	std::string get_name() const { return "StandardPivots"; }

	/// @copydoc PivotPicker::pick
	Loop pick(Pose const & pose, Loop const & loop);

private:
	Size counter_;

};

}
}
}

#endif

