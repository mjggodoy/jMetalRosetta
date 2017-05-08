// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/hbs/HbsPatcher.hh
/// @brief
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_hbs_HbsPatcher_hh
#define INCLUDED_protocols_simple_moves_hbs_HbsPatcher_hh
// Unit Headers
#include <protocols/simple_moves/hbs/HbsPatcher.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace hbs {

/// @details
class HbsPatcher : public protocols::moves::Mover {

public:

	/// @brief
	HbsPatcher( core::Size hbs_pre_position );

	virtual ~HbsPatcher();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	core::Size const hbs_pre_pos_;
	core::Size const hbs_post_pos_;
};//end HbsPatcher


}//namespace hbs
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_hbs_HbsPatcher_hh
