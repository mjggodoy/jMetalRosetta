// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/oop/OopRandomPuckMover.hh
/// @brief
/// @author
#ifndef INCLUDED_protocols_simple_moves_oop_OopRandomPuckMover_hh
#define INCLUDED_protocols_simple_moves_oop_OopRandomPuckMover_hh
// Unit Headers
#include <protocols/simple_moves/oop/OopRandomPuckMover.fwd.hh>
#include <protocols/simple_moves/oop/OopMover.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace oop {

/// @details
class OopRandomPuckMover : public protocols::moves::Mover {

public:

	/// @brief
	OopRandomPuckMover(
	);

	OopRandomPuckMover( utility::vector1< core::Size > oop_seq_positions );

	virtual ~OopRandomPuckMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	utility::vector1< core::Size > const oop_seq_positions_;
	utility::vector1< std::string > available_moves_;

};//end OopRandomPuckMover


}//namespace oop
}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_oop_OopRandomPuckMover_hh
