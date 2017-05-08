// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/CompositionMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_CompositionMover_hh
#define INCLUDED_protocols_moves_CompositionMover_hh

#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/CompositionMover.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

class CompositionMover : public Mover {

public:
	CompositionMover();

	void add_mover( MoverOP m );
	void clear();

	using Mover::apply;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	utility::vector1< MoverOP > get_movers();

	void apply( core::pose::Pose & pose, Size const i, Size const j);

private:

	void apply_mover_if_defined( core::pose::Pose & pose, Size const n );

private:
	utility::vector1< MoverOP > movers_;
};

} // moves
} // protocols


#endif
