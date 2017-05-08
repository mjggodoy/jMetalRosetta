// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/LoopCycle.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_LoopCycle_HH
#define INCLUDED_core_scoring_loop_graph_LoopCycle_HH

#include <core/scoring/loop_graph/LoopCycle.fwd.hh>
#include <core/scoring/loop_graph/Loop.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace loop_graph {

class LoopCycle: public utility::pointer::ReferenceCount {

public:

	//constructor
	LoopCycle();

	LoopCycle( utility::vector1< Loop > const & loops );

	//destructor
	~LoopCycle();

public:

	Loop const & loop( Size const n ) const;

	utility::vector1< Loop > const & loops() const { return loops_; }

	Size size() const { return loops_.size(); }

	Size find_index_for_loop_landing_at_domain( Size const & takeoff_domain ) const;

	friend
	/// @brief Test IO operator for debug and Python bindings
	std::ostream & operator << ( std::ostream & os, LoopCycle const & loop_cycle);

	/// @brief a and b are the same up to circular permutation
	friend
	inline
	bool
	operator == (
		LoopCycle const & a,
		LoopCycle const & b )
	{
		if ( a.loops_.size() != b.loops_.size() ) return false;
		int const N( a.loops_.size() );
		for ( int offset = 0; offset < N; offset++ ) {
			bool match = true;
			for ( int n = 1; n <= N; n++ ) {
				int n_offset = 1 + ( ( n + offset - 1 ) % N );
				if ( !( a.loops_[ n ] == b.loops_[ n_offset ] ) ) {
					match = false;
					break;
				}
			} // n loops
			if ( match ) return true;
		}
		return false;
	}

private:

	utility::vector1< Loop > loops_;

};


} //loop_graph
} //scoring
} //core

#endif
