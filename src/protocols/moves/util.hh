// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_moves_util_HH
#define INCLUDED_protocols_moves_util_HH

// Project headers
#include <core/types.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace moves {

/// @brief Searches <movers> for the named mover, returning it if it exists,
/// otherwise halts execution with an error message.
protocols::moves::MoverOP
find_mover_or_die(
	const std::string& mover_name,
	utility::tag::TagCOP tag,
	const protocols::moves::Movers_map& movers
);

/// @brief Searches <filters> for the named filter, returning it if it exists,
/// otherwise halts execution with an error message.
protocols::filters::FilterOP
find_filter_or_die(
	const std::string& filter_name,
	utility::tag::TagCOP tag,
	const protocols::filters::Filters_map& filters
);

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_UTIL_HH_
