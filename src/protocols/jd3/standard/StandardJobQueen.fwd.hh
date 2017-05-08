// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardJobQueen.fwd.hh
/// @brief  class declaration for StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_standard_StandardJobQueen_FWD_HH
#define INCLUDED_protocols_jd3_standard_StandardJobQueen_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace standard {

class StandardJobQueen;
class PreliminaryLarvalJob;

typedef utility::pointer::shared_ptr< StandardJobQueen > StandardJobQueenOP;
typedef utility::pointer::shared_ptr< StandardJobQueen const > StandardJobQueenCOP;

typedef utility::pointer::shared_ptr< PreliminaryLarvalJob > PreliminaryLarvalJobOP;
typedef utility::pointer::shared_ptr< PreliminaryLarvalJob const > PreliminaryLarvalJobCOP;

} // namespace standard
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_standard_StandardJobQueen_FWD_HH
