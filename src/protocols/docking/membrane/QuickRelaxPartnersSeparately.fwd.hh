// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Mover partners apart and relax them separately
/// @details Run quick relax on separated partners; this emulates unbound docking
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_fwd_hh
#define INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace docking {
namespace membrane {

class QuickRelaxPartnersSeparately;
typedef utility::pointer::shared_ptr< QuickRelaxPartnersSeparately > QuickRelaxPartnersSeparatelyOP;
typedef utility::pointer::shared_ptr< QuickRelaxPartnersSeparately const > QuickRelaxPartnersSeparatelyCOP;

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_fwd_hh
