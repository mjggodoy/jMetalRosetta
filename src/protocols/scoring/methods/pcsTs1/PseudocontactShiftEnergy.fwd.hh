// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////
///
/// @file protocols/scoring/PseudocontactShiftEnergy.fwd.hh
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcsTs1_PseudocontactShiftEnergy_fwd_hh
#define INCLUDED_protocols_scoring_methods_pcsTs1_PseudocontactShiftEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs1 {

class PCS_Energy_Ts1;
class PCS_Energy_parameters_manager_Ts1;

typedef utility::pointer::shared_ptr< PCS_Energy_Ts1 > PCS_Energy_Ts1OP;
typedef utility::pointer::shared_ptr< PCS_Energy_Ts1 const > PCS_Energy_Ts1COP;
}
}
}
}
#endif // INCLUDED_protocols_scoring_methods_PseudocontactShiftEnergy_FWD_HH
