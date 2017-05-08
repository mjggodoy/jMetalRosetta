// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/mainchain_potential/MainchainScoreTable.hh
/// @brief  Forward declarations for a general class for storing a torsional potential for mainchain resiudes.
/// @details Can be used by terms like rama, rama_prepro, p_aa_pp.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_chemical_mainchain_potential_MainchainTorsionPotential_fwd_hh
#define INCLUDED_core_chemical_mainchain_potential_MainchainTorsionPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace mainchain_potential {

class MainchainScoreTable;

typedef utility::pointer::shared_ptr< MainchainScoreTable > MainchainScoreTableOP;
typedef utility::pointer::shared_ptr< MainchainScoreTable const > MainchainScoreTableCOP;

}
}
}

#endif
