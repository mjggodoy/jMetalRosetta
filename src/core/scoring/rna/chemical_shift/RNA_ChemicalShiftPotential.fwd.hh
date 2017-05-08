// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.fwd.hh
/// @brief  RNA_ChemicalShiftPotential potential class forward delcaration
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)

#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftPotential_FWD_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftPotential_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

class RNA_ChemicalShiftPotential;
typedef utility::pointer::shared_ptr< RNA_ChemicalShiftPotential > RNA_ChemicalShiftPotentialOP;
typedef utility::pointer::shared_ptr< RNA_ChemicalShiftPotential const > RNA_ChemicalShiftPotentialCOP;

} //chemical_shift
} //rna
} //scoring
} //core

#endif
