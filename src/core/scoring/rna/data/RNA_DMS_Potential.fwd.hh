// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/data/RNA_DMS_Potential.fwd.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_data_RNA_DMS_Potential_FWD_HH
#define INCLUDED_core_scoring_rna_data_RNA_DMS_Potential_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {
namespace data {

class RNA_DMS_Potential;
typedef utility::pointer::shared_ptr< RNA_DMS_Potential > RNA_DMS_PotentialOP;
typedef utility::pointer::shared_ptr< RNA_DMS_Potential const > RNA_DMS_PotentialCOP;

} //data
} //rna
} //scoring
} //core

#endif
