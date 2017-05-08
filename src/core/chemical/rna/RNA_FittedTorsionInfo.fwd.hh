// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/rna/RNA_Info.fwd.hh
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)

#ifndef INCLUDED_core_chemical_rna_RNA_FittedTorsionInfo_fwd_HH
#define INCLUDED_core_chemical_rna_RNA_FittedTorsionInfo_fwd_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace rna {

class RNA_FittedTorsionInfo;
typedef utility::pointer::shared_ptr< RNA_FittedTorsionInfo > RNA_FittedTorsionInfoOP;

class GaussianParameter;

}
}
}

#endif
