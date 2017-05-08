// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Align_RotamerEvaluator.fwd.hh
/// @brief forward declaration for Align_RotamerEvaluator class.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_Align_RotamerEvaluator_fwd_hh
#define INCLUDED_protocols_comparative_modeling_Align_RotamerEvaluator_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace comparative_modeling {

class Align_RotamerEvaluator;
typedef utility::pointer::shared_ptr< Align_RotamerEvaluator > Align_RotamerEvaluatorOP;
typedef utility::pointer::shared_ptr< Align_RotamerEvaluator const > Align_RotamerEvaluatorCOP;

} // comparative_modeling
} // protocols

#endif
