// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/Pose.fwd.hh
/// @brief  Pose forward declarations header
/// @author Phil Bradley

#ifndef INCLUDED_protocols_stepwise_rna_RNA_AnalyticLoopCloser_FWD_HH
#define INCLUDED_protocols_stepwise_rna_RNA_AnalyticLoopCloser_FWD_HH


//Auto Headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

class RNA_AnalyticLoopCloser;
typedef utility::pointer::shared_ptr< RNA_AnalyticLoopCloser > RNA_AnalyticLoopCloserOP;

} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
