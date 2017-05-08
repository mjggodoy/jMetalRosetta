// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/SeqprofConsensusOperation.fwd.hh
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2011

#ifndef INCLUDED_protocols_toolbox_task_operations_SeqprofConsensusOperation_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_SeqprofConsensusOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class SeqprofConsensusOperation;
class RestrictConservedLowDdgOperation;

typedef utility::pointer::shared_ptr< SeqprofConsensusOperation > SeqprofConsensusOperationOP;
typedef utility::pointer::shared_ptr< RestrictConservedLowDdgOperation > RestrictConservedLowDdgOperationOP;

} // task_operations
} // toolbox
} // protocols


#endif
