// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/ResidueProbDesignOperationCreator.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_ResidueProbDesignOperationCreator_hh
#define INCLUDED_protocols_antibody_design_ResidueProbDesignOperationCreator_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>


namespace protocols {
namespace toolbox {
namespace task_operations {

// Undefined, commenting out to fix PyRosetta build
/*
class ResidueProbDesignOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
virtual std::string keyname() const { return "ResidueProbDesignOperation"; }
};
*/

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif //INCLUDED_protocols_antibody_design_ResidueProbDesignOperationCreator_hh
