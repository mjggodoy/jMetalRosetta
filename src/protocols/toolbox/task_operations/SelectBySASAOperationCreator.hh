// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_SelectBySASAOperationCreator_hh
#define INCLUDED_protocols_toolbox_task_operations_SelectBySASAOperationCreator_hh


#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

class SelectBySASAOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
#endif

