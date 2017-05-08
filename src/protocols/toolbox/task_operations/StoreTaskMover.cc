// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/StoreTaskMover.cc
/// @brief The StoreTaskMover allows you to create a PackerTask at some point during a
/// RosettaScripts run and save it for access later during the same run. Can be useful
/// for mutating/analyzing a particular set of residues using many different movers/filters.
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/StoreTaskMover.hh>
#include <protocols/toolbox/task_operations/StoreTaskMoverCreator.hh>

//project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <protocols/toolbox/task_operations/STMStoredTask.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.StoreTaskMover" );

namespace protocols {
namespace toolbox {
namespace task_operations {

// @brief default constructor
StoreTaskMover::StoreTaskMover() {}

// @brief destructor
StoreTaskMover::~StoreTaskMover() {}

void
StoreTaskMover::apply( core::pose::Pose & pose )
{

	// Create the PackerTask using the TaskFactory created in parse_my_tag
	runtime_assert( task_factory_ != 0 );
	core::pack::task::PackerTaskOP  task = task_factory_->create_task_and_apply_taskoperations( pose );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	}
	// If the pose doesn't have STM_STORED_TASK data, put a blank STMStoredTask in there.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		protocols::toolbox::task_operations::STMStoredTaskOP blank_tasks( new protocols::toolbox::task_operations::STMStoredTask() );
		pose.data().set( core::pose::datacache::CacheableDataType::STM_STORED_TASKS, blank_tasks );
	}
	// Grab a reference to the data
	protocols::toolbox::task_operations::STMStoredTask & stored_tasks = *( utility::pointer::static_pointer_cast< protocols::toolbox::task_operations::STMStoredTask > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) );
	// If you haven't set overwrite to true and your task name already exists, fail. Otherwise, put the task you've made into the data cache.
	if ( overwrite_ || !stored_tasks.has_task(task_name_) ) {
		stored_tasks.set_task( task, task_name_ );
	} else {
		utility_exit_with_message("A stored task with the name " + task_name_ + " already exists; you must set overwrite flag to true to overwrite." );
	}
}

void
StoreTaskMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data_map );
	task_name_ = tag->getOption< std::string >( "task_name" );
	overwrite_ = tag->getOption< bool >( "overwrite", false );

}

// @brief Identification
// XRW TEMP std::string StoreTaskMoverCreator::keyname() const { return StoreTaskMover::mover_name(); }
// XRW TEMP std::string StoreTaskMover::mover_name() { return "StoreTaskMover"; }
// XRW TEMP std::string StoreTaskMover::get_name() const { return "StoreTaskMover"; }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StoreTaskMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new StoreTaskMover );
// XRW TEMP }

protocols::moves::MoverOP
StoreTaskMover::clone() const {
	return protocols::moves::MoverOP( new StoreTaskMover( *this ) );
}

protocols::moves::MoverOP
StoreTaskMover::fresh_instance() const {
	return protocols::moves::MoverOP( new StoreTaskMover );
}

std::string StoreTaskMover::get_name() const {
	return mover_name();
}

std::string StoreTaskMover::mover_name() {
	return "StoreTaskMover";
}

void StoreTaskMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "task_name", xs_string, "XRW_TODO" )
		+ XMLSchemaAttribute::attribute_w_default( "overwrite", xsct_rosetta_bool, "XRW TO DO", "false" );

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string StoreTaskMoverCreator::keyname() const {
	return StoreTaskMover::mover_name();
}

protocols::moves::MoverOP
StoreTaskMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StoreTaskMover );
}

void StoreTaskMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StoreTaskMover::provide_xml_schema( xsd );
}


} // task_operations
} // toolbox
} // protocols

