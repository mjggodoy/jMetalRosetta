// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperationCreators.hh
/// @brief  Declaration for the class that connects ResLvlTaskOperations with the ResLvlTaskOperationFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperationCreators_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperationCreators_hh

#include <core/pack/task/operation/ResLvlTaskOperationCreator.hh>

#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictToRepackingRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class RestrictAbsentCanonicalAASRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class DisallowIfNonnativeRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class PreventRepackingRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class AddBehaviorRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class IncludeCurrentRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class PreserveCBetaRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ExtraChiCutoffRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ExtraRotamersGenericRLTCreator : public ResLvlTaskOperationCreator {
public:
	virtual ResLvlTaskOperationOP create_res_level_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
