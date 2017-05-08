// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResFilterCreators.hh
/// @brief  Declaration for the class that connects ResFilters with the ResFilterFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilterCreators_hh
#define INCLUDED_core_pack_task_operation_ResFilterCreators_hh

#include <core/pack/task/operation/ResFilterCreator.hh>

#include <core/pack/task/operation/ResFilter.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class AnyResFilterCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class AllResFilterCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class NoResFilterCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueTypeFilterCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueHasPropertyCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueLacksPropertyCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResiduePDBInfoHasLabelCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResiduePDBInfoLacksLabelCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueName3IsCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueName3IsntCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueIndexIsCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResidueIndexIsntCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResiduePDBIndexIsCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ResiduePDBIndexIsntCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ChainIsCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

class ChainIsntCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};


} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
