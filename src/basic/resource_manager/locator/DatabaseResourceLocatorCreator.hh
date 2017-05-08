// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/DatabaseResourceLocatorCreator.hh
/// @brief  Header for DatabaseResourceLocator Creator for the load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_basic_resource_manager_locator_DatabaseResourceLocatorCreator_hh
#define INCLUDED_basic_resource_manager_locator_DatabaseResourceLocatorCreator_hh

// Unit Headers
#include <basic/resource_manager/locator/DatabaseResourceLocator.fwd.hh>
#include <basic/resource_manager/ResourceLocatorCreator.hh>

#include <utility/vector1.hh>
#include <string>

namespace basic {
namespace resource_manager {
namespace locator {

/// @brief creator for the DatabaseResourceLocator class
class DatabaseResourceLocatorCreator : public ResourceLocatorCreator
{
public:
	DatabaseResourceLocatorCreator();
	virtual ~DatabaseResourceLocatorCreator();

	virtual ResourceLocatorOP create_resource_locator() const;
	virtual std::string locator_type() const;
};

} //namespace
} //namespace
} //namespace

#endif
