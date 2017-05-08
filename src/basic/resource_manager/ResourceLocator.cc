// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLocater.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <utility>

namespace basic {
namespace resource_manager {

ResourceStream::~ResourceStream() = default;

ResourceLocator::ResourceLocator() :
	locator_tag_("")
{}

ResourceLocator::ResourceLocator(
	std::string const & locator_tag
) :
	locator_tag_(locator_tag)
{}

ResourceLocator::ResourceLocator(
	ResourceLocator const & src
) :
	utility::pointer::ReferenceCount(),
	locator_tag_( src.locator_tag() )
{}

ResourceLocator::~ResourceLocator() = default;

void
ResourceLocator::locator_tag(
	std::string const & locator_tag
) {
	locator_tag_ = locator_tag;
}

std::string
ResourceLocator::locator_tag() const {
	return locator_tag_;
}


} // namespace resource_manager
} // namespace basic
