// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DbDataType.fwd.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_basic_database_schema_generator_DbDataType_FWD_HH
#define INCLUDED_basic_database_schema_generator_DbDataType_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace basic {
namespace database {
namespace schema_generator {

class DbDataType;
typedef utility::vector1< DbDataType > DbDataTypes;

typedef utility::pointer::shared_ptr< DbDataType > DbDataTypeOP;
typedef utility::pointer::shared_ptr< DbDataType const > DbDataTypeCOP;

} //namesapce
} //namespace
} //namespace

#endif
