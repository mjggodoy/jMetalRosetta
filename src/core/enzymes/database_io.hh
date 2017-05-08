// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/database_io.hh
/// @brief   Database input/output function declarations for enzyme data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_enzymes_database_io_HH
#define INCLUDED_core_enzymes_database_io_HH

// Unit header
#include <core/enzymes/EnzymeData.hh>

// C++ headers
#include <string>


namespace core {
namespace enzymes {

EnzymeData read_enzyme_data_from_file( std::string const & filename );

}  // namespace enzymes
}  // namespace core

#endif  // INCLUDED_core_enzymes_database_io_HH
