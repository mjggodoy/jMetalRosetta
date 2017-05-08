// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh
/// @brief  Residue Set Factory class for symmetric packing
/// @author Ingemar Andre

// Unit header
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh>

// Package headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/exit.hh>

// STL Headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

RotamerSetOP
SymmetricRotamerSetFactory::create_rotamer_set( conformation::Residue const & /*res*/ )
{
	return RotamerSetOP( new SymmetricRotamerSet_() );
}

}
}
}
}
