// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueVector.hh
/// @brief  The ResidueVector class identifies a subset of residues from a Pose
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueVector_HH
#define INCLUDED_core_select_residue_selector_ResidueVector_HH

// Unit headers
#include <core/select/residue_selector/ResidueVector.fwd.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Basic headers

// Utility Headers
#include <utility/vector1.hh>

// C++ headers

namespace core {
namespace select {
namespace residue_selector {

class ResidueVector : public utility::vector1< core::Size > {
public:
	/// @brief Default constructor
	ResidueVector();

	/// @brief Constructor from ResidueSubset.
	ResidueVector( ResidueSubset const & subset );
	ResidueVector( utility::vector1< core::Size > const & vec );

	/// @brief Destructor.
	virtual ~ResidueVector();

	void from_subset( ResidueSubset const & subset );

private:
};


} //namespace residue_selector
} //namespace select
} //namespace core

#endif
