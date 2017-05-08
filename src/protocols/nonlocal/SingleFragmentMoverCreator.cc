// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/SingleFragmentMoverCreator.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/SingleFragmentMoverCreator.hh>

// C/C++ headers
#include <string>

// Project headers
#include <protocols/moves/Mover.fwd.hh>

// Package headers
#include <protocols/nonlocal/SingleFragmentMover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

// XRW TEMP protocols::moves::MoverOP SingleFragmentMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SingleFragmentMover() );
// XRW TEMP }

// XRW TEMP std::string SingleFragmentMoverCreator::keyname() const {
// XRW TEMP  return mover_name();
// XRW TEMP }

// XRW TEMP std::string SingleFragmentMoverCreator::mover_name() {
// XRW TEMP  return "SingleFragmentMover";
// XRW TEMP }

}  // namespace nonlocal
}  // namespace protocols
