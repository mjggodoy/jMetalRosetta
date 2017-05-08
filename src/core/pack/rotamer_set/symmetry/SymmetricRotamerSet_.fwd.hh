// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.fwd.hh
/// @brief  rotamer set implementation class for symmetric packing forward declarations
/// @author Ingemar Andre


#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSet__fwd_hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSet__fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

class SymmetricRotamerSet_;

typedef utility::pointer::shared_ptr< SymmetricRotamerSet_ > SymmetricRotamerSet_OP;
typedef utility::pointer::shared_ptr< SymmetricRotamerSet_ const > SymmetricRotamerSet_COP;

} // namespace symmetry
} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSet__fwd_HH
