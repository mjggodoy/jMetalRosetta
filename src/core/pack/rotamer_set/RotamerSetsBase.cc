// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSets.fwd.hh
/// @brief  Residue Sets Base class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/rotamer_set/RotamerSetsBase.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

RotamerSetsBase::RotamerSetsBase() {}
RotamerSetsBase::~RotamerSetsBase() {}


} // namespace rotamer_set
} // namespace pack
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerSetsBase::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerSetsBase::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::RotamerSetsBase );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::RotamerSetsBase )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSetsBase )
#endif // SERIALIZATION
