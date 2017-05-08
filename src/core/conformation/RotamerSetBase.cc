// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/RotamerSetBase.cc
/// @brief  RotamerSetBase class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


// Unit headers
#include <core/conformation/RotamerSetBase.hh>

// Package headers
#include <core/conformation/RotamerSetCacheableDataType.hh>

// Utility headers
#include <basic/datacache/BasicDataCache.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {

RotamerSetBase::RotamerSetBase() :
	parent(),
	data_cache_( BasicDataCacheOP( new BasicDataCache( RotamerSetCacheableDataType::num_cacheable_data_types) ) )
{}

RotamerSetBase::~RotamerSetBase() = default;

/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
RotamerSetBase::BasicDataCache &
RotamerSetBase::data()
{
	return * data_cache_;
}

/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
RotamerSetBase::BasicDataCache const &
RotamerSetBase::data() const
{
	return * data_cache_;
}


} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::RotamerSetBase::save( Archive & arc ) const {
	arc( CEREAL_NVP( data_cache_ ) ); // BasicDataCacheOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::RotamerSetBase::load( Archive & arc ) {
	arc( data_cache_ ); // BasicDataCacheOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::RotamerSetBase );
CEREAL_REGISTER_TYPE( core::conformation::RotamerSetBase )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_RotamerSetBase )
#endif // SERIALIZATION
