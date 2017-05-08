// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/WriteableCacheableData.hh
/// @brief
/// @author Justin Porter


#ifndef INCLUDED_basic_datacache_WriteableCacheableData_hh
#define INCLUDED_basic_datacache_WriteableCacheableData_hh

// unit headers
#include <basic/datacache/WriteableCacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>

#include <istream>
#include <ostream>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {


/// @brief base class for data storable within a DataCache
class WriteableCacheableData : public CacheableData {

public:

	~WriteableCacheableData() override;

	virtual
	void write( std::ostream &out ) const = 0;

	virtual
	std::string datatype() const = 0;

	WriteableCacheableDataOP shared_from_this();

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_WriteableCacheableData )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_WriteableCacheableData_HH */
