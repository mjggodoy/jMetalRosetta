// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/ring_conformer_io.hh
/// @brief   Database input/output function declarations for ring-conformer-specific data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_ring_conformer_io_HH
#define INCLUDED_core_chemical_rings_ring_conformer_io_HH

// Unit header
// no fwd.hh exists
namespace core { namespace chemical { namespace rings { struct RingConformer; } } }
//#include <core/chemical/rings/RingConformer.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <iosfwd>


namespace core {
namespace chemical {
namespace rings {

/// @brief  Return a list of ring conformers, read from a database file.
utility::vector1< RingConformer > read_conformers_from_database_file_for_ring_size( std::string const & filename,
	core::Size ring_size );

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_ring_conformer_io_HH
