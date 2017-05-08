// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key4Tuple.srlz.cc
/// @brief  Serlialization routines for Key4Tuples
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/keys/Key4Tuple.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Boost headers
#include <boost/preprocessor/punctuation/comma.hpp>

namespace utility {
namespace keys {

	template < class Archive > void save( Archive & arc, Key4Tuple< platform::Size, platform::Size, platform::Size, platform::Size > const & k4t )
{ save_key4tuple( arc, k4t ); }

template < class Archive > void load( Archive & arc, Key4Tuple< platform::Size, platform::Size, platform::Size, platform::Size > & k4t )
{ load_key4tuple( arc, k4t ); }

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( Key4Tuple< platform::Size BOOST_PP_COMMA() platform::Size BOOST_PP_COMMA() platform::Size BOOST_PP_COMMA() platform::Size > );

}
}

#endif
