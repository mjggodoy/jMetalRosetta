// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/Quaternion.fwd.hh
/// @brief  numeric::Quaternion forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_Quaternion_fwd_hh
#define INCLUDED_numeric_Quaternion_fwd_hh


namespace numeric {


// Forward
template< typename > class Quaternion;


// Types
typedef  Quaternion< float >        Quaternion_float;
typedef  Quaternion< double >       Quaternion_double;
typedef  Quaternion< long double >  Quaternion_longdouble;


} // namespace numeric


#endif // INCLUDED_numeric_Quaternion_FWD_HH

