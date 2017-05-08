// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/AveragePathLengthFilter.fwd.hh
/// @brief Filter on the average covalent path length between residues, including disulfides
/// @author Gabriel Rocklin (grocklin@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_AveragePathLengthFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_AveragePathLengthFilter_fwd_hh
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_filters {

class AveragePathLengthFilter;

typedef utility::pointer::shared_ptr< AveragePathLengthFilter > AveragePathLengthFilterOP;
typedef utility::pointer::shared_ptr< AveragePathLengthFilter const > AveragePathLengthFilterCOP;

}
}

#endif


