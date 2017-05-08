// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/InterfaceBindingEnergyDensityFilter.fwd.hh
/// @brief Forward declaration and typedefs for InterfaceBindingEnergyDensityFilter
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_InterfaceBindingEnergyDensityFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_InterfaceBindingEnergyDensityFilter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_filters {

class InterfaceBindingEnergyDensityFilter;

typedef utility::pointer::shared_ptr< InterfaceBindingEnergyDensityFilter > InterfaceBindingEnergyDensityFilterOP;
typedef utility::pointer::shared_ptr< InterfaceBindingEnergyDensityFilter const > InterfaceBindingEnergyDensityFilterCOP;

}
}

#endif
