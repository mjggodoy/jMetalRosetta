//-*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/SSElementBisectddGFilter.fwd.hh
/// @brief  Filter for the interconnection between secondary elements
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_SSElementBisectddGFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_SSElementBisectddGFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace simple_filters {

// Forward
class SSElementBisectddGFilter;

// Types
typedef utility::pointer::shared_ptr< SSElementBisectddGFilter >  SSElementBisectddGFilterOP;
typedef utility::pointer::shared_ptr< SSElementBisectddGFilter const >  SSElementBisectddGFilterCOP;

} // namespace protocols
} // namespace filters

#endif
