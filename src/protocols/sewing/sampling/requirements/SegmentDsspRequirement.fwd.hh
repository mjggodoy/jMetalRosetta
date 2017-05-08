// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SegmentDsspRequirement.fwd.hh
///
/// @brief
/// @author Frank Teets



#ifndef INCLUDED_protocols_sewing_sampling_requirements_SegmentDsspRequirement_FWD_HH
#define INCLUDED_protocols_sewing_sampling_requirements_SegmentDsspRequirement_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

class SegmentDsspRequirement;
typedef utility::pointer::shared_ptr< SegmentDsspRequirement > SegmentDsspRequirementOP;
typedef utility::pointer::shared_ptr< SegmentDsspRequirement const > SegmentDsspRequirementCOP;

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
