// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/PoseMetricCalculators/RotamerBoltzCalculator.hh
/// @brief Calculates Rotamer occupancy of each rotameric state in a given set of residues.
/// @author Hetu Kamisetty


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator_fwd_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class RotamerBoltzCalculator;

typedef utility::pointer::shared_ptr< RotamerBoltzCalculator> RotamerBoltzCalculatorOP;
typedef utility::pointer::shared_ptr< RotamerBoltzCalculator const > RotamerBoltzCalculatorCOP;

} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
