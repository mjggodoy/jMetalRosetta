// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_simple_calculators_SasaCalculatorLegacy_HH
#define INCLUDED_core_pose_metrics_simple_calculators_SasaCalculatorLegacy_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {

class SasaCalculatorLegacy : public core::pose::metrics::StructureDependentCalculator {

public:

	SasaCalculatorLegacy( core::Real probe_r = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ) : probe_radius_(probe_r) {}

	core::pose::metrics::PoseMetricCalculatorOP clone() const { return core::pose::metrics::PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() ); };

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::Real total_sasa_;
	core::id::AtomID_Map< core::Real > atom_sasa_;
	utility::vector1< core::Real > residue_sasa_;
	//utility::vector1< core::Real > residue_hsasa_; //Hydrophobic
	core::Real probe_radius_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace simple_calculators
} // namespace metrics
} // namespace pose
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_metrics_simple_calculators_SasaCalculatorLegacy )
#endif // SERIALIZATION


#endif //INCLUDED_core_pose_metrics_simple_calculators_SasaCalculatorLegacy_HH
