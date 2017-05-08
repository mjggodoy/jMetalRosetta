// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/PackstatCalculator.cc
/// @brief  packstat calculator class
/// @author Florian Richter

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>


#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

static THREAD_LOCAL basic::Tracer TR( "protocols/toolbox/PoseMetricCalculators/PackstatCalculator" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {


PackstatCalculator::PackstatCalculator(
	core::Size oversample,
	bool remove_nonprotein_res
) : total_packstat_(0),
	special_region_packstat_(0),
	oversample_(oversample),
	remove_nonprotein_res_(remove_nonprotein_res)
{
	special_region_.clear();
	residue_packstat_.clear();
}


PackstatCalculator::PackstatCalculator(
	std::set< core::Size > const & special_region,
	core::Size oversample,
	bool remove_nonprotein_res
) : total_packstat_(0),
	special_region_packstat_(0),
	oversample_(oversample),
	remove_nonprotein_res_(remove_nonprotein_res),
	special_region_( special_region )
{
	residue_packstat_.clear();
}


void
PackstatCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{

	if ( key == "total_packstat" ) {
		basic::check_cast( valptr, &total_packstat_, "total_packstat expects to return a real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_packstat_ );

	} else if ( key == "special_region_packstat" ) {
		basic::check_cast( valptr, &special_region_packstat_, "special_region_packstat expects to return a real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( special_region_packstat_ );

	} else if ( key == "residue_packstat" ) {
		basic::check_cast( valptr, &residue_packstat_, "residue_packstat expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_packstat_ );

	} else {
		basic::Error() << "PackstatCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


std::string
PackstatCalculator::print( std::string const & key ) const
{

	if ( key == "total_packstat" ) {
		return utility::to_string( total_packstat_ );
	} else if ( key == "special_region_packstat" ) {
		return utility::to_string( special_region_packstat_ );
	} else if ( key == "residue_packstat" ) {
		return utility::to_string( residue_packstat_ );
	}

	basic::Error() << "PackstatCalculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

} //print


/// @brief this function doesn't actually recompute anything by itself, but calls the
/// @brief stuff in the packstat code
void
PackstatCalculator::recompute( Pose const & this_pose )
{
	using namespace core::scoring::packstat;

	if ( remove_nonprotein_res_ ) {

		bool has_nonprot_res(false);

		for ( core::Size i = 1; i <= this_pose.size(); ++i ) {
			if ( ! this_pose.residue_type(i).is_protein() ) {
				has_nonprot_res = true;
				break;
			}
		}

		if ( has_nonprot_res ) {

			PoseOP pureprotpose( new Pose( this_pose ) );

			pose_manipulation::remove_non_protein_residues( *pureprotpose );

			total_packstat_ = compute_packing_score( *pureprotpose, oversample_ );
			residue_packstat_ = compute_residue_packing_scores( *pureprotpose, oversample_ );
			runtime_assert( pureprotpose->size() == residue_packstat_.size() );
		} else {
			total_packstat_ = compute_packing_score( this_pose, oversample_ );
			residue_packstat_ = compute_residue_packing_scores( this_pose, oversample_ );
			runtime_assert( this_pose.size() == residue_packstat_.size() );
		}
	} else {
		total_packstat_ = compute_packing_score( this_pose, oversample_ );
		residue_packstat_ = compute_residue_packing_scores( this_pose, oversample_ );
		runtime_assert( this_pose.size() == residue_packstat_.size() );
	}


	core::Real respackstat_sum(0.0);
	core::Real special_region_sum(0.0);
	for ( Size i = 1; i <= residue_packstat_.size(); ++i ) {

		respackstat_sum = respackstat_sum + residue_packstat_[i];
		if ( special_region_.find( i ) != special_region_.end() ) special_region_sum = special_region_sum + residue_packstat_[i];

	}

	//runtime_assert( total_packstat_ == ( respackstat_sum / residue_packstat_.size() ) );
	//std::cerr << "total packstat is " << total_packstat_ << ", average residue packstat is " << respackstat_sum / residue_packstat_.size() << std::endl;

	special_region_packstat_ = ( special_region_sum / special_region_.size() );

} //recompute


} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::PackstatCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( total_packstat_ ) ); // core::Real
	arc( CEREAL_NVP( special_region_packstat_ ) ); // core::Real
	arc( CEREAL_NVP( residue_packstat_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( oversample_ ) ); // core::Size
	arc( CEREAL_NVP( remove_nonprotein_res_ ) ); // _Bool
	arc( CEREAL_NVP( special_region_ ) ); // std::set<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::PackstatCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( total_packstat_ ); // core::Real
	arc( special_region_packstat_ ); // core::Real
	arc( residue_packstat_ ); // utility::vector1<core::Real>
	arc( oversample_ ); // core::Size
	arc( remove_nonprotein_res_ ); // _Bool
	arc( special_region_ ); // std::set<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::PackstatCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::PackstatCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_PackstatCalculator )
#endif // SERIALIZATION
