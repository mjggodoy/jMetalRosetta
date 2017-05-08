// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ResidualDipolarCouplingRigidSegments.cc
/// @brief  Uses NMR RDC for scoring
/// @author Nikolas Sgourakis
/// @author Oliver Lange
/// @author Srivatsan Raman

//Unit headers

#include <protocols/scoring/ResidualDipolarCouplingRigidSegments.hh>

// Package headers

// Project headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <basic/options/option.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>


// ObjexxFCL headers

//C++ headers
#include <iostream>
#include <string>
//#include <iostream.h>

/// Utility headers
#include <utility/vector1.hh>

// option key includes
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring.ResidualDipolarCouplingRigidSegments" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {
using namespace core;
using namespace scoring;

//typedef utility::vector1<core::scoring::RDC> RDC_lines;
typedef utility::vector1< core::scoring::ResidualDipolarCoupling::RDC_lines > RDC_lines_collection;
//typedef core::Real Tensor[3][3];
//typedef core::Real rvec[3];


extern void store_RDC_segments_in_pose(ResidualDipolarCouplingRigidSegmentsOP rdcrs_info,
	core::pose::Pose& pose) {
	//using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA;
	pose.data().set(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA, rdcrs_info);
}

extern ResidualDipolarCouplingRigidSegmentsCOP retrieve_RDC_segments_from_pose(
	core::pose::Pose const& pose) {
	//using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA;
	if ( pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA) ) {
		return utility::pointer::static_pointer_cast< protocols::scoring::ResidualDipolarCouplingRigidSegments const > ( pose.data().get_const_ptr(
			core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA) );
	};
	return nullptr;
}

extern ResidualDipolarCouplingRigidSegmentsOP retrieve_RDC_segments_from_pose(core::pose::Pose& pose) {
	//using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA;
	if ( pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA) ) {
		return utility::pointer::static_pointer_cast< protocols::scoring::ResidualDipolarCouplingRigidSegments > ( pose.data().get_ptr(
			core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA) );
	};
	return nullptr;
}


Real ResidualDipolarCouplingRigidSegments::compute_pairwise_score() const{
	Real score(0);
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::rdc::segment_scoring_mode ].user() ) {
		tr.Warning << "no al.tensor scoring method specified" << std::endl;
		return score;
	}

	if ( option[ OptionKeys::rdc::segment_scoring_mode ]() == "pairwise"  ) {
		core::Size n_of_exps = rdc_segments_[1]->get_n_alignments();
		//Loop over all experiments caution: indexing starts at 0
		for ( Size exp=0; exp < n_of_exps; ++exp ) {
			for ( auto it1=rdc_segments_.begin(); it1 != rdc_segments_.end(); ++it1 ) {
				for ( auto it2=it1; (it2+1) != rdc_segments_.end(); ++it2 ) {
					utility::vector0< ResidualDipolarCoupling::Tensor > & S1 = (*it1)->tensor();
					utility::vector0< ResidualDipolarCoupling::Tensor > & S2 = (*it2)->tensor();
					Real dot = (S1[exp][0][0] * S2[exp][0][0]) + (S1[exp][0][1] * S2[exp][0][1]) +
						(S1[exp][0][2] * S2[exp][0][2]) + (S1[exp][1][1] * S2[exp][1][1]) + (S1[exp][1][2] * S2[exp][1][2]);
					score += dot;

					/*Real Delta_FA(0);
					Delta_FA = std::abs ( (*it1)->get_fractional_anisotropy(exp) - (*it2)->get_fractional_anisotropy(exp) );
					score += Delta_FA;*/
				}
			}
		}
		return score;
	}

	if ( option[ OptionKeys::rdc::segment_scoring_mode ]() == "fixed_sum"  ) {

		core::Size n_of_exps = rdc_segments_[1]->get_n_alignments();
		core::Size ndata_segment (0), ndata(0);
		core::Real score_temp;
		core::Size id;

		for ( Size exp=0; exp < n_of_exps; ++exp ) {
			ndata = 0;
			score_temp = 0;
			for ( auto const & rdc_segment : rdc_segments_ ) {
				Real trace(0);
				ndata_segment = 0;
				//loop over RDCs in segment and find number of data points for this experiment
				RDC_lines myrdcs = rdc_segment->get_RDC_data();
				for ( RDC_lines::const_iterator it = myrdcs.begin(); it!=myrdcs.end(); ++it ) {
					id = it->expid();
					if ( id == exp ) { ++ndata_segment; }
				}

				trace = rdc_segment->get_al_tensor_trace(exp);
				score_temp +=  ( 3 - trace ) * ndata_segment;
				ndata += ndata_segment;
			}

			score += score_temp / ndata;
		} //close loop over exps


		return score;
	}

	if ( option[ OptionKeys::rdc::segment_scoring_mode ]() == "fixed_axis_z"  ) {

		core::Size n_of_exps = rdc_segments_[1]->get_n_alignments();
		core::Size ndata_segment (0), ndata(0);
		core::Real score_temp;
		core::Size id;

		for ( Size exp=0; exp < n_of_exps; ++exp ) {
			ndata = 0;
			score_temp = 0;
			for ( auto const & rdc_segment : rdc_segments_ ) {
				Real maxz(0);
				ndata_segment = 0;
				//loop over RDCs in segment and find number of data points for this experiment
				RDC_lines myrdcs = rdc_segment->get_RDC_data();
				for ( RDC_lines::const_iterator it = myrdcs.begin(); it!=myrdcs.end(); ++it ) {
					id = it->expid();
					if ( id == exp ) { ++ndata_segment; }
				}

				maxz = rdc_segment->get_al_tensor_max_z(exp);
				score_temp +=  ( 1 - maxz ) * ndata_segment;
				ndata += ndata_segment;
			}

			score += score_temp / ndata;
		} //close loop over exps


		return score;
	}
	return 0;
}

Real ResidualDipolarCouplingRigidSegments::compute_total_score(core::pose::Pose const& pos) const{
	Real score(0);
	Real n_rdcs(0);
	Real total_lines(0);
	for ( auto const & rdc_segment : rdc_segments_ ) {
		n_rdcs = rdc_segment->get_RDC_data().size();
		total_lines += n_rdcs;
		score += ( rdc_segment->compute_dipscore(pos) ) * n_rdcs;
	}
	return score / total_lines;
}

Size ResidualDipolarCouplingRigidSegments::find_effective_plane(core::scoring::RDC const& line) const{
	return line.res1() > line.res2() ? line.res1() : line.res2();
}

Size ResidualDipolarCouplingRigidSegments::find_segid_from_RDC_line(core::scoring::RDC const& line) const {
	Size eff_plane(find_effective_plane(line));
	Size segid (segment_definitions_.loop_index_of_residue(eff_plane));
	//    cout <<segid << std::endl;
	return segid;
}

void ResidualDipolarCouplingRigidSegments::read_RDC_segment_file_from_cmdline(){
	read_RDC_segment_file( basic::options::option[ basic::options::OptionKeys::rdc::segment_file ]() );
}


void ResidualDipolarCouplingRigidSegments::read_RDC_segment_file(std::string const& filename){
	//make an ifstream
	std::ifstream is( filename.c_str() );

	if ( !is.good() ) {
		utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + filename + "'" );
	}

	loops::PoseNumberedLoopFileReader reader;
	reader.hijack_loop_reading_code_set_loop_line_begin_token( "RDC_SEGMENT" );
	loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(
		is, filename, false /*no strict checking */ );
	segment_definitions_ = loops::Loops( loops );
}

core::scoring::ResidualDipolarCoupling::RDC_lines ResidualDipolarCouplingRigidSegments::read_RDCs_from_cmdline()const {
	core::scoring::ResidualDipolarCoupling all_rdcs;
	return all_rdcs.get_RDC_data();
}

void ResidualDipolarCouplingRigidSegments::sort_into_segments(RDC_lines all_rdcs) {
	//  std::cout <<"calling sorter" <<std::endl;
	RDC_lines_collection rdc_segm_data( segment_definitions_.num_loop() );

	//process all lines and assign into segments

	for ( auto & all_rdc : all_rdcs ) {
		// std::cout << line_it->res1() <<std::endl;
		Size segid( find_segid_from_RDC_line( all_rdc ) );
		if ( segid >0 ) {
			rdc_segm_data[ segid ].push_back( all_rdc );
		}
	}


	//create new RDC objects for all Segments
	for ( auto & it : rdc_segm_data ) {
		rdc_segments_.push_back( core::scoring::ResidualDipolarCouplingOP( new core::scoring::ResidualDipolarCoupling( it ) ) );
	}
}

void ResidualDipolarCouplingRigidSegments::show(std::ostream& out) const {
	Size ct (0);
	for ( auto const & rdc_segment : rdc_segments_ ) {
		out <<"SEGMENT " << ++ct <<std::endl;
		rdc_segment->show(out);
		out << std::endl;
	}
}

std::ostream& operator<<(std::ostream& out, ResidualDipolarCouplingRigidSegments const& rdcrs) {
	rdcrs.show(out);
	return out;
}
} //namespace protocol
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::scoring::ResidualDipolarCouplingRigidSegments::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( rdc_segments_ ) ); // RDC_Segments
	arc( CEREAL_NVP( segment_definitions_ ) ); // protocols::loops::Loops
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::scoring::ResidualDipolarCouplingRigidSegments::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( rdc_segments_ ); // RDC_Segments
	arc( segment_definitions_ ); // protocols::loops::Loops
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::scoring::ResidualDipolarCouplingRigidSegments );
CEREAL_REGISTER_TYPE( protocols::scoring::ResidualDipolarCouplingRigidSegments )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_scoring_ResidualDipolarCouplingRigidSegments )
#endif // SERIALIZATION
