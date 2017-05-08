// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsDataCenterManagerSingleton.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

// Unit headers
#include <protocols/scoring/methods/pcs2/PcsDataCenterManagerSingleton.hh>

// Package headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsInputCenterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenter.fwd.hh>
// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static THREAD_LOCAL basic::Tracer TR_PcsDataCenterManagerSingleton( "protocols.scoring.methods.pcs.PcsDataCenterManagerSingleton" );

PcsDataCenterManagerSingleton::PcsDataCenterManagerSingleton() :
	PcsDataCenterManagerSingleton( *PcsEnergyParameterManager::get_instance() )
{}

PcsDataCenterManagerSingleton::PcsDataCenterManagerSingleton(PcsEnergyParameterManager & pcs_e_p_m){

	core::Size i_multi_data;
	core::Size n_multi_data;

	n_multi_data = pcs_e_p_m.get_n_multi_data();

	// using namespace basic::options;
	// using namespace basic::options::OptionKeys;

	for ( i_multi_data = 1; i_multi_data <= n_multi_data; ++i_multi_data ) {

		utility::vector1<std::string> vec_filename;
		utility::vector1<core::Real> vec_weight;
		vec_filename = pcs_e_p_m.get_PcsEnergyParameter_for(i_multi_data).get_vector_filename();
		vec_weight = pcs_e_p_m.get_PcsEnergyParameter_for(i_multi_data).get_vector_weight();

		core::Size start(pcs_e_p_m.get_PcsEnergyParameter_for(i_multi_data).get_include_only_start());
		core::Size end(pcs_e_p_m.get_PcsEnergyParameter_for(i_multi_data).get_include_only_end());
		core::Real individual_scale(pcs_e_p_m.get_PcsEnergyParameter_for(i_multi_data).get_individual_scale());

		if ( vec_filename.size() == 0 ) {
			utility_exit_with_message("Missing input file for PCS. Review your setup file");
		}
		PcsInputCenter pcs_i_c = PcsInputCenterManager::get_instance()->get_PcsInputCenter_for(vec_filename, vec_weight);
		PcsDataCenterOP pcs_d_c_OP;
		pcs_d_c_OP = PcsDataCenterOP( new PcsDataCenter(pcs_i_c, start, end, individual_scale) );
		(*this).get_PCS_data_all().push_back(*pcs_d_c_OP);
	}
}

utility::vector1<PcsDataCenter> &
PcsDataCenterManagerSingleton::get_PCS_data_all() {
	return (PCS_data_all_);
}

std::ostream &
operator<<(std::ostream& out, const PcsDataCenterManagerSingleton & m){
	core::Size i;

	out << "n paramagnetic center: " << m.get_n_multi_data() << std::endl;
	for ( i = 1 ; i <= m.get_n_multi_data(); ++i ) {
		out << m.PCS_data_all_[i] << std::endl;
	}
	return out;
}

core::Size
PcsDataCenterManagerSingleton::get_n_multi_data() const{
	return (PCS_data_all_.size());
}

}//namespcacs PCS
}//namespace methods
}//namespace scoring
}//namespace protocols
