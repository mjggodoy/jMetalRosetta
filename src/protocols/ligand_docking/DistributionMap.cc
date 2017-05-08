// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/distributions.cc
/// @brief  enumerate some distributions and map them to strings
/// @author Gordon Lemmon

// Unit headers
#include <protocols/ligand_docking/DistributionMap.hh>

// Utility headers
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace ligand_docking {

Distribution DistributionMap::operator[](std::string distribution){
	return distribution_map_[distribution];
}

DistributionMap::DistributionMap(){
	distribution_map_["uniform"]= Uniform;
	distribution_map_["gaussian"]= Gaussian;
}// private constructor

Distribution get_distribution(std::string distribution_str){
	DistributionMap* distribution_map= DistributionMap::get_instance();
	return (*distribution_map)[distribution_str];
}

} //namespace ligand_docking
} //namespace protocols
