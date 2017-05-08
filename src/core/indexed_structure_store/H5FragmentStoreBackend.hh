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
/// @author Alex Ford <fordas@uw.edu>

#ifndef INCLUDED_core_indexed_structure_store_H5FragmentStoreBackend_hh
#define INCLUDED_core_indexed_structure_store_H5FragmentStoreBackend_hh

#ifdef USEHDF5

// Utility Headers
#include <vector>
#include <string>

#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/indexed_structure_store/H5FragmentStoreBackend.fwd.hh>
#include <core/indexed_structure_store/FragmentStore.fwd.hh>


#include "H5Cpp.h"

namespace core
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from data store and manages retrieval on indices on store.
class H5FragmentStoreBackend
{
public:
	// Structure database contains structure data and structure geometry indices for a collection
	// of target structures. The data is stored in the following tables:
	//    '/fragments' =
  //    TODO alexford document
	//
	// Opens file handle
	H5FragmentStoreBackend(std::string target_filename);

	// @brief Retrieves fragment store from backend.
	FragmentStoreOP get_fragment_store(std::string store_name);
	void append_to_fragment_store(FragmentStoreOP fragment_store, std::string store_name, std::string group_field, std::string group_type);


	//H5 DataTypes for store datatypes.
	static H5::DataType FragmentThresholdDistanceEntryDatatype();
	static H5::DataType FragmentCoordinateEntryDatatype(FragmentSpecification fragment_spec);
	static H5::DataType FragmentInt64GroupEntryDatatype(std::string group_field);
	static H5::DataType FragmentRealGroupEntryDatatype(std::string group_field);
	static H5::DataType FragmentReal1PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec);
	static H5::DataType FragmentString1PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec);
	static H5::DataType FragmentString5PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec);


private:
	H5::H5File target_file_;
	std::string target_filename_;
};

}
}
#endif

#endif
