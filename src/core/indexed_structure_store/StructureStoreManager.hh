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

#ifndef INCLUDED_core_indexed_structure_store_StructureStoreManager_hh
#define INCLUDED_core_indexed_structure_store_StructureStoreManager_hh

// Project headers
#include <core/types.hh>
#include <core/indexed_structure_store/StructureStoreManager.fwd.hh>
#include <core/indexed_structure_store/FragmentLookup.fwd.hh>
#include <core/indexed_structure_store/FragmentStore.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <numeric/types.hh>
#include <string>
#include <map>
#include <utility/vector1.hh>

namespace core
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from data store and manages retrieval on indices on store.
class StructureStoreManager : public utility::SingletonBase< StructureStoreManager >
{
public:
	friend class utility::SingletonBase< StructureStoreManager >;

	typedef std::map< std::string, std::map<Size, FragmentStoreOP> > GroupedFragmentMap;
	typedef std::map< std::string, FragmentLookupOP > FragmentLookupMap; //for later
public:
	// @brief Load fragment lookup from the default structure store. Default store is defined
	// via indexed_structure_store:fragment_store option.
	FragmentLookupOP load_fragment_lookup(std::string lookup_name);

	// @brief Load fragment lookup from the provided structure store.
	//
	// store_path - Target store path. HDF5 file path (extension is optional) if
	//     HDF5 support is enabled. Root directory of binary store if not HDF5
	//     support.
	// lookup_name - Lookup name within store. Fragment lookups within store are
	//     under <store_path>/fragments/<lookup_name>
	FragmentLookupOP load_fragment_lookup(std::string lookup_name, std::string store_path);
	FragmentStoreOP load_fragment_store(std::string lookup_name, std::string store_path, utility::vector1<std::string> fields_to_load, utility::vector1<std::string> fields_to_load_types);
	std::map<numeric::Size,FragmentStoreOP> load_grouped_fragment_store(std::string group_field, std::string lookup_name, std::string store_path, utility::vector1<std::string> fields_to_load, utility::vector1<std::string> fields_to_load_types);



	std::string resolve_store_path(std::string store_path);

private:

	StructureStoreManager();
	GroupedFragmentMap grouped_fragment_map_;

};

}
}
#endif
