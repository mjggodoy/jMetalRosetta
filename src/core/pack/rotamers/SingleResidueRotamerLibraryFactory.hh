// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/residue_selector/SingleResidueRotamerLibraryFactory.hh
/// @brief  Class for instantiating arbitrary SingleResidueRotamerLibrarys from a string --> SingleResidueRotamerLibraryCreator map
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_pack_rotamers_SingleResidueRotamerLibraryFactory_HH
#define INCLUDED_core_pack_rotamers_SingleResidueRotamerLibraryFactory_HH

// Package headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryCreator.fwd.hh>

// Program headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Basic headers

// Utility headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <map>
#include <string>

#ifdef MULTI_THREADED
#include <mutex>
#include <thread>
#endif

namespace core {
namespace pack {
namespace rotamers {

class SingleResidueRotamerLibraryFactory : public utility::SingletonBase< SingleResidueRotamerLibraryFactory > {
public:
	friend class utility::SingletonBase< SingleResidueRotamerLibraryFactory >;

	void factory_register( SingleResidueRotamerLibraryCreatorOP creator );

	bool has_type( std::string const & ) const;

	std::string type_for_residuetype( core::chemical::ResidueType const & restype ) const;

	/// @brief Get the SingleResidueRotamerLibrary coresponding to the given ResidueType
	/// If forcebasic is true, a SingleBasicRotamerLibrary will be returned instead of a null pointer
	core::pack::rotamers::SingleResidueRotamerLibraryCOP
	get( core::chemical::ResidueType const & , bool forcebasic = false) const;

	/// @brief Get the SingleResidueRotamerLibrary coresponding to the given ResidueType
	core::pack::rotamers::SingleResidueRotamerLibraryCOP
	get( core::chemical::ResidueType const &, core::conformation::Residue const & ) const;

private:
	std::string get_cachetag( core::chemical::ResidueType const & restype ) const;

private:
	SingleResidueRotamerLibraryFactory();
	SingleResidueRotamerLibraryFactory( SingleResidueRotamerLibraryFactory const & ); // unimplemented
	SingleResidueRotamerLibraryFactory const & operator = ( SingleResidueRotamerLibraryFactory const & ); // unimplemented

private:
	typedef std::map< std::string, SingleResidueRotamerLibraryCreatorOP > CreatorMap;
	CreatorMap creator_map_;

#ifdef MULTI_THREADED
	/// @brief The mutex for the cache_ - aquire the mutex prior to reading/writing to the cache.
	/// In practice, you shouldn't be touching the cache outside the get(ResidueType) method.
	static std::mutex cache_mutex_;
#endif
	mutable std::map< std::pair< std::string, std::string >, core::pack::rotamers::SingleResidueRotamerLibraryCOP > cache_;
};

} //namespace rotamers
} //namespace pack
} //namespace core


#endif
