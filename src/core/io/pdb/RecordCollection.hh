// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/io/pdb/RecordCollection.hh
/// @brief   Declarations and simple accessor/mutator definitions for RecordCollection.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_RecordCollection_HH
#define INCLUDED_core_io_pdb_RecordCollection_HH

// Unit headers
#include <core/io/pdb/RecordCollection.fwd.hh>
#include <core/io/pdb/RecordType.hh>

// Utility header
#include <utility/SingletonBase.hh>

// C++ header
#include <string>
#include <map>


namespace core {
namespace io {
namespace pdb {

/// @details  This class is a singleton and manages the definitions of PDB Records, which should only be read from the
/// database one time.
class RecordCollection : public utility::SingletonBase< RecordCollection > {
	friend class utility::SingletonBase< RecordCollection >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief  Is the given string a valid 6-letter PDB record type?
	static bool is_valid_record_type( std::string const & type );

	/// @brief  Get the corresponding PDB record from the corresponding record type.
	static Record record_from_record_type( RecordType const & type );

	/// @brief  Get the corresponding PDB record from the corresponding record type string.
	static Record record_from_record_type( std::string const & type );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	RecordCollection();

private:  // Private data /////////////////////////////////////////////////////
	std::map< std::string, RecordType > string_to_record_type_map_;
	RecordDef record_definitions_;
};

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_RecordCollection_HH
