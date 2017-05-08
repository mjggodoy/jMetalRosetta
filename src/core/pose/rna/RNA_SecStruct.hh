// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/RNA_SecStruct.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_RNA_SecStruct_HH
#define INCLUDED_core_pose_rna_RNA_SecStruct_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/rna/RNA_SecStruct.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace core {
namespace pose {
namespace rna {

class RNA_SecStruct: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_SecStruct( std::string const & secstruct,
		std::string const & secstruct_file = "",
		std::string const & sequence = "" );

	RNA_SecStruct(){} // everything empty.

	//destructor
	~RNA_SecStruct();

public:

	utility::vector1< std::pair< core::Size, core::Size > > base_pairs() const { return base_pairs_; }

	utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > stems() const { return stems_; }

	void
	check_compatible_with_sequence( std::string const & sequence_in,
		bool const check_complementarity = true ) const;

	core::Size size() const { return secstruct_.size(); }

	std::string secstruct() const { return secstruct_; }

	void remove_pair( std::pair<core::Size, core::Size > pair );

	bool blank() const;

	// @brief in any base pair?
	bool in_helix( core::Size const & i ) const;

	// @brief in same helix -- not necessarily paired though. great for checking helix boundaries.
	bool in_same_helix( core::Size const & i, core::Size const & j ) const;

private:

	void
	set_secstruct( std::string const & secstruct );

	void
	set_basepairs_from_secstruct();

	bool
	check_balanced_secstruct() const;

	void
	blank_secstruct( std::string const & sequence_in );

	void
	read_secstruct_from_file( std::string const & filename );

	void
	figure_out_stems();

private:

	std::string secstruct_;
	utility::vector1< std::pair< core::Size, core::Size > > base_pairs_;
	utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > stems_;

	utility::vector1< Size > spacer_positions_;

	std::map< char, utility::vector1< char > > rna_complement_;

};

} //secstruct
} //farna
} //protocols

#endif
