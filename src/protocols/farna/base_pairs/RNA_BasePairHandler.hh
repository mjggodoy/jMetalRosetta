// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/base_pairs/RNA_BasePairHandler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_BasePairHandler_HH
#define INCLUDED_protocols_farna_RNA_BasePairHandler_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/farna/base_pairs/RNA_BasePairHandler.fwd.hh>
#include <protocols/farna/base_pairs/BasePairStep.fwd.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/BasePair.hh>
#include <map>

namespace protocols {
namespace farna {
namespace base_pairs {

class RNA_BasePairHandler: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_BasePairHandler( core::pose::Pose const & pose );

	RNA_BasePairHandler( protocols::farna::setup::RNA_DeNovoParameters const & rna_params );

	//destructor
	~RNA_BasePairHandler();

public:

	bool
	check_base_pairs( core::pose::Pose & pose, protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const;

	std::map< Size, Size >
	connections() const;

	void
	setup_base_pair_constraints( core::pose::Pose & pose,
		protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
		core::Real const suppress_bp_constraint = 1.0 ) const;

	core::pose::rna::RNA_BasePairList
	rna_pairing_list() const { return rna_pairing_list_; }

	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > >  const &
	chain_connections() const { return chain_connections_; }

	utility::vector1< Size >
	get_stem_residues(  core::pose::Pose const & pose ) const;

	utility::vector1< BasePairStep >
	get_canonical_base_pair_steps() const;

	utility::vector1< BasePairStep >
	get_noncanonical_base_pair_steps() const;

	utility::vector1< BasePairStep >
	get_base_pair_steps( bool const just_canonical ) const;

private:

	void
	figure_out_partner( std::map< Size, Size > & partner, bool const force_canonical ) const;

private:

	core::pose::rna::RNA_BasePairList rna_pairing_list_;
	utility::vector1< Size > cutpoints_open_;
	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;

};

} //base_pairs
} //farna
} //protocols

#endif
