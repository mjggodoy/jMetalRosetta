// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SWA_ProtocolUtil.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_MainChainTorsionSet_HH
#define INCLUDED_protocols_stepwise_protein_MainChainTorsionSet_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

//////////////////////////////////
class MainChainTorsionSet: public utility::pointer::ReferenceCount {


public:

	MainChainTorsionSet( core::Real const & phi, core::Real const & psi, core::Real const & omega );

	MainChainTorsionSet( core::Real const & phi, core::Real const & psi );

	virtual ~MainChainTorsionSet();

	MainChainTorsionSet &
	operator=( MainChainTorsionSet const & src );

	core::Real phi() const;
	core::Real psi() const;
	core::Real omega() const;

private:
	core::Real phi_;
	core::Real psi_;
	core::Real omega_;
};

typedef utility::vector1< MainChainTorsionSet > MainChainTorsionSetList ;


} //protein
} //modeler
} //stepwise
} //protocols

#endif
