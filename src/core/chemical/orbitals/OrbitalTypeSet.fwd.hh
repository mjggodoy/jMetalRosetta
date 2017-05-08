// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeSet_fwd_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeSet_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace chemical {
namespace orbitals {

class OrbitalTypeSet;

typedef  utility::pointer::weak_ptr< OrbitalTypeSet > OrbitalTypeSetAP;
typedef  utility::pointer::weak_ptr< OrbitalTypeSet const > OrbitalTypeSetCAP;
typedef  utility::pointer::shared_ptr< OrbitalTypeSet > OrbitalTypeSetOP;
typedef  utility::pointer::shared_ptr< OrbitalTypeSet const > OrbitalTypeSetCOP;

}
}
}


#endif /* ORBITALTYPESET_FWD_HH_ */
