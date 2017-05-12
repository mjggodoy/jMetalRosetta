// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyDesignProtocol.fwd.hh
/// @brief Handles all loading of CDR, Framework, and cluster/dmap info from database file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)



#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_fwd_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace antibody {
namespace design {

// Forward
class AntibodyDesignProtocol;

typedef utility::pointer::shared_ptr< AntibodyDesignProtocol > AntibodyDesignProtocolOP;
typedef utility::pointer::shared_ptr< AntibodyDesignProtocol const > AntibodyDesignProtocolCOP;


}
} //namespace antibody_design
} //namespace protocols

#endif //INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_fwd_hh
