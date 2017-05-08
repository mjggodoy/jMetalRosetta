// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/MakeAsymmetricStructureDataMover.fwd.hh
/// @brief Converts a StructureData for a symmetric pose into an asymmetric representation
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_denovo_design_movers_MakeAsymmetricStructureDataMover_fwd_hh
#define INCLUDED_protocols_denovo_design_movers_MakeAsymmetricStructureDataMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace denovo_design {
namespace movers {

class MakeAsymmetricStructureDataMover;

typedef utility::pointer::shared_ptr< MakeAsymmetricStructureDataMover > MakeAsymmetricStructureDataMoverOP;
typedef utility::pointer::shared_ptr< MakeAsymmetricStructureDataMover const > MakeAsymmetricStructureDataMoverCOP;



} //protocols
} //denovo_design
} //movers


#endif //INCLUDED_protocols_denovo_design_movers_MakeAsymmetricStructureDataMover_fwd_hh





