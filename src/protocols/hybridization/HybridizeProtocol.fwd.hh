// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file HybridizeProtocol.fwd.hh
/// @brief
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_HybridizeProtocol_fwd_hh
#define INCLUDED_protocols_hybridization_HybridizeProtocol_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace hybridization {

class HybridizeProtocol;
typedef utility::pointer::shared_ptr< HybridizeProtocol > HybridizeProtocolOP;
typedef utility::pointer::shared_ptr< HybridizeProtocol const > HybridizeProtocolCOP;


} // hybridization
} // protocols

#endif
