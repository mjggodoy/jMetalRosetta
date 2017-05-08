// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
///           maintains list of ToplogyClaimers
///           maintains SmallFragWeights -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by SmallFragWeights
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_weights_SmoothFragWeight_hh
#define INCLUDED_protocols_topology_broker_weights_SmoothFragWeight_hh

// Unit Headers
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Package Headers
#include <protocols/abinitio/FragmentSampler.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace weights {


class SmoothFragWeight : public AbinitioMoverWeight {
public:
	SmoothFragWeight( core::Real weight = 1.0 ) : weight_( weight ) {};
	virtual core::Real weight( core::Size stageID, core::Real progress /* progress within stage */ ) const {
		if ( stageID < abinitio::STAGE_4 ) return 0.0;
		if ( stageID == abinitio::STAGE_4 && progress >0.5 ) return weight_;
		return 0.0;
	};
private:
	core::Real weight_;
}; //class SmoothFragWeight

// Types
typedef  utility::pointer::shared_ptr< SmoothFragWeight >  SmoothFragWeightOP;
typedef  utility::pointer::shared_ptr< SmoothFragWeight const >  SmoothFragWeightCOP;

typedef  utility::pointer::weak_ptr< SmoothFragWeight >  SmoothFragWeightAP;
typedef  utility::pointer::weak_ptr< SmoothFragWeight const >  SmoothFragWeightCAP;

}
}
}

#endif
