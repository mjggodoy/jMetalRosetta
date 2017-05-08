// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/DofUnlock.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/DofUnlock.hh>

// Package headers
#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/pose/Pose.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.DofUnlock", basic::t_info );

namespace protocols {
namespace environment {

using core::environment::DofPassport;
using core::environment::DofPassportCOP;

DofUnlock::DofUnlock( core::conformation::Conformation& conf,
	DofPassportCOP pass ) :
	conformation_( conf ),
	pass_( pass )
{
	if ( pass ) {
		conformation_.push_passport( pass );
	}
}

DofUnlock::~DofUnlock(){
	if ( pass_ && conformation_.is_protected() ) {
		// should always be valid since pass_ exists.
		DofPassportCOP pass_out( conformation_.pop_passport() );
		if ( pass_out->mover() != pass_->mover() ||
				pass_out->env_id() != pass_out->env_id() ) {
			std::ostringstream ss;
			ss << "DofUnlock popped a passport belonging to mover " << pass_out->mover()
				<< " and environment id " << pass_out->env_id() << " when it expected a passport from "
				<< pass_->mover() << " and environment id " << pass_->env_id() << "."
				<< "Something has gone horribly wrong, probably a result of strange pose- or conformation-copying behavior." << std::endl;
			tr.Error << "[ERROR]" << ss.str() << std::endl;
			// This has to be a utility_exit_with_message(), as throwing directly results in compiler errors
			// *Destructors shouldn't throw*
			utility_exit_with_message( ss.str() );
		}
	}
}

} // environment
} // protocols
