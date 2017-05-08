// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SymmetryClaimer.hh
/// @brief  Basic Claimer for making Symmetry Claims
/// @author Justin Porter, Tatjana Braun


#ifndef INCLUDED_protocols_topology_broker_SymmetryClaimer_hh
#define INCLUDED_protocols_topology_broker_SymmetryClaimer_hh


// Unit Headers
#include <protocols/topology_broker/SymmetryClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>


// Project Headers
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers


//// C++ headers

#include <string>
#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {


class SymmetryClaimer : public virtual TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	SymmetryClaimer(); //for factory
	~SymmetryClaimer() override = default;

	SymmetryClaimer( SymmetryClaimer const & src );

	SymmetryClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<SymmetryClaimer>( TopologyClaimer::shared_from_this() ); }


	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new SymmetryClaimer( *this ) );
	}

	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "SymmetryClaimer";
	}

	void generate_symmetry_claims( claims::SymmetryClaims& ) override;

	virtual void symmetry_duplicate( claims::DofClaims&, core::pose::Pose& );


protected:
	bool read_tag( std::string tag, std::istream& ) override;

	//    virtual void symmdup_sequence_claims( claims::DofClaims& );
	//    virtual void symmdup_pose_sequence( core::pose::Pose& );
	//    virtual void build_symm_virtual_residues( claims::DofClaims&, core::pose::Pose& );

private:
	core::conformation::symmetry::SymmDataOP symm_data_;
	Size asymmetric_res_;


};

}
}

#endif
