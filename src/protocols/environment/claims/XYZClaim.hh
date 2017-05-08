// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_environment_claims_XYZClaim_hh
#define INCLUDED_protocols_environment_claims_XYZClaim_hh


// Unit Headers
#include <protocols/environment/claims/XYZClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
#include <string>
#include <sstream>
#include <utility>

namespace protocols {
namespace environment {
namespace claims {

class XYZClaim : public EnvClaim {
	typedef core::environment::FoldTreeSketch FoldTreeSketch;
	typedef EnvClaim Parent;
	typedef core::environment::LocalPosition LocalPosition;
	typedef core::environment::LocalPositions LocalPositions;
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

public:

	XYZClaim( ClientMoverOP owner,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap const& );

	// Initializer for a single backbone angle
	XYZClaim( ClientMoverOP owner,
		LocalPosition const& local_pos );

	// Initializer for a contiguous range of residues.
	// TODO: replace this function in RigidChunkCM with the constructor using ResidueSelectorsand remove
	XYZClaim( ClientMoverOP owner,
		std::string const& label,
		std::pair< core::Size, core::Size > const& range );

	XYZClaim( ClientMoverOP owner,
		ResidueSelectorCOP selector );

	virtual void yield_elements( core::pose::Pose const&,
		DOFElements& elements ) const;

	ControlStrength const& ctrl_strength() const;

	void strength( ControlStrength const& control, ControlStrength const& initialization );

	ControlStrength const& init_strength() const;

	ResidueSelectorCOP selector() const { return selector_; }

	void set_relative( bool s ) { bRelative_ = s; }

	bool relative() const { return bRelative_; }


	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string class_name();

	virtual EnvClaimOP clone() const;

	virtual std::string type() const;

	virtual void show( std::ostream& os ) const;

protected:
	virtual
	DOFElement wrap_dof_id( core::id::DOF_ID const& id ) const;

	void build_bond_length_elements( core::Size seqpos,
		ProtectedConformationCOP const&,
		DOFElements& elements ) const;

	void build_bond_angle_elements( core::Size seqpos,
		ProtectedConformationCOP const&,
		DOFElements& elements ) const;

	void build_bond_torsion_elements( core::Size seqpos,
		ProtectedConformationCOP const&,
		DOFElements& elements ) const;

private:
	ResidueSelectorCOP selector_;
	ControlStrength c_str_;
	ControlStrength i_str_;
	bool bRelative_;

}; //class XYZClaim


}
}
}

#endif
