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
/// @author Justin R. Porter


#ifndef INCLUDED_protocols_environment_claims_JumpClaim_hh
#define INCLUDED_protocols_environment_claims_JumpClaim_hh


// Unit Headers
#include <protocols/environment/claims/JumpClaim.fwd.hh>
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <utility/tag/Tag.fwd.hh>

//// C++ headers

// option key includes


namespace protocols {
namespace environment {
namespace claims {

class JumpClaim : public EnvClaim {
	typedef EnvClaim Parent;
	typedef core::environment::LocalPosition LocalPosition;
	typedef core::select::residue_selector::ResidueSelector ResidueSelector;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

public:

	JumpClaim( ClientMoverOP owner,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap const& );

	JumpClaim( ClientMoverOP owner,
		std::string const& label,
		LocalPosition const& jpos1,
		LocalPosition const& jpos2,
		LocalPosition const& cutp = core::environment::NO_POSITION );

	virtual void yield_elements( core::environment::FoldTreeSketch const& fts, ResidueElements& elements ) const;

	virtual void yield_elements( core::environment::FoldTreeSketch const& fts, JumpElements& elements ) const;

	virtual void yield_elements( core::environment::FoldTreeSketch const& fts, CutElements& elements ) const;

	virtual void yield_elements( core::pose::Pose const&, DOFElements& elements ) const;

	/// @brief set the two atom names to use as the start and beginning of the jump. Both must be set simultaneously
	///        because the FoldTree requires this.
	void set_atoms( std::string const& a1, std::string const& a2 );

	void strength( ControlStrength const& cstr, ControlStrength const& istr );

	bool physical() const { return physical_cut_; }

	void physical( bool setting ) { physical_cut_ = setting; }

	void create_vrt_if_necessary( bool setting );

	/// @brief configure this JumpClaim to create a virtual residue to be the jump point if label of the
	///        proposed jump point (i.e. pos1_.label()) does not already exist.
	/// @param set this behavior for pos1()
	/// @param set this behavior for pos2()
	void create_vrt_if_necessary( bool setting_p1,
		bool setting_p2 );

	void cut( LocalPosition const& );

	void stubs_intra_residue( bool setting ){ stubs_intra_residue_ = setting; }

	bool stubs_intra_residue() const { return stubs_intra_residue_; }

	std::string const& label() const;

	LocalPosition const& pos1() const;

	LocalPosition const& pos2() const;

	std::pair< std::string, std::string > const& atoms() const { return atoms_; }


	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string class_name();



	virtual EnvClaimOP clone() const;

	virtual std::string type() const;

	virtual void show( std::ostream& os ) const;

protected:

	//protected and not public because it can return garbage if the value isn't set.
	LocalPosition const& cut() const;

	virtual void validate() const;

private:

	std::string label_;

	LocalPosition pos1_;
	LocalPosition pos2_;

	LocalPosition cut_;

	std::pair< std::string, std::string > atoms_;

	bool physical_cut_;
	bool create_vrt_p1_, create_vrt_p2_;

	bool stubs_intra_residue_;

	ControlStrength c_str_;
	ControlStrength i_str_;

}; //JumpClaim

}
}
}

#endif
