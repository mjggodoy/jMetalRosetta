// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_JumpSRFD_HH
#define INCLUDED_core_fragment_JumpSRFD_HH

// Unit Headers
#include <core/fragment/JumpSRFD.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>


#include <core/kinematics/RT.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

class UpJumpSRFD : public SingleResidueFragData {
	typedef SingleResidueFragData Parent;
public:

	UpJumpSRFD( char sequence = 'X'  )
	: SingleResidueFragData( sequence)
	{};

	/// @brief clone
	SingleResidueFragDataOP clone() const override {
		return SingleResidueFragDataOP( new UpJumpSRFD( *this ) );
	};

	/// @brief create a new instance of this object

	SingleResidueFragDataOP create() const override {
		return SingleResidueFragDataOP( new UpJumpSRFD() );
	}

	bool apply( pose::Pose&, Size const, Frame const& ) const override {
		return true;
	}

	bool apply( kinematics::MoveMap const &, pose::Pose &, Size const, Frame const & ) const override {
		return true;
	}

	bool apply_ss( std::string&, Size, Frame const& ) const override {
		return true;
	};

	bool steal( pose::Pose const&, Size, Frame const& ) override {
		return true;
	}

	bool is_compatible( SingleResidueFragData const& ) const override;

	bool is_applicable( kinematics::MoveMap const&, Size ) const override {
		return true;
	}


	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "UpJump";
	}


	// for UpJumpSRFD these functions should never be called
	// might obsolet these completely and only use SRFDs via the func (pose, ipos, Frame ) variant
	bool apply(pose::Pose &, Size) const override {
		debug_assert( 0 );
		return false;
	}

	/// @brief for UpJumpSRFD this function should never be called, instead use Frame version
	/// @return always false
	/// @warning will trigger a false runtime assert
	bool apply( kinematics::MoveMap const &, pose::Pose &, Size const ) const override {
		runtime_assert( 0 );
		return false;
	}

	// insert fragment_data sec-struct into ss-string at position seq_pos
	bool apply_ss( std::string&, Size) const override {
		debug_assert( 0 );
		return false;
	}

	// insert fragment_data into pose at position seq_pos
	bool steal(pose::Pose const&, Size) override  {
		debug_assert( 0 );
		return false;
	}

};


class DownJumpSRFD : public SingleResidueFragData {
	typedef SingleResidueFragData Parent;
public:
	typedef utility::vector1<std::string> AtomList;

	DownJumpSRFD( char sequence = 'X'  )
	: SingleResidueFragData( sequence)
	{
		set_standard_stub_atoms();
	};

	DownJumpSRFD( kinematics::RT rt, AtomList downstream_stub_atoms, AtomList upstream_stub_atoms, char sequence ) :
		SingleResidueFragData( sequence ),
		rt_( rt ),
		downstream_stub_atoms_( downstream_stub_atoms ),
		upstream_stub_atoms_( upstream_stub_atoms )
	{
		debug_assert( downstream_stub_atoms.size() == 3 || downstream_stub_atoms.size() == 4 );
		debug_assert( upstream_stub_atoms.size() == 3 || upstream_stub_atoms.size() == 4 );
	};

	/// @brief clone
	SingleResidueFragDataOP clone() const override {
		return SingleResidueFragDataOP( new DownJumpSRFD( *this ) );
	};

	/// @brief create a new instance of this object

	SingleResidueFragDataOP create() const override {
		return SingleResidueFragDataOP( new DownJumpSRFD() );
	}

	/// @brief set value of jump
	void
	set_jump( kinematics::RT setting ) {
		rt_ = setting;
	};

	void
	set_stub_atoms( AtomList downstream, AtomList upstream );

	//@brief set stub atoms to "CA" as center and "N", "CA", "C" as stub atoms
	void
	set_standard_stub_atoms();

	/// @brief the seq_pos is hijacked for the rt_nr

	bool apply( pose::Pose&, Size const intra_frame_pos, Frame const& ) const override;
	bool apply( kinematics::MoveMap const &, pose::Pose &, Size const intra_frame_pos, Frame const & ) const override;
	bool apply_ss( std::string&, Size, Frame const& ) const override { return true; }; //does nothing as jump has no ss

	bool steal( pose::Pose const&, Size pos, Frame const& ) override;
	bool is_compatible( SingleResidueFragData const& ) const override;
	bool is_applicable( kinematics::MoveMap const&, Size pos, Frame const& ) const override;


	// for DownJumpSRFD these functions should never be called
	// might obsolet these completely and only use SRFDs via the func (pose, ipos, Frame ) variant
	bool apply(pose::Pose &, Size) const override {
		debug_assert( 0 );
		return false;
	}

	/// @brief for DownJumpSRFD this function should never be called, instead use Frame version
	/// @return always false
	/// @warning will trigger a false runtime assert
	bool apply( kinematics::MoveMap const &, pose::Pose &, Size const ) const override;

	// insert fragment_data sec-struct into ss-string at position seq_pos
	bool apply_ss( std::string&, Size) const override {
		debug_assert( 0 );
		return false;
	}

	// insert fragment_data into pose at position seq_pos
	bool steal(pose::Pose const&, Size) override  {
		debug_assert( 0 );
		return false;
	}


	bool
	is_applicable( kinematics::MoveMap const&, Size ) const override {
		debug_assert( 0 );
		return false;
	}

	//  void get_stubs(
	//     pose::Pose const& pose,
	//   Size downstream_res_nr,
	//   id::StubID &up_stub,
	//   id::StubID &down_stub
	//  ) const;


	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "DownJump";
	}


	void show( std::ostream &out ) const override;

	virtual
	void read( std::istream &in );

private:
	kinematics::RT rt_;
	AtomList downstream_stub_atoms_;
	AtomList upstream_stub_atoms_;
};

}
}

#endif
