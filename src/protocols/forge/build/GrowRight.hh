// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/GrowRight.hh
/// @brief instruction to create a c-side extension
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_GrowRight_hh
#define INCLUDED_protocols_forge_build_GrowRight_hh

// unit headers
#include <protocols/forge/build/GrowRight.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief instruction to create a c-side extension
/// @remarks Use this for c-side insertions but typically not c-terminal
///  extensions unless necessary.  It does not automatically cover the
///  additional residue on the left endpoint that needs to move during
///  c-terminal extensions due to invalid psi torsion. For that case,
///  use the SegmentRebuild class replacing the c-terminal residue with
///  desired length+1.
class GrowRight : public BuildInstruction {


private: // typedefs


	typedef BuildInstruction Super;


public: // typedefs


	typedef Super::Size Size;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::MoveMap MoveMap;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;


public: // construct/destruct


	/// @brief default constructor
	GrowRight();


	/// @brief constructor
	/// @param[in] pos grow a c-side extension after this position
	/// @param[in] ss the secondary structure desired, also defines length of extension
	/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @remarks length of the *one-letter* aa must equal the length of ss
	GrowRight(
		Size const pos,
		String const & ss,
		String const & aa = String(),
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);


	/// @brief copy constructor
	GrowRight( GrowRight const & rval );


	/// @brief default destructor
	virtual
	~GrowRight();


public: // assignment


	/// @brief copy assignment
	GrowRight & operator =( GrowRight const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // accessors


	/// @brief grow right from this anchor position
	/// @remarks This can change if listening to Conformation LengthEvents.
	///  Use original_interval() to get the original anchor position.
	inline
	Size pos() const {
		return pos_;
	}


	/// @brief get secondary structure string
	inline
	String const & ss() const {
		return ss_;
	}


	/// @brief get annotated amino acid string
	inline
	String const & aa() const {
		return aa_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return false, stores invalid interval
	inline
	virtual
	bool original_interval_valid() const {
		return false;
	}


	/// @brief a copy of the working range of residues specifying the modified region
	/// @remarks This can change if listening to Conformation LengthEvents
	inline
	virtual
	Interval interval() const {
		return Interval( pos_ + 1, pos_ + ss_.length() );
	}


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return An empty set -- no positions are pre-existing.
	virtual
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return A set of positions spanning the entire region -- all positions
	///  are new.
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return An empty set -- no positions are defined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @return A set of positions spanning the entire region -- all positions
	///  are undefined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions undefined_positions() const;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
	///  at the MoveMapTorsionID level
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	MoveMap movemap() const ;


public: // virtual Conformation observer interface


	/// @brief update indexing on residue append
	virtual
	void on_residue_append( LengthEvent const & event );


	/// @brief update indexing on residue prepend
	virtual
	void on_residue_prepend( LengthEvent const & event );


	/// @brief update indexing on residue delete
	virtual
	void on_residue_delete( LengthEvent const & event );


public: // original positions


	/// @brief return the set of positions within the original interval that
	///  will be kept in this BuildInstruction
	/// @return An empty set -- no positions are kept.
	virtual
	Positions original_kept_positions() const;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	/// @return An empty set -- no positions are deleted.
	virtual
	Positions original_deleted_positions() const;


public: // instruction comparison


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no fixed positions
	virtual
	Positions original_fixed_positions() const;


	/// @brief return set of any mutable positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no mutable positions
	virtual
	Positions original_mutable_positions() const;


public: // virtual object descriptor


	/// @brief does this object create undefined backbone in the modified region?
	/// @return true
	inline
	virtual
	bool creates_undefined_backbone() const {
		return true;
	}


protected: // virtual Pose modification methods


	/// @brief are dependencies satisfied so that modify_impl() can complete
	///  successfully?
	/// @return always True, this BuildInstruction has no dependencies
	inline
	virtual
	bool dependencies_satisfied() const {
		return true;
	}


	/// @brief do the actual work of modifying the Pose
	virtual
	void modify_impl( Pose & pose );


protected: // virtual mutators


	/// @brief do the actual reset of intervals, positions, etc to initial state
	virtual
	void reset_accounting_impl();


private: // data


	/// @brief make a c-side extension after this position
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Size pos_;


	/// @brief secondary structure string, also defines length of extension
	String ss_;


	/// @brief annotated amino acid string, length of the one-letter version
	///  must be equal to length of ss
	String aa_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_GrowRight_HH */
