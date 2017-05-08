// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rotamers/SingleLigandRotamerLibrary.hh
/// @brief  SingleLigandRotamerLibrary class
/// @author Ian W. Davis

#ifndef INCLUDED_core_pack_rotamers_SingleLigandRotamerLibrary_hh
#define INCLUDED_core_pack_rotamers_SingleLigandRotamerLibrary_hh

// Unit Headers
#include <core/pack/rotamers/SingleLigandRotamerLibrary.fwd.hh>

// Package Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Project Headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>

#include <map>

namespace core {
namespace pack {
namespace rotamers {

typedef std::map< std::string, core::Vector > NamePosMap;

/// @brief A fixed library of conformations for some residue type (doesn't have to be a ligand).
/// @details Reads residue conformations in PDB format separated by mandatory TER records.
/// "Included" from a residue .params file with the PDB_ROTAMERS keyword.
class SingleLigandRotamerLibrary : public SingleResidueRotamerLibrary
{
	//typedef utility::vector1< core::Size > Fragment; // a list of atom indices (variable len)
	//typedef utility::vector1< Fragment > Fragments;
	//typedef utility::vector1< core::Size > Automorphism; // a remapping of atom indices (fixed len)
	//typedef utility::vector1< Automorphism > Automorphisms;

public:

	SingleLigandRotamerLibrary();

	virtual ~SingleLigandRotamerLibrary();

	/// @brief Reads conformers from PDB-format file (must be separated by TER records!)
	virtual
	void
	init_from_file(
		std::string const & filename,
		chemical::ResidueType const & restype
	);

	/// @brief Reads conformers from a vector of name:coordinate maps
	void
	init_from_vector(
		utility::vector1< NamePosMap > const & coordinates
	);

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const;

	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		bool curr_rotamer_only,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const;

	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const &,// rsd,
		pose::Pose const & /*pose*/,
		dunbrack::RotamerLibraryScratchSpace &,// scratch,
		numeric::random::RandomGenerator &,// RG,
		dunbrack::ChiVector &,// new_chi_angles,
		bool //perturb_from_rotamer_center
	) const {} // stubbed out for the moment.

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		utility::graph::GraphCOP,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	void
	write_to_file( utility::io::ozstream &out ) const;

	Real
	get_reference_energy() const{
		return ref_energy_;
	}

	void set_reference_energy( Real ref_E_in){
		ref_energy_ = ref_E_in;
	}

	/// @brief Build a set of rotamers for the given ResidueType
	void
	build_base_rotamers( chemical::ResidueType const & restype, RotamerVector & base_rotamers ) const;

private:

	/// @brief Fills in missing hydrogens/virtual atoms from library load
	///
	/// @details "missing" is the vector for atoms missing for the current "rsd"
	//           "missed" is the vector annotating which atoms have already been filled (for diagnostic output tracking)
	void
	fill_missing_atoms( utility::vector1< bool > missing, conformation::ResidueOP rsd, utility::vector1< bool > & missed ) const;

	// Breaking the ligand into rigid fragments that would supply (putative) pharamacophores
	// to superimpose on was a nice idea, but it breaks the packer assumption that nbr_atom doesn't move.

	//void find_fragments(chemical::ResidueTypeCOP restype);
	//void list_automorphisms(chemical::ResidueTypeCOP restype);
	//void unique_auto_for_frags();

	//void superimpose(
	// conformation::Residue const & existing,
	// conformation::Residue & conformer,
	// Fragment const & frag,
	// Automorphism const & morph
	//) const;

private:

	utility::vector1< NamePosMap > atom_positions_;

	// A baseline reference energy applied to *all* conformers in this library,
	// like amino acid reference energies.
	Real ref_energy_;

	// Breaking the ligand into rigid fragments that would supply (putative) pharamacophores
	// to superimpose on was a nice idea, but it breaks the packer assumption that nbr_atom doesn't move.

	/// A ligand with N torsions decomposes into N+1 rigid fragments,
	/// or fewer if some consist of 2 atoms or less.
	/// Each fragement can be used for superimposing on, possibly in multiple ways.
	//Fragments rigid_frags_;
	/// All the automorphisms for the ligand as a whole, shared among frags (below).
	//Automorphisms automorphs_;
	/// A subset of automorphs_ that reflects the *unique* automorphisms within
	/// a particular rigid fragment.  These represent different ways of
	/// superimposing two structures using the atoms from that fragment.
	//utility::vector1< utility::vector1 < Automorphism * > > frag_automorphs_;
	/// Total number of automorphisms in frag_automorphs_ -- number of possible ways to superimpose
	//Size total_superpos_;

}; // SingleLigandRotamerLibrary


} // namespace rotamers
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_rotamers_SingleLigandRotamerLibrary_hh
