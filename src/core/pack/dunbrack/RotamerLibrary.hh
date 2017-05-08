// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dunbrack/RotamerLibrary.hh
/// @brief  Rotamer Library classes and utility functions
/// @author Phil Bradley


#ifndef INCLUDED_core_pack_dunbrack_RotamerLibrary_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibrary_hh

// Unit headers
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
//#include <utility/graph/Graph.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
//#include <utility/graph/Graph.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#endif

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace dunbrack {

void
rotamer_from_chi(
	conformation::Residue const & rsd,
	RotVector & rot
);

void
rotamer_from_chi(
	chemical::ResidueType const & rsd_type,
	ChiVector const & chi,
	RotVector & rot
);

/// @brief Do not bother calling this function if you're interested in staying
/// up to date with changes in the rotamer library.  Call rotamer_from_chi instead.
/// It will eventually call this function iff dun10 is not set to true.
void
rotamer_from_chi_02(
	ChiVector const & chi,
	//ChiVector chi, //pass in by value and locally modified
	chemical::AA const res, // pass an AA, not an AA &
	RotVector & rot
);

void
rotamer_from_chi_02(
	Real4 const & chi,
	chemical::AA const res,
	Size nchi,
	Size4 & rot
);

/// @brief Find the difference in angles between two chi values using hard-coded
/// symmetry information for the symmetric amino acids.  Disappears for 2010 library.
Real
subtract_chi_angles(
	Real chi1,
	Real chi2,
	chemical::AA const & aa,
	int chino
);

///////////////////////////////////////////////////////////////////////////////
/// @brief A class to manage the Dunbrack Rotamer Libraries
/// @details This class no longer manages arbitrary rotamer libraries.
/// Use core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( restype );
/// to get rotamer libraries for an arbitrary ResidueType.

class RotamerLibrary : public utility::SingletonBase< RotamerLibrary >
{
public:
	friend class utility::SingletonBase< RotamerLibrary >;

	typedef chemical::AA AA;
	typedef conformation::Residue Residue;
	typedef chemical::ResidueType ResidueType;

public:
	virtual ~RotamerLibrary();

	Real
	rotamer_energy(
		Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	Real
	best_rotamer_energy(
		Residue const & rsd,
		bool curr_rotamer_only,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// to do:
	Real
	rotamer_energy_deriv(
		Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	// // get_rsd_library() is no longer valid! Use core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( restype ); instead.
	// SingleResidueRotamerLibraryCOP
	// get_rsd_library( chemical::ResidueType const & rsd_type ) const;

	rotamers::SingleResidueRotamerLibraryCOP
	get_library_by_aa( chemical::AA const & aa ) const;

	/// @brief Reload the Dunbrack Rotamer libraries from ASCII, and make sure that they match the ones loaded from binary.
	/// Return true if the binary file is valid, false if the binary is invalid.
	/// NOTE WELL: This is *not* a const function, as reloading from ASCII modifies internals.
	/// It's also *VERY* thread unsafe. Never call this function from a multithreaded context.
	bool
	validate_dunbrack_binary();

public: // public such that people can see if we're reading from the binary or ASCII, for debugging purposes

	bool decide_read_from_binary() const;

private:

	bool decide_read_from_binary_02() const;
	bool decide_read_from_binary_10() const;

	bool binary_is_up_to_date( utility::io::izstream & binlib ) const;

	bool binary_is_up_to_date_02( utility::io::izstream & binlib ) const;
	bool binary_is_up_to_date_10( utility::io::izstream & binlib ) const;

	bool decide_write_binary() const;

	bool decide_write_binary_02() const;
	bool decide_write_binary_10() const;

	std::string get_library_name_02() const;
	std::string get_library_name_10() const;

public: // public such that people can see which binary we're reading from.

	std::string get_binary_name() const;

private:

	std::string get_binary_name_02(bool for_writing = false) const;
	std::string get_binary_name_10(bool for_writing = false) const;

	Size current_binary_format_version_id_02() const;
	Size current_binary_format_version_id_10() const;

	/// @brief Main interface for reading in dunbrack libraries.
	/// Version option checks are handled inside.
	void create_fa_dunbrack_libraries();

	void create_fa_dunbrack_libraries_from_ASCII();

	void create_fa_dunbrack_libraries_from_binary();

	void create_fa_dunbrack_libraries_02_from_ASCII();
	void create_fa_dunbrack_libraries_02_from_binary();

	void create_fa_dunbrack_libraries_10_from_ASCII();
	void create_fa_dunbrack_libraries_10_from_binary();

	/// @brief Add a fullatom canonical AA Dunbrack library
	/// @details This, along with the create_* functions, are
	/// private and should only be called during singleton construction.
	/// If you need to make it public, you need to add a mutex
	void
	add_residue_library(
		AA const & aa,
		rotamers::SingleResidueRotamerLibraryCOP rot_lib
	);

	void write_binary_fa_dunbrack_libraries() const;

	void write_binary_fa_dunbrack_libraries_02() const;
	void write_binary_fa_dunbrack_libraries_10() const;

	std::string random_tempname( std::string const & same_dir_as, std::string const & prefix ) const;

	/// @brief Instantiate the appropriate RSRDL< T > library given the n_chi input and initialize
	/// the library from the input stream.  Return the string that the library returns after reading
	/// the from the stream (which may be the empty string, or the three letter code for the next AA
	/// in the input stream).  If the first three letter code for this amino acid has already been
	/// retrieved from the input stream, indicate so with the "first_three_letter_code_alread_read"
	/// boolean.
	SingleResidueDunbrackLibraryOP
	create_rotameric_dunlib(
		chemical::AA aa,
		Size const n_chi,
		Size const n_bb,
		utility::io::izstream & library,
		bool dun02,
		std::string & next_aa_in_library,
		bool first_three_letter_code_already_read
	) const;

	SingleResidueDunbrackLibraryOP
	create_rotameric_dunlib(
		chemical::AA aa,
		Size const n_chi,
		Size const n_bb,
		bool dun02
	) const;

	SingleResidueDunbrackLibraryOP
	create_rotameric_dunlib(
		chemical::AA aa,
		Size const n_chi,
		bool dun02
	) const;

	/// @brief Instantiate the appropriate SRSRDL< T > library given the nchi input.  nchi indicates
	/// the number of rotameric chi -- there are nchi+1 chi described by the library, the last of which
	/// is non-rotameric.  (E.g. Tyrosine has an nchi of 1, its chi2 is non-rotameric and its chi3 is not
	/// described by the library.)  The input streams are either read or not read depending on the
	/// boolean variables.
	SingleResidueDunbrackLibraryOP
	create_semi_rotameric_dunlib(
		chemical::AA aa,
		Size const nchi,
		Size const nbb,
		bool const use_bbind_rnchi_scoring,
		bool const use_bbind_rnchi_sampling,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle,
		utility::io::izstream & rotamer_definitions,
		utility::io::izstream & regular_library,
		utility::io::izstream & continuous_minimization_bbdep
		// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
	) const;

	SingleResidueDunbrackLibraryOP
	create_semi_rotameric_dunlib(
		chemical::AA aa,
		Size const nchi,
		bool const use_bbind_rnchi_scoring,
		bool const use_bbind_rnchi_sampling,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle
	) const;


	SingleResidueDunbrackLibraryOP
	create_semi_rotameric_dunlib(
		chemical::AA aa,
		Size const nchi,
		Size const nbb,
		bool const use_bbind_rnchi_scoring,
		bool const use_bbind_rnchi_sampling,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle
	) const;

	SingleResidueDunbrackLibraryOP
	create_srdl(
		chemical::AA aa_in
	) const;

public:
	static
	void
	initialize_dun10_aa_parameters(
		utility::vector1< chemical::AA > & rotameric_amino_acids,
		utility::vector1< Size > & rotameric_n_chi,
		utility::vector1< Size > & rotameric_n_bb,
		utility::vector1< chemical::AA > & sraa,
		utility::vector1< Size > & srnchi,
		utility::vector1< Size > & srnbb,
		utility::vector1< bool > & scind,
		utility::vector1< bool > & sampind,
		utility::vector1< bool > & sym,
		utility::vector1< Real > & astr
	);

	static
	void
	initialize_dun02_aa_parameters(
		utility::vector1< chemical::AA > & rotameric_amino_acids,
		utility::vector1< Size > & rotameric_n_chi,
		utility::vector1< Size > & rotameric_n_bb
	);

	static
	void
	initialize_dun02_aa_parameters(
		utility::vector1< Size > & nchi_for_aa,
		utility::vector1< Size > & nbb_for_aa,
		Size & n_rotameric_aas
	);

private:
	void
	write_to_binary( utility::io::ozstream & out ) const;

	/// @brief Initialize the RotamerLibrary from the stored binary
	/// @details Invoked through create_singleron_instance()
	/// Not threadsafe in itself. Do not call directly.
	void
	read_from_binary( utility::io::izstream & in );


	RotamerLibrary();
	RotamerLibrary( RotamerLibrary const & ); // unimplemented
	RotamerLibrary const & operator = ( RotamerLibrary const & ); // unimplemented

private:

	// NOTE: If you make the data in this singleton mutable,
	// you will need to add a ReadWriteMutex to make it threadsafe.

	utility::vector1< rotamers::SingleResidueRotamerLibraryCOP > aa_libraries_;
	typedef utility::vector1< rotamers::SingleResidueRotamerLibraryCOP > library_iterator;

};

} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_pack_dunbrack_RotamerLibrary_HH
