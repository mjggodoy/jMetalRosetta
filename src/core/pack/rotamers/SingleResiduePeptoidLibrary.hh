// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResiduePeptoidLibrary.hh
/// @brief  SingleResiduePeptoidLibrary class
/// @brief  Similar to SingleResidueDunbrackLibrary class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_hh
#define INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_hh

// Unit Headers
#include <core/pack/rotamers/SingleResiduePeptoidLibrary.fwd.hh>

// Package Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/chemical/AA.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>

// Numeric Headers
#include <numeric/numeric.functions.hh>

namespace core {
namespace pack {
namespace rotamers {

class SingleResiduePeptoidLibrary : public SingleResidueRotamerLibrary {
public:
	typedef chemical::AA AA;

public:
	/// constants

	/// A good "omega" value to use for N-term residues
	static Real const NEUTRAL_OMG;

	/// A good "phi" value to use for N-term residues
	static Real const NEUTRAL_PHI;

	/// A good "psi" value to use for C-term residues
	static Real const NEUTRAL_PSI;

	// Indicies of backbone dihedrals for peptides and peptoids
	static Size const RSD_PHI_INDEX = 1;
	static Size const RSD_PSI_INDEX = 2;
	static Size const RSD_OMG_INDEX = 3;

public:

	/// c-tor
	SingleResiduePeptoidLibrary(
		Size const n_rotameric_chi
	);

	virtual ~SingleResiduePeptoidLibrary();

public:

	void read_options();

	virtual void write_to_binary( utility::io::ozstream & out ) const;
	virtual void read_from_binary( utility::io::izstream & in );

	/// DOUG DOUG DOUG The definitions of these need to change (and have been changed)
	/// @brief Return all of the rotamer sample data given a particular phi/psi.
	/// For N-terminus residues, hand in the omega value SingleResiduePeptoidLibrary::NEUTRAL_OMG and
	/// phi value SingleResiduePeptoidLibrary::NEUTRAL_PHI and for C-terminus residues, hand in the
	/// psi value SingleResiduePeptoidLibrary::NEUTRAL_PSI. The returned samples should be in
	/// semi-decrasing order by probability; semi, because the rotamers are constructed in sorted order
	/// by their probability in the lower phi-psi bin that the input phi/psi perscribes.
	virtual
	utility::vector1< dunbrack::DunbrackRotamerSampleData >
	get_all_rotamer_samples(
		Real omg,
		Real phi,
		Real psi
	) const = 0;

	/// @brief Return the probability for a particular rotamer where rotamers are
	/// indexed in order of decreasing probability (or something very close to
	/// decreasing probability).
	virtual
	Real
	get_probability_for_rotamer(
		Real omg,
		Real phi,
		Real psi,
		Size rot_ind
	) const = 0;

	virtual
	dunbrack::DunbrackRotamerSampleData
	get_rotamer(
		Real omg,
		Real phi,
		Real psi,
		Size rot_ind
	) const = 0;

	// DOUG DOUG DOUG I suspect this will be needed but will be difficult to get since conf::res doesn't have preceding omega
	// DOUG DOUG DOUG get_omg_from_rsd, get_phi_from_rsd, get_psi_from_rsd functions are used in...
	// RSRPL::assign_random_rotamer
	// RSRPL::find_another_representative_for_unlikely_rotamer
	// RSRPL::best_rotamer_energy
	// RSRPL::interpolate_rotamers(one version doesn't)
	// RSRPL::fill_rotamer_vector(only to get phi and psi but it also gets the pose)
	// DOUG DOUG DOUG maybe just write new ones that take poses, they do not seem to be used outside of this area of code, assign_random_rotamer is though
	// virtual
	// Real
	// get_omg_from_rsd(
	//  conformation::Residue const & rsd
	// ) const = 0;

	// virtual
	// Real
	// get_phi_from_rsd(
	//  conformation::Residue const & rsd
	// ) const = 0;

	// virtual
	// Real
	// get_psi_from_rsd(
	//  conformation::Residue const & rsd
	// ) const = 0;

	virtual
	Real
	get_omg_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const = 0;

	virtual
	Real
	get_phi_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const = 0;

	virtual
	Real
	get_psi_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const = 0;

public:
	/// Virtual functions the derived classes must implement

	/// @brief Derived classes should invoke base class function as well.
	virtual Size memory_usage_in_bytes() const;

	/// @brief The number of chi represented by the library.
	virtual
	Size nchi() const = 0;

	/// @brief The number of rotamer bins represented by the library.
	virtual
	Size n_rotamer_bins() const = 0;

	/// @brief Tell the base class the number of chi bins for each rotameric
	/// chi dimension
	void
	set_n_chi_bins( utility::vector1< Size > const & );

protected:
	/// Read access for the derived class
	bool dun02() const { return dun02_; }

	/// Worker functions available to the derived classes

	virtual
	Size memory_usage_static() const = 0;

	virtual
	Size memory_usage_dynamic() const;

	/// @brief Read access to the n_chi_bins_ vector
	utility::vector1< Size > const &
	n_chi_bins() const {
		return n_chi_bins_;
	}

	/// @brief The base class needs to be informed about which rotamer wells
	/// exist in order to create the rotwell to packed rot conversion data.
	/// set_chi_nbins must be called first.
	void
	mark_rotwell_exists( utility::vector1< Size > const & rotwell );

	/// @brief After the derived class has marked all the rotwells that do exist,
	/// the base class will create the rotwell to packerot conversion data.
	void
	declare_all_existing_rotwells_encountered();

	/// @brief The number of existing rotamers
	Size
	n_packed_rots() const {
		return n_packed_rots_;
	}

	/// @brief The number of possible rotamers -- product of the chi_nbins_ array
	Size
	n_possible_rots() const {
		return n_possible_rots_;
	}

public:

	/// @brief Convert a vector of chi angles (degrees) into a integer vector of rotamer wells.
	/// Derived class should be consistent, but may be arbitrary in how wells divide angle space.
	virtual
	void
	get_rotamer_from_chi(
		dunbrack::ChiVector const & chi,
		dunbrack::RotVector & rot ) const = 0;

public:
	/// Conversion functions

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// (non-compact) "rotamer number"
	Size
	rotwell_2_rotno( utility::vector1< Size > const & rotwell ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// (non-compact) "rotamer number"
	Size
	rotwell_2_rotno( dunbrack::Size4 const & rotwell ) const;

	/// @brief Convert from the rotamer number to the compacted
	/// "packed rotamer number".  Returns 0 if rotno has no corresponding packed rotno.
	Size
	rotno_2_packed_rotno( Size const rotno ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
	Size
	rotwell_2_packed_rotno( utility::vector1< Size > const & rotwell ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
	Size
	rotwell_2_packed_rotno( dunbrack::Size4 const & rotwell ) const;

	/// @brief Convert from the packed rotamer number to the rotamer well
	void
	packed_rotno_2_rotwell( Size const packed_rotno, utility::vector1< Size > & rotwell ) const;

	void
	packed_rotno_2_rotwell(
		Size const packed_rotno,
		dunbrack::Size4 & rotwell
	) const;

	utility::vector1< Size > const &
	packed_rotno_2_rotwell( Size const packed_rotno ) const;

	/// @brief Convert from the rotamer number to the rotamer well
	void
	rotno_2_rotwell( Size const rotno, utility::vector1< Size > & rotwell ) const;

	/// @brief, Turns out, when non-rotameric chi are taken out of the picture,
	/// all remaining chi are binned the same way, except proline. Valid only for
	/// Dun08 library.
	inline
	Size
	bin_rotameric_chi(
		Real chi,
		Size which_chi
	) const {
		debug_assert( ! dun02_ );
		debug_assert( -180.0 <= chi && chi <= 180.0 );

		if ( aa_ == chemical::aa_pro ) {
			if ( which_chi == 1 ) {
				if ( chi > 0 ) { return 1; }
				else { return 2; }
			} else {
				return 1;
			}
		}

		if ( ( chi >= 0.0 ) && ( chi <= 120.0 ) ) { return 1; }
		else if ( std::abs(chi) >= 120.0 ) { return 2; }
		else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { return 3; }
	}

	/// @brief This is not the right place for this code, but the numeric interpolation library
	/// uselessly indexes by 0 and the basic functions aren't inlined...
	inline
	void bin_angle(
		Real const angle_start,
		Real const angle_step,
		Real const ASSERT_ONLY( angle_range ),
		Size const nbins,
		Real const ang,
		Size & bin_lower,
		Size & bin_upper,
		Real & angle_alpha
	) const {
		/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
		/// though it is supposed to return values in the range [-180, 180).
		debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
		debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

		Real real_bin_lower = ( ang - angle_start ) / angle_step;
		Size bin_prev = static_cast< Size > ( real_bin_lower );
		bin_lower = 1 + numeric::mod( bin_prev, nbins );
		bin_upper = numeric::mod( bin_lower, nbins ) + 1;
		angle_alpha = ( (ang - angle_start ) - ( bin_prev * angle_step ) ) / angle_step;
	}


public:
	/// @brief The amino acid this library is representing
	AA
	aa() const {
		return aa_;
	}

	/// @brief When creating rotamer, what position in the CDF should one build until?
	/// Unlikely rotamers ( < 0.5 %) are numerous, but are very infrequently useful.
	Real
	probability_to_accumulate_while_building_rotamers( bool buried ) const;

	/// @brief setters for accumulation probability cutoff (to support externally-controlled option dependence)
	void prob_to_accumulate( Real, Real );
	void prob_to_accumulate_buried( Real );
	void prob_to_accumulate_nonburied( Real );

	///DOUG DOUG DOUG This can be removed but might need to be converted to strings. However it seems that in our case they could be member variables since we really are having SINGLE RESIDUE rot libs and these provide look up for each.
	/// @brief Hard coded specifics about the amino acids
	static
	void
	n_rotamer_bins_for_aa(
		chemical::AA const aa,
		dunbrack::RotVector & rot
	);

	/// @brief Reports information about the *rotameric* chi only; no details
	/// about the non rotameric chi.
	static
	void
	n_rotameric_bins_for_aa(
		chemical::AA const aa,
		dunbrack::RotVector & rot,
		bool dun02
	);

	/// @brief Hard coded rotamer well info for the 2002 library.
	static
	void
	n_rotamer_bins_for_aa_02(
		chemical::AA const aa,
		dunbrack::RotVector & rot
	);

private:

	/// @brief This function forces the instantiation of virtual templated methods in the derived classes.
	/// Functions like this one are necessary when combining polymorphism and templates.  Though
	/// these functions must be compiled, they need never be called. Do not call this function.
	void hokey_template_workaround();


private:
	/// data
	bool const dun02_; // Are we using the 2002 definitions for rotamer wells?
	AA const aa_;
	Size const n_rotameric_chi_;
	utility::vector1< Size > n_chi_bins_;
	utility::vector1< Size > n_chi_products_; // n_chi_products_[ i ] = prod( j in i+1 to nchi, n_chi_bins_[ j ] );

	Size n_packed_rots_;
	Size n_possible_rots_; // prod( i in 1 to nchi, n_chi_bins_[ i ] );

	Real prob_to_accumulate_buried_, prob_to_accumulate_nonburied_;

	utility::vector1< bool > rotwell_exists_;
	bool packed_rotno_conversion_data_current_;

	utility::vector1< Size > rotno_2_packed_rotno_;
	utility::vector1< Size > packed_rotno_2_rotno_;
	utility::vector1< utility::vector1< Size > > packed_rotno_2_rotwell_;

};


} // rotamers
} // pack
} // core

#endif // INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_HH

