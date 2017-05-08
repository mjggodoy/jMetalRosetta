// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   core/chemical/rna/RNA_FittedTorsionInfo.hh
/// @brief  Statistically derived torsion information for RNA
/// @author Rhiju Das

#ifndef INCLUDED_core_chemical_rna_RNA_FittedTorsionInfo_HH
#define INCLUDED_core_chemical_rna_RNA_FittedTorsionInfo_HH

//Unit header
#include <core/chemical/rna/RNA_FittedTorsionInfo.fwd.hh>

#include <core/types.hh>
//#include <core/chemical/rna/util.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++

namespace core {
namespace chemical {
namespace rna {

/////////////////////////////////////////
class GaussianParameter {
public:
	Real amplitude, center, width;

	GaussianParameter( Real const amplitude_in, Real const center_in, Real const width_in ):
		amplitude( amplitude_in ),
		center   ( center_in ),
		width    ( width_in )
	{}
};
/////////////////////////////////

class RNA_FittedTorsionInfo : public utility::pointer::ReferenceCount {

public:
	typedef utility::vector1< GaussianParameter > GaussianParameterSet;

	RNA_FittedTorsionInfo();

	virtual ~RNA_FittedTorsionInfo();

	GaussianParameterSet gaussian_parameter_set_alpha() const{ return gaussian_parameter_set_alpha_; }
	GaussianParameterSet gaussian_parameter_set_beta() const{ return gaussian_parameter_set_beta_; }
	GaussianParameterSet gaussian_parameter_set_gamma() const{ return gaussian_parameter_set_gamma_; }
	GaussianParameterSet gaussian_parameter_set_delta_north() const{ return gaussian_parameter_set_delta_north_; }
	GaussianParameterSet gaussian_parameter_set_delta_south() const{ return gaussian_parameter_set_delta_south_; }
	GaussianParameterSet gaussian_parameter_set_epsilon_north() const{ return gaussian_parameter_set_epsilon_north_; }
	GaussianParameterSet gaussian_parameter_set_epsilon_south() const{ return gaussian_parameter_set_epsilon_south_; }
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_sc_minus() const{ return gaussian_parameter_set_zeta_alpha_sc_minus_; }
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_sc_plus() const{ return gaussian_parameter_set_zeta_alpha_sc_plus_; }
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_ap() const{ return gaussian_parameter_set_zeta_alpha_ap_; }
	GaussianParameterSet gaussian_parameter_set_chi_north() const{ return gaussian_parameter_set_chi_north_; }
	GaussianParameterSet gaussian_parameter_set_chi_south() const{ return gaussian_parameter_set_chi_south_; }
	GaussianParameterSet gaussian_parameter_set_nu2_north() const{ return gaussian_parameter_set_nu2_north_; }
	GaussianParameterSet gaussian_parameter_set_nu2_south() const{ return gaussian_parameter_set_nu2_south_; }
	GaussianParameterSet gaussian_parameter_set_nu1_north() const{ return gaussian_parameter_set_nu1_north_; }
	GaussianParameterSet gaussian_parameter_set_nu1_south() const{ return gaussian_parameter_set_nu1_south_; }

	Real delta_north() const{ return ideal_delta_north_; }
	Real nu2_north() const{ return ideal_nu2_north_; }
	Real nu1_north() const{ return ideal_nu1_north_; }

	Real delta_south() const{ return ideal_delta_south_; }
	Real nu2_south() const{ return ideal_nu2_south_; }
	Real nu1_south() const{ return ideal_nu1_south_; }

	Real alpha_aform() const{ return gaussian_parameter_set_alpha_[1].center; }
	Real beta_aform() const{ return gaussian_parameter_set_beta_[1].center; }
	Real gamma_aform() const{ return gaussian_parameter_set_gamma_[1].center; }
	Real epsilon_aform() const{ return gaussian_parameter_set_epsilon_north_[1].center; }
	Real zeta_aform() const{ return gaussian_parameter_set_zeta_alpha_sc_minus_[1].center; }

	Real epsilon_north() const{ return gaussian_parameter_set_epsilon_north_[1].center; }
	Real epsilon_south() const{ return gaussian_parameter_set_epsilon_south_[1].center; }

	Real chi_north_anti() const{ return gaussian_parameter_set_chi_north_[1].center; }
	Real chi_north_syn() const{ return gaussian_parameter_set_chi_north_[2].center; }
	Real chi_south_anti() const{ return gaussian_parameter_set_chi_south_[1].center; }
	Real chi_south_syn() const{ return gaussian_parameter_set_chi_south_[2].center; }

	Real delta_cutoff() const { return delta_cutoff_; }

private:


	void
	init_rna_torsion_gaussian_parameters();

	bool rna_tight_torsions_;

	GaussianParameterSet gaussian_parameter_set_alpha_;
	GaussianParameterSet gaussian_parameter_set_beta_;
	GaussianParameterSet gaussian_parameter_set_gamma_;
	GaussianParameterSet gaussian_parameter_set_delta_north_;
	GaussianParameterSet gaussian_parameter_set_delta_south_;
	GaussianParameterSet gaussian_parameter_set_epsilon_north_;
	GaussianParameterSet gaussian_parameter_set_epsilon_south_;
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_sc_minus_;
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_sc_plus_;
	GaussianParameterSet gaussian_parameter_set_zeta_alpha_ap_;
	GaussianParameterSet gaussian_parameter_set_chi_north_;
	GaussianParameterSet gaussian_parameter_set_chi_south_;
	GaussianParameterSet gaussian_parameter_set_nu2_north_;
	GaussianParameterSet gaussian_parameter_set_nu2_south_;
	GaussianParameterSet gaussian_parameter_set_nu1_north_;
	GaussianParameterSet gaussian_parameter_set_nu1_south_;

	Real const ideal_delta_north_;
	Real const ideal_nu2_north_;
	Real const ideal_nu1_north_;

	Real const ideal_delta_south_;
	Real const ideal_nu2_south_;
	Real const ideal_nu1_south_;

	Real const delta_cutoff_;

};

}
}
}

#endif
