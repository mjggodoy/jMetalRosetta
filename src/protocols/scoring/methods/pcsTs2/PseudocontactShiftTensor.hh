// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @brief Hold chi- tensor information for the Pseudocontact Shift calculation
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890 JBNMR 2008  41:179-189 schmitz et all will explains the tensor convention used
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcsTs2_PseudocontactShiftTensor_hh
#define INCLUDED_protocols_scoring_methods_pcsTs2_PseudocontactShiftTensor_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <string>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs2 {

class PCS_tensor_Ts2 {
private:

	core::Real a_; //alpha
	core::Real b_; //beta
	core::Real g_; //gamma
	core::Real ax_; //axial component
	core::Real rh_; //rhombic component
	//ideally, alpha beta gamma ax and rh will follow the UTR convention
	//see numbat paper JBNMR 2008  41:179-189 schmitz et all

	core::Real chi_xx_;
	core::Real chi_yy_;
	core::Real chi_xy_;
	core::Real chi_xz_;
	core::Real chi_yz_;

	std::string label_; //most likely the filename of the pcs exepriment
	//Used to identify which data has been used in case of multiple lanthanides

private:
	void
	set_abgar(); //Calculate and set the a_, b_, g_, ax_ rh_ members
	//from the chi_xx_ etc...

public:
	PCS_tensor_Ts2();

	~PCS_tensor_Ts2();

	PCS_tensor_Ts2(PCS_tensor_Ts2 const & other);

	PCS_tensor_Ts2 &
	operator=(PCS_tensor_Ts2 const & other);

	PCS_tensor_Ts2(  core::Real const chi_xx,
		core::Real const chi_yy,
		core::Real const chi_xy,
		core::Real const chi_xz,
		core::Real const chi_yz,
		std::string const label);

	core::Real
	delta_X_xx() const;

	core::Real
	delta_X_yy() const;

	core::Real
	delta_X_zz() const;

	core::Real
	delta_chi_zz() const;

	core::Real
	chi_xx() const;

	core::Real
	chi_yy() const;

	core::Real
	chi_xy() const;

	core::Real
	chi_xz() const;

	core::Real
	chi_yz() const;

	void
	copy_from_ref(PCS_tensor_Ts2 & other);

	std::string const &
	get_label() const;

	void
	reset_tensor(core::Real const chi_xx,
		core::Real const chi_yy,
		core::Real const chi_xy,
		core::Real const chi_xz,
		core::Real const chi_yz);

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_tensor_Ts2 &PCS_t);
};


}//namespace pcsTs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
