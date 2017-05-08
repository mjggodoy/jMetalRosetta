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
/// @file GridSearchIterator.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcsTs3_GridSearchIterator_hh
#define INCLUDED_protocols_scoring_methods_pcsTs3_GridSearchIterator_hh


// Unit headers

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers

// Objexx headers

// C++ headers


namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs3 {

class GridSearchIterator_Ts3 : public utility::pointer::ReferenceCount {
public:
	GridSearchIterator_Ts3();

	~GridSearchIterator_Ts3();

	GridSearchIterator_Ts3(GridSearchIterator_Ts3 const & other);

	GridSearchIterator_Ts3 &
	operator=(GridSearchIterator_Ts3 const & other);

	GridSearchIterator_Ts3(numeric::xyzVector< core::Real > const coo1,
		numeric::xyzVector< core::Real > const coo2,
		core::Real const k,
		core::Real const edge_size,
		core::Real const step_size,
		core::Real const small_cutoff,
		core::Real const large_cutoff,
		core::Real const cone_angle);

	bool
	next(core::Real &x,
		core::Real &y,
		core::Real &z);

	bool
	next_center(core::Real &x,
		core::Real &y,
		core::Real &z);

	void
	reset();

private:
	core::Real const x_center_;
	core::Real const y_center_;
	core::Real const z_center_;
	core::Real x_current_;
	core::Real y_current_;
	core::Real z_current_;
	core::Real const step_;
	core::Real const edge_;
	core::Real step_x_;
	core::Real step_y_;
	core::Real step_z_;
	bool next_to_return_;
	core::Real const delta_;
	core::Real const small_cutoff_square_;
	core::Real const large_cutoff_square_;
	core::Real const x_vector_;
	core::Real const y_vector_;
	core::Real const z_vector_;
	core::Real const norme_vector_;
	core::Real const cone_angle_cos_;
};

} //namespace pcsTs3
} //namespace methods
} //namespace scoring
} //namespace protocols

#endif
