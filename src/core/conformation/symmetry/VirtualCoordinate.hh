// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Symmetry data container
/// @file   core/conformation/symmetry/SymmData.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_VirtualCoordinate_hh
#define INCLUDED_core_conformation_symmetry_VirtualCoordinate_hh

// Utility headers
#include <core/conformation/symmetry/VirtualCoordinate.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iosfwd>

namespace core {
namespace conformation {
namespace symmetry {

class VirtualCoordinate {

public:

	/// @brief Default constructor.
	///
	VirtualCoordinate();

	/// @brief copy constructor
	VirtualCoordinate( VirtualCoordinate const & src );

	/// @brief Non-mirror constructor
	///
	VirtualCoordinate(
		numeric::xyzVector< core::Real> axis_x,
		numeric::xyzVector< core::Real> axis_y,
		numeric::xyzVector< core::Real> axis_origin
	);

	/// @brief Mirror constructor
	///
	VirtualCoordinate(
		numeric::xyzVector< core::Real> axis_x,
		numeric::xyzVector< core::Real> axis_y,
		numeric::xyzVector< core::Real> axis_origin,
		bool mirror_z
	);

	VirtualCoordinate &
	operator=( VirtualCoordinate const & src );

	~VirtualCoordinate();

	// @details accessor functions
	numeric::xyzVector< core::Real> &
	get_x();

	numeric::xyzVector< core::Real> &
	get_y();

	numeric::xyzVector< core::Real> &
	get_origin();

	numeric::xyzVector< core::Real> const &
	get_x() const;

	numeric::xyzVector< core::Real> const &
	get_y() const;

	numeric::xyzVector< core::Real> const &
	get_origin() const;

	bool
	get_mirror_z() const;

	void
	set_mirror_z( bool val );

	void
	add_coordinate_from_string(
		utility::vector1< std::string > coords,
		core::Size coord_start=2
	);

	friend
	bool
	operator==(VirtualCoordinate const & a, VirtualCoordinate const & b);

	friend
	bool
	operator!=(VirtualCoordinate const & a, VirtualCoordinate const & b);


private:

	numeric::xyzVector< core::Real> axis_x_; // store unit vector for X
	numeric::xyzVector< core::Real> axis_y_; // store unit vector for Y
	numeric::xyzVector< core::Real> axis_origin_; // store origin for coordinate system
	bool mirror_Z_;
};

} // symmetry
} // conformation
} // core
#endif
