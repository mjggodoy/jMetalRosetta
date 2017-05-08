// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/simple_moves/RingConformationMoverCreator.hh
/// @brief   Method declarations for RingConformationMoverCreator.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_simple_moves_RingConformationMoverCreator_HH
#define INCLUDED_protocols_simple_moves_RingConformationMoverCreator_HH

// Project headers
#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace simple_moves {

/// @brief  MoverCreator allowing the MoverFactory to create a RingConformationMover
class RingConformationMoverCreator: public protocols::moves::MoverCreator {
public:
	/// @brief  Return an up-casted owning pointer (MoverOP) to the mover.
	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;

	/// @brief  Return the string identifier for the associated Mover (RingConformationMover).
	// XRW TEMP  std::string keyname() const override;

	/// @brief  Static method that returns the keyname for performance reasons.
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_simple_moves_protocols_RingConformationMoverCreator_HH
