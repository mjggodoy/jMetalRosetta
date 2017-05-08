// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/generalized_kinematic_closure/GeneralizedKICCreator.hh
/// @brief This class will create instances of Mover GeneralizedKIC for the MoverFactory
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKICCreator_hh
#define INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKICCreator_hh

//#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace generalized_kinematic_closure {

class GeneralizedKICCreator : public protocols::moves::MoverCreator {
public:
	// XRW TEMP  moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace generalized_kinematic_closure
} //namespace protocols

#endif //INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKICCreator_hh

