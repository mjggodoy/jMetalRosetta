// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/hybridization/BackboneTorsionPerturbation.hh
/// @brief  Declaration of the MoverCreator class for the BackboneTorsionPerturbation
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_BackboneTorsionPerturbationCreator_hh
#define INCLUDED_protocols_hybridization_BackboneTorsionPerturbationCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace hybridization {

class BackboneTorsionPerturbationCreator : public moves::MoverCreator
{
public:
	// XRW TEMP  moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	// XRW TEMP  static  std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}

#endif
