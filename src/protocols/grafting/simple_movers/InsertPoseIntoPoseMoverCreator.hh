// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/grafting/simple_movers/DeleteRegionMover.hh
/// @brief
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_simple_movers_INSERTPOSEINTOPOSEMOVERCREATOR_HH
#define INCLUDED_protocols_grafting_simple_movers_INSERTPOSEINTOPOSEMOVERCREATOR_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

class InsertPoseIntoPoseMoverCreator : public protocols::moves::MoverCreator {
public:
	// XRW TEMP  virtual protocols::moves::MoverOP create_mover() const;
	// XRW TEMP  virtual std::string keyname() const;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;



};

}
}
}

#endif //#ifndef INCLUDED_protocols/antibody_design_INSERTPOSEINTOPOSEMOVERCREATOR_HH

