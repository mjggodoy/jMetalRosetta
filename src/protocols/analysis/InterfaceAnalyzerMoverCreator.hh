// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/InterfaceAnalyzerMoverCreator.hh
/// @brief This class will create instances of Mover InterfaceAnalyzerMover for the MoverFactory
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_analysis_InterfaceAnalyzerMoverCreator_HH
#define INCLUDED_protocols_analysis_InterfaceAnalyzerMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace analysis {

class InterfaceAnalyzerMoverCreator : public protocols::moves::MoverCreator {
public:
	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace analysis
} //namespace protocols

#endif //INCLUDED_protocols_analysis_InterfaceAnalyzerMoverCreator_HH

