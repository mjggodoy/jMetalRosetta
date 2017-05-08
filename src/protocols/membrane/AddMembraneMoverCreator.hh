// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneMover.hh
///
/// @brief      Initialize the RosettaMP Framework by adding membrane representations to the pose
/// @details Given a pose, initialize and configure with the RosettaMP framework by taking the
///    following steps:
///     (1) Add a membrane residue to the pose (type MEM)
///      (a) Append residue to the pose or accept a new one
///      (b) Update PDB info to acknowledge the membrane residue
///      (c) Set the MEM residue at the root of the foldtree
///     (2) Initialize transmembrane spanning topology (Object: SpanningTopology)
///     (3) Initialize the MembraneInfo object either with or without the LipidAccInfo object.
///     (4) Set the membrane starting position (either default or based on user input)
///
///    This object does a massive amount of reading from CMD, RosettaScripts or Constructor. If you add
///    a new piece of data - you must modify MembraneInfo, all of these channels for input AND the apply function!
///    If and only if AddMembraneMover is applied to the pose, pose.conformation().is_membrane() MUST return true.
///
///    Last Updated: 7/23/15
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh
#define INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class AddMembraneMoverCreator : public protocols::moves::MoverCreator {

public:

	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh
