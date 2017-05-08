// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pose_selectors/ClusterPoseSelectorCreator.hh
/// @brief  Creator for ClusterPoseSelector
/// @author Luki Goldschmidt <luki@uw.edu>

#ifndef INCLUDED_protocols_pose_selectors_ClusterPoseSelectorCreator_hh
#define INCLUDED_protocols_pose_selectors_ClusterPoseSelectorCreator_hh

#include <protocols/rosetta_scripts/PoseSelectorCreator.hh>
#include <string>

namespace protocols {
namespace pose_selectors {

class ClusterPoseSelectorCreator : public protocols::rosetta_scripts::PoseSelectorCreator {
public:
	protocols::rosetta_scripts::PoseSelectorOP create_selector() const override;
	std::string keyname() const override { return "ClusterPoseSelector"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}

#endif

