// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PDBPoseInputterCreator.hh
/// @brief  Declaration of the %PDBPoseInputterCreator class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputterCreator_HH
#define INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputterCreator_HH

//unit headers
#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>

// Package headers
#include <protocols/jd3/pose_inputters/PoseInputter.fwd.hh>

// utility headers
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

class PDBPoseInputterCreator : public PoseInputterCreator
{
public:
	virtual PoseInputterOP create_inputter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	virtual void list_options_read( utility::options::OptionKeyList & read_options ) const;
};

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_inputters_PoseInputter_HH
