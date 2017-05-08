// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   SavePoseMover.hh
///
/// @brief
/// @author Florian Richter

#ifndef INCLUDED_protocols_rosetta_scripts_SavePoseMover_hh
#define INCLUDED_protocols_rosetta_scripts_SavePoseMover_hh

//unit headers
#include <protocols/rosetta_scripts/SavePoseMover.fwd.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace rosetta_scripts {

/// @brief mover that can be used to save or restore a pose at an arbitrary
/// point during a rosetta scripts protocol. other movers or filters can be
/// set up to access poses saved by this mover during their apply calls.
class SavePoseMover : public moves::Mover
{
public:

	SavePoseMover();
	~SavePoseMover() override;

	void apply( core::pose::Pose & pose  ) override;
	std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	///@brief Returns a cloned PoseOP of reference_pose_
	core::pose::PoseOP get_cloned_saved_pose() const;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::pose::PoseOP reference_pose_;
	bool restore_pose_; //determines whether this mover saves or restores a pose
};


} // rosetta_scripts
} // protocols

#endif
