// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/copy_dofs/CopyDofs.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_copy_dofs_CopyDofs_HH
#define INCLUDED_core_pose_copy_dofs_CopyDofs_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/copydofs/CopyDofs.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

namespace core {
namespace pose {
namespace copydofs {

extern core::Size const FIXED_DOMAIN;

class CopyDofs: public utility::pointer::ReferenceCount {

public:

	//constructor
	CopyDofs( pose::MiniPose const & template_pose,
		std::map< id::AtomID, id::AtomID > const & atom_id_map,
		std::map< id::AtomID, Size > const & atom_id_domain_map );

	//constructor
	CopyDofs( pose::MiniPose const & template_pose,
		std::map< id::AtomID, id::AtomID > const & atom_id_map );

	//destructor
	~CopyDofs();

public:

	void apply( pose::Pose & pose );

	void figure_out_dofs( pose::Pose & pose );

	CopyDofsInfo copy_dofs_info() const { return copy_dofs_info_; }

private:

	bool
	get_scratch_atom_id( id::AtomID & other_scratch_atom_id,
		std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
		core::kinematics::tree::AtomCOP other_atom );

	bool
	check_domain_map( std::map< id::AtomID, Size > const & atom_id_domain_map,
		id::AtomID const & atom_id1,
		id::AtomID const & atom_id2 ) const;

	bool
	check_domain_map( std::map< id::AtomID, Size > const & atom_id_domain_map,
		utility::vector1< id::AtomID > const & atom_ids1,
		utility::vector1< id::AtomID > const & atom_ids2 ) const;

private:

	pose::MiniPose const & scratch_pose_; // template_pose
	std::map< id::AtomID, id::AtomID > const & atom_id_map_;
	std::map< id::AtomID, Size > atom_id_domain_map_; // blank by default.
	CopyDofsInfo  copy_dofs_info_;
};


} //copy_dofs
} //pose
} //core

#endif
