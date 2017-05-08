// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CopyDofMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/id/AtomID.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <utility>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.CopyDofMover" );

/////////////////////////////////////////////////////////////////////////////////////
//
// copy_dofs() is used a lot in the Das Lab -- identifies atoms in pose that
//  have same names as in template pose, and then smartly updates internal coordinates
//  to match.
//
// This Mover could be made even more powerful if it cached the DOF_IDs and values that
//  need to be applied based on an 'example' pose.
//
//  --rhiju, feb. 2014

namespace protocols {
namespace simple_moves {

using namespace core;

//Constructor
CopyDofMover::CopyDofMover( pose::Pose const & template_pose, std::map< Size, Size > res_map ):
	template_pose_( template_pose ),
	template_mini_pose_( template_pose ),
	res_map_(std::move( res_map )),
	backbone_only_( false ),
	side_chain_only_( false ),
	ignore_virtual_( false ),
	use_hash_( true ),
	pose_string_( "" )
{}

//Destructor
CopyDofMover::~CopyDofMover() = default;

////////////////////////////////////////////////////////////
void
CopyDofMover::apply( core::pose::Pose & pose )
{
	using namespace core::pose::copydofs;

	pose_string_ = pose_string( pose );
	if ( use_hash_ && check_for_precomputed_copy_dofs_info( pose ) ) {
		// use precomputed copy dofs info
		core::pose::copydofs::apply_dofs( pose, copy_dofs_info_[ pose_string_ ] );
	} else {
		std::map < core::id::AtomID , core::id::AtomID > atom_id_map;
		setup_atom_id_map_match_atom_names( atom_id_map, res_map_, pose, template_pose_, backbone_only_, side_chain_only_, ignore_virtual_ );
		std::map< id::AtomID, Size > atom_id_domain_map = blank_atom_id_domain_map( pose );
		CopyDofs copy_dofs( template_mini_pose_, atom_id_map, atom_id_domain_map );
		copy_dofs.apply( pose );
		if ( use_hash_ ) {
			copy_dofs_info_[ pose_string_ ]  = copy_dofs.copy_dofs_info();
		}
	}

}

////////////////////////////////////////////////////////////
// might be a useful function for core/pose/util.hh
std::string
CopyDofMover::pose_string(
	pose::Pose const & pose,
	utility::vector1< Size > const & res_list ) const
{
	std::string pose_string = "";
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		pose_string += pose.residue( res_list[ n ] ).name();
	}
	pose_string += pose.fold_tree().to_string();
	return pose_string;
}

////////////////////////////////////////////////////////////
std::string
CopyDofMover::pose_string(
	pose::Pose const & pose
) const
{
	utility::vector1< Size > res_list;
	for ( auto const & it : res_map_ ) {
		res_list.push_back( it.first );
	}
	return pose_string( pose, res_list );
}
////////////////////////////////////////////////////////////
bool
CopyDofMover::check_for_precomputed_copy_dofs_info( pose::Pose const & pose ) const
{
	return ( copy_dofs_info_.find( pose_string( pose ) ) != copy_dofs_info_.end() );
}

////////////////////////////////////////////////////////////
pose::copydofs::CopyDofsInfo const &
CopyDofMover::copy_dofs_info( pose::Pose const & pose ) const
{
	runtime_assert( check_for_precomputed_copy_dofs_info( pose ) );
	return copy_dofs_info_.find( pose_string( pose ) )->second;
}

////////////////////////////////////////////////////////////
void
CopyDofMover::set_copy_dofs_info(
	pose::Pose const & pose,
	pose::copydofs::CopyDofsInfo const & copy_dofs_info ) {
	copy_dofs_info_[ pose_string( pose ) ] = copy_dofs_info;
}


} //simple_moves
} //protocols
