// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_setup_FullModelInfoSetupFromCommandLine_HH
#define INCLUDED_protocols_stepwise_setup_FullModelInfoSetupFromCommandLine_HH

#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
#include <map>

namespace protocols {
namespace stepwise {
namespace setup {

core::pose::PoseOP
get_pdb_with_full_model_info( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set );

core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file );

core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set );

void
get_other_poses(  utility::vector1< core::pose::PoseOP > & other_poses,
	utility::vector1< std::string > const & other_files,
	core::chemical::ResidueTypeSetCAP rsd_set );

void
initialize_native_and_align_pose( core::pose::PoseOP & native_pose,
	core::pose::PoseOP & align_pose,
	core::chemical::ResidueTypeSetCAP rsd_set,
	core::pose::PoseCOP start_pose );

core::pose::PoseOP
initialize_pose_and_other_poses_from_command_line( core::chemical::ResidueTypeSetCAP rsd_set );

void
cleanup( core::pose::Pose & pose,
	bool const force_cut_at_rna_chainbreak = false );

void
get_extra_cutpoints_from_names( core::Size const nres,
	utility::vector1< core::Size > & cutpoints,
	std::map< core::Size, std::string > const & non_standard_residue_map );

core::pose::full_model_info::FullModelParametersOP
get_sequence_information( std::string const & fasta_file,
	utility::vector1< core::Size > & cutpoint_open_in_full_model );

void
fill_full_model_info_from_command_line( core::pose::Pose & pose );

void
fill_full_model_info_from_command_line( utility::vector1< core::pose::PoseOP > & pose_ops );

void
fill_full_model_info_from_command_line( core::pose::Pose & pose, utility::vector1< core::pose::PoseOP > & other_pose_ops );

void
fill_full_model_info_from_command_line( utility::vector1< core::pose::Pose * > & pose_pointers );

void
setup_fold_trees(
	utility::vector1< core::pose::Pose * > & pose_pointers,
	utility::vector1< core::Size > & cutpoint_open_in_full_model /* can be updated here*/,
	utility::vector1< core::Size > & fixed_domain_map /* can be updated here*/,
	utility::vector1< core::Size > const & cutpoint_closed,
	utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & extra_minimize_jump_res,
	utility::vector1< core::Size > const & sample_res,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & jump_res,
	utility::vector1< core::Size > const & preferred_root_res,
	utility::vector1< core::Size > const & virtual_sugar_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists );
void
update_pose_fold_tree( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_min_res,
	utility::vector1< core::Size > const & sample_res,
	utility::vector1< core::Size > const & jump_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters );

void
define_chains( core::pose::Pose const & pose,
	utility::vector1< utility::vector1< core::Size > > & all_res_in_chain,
	utility::vector1< utility::vector1< core::Size > > & all_fixed_res_in_chain,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_min_res );

void
setup_user_defined_jumps(
	utility::vector1< core::Size > const & jump_res,
	utility::vector1< core::Size > & jump_partners1,
	utility::vector1< core::Size > & jump_partners2,
	utility::vector1< std::pair< core::Size, core::Size > > & chain_connections,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain );

core::Size
get_chain( core::Size const i, utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain );

void
setup_jumps( core::pose::Pose const & pose,
	utility::vector1< core::Size > & jump_partners1,
	utility::vector1< core::Size > & jump_partners2,
	utility::vector1< std::pair< core::Size, core::Size > > & chain_connections,
	utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain,
	std::pair< utility::vector1< int >, utility::vector1< char > > const & resnum_and_chain_in_pose );

core::kinematics::FoldTree
get_tree( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & cuts,
	utility::vector1< core::Size > const & jump_partners1,
	utility::vector1< core::Size > const & jump_partners2 );

core::kinematics::FoldTree
get_tree( core::Size const nres,
	utility::vector1< core::Size > const & cuts,
	utility::vector1< core::Size > const & jump_partners1,
	utility::vector1< core::Size > const & jump_partners2,
	utility::vector1< std::string > const & jump_atoms1,
	utility::vector1< std::string > const & jump_atoms2 );

void
update_fixed_domain_from_extra_minimize_jump_res( utility::vector1< core::Size > & fixed_domain,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_minimize_jump_res );
void
add_cutpoint_closed( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & cutpoint_closed );

void
put_in_cutpoint( core::pose::Pose & pose, core::Size const i );

void
add_virtual_sugar_res( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & virtual_sugar_res );

utility::vector1< core::Size >
figure_out_working_res( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & sample_res );

utility::vector1< core::Size >
figure_out_sample_res( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & working_res );

void
check_working_res( utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & sample_res );

void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map );

void
figure_out_motif_mode( utility::vector1< core::Size > & extra_min_res,
	utility::vector1< core::Size > & terminal_res,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & cutpoint_open_in_full_model );

void
add_block_stack_variants( utility::vector1< core::pose::Pose * > const & pose_pointers,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists,
	utility::vector1< core::Size > const & block_stack_above_res,
	utility::vector1< core::Size > const & block_stack_below_res );

void
update_jump_res( utility::vector1< core::Size > & jump_res,
	utility::vector1< core::Size > const & extra_minimize_jump_res );

void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map );

utility::vector1< core::Size >
figure_out_fixed_domain_map( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & extra_minimize_res );

utility::vector1< core::Size >
figure_out_dock_domain_map( utility::vector1< core::Size > & cutpoint_open_in_full_model,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & sample_res,
	core::Size const nres );

void
reorder_pose( core::pose::Pose & pose, utility::vector1< core::Size > & res_list );

bool
just_modeling_RNA( utility::vector1< std::string > const & fasta_files );

void
setup_for_density_scoring( core::pose::Pose & pose );

} //setup
} //stepwise
} //protocols

#endif
