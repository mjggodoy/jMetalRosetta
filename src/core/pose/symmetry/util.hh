// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pose/symmetry/util.hh
/// @brief utility functions for handling with symmetric conformations
/// @author Ingemar Andre
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_core_pose_symmetry_util_hh
#define INCLUDED_core_pose_symmetry_util_hh


// Unit headers
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.hh>
#include <string>

namespace core {
namespace pose {
namespace symmetry {

bool
is_symmetric( pose::Pose const & pose );

bool
is_mirror_symmetric( pose::Pose const & pose );

bool
is_symmetric( scoring::ScoreFunction const & scorefxn );

bool
is_symmetric( scoring::Energies const & energies );

// This is a stopgap measure to refactoring SymmetricScoreFunction to
// be interalized into the ScoreFunciton class to shield this logic
// from the user.
void
make_score_function_consistent_with_symmetric_state_of_pose(
	pose::Pose const & pose,
	scoring::ScoreFunctionOP & scorefxn
);

conformation::symmetry::SymmetryInfoCOP
symmetry_info( pose::Pose const & pose );

void
make_symmetric_pose(
	pose::Pose & pose,
	conformation::symmetry::SymmetryInfo symmetry_info
);

void
make_symmetric_pose(
	pose::Pose & pose,
	conformation::symmetry::SymmData & symmdata
);

void
make_symmetric_pose(
	pose::Pose & pose,
	std::string symmdef_file=""
);

void
make_asymmetric_pose(
	pose::Pose & pose
);

/// @brief extract the asu from a pose... unlike previous function symmetric clones are thrown out
/// @param[in]  pose_in            Symmetric input pose containing the asymmetric subunit of interest
/// @param[out] pose_out           Asymmetric subunit will be placed into this object
/// @param[in]  with_virtual_atoms If true, virtual atoms related to symmetry will be kept with the asymmetric subunit.
///                                If false, virtual atoms will be removed (default=true)
/// @param[in]  with_unknown_aa    If false, amino acids with type aa_unk will be ignored.  If true, amino acids with
///                                type aa_unk will be extracted (default=false)
void
extract_asymmetric_unit(
	core::pose::Pose const & pose_in,
	core::pose::Pose & pose_out,
	bool const with_virtual_atoms = true,
	bool const with_unknown_aa = false
);

core::pose::Pose
get_asymmetric_pose_copy_from_symmetric_pose(
	pose::Pose const & pose
);

// @details make symmetric PDBIinfo
void
make_symmetric_pdb_info(
	pose::Pose const & pose,
	pose::PDBInfoOP pdb_info_src,
	pose::PDBInfoOP pdb_info_target
);

// @details extract the pdbInfo from the asymmetric unit
void
extract_asymmetric_unit_pdb_info(
	pose::Pose const & pose,
	pose::PDBInfoCOP pdb_info_src,
	pose::PDBInfoOP pdb_info_target
);

void
make_symmetric_movemap(
	pose::Pose const & pose,
	kinematics::MoveMap & movemap
);

int
find_symmetric_basejump_anchor( pose::Pose & pose );

int
find_symmetric_basejump_anchor( pose::Pose & pose );

void
find_new_symmetric_jump_residues( core::pose::Pose & pose );

void
rotate_anchor_to_x_axis( core::pose::Pose & pose );

// find symm axis
numeric::xyzVector< core::Real >
get_symm_axis( core::pose::Pose & pose );

/// @brief    Converts an asymmetric foldtree (f) with virtual root into a
///           symmetric foldtree compatible with symmetric pose (p)
/// @param    p - A symmetric pose
/// @param    f - An asymmetric foldtree. This foldtree MUST have a virtual root
/// @details  This function does not require the symm data
void
symmetrize_fold_tree( core::pose::Pose const &p, kinematics::FoldTree &f );

void
set_asymm_unit_fold_tree( core::pose::Pose &p, kinematics::FoldTree const &f);

// symmetry-aware version of FoldTree::partition_by_jump().  Accepts multiple jumps.
void
partition_by_symm_jumps(
	utility::vector1< int > jump_numbers,
	core::kinematics::FoldTree const & ft,
	conformation::symmetry::SymmetryInfoCOP symm_info,
	ObjexxFCL::FArray1D_bool &partner1 );

// For single component systems
core::pose::Pose
get_buildingblock_and_neighbor_subs(core::pose::Pose const & pose_in, utility::vector1<Size> intra_subs);

// For Multicomponent systems
core::pose::Pose
get_subpose(core::pose::Pose const & pose, utility::vector1<std::string> subs);

utility::vector1<Size>
get_resis(core::pose::Pose const & pose, utility::vector1<std::string> subs);

utility::vector1<std::string>
get_full_intracomponent_and_neighbor_subs(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_full_intracomponent_and_neighbor_subpose(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_full_intracomponent_and_neighbor_resis(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_full_intracomponent_subpose(core::pose::Pose const &pose, std::string sym_dof_name);

utility::vector1<Size>
get_full_intracomponent_resis(core::pose::Pose const &pose, std::string sym_dof_name);

utility::vector1<std::string>
get_full_intracomponent_neighbor_subs(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_full_intracomponent_neighbor_subpose(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_full_intracomponent_neighbor_resis(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

bool
intracomponent_contact(core::pose::Pose const &pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<std::string>
get_intracomponent_and_neighbor_subs(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_intracomponent_and_neighbor_subpose(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_intracomponent_and_neighbor_resis(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<std::string>
get_intracomponent_subs(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_intracomponent_subpose(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_intracomponent_resis(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<std::string>
get_neighbor_subs(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_neighbor_subpose(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_neighbor_resis(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<std::string>
get_intracomponent_and_intraneighbor_subs(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_intracomponent_and_intraneighbor_subpose(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_intracomponent_and_intraneighbor_resis(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<std::string>
get_intracomponent_and_interneighbor_subs(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

core::pose::Pose
get_intracomponent_and_interneighbor_subpose(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

utility::vector1<Size>
get_intracomponent_and_interneighbor_resis(core::pose::Pose const & pose, std::string sym_dof_name, Real contact_dist=10.0);

// make a residue mask (like that used to restrict residues to repack) symmetric
void
make_residue_mask_symmetric( core::pose::Pose const &p, utility::vector1< bool > & msk );

core::kinematics::FoldTree
sealed_symmetric_fold_tree( core::pose::Pose  & pose );

int
get_sym_aware_jump_num( core::pose::Pose const & pose, core::Size jump_num );

utility::vector1<std::string>
sym_dof_names(core::pose::Pose const & pose);

int
sym_dof_jump_num(core::pose::Pose const & pose, std::string const & jname);

std::string
jump_num_sym_dof(core::pose::Pose const & pose, Size const & jnum);

utility::vector1<Size>
get_symdof_subunits(core::pose::Pose const & pose, std::string const & jname);

utility::vector1<char> symmetric_components(core::pose::Pose const & pose);
bool is_multicomponent (core::pose::Pose const & pose);
bool is_singlecomponent(core::pose::Pose const & pose);
Size get_component_lower_bound(core::pose::Pose const & pose, char c);
Size get_component_upper_bound(core::pose::Pose const & pose, char c);
char get_component_of_residue(core::pose::Pose const & pose, Size ir);
char get_subunit_name_to_component(core::pose::Pose const & pose, std::string const & vname);
utility::vector1<char> const & get_jump_name_to_components(core::pose::Pose const & pose, std::string const & jname);
utility::vector1<Size> const & get_jump_name_to_subunits(core::pose::Pose const & pose, std::string const & jname);
std::string get_resnum_to_subunit_component(core::pose::Pose const & pose, Size const & resnum);
utility::vector1<std::string> get_full_intracomponent_subs(core::pose::Pose const & pose, std::string const & jname);

} // symmetry
} // pose
} // core


#endif
