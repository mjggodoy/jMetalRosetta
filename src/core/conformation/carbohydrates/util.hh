// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/carbohydrates/util.hh
/// @brief Utility functions that DO NOT require a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_conformation_carbohydrates_util_hh
#define INCLUDED_core_conformation_carbohydrates_util_hh

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>

namespace core {
namespace conformation {
namespace carbohydrates {

/// @brief  Use a saccharide residue's connections to find the residue from which it follows or branches.
///         Returns 0 if it has no parent.
core::uint
find_seqpos_of_saccharides_parent_residue( conformation::Residue const & residue );

/// @brief  Use the mainchain acceptor to find the mainchain child.
///         Typically, this is N+1, but this is independant of the residue number.
///         Returns 0 if it has no child.
core::uint
find_seqpos_of_saccharides_mainchain_child( conformation::Residue const & residue );


/// @brief  Use a saccharide residue's connections to find the residue following it from a given linkage position.
core::uint
find_seqpos_of_saccharides_child_residue_at( conformation::Residue const & residue,
	core::uint linkage_position );
	
// Use a saccharide residue's connections to find its linkage number on the previous residue.
/// @return an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
/// parent monosaccharide residue; e.g., 4 specifies O4; n = 0 specifies that the residue at <seqpos> is a lower
/// terminus or connected to a non-sugar.
core::uint
get_linkage_position_of_saccharide_residue( conformation::Residue const & rsd,
	conformation::Residue const & parent_rsd);
	
	
/// @brief Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic
/// carbon.
bool
has_exocyclic_glycosidic_linkage( conformation::Conformation const & conf, uint const seqpos );

/// @brief Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic
/// carbon.
bool
has_exocyclic_glycosidic_linkage( conformation::Residue const & rsd, conformation::Residue const & parent_rsd );



////////////////////////    Branch Information    ///////////////////////////////////
///
///
///



/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
get_branching_residues(
	conformation::Conformation const & conf,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips );
	
	
/// @brief  Find all children residues, list of residues, and any found tips from a given residue not including parent
///
///  Children Residues:  Filled in list of children residues found if not tips.
///  List Of  Residues:  All the residue nums found.
///  Tips:  All 'ends' of of children found.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
fill_upstream_children_res_and_tips(
	conformation::Conformation const & conf,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips );
	
///@brief Get the size of the glycan tree given the first carbohydrate residue in the tree.
/// On-the-fly calculation
core::Size
get_glycan_tree_size(conformation::Conformation const & conf, core::Size const first_glycan_resnum);

///@brief Get the largest glycan tree size int he pose.
core::Size
get_largest_glycan_tree_size(conformation::Conformation const & conf);

///@brief Get the residue distance from the position to the root/end of the glycan.
/// On-the-fly calculation
core::Size
get_distance_to_root(conformation::Conformation const & conf, core::Size const position);

///@brief Get which residues denote starting a glycan.
///  These are the first residue of the glycan tree, and the tree can be
///  branching from protein  or not.
///
utility::vector1< bool >
get_glycan_start_points(conformation::Conformation const & conf);

/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Convenience function. Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_residues_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position);
	
/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Convenience function. Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position);
	
/// @brief Get residues further down the branch from this residue.  starting_position ->
///  Returns pair of all_upstream_residues, tips.
///  Tips are the ends of linear glycan branches.
std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position,
	bool include_starting_position = false);
	
/// @brief Get the carbohydrate residue connecting the protein branch point.
/// Returns 0 if branch point is not connected to carbohydrate downstream.
///
core::Size
get_glycan_connecting_protein_branch_point(
	conformation::Conformation const & conf,
	core::Size const protein_branch_point_resnum);
	
///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
/// Returns 0 if that glycan_position is not part of the glycan we are interested in or not in pose.
///
core::Size
get_resnum_from_glycan_position(
	conformation::Conformation const & conf,
	core::Size const first_glycan_resnum,
	core::Size const glycan_position);
	
///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
/// Returns 0 if that glycan_position is not part of the glycan we are interested in or not in pose.
///
core::Size
get_glycan_position_from_resnum(
	conformation::Conformation const & conf,
	core::Size const first_glycan_resnum,
	core::Size const resnum);
	
	
} //core
} //conformation
} //carbohydrates


#endif	//core/conformation/carbohydrates_util_hh

