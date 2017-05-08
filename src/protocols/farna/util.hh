// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_ProtocolUtil_HH
#define INCLUDED_protocols_rna_RNA_ProtocolUtil_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh> //courtesy to other namespaces that want clear_rna_secstruct_legacy_info()
#include <protocols/farna/base_pairs/BasePairStep.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <map>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


namespace protocols {
namespace farna {

typedef  numeric::xyzVector< core::Length >  Vector;

void
figure_out_secstruct( core::pose::Pose & pose );

void
get_base_pairing_info( core::pose::Pose const & pose,
	core::Size const & seqpos,
	char & secstruct,
	ObjexxFCL::FArray1D <bool> & edge_is_base_pairing );

void
get_base_pairing_list( core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size> > & base_pairing_list );

void
create_rna_vall_torsions( core::pose::Pose & pose,
	utility::io::ozstream & torsions_out,
	utility::vector1 <core::Size> const & exclude_res_list );

void
create_rna_vall_torsions( core::pose::Pose & pose,
	std::string const & outfile,
	utility::vector1 <core::Size> const & exclude_res_list );

void
assert_phosphate_nomenclature_matches_mini( core::pose::Pose const & pose);

void
ensure_phosphate_nomenclature_matches_mini( core::pose::Pose & pose );

void
export_packer_results(
	utility::vector1< std::pair< core::Real, std::string > > & results,
	utility::vector1< core::pose::PoseOP > pose_list,
	core::scoring::ScoreFunctionOP & scorefxn,
	std::string const & outfile,
	bool const dumo = false );

void
check_base_pair( core::pose::Pose & pose, ObjexxFCL::FArray1D_int & struct_type );

void
setup_base_pair_constraints(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size > > const &  pairings,
	core::Real const scale_factor = 1.0,
	bool const use_flat_harmonic = false );

void
setup_coarse_chainbreak_constraints( core::pose::Pose & pose, core::Size const & n );

std::string const
convert_based_on_match_type( std::string const & RNA_string, core::Size const type );

bool
compare_RNA_char( char const char1, char const char2 );

bool
compare_RNA_secstruct( char const char1, char const char2 );

Vector
get_sugar_centroid( core::conformation::Residue const & rsd );

void
make_extended_coarse_pose( core::pose::Pose & coarse_pose, std::string const & full_sequence );

void
make_coarse_pose( core::pose::Pose const & pose, core::pose::Pose & coarse_pose );

void
print_internal_coords( core::pose::Pose const & pose );

void
remove_cutpoint_closed( core::pose::Pose & pose, core::Size const i );

void
remove_cutpoints_closed( core::pose::Pose & pose );

bool
possible_root( core::kinematics::FoldTree const & f, core::Size const & n );

inline bool is_num_in_list ( core::Size const i,
	utility::vector1 <core::Size> const & list )
{
	return std::find(list.begin(), list.end(), i)!=list.end();
}

utility::vector1< core::Size >
get_rigid_body_jumps( core::pose::Pose const & pose );

bool
let_rigid_body_jumps_move( core::kinematics::MoveMap & movemap,
	core::pose::Pose const & pose,
	bool const move_first_rigid_body  = false );

core::Size
get_anchor_rsd( core::pose::Pose const & pose );

void
translate_virtual_anchor_to_first_rigid_body( core::pose::Pose & pose );

bool
involved_in_phosphate_torsion( std::string atomname );

void
set_output_res_and_chain( core::pose::Pose & extended_pose,
	std::pair< utility::vector1< int >, utility::vector1< char > > const & output_resnum_and_chain );

void
figure_out_base_pair_partner( core::pose::Pose & pose, std::map< core::Size, core::Size > & partner,
	bool const strict = true );

utility::vector1< core::pose::rna::BasePair >
classify_base_pairs_lores( core::pose::Pose const & pose );

void
process_input_file( std::string const & silent_file,
	utility::vector1< core::pose::PoseOP > & pose_list,
	bool is_pdb = false,
	bool coarse_rna  = false );

void
print_hbonds( core::pose::Pose & pose );

bool
moveable_jump( core::id::AtomID const & jump_atom_id1,
	core::id::AtomID const & jump_atom_id2,
	protocols::toolbox::AtomLevelDomainMap const & atom_level_domain_map);

bool
moveable_jump( core::Size const jump_pos1,
	core::Size const jump_pos2,
	protocols::toolbox::AtomLevelDomainMap const & atom_level_domain_map);

void
fill_in_default_jump_atoms( core::kinematics::FoldTree & f, core::pose::Pose const & pose );

core::Size
virtualize_bulges( core::pose::Pose & input_pose,
	utility::vector1< core::Size > const & in_allow_bulge_res_list,
	core::scoring::ScoreFunctionCOP const & scorefxn,
	std::string const & tag,
	bool const allow_pre_virtualize,
	bool const allow_consecutive_bulges,
	bool const verbose );

core::scoring::ScoreFunctionOP
get_rna_hires_scorefxn();

utility::vector1< core::Size >
get_moving_res( core::pose::Pose const & pose,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

bool
base_pair_step_moving( base_pairs::BasePairStep const & base_pair_step,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::pose::Pose const & pose );

bool
base_pair_moving( core::pose::rna::BasePair const & base_pair,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::pose::Pose const & pose );

void delete_non_protein_from_pose( core::pose::Pose & pose );

} //farna
} //protocols

#endif
