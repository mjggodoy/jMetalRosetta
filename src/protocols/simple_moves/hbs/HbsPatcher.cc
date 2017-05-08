// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/hbs/HbsPatcher.cc
/// @brief HbsPatcher methods implemented
/// @author Andy Watkins, amw579@nyu.edu

// Unit Headers
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.hbs.HbsPatcher" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace hbs {

void HbsPatcher::apply( core::pose::Pose & pose )
{
	TR<< "patching residues" <<std::endl;

	//awatkins: an hbs pre position cannot be last position
	runtime_assert_msg ( hbs_pre_pos_ != pose.size(), "beginning of hbs cannot be last residue" );
	// I believe that this should be terminal, but since we're manually cutting off later residues
	// I bet we can't require it that way, so we have to trust hbs creator not to be dumb.
	//awatkins: an hbs post position cannot be first position
	runtime_assert ( hbs_post_pos_ != 1 );

	chemical::ResidueTypeSetCOP restype_set = pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );

	std::string const pre_base_name( core::chemical::residue_type_base_name( pose.residue_type( hbs_pre_pos_ ) ) );
	std::string const post_base_name( core::chemical::residue_type_base_name( pose.residue_type( hbs_post_pos_ ) ) );
	TR << "pre restype basename: " << pre_base_name << std::endl;
	TR << "post restype basename: " << post_base_name << std::endl;

	//awatkins: check for proline
	if ( pre_base_name == "PRO" || pre_base_name == "DPRO" ||
			post_base_name == "PRO" || post_base_name == "DPRO" ) {
		utility_exit_with_message("Cannot patch proline");
	}
	if ( pose.residue(hbs_pre_pos_).has_variant_type(chemical::HBS_POST) == 1 ) {
		utility_exit_with_message("Cannot patch HBS_PRE on an HBS_POST");
	}
	if ( pose.residue(hbs_post_pos_).has_variant_type(chemical::HBS_PRE) == 1 ) {
		utility_exit_with_message("Cannot patch HBS_POST on an HBS_PRE");
	}

	//awatkins: check if already patched
	if ( !pose.residue( hbs_pre_pos_ ).has_variant_type( chemical::HBS_PRE ) ) {
		TR<< "patching pre" <<std::endl;

		//awatkins: get base residue type
		chemical::ResidueType const & pre_base_type = pose.residue(hbs_pre_pos_).type();
		TR<< pre_base_type.name() << std::endl;

		//awatkins: add variant

		std::string const base_name( core::chemical::residue_type_base_name( pre_base_type ) );

		// the desired set of variant types:
		utility::vector1< std::string > target_variants( pre_base_type.properties().get_list_of_variants() );
		if ( !pre_base_type.has_variant_type( chemical::HBS_PRE ) ) {
			target_variants.push_back( "HBS_PRE" );
			target_variants.push_back( "LOWER_TERMINUS_VARIANT" );
		}
		if ( pre_base_type.has_variant_type( chemical::ACETYLATED_NTERMINUS_VARIANT ) ) {
			target_variants.erase( std::find( target_variants.begin(), target_variants.end(), "ACETYLATED_NTERMINUS_VARIANT" ) );
		}
		if ( pre_base_type.has_variant_type( chemical::METHYLATED_NTERM_VARIANT ) ) {
			target_variants.erase( std::find( target_variants.begin(), target_variants.end(), "METHYLATED_NTERM_VARIANT" ) );
		}
		if ( pre_base_type.has_variant_type( chemical::N_ACETYLATION ) ) {
			target_variants.erase( std::find( target_variants.begin(), target_variants.end(), "N_ACETYLATION" ) );
		}

		ResidueTypeCOP rsd = ResidueTypeFinder( *restype_set ).residue_base_name( base_name ).variants( target_variants ).get_representative_type();
		conformation::Residue replace_res_pre( *rsd, true );
		replace_res_pre.set_all_chi(pose.residue(hbs_pre_pos_).chi());
		replace_res_pre.residue_connection_partner( 2, hbs_post_pos_, 3 );
		//replace_res_pre.mainchain_torsions(pose.residue(hbs_pre_pos_).mainchain_torsions());
		//replace_res_pre.update_residue_connection_mapping();
		pose.replace_residue( hbs_pre_pos_, replace_res_pre, true );
		conformation::idealize_position( hbs_pre_pos_, pose.conformation() );
		TR<< replace_res_pre.name() << std::endl;

		//conformation::Residue replace_res_pre( *new_type, true );
		//pose.dump_pdb( "rosetta_out_hbs_post_patch.pdb" );

	}// if pre
	if ( !pose.residue(hbs_post_pos_).has_variant_type( chemical::HBS_POST ) ) {
		TR<< "patching post" <<std::endl;
		//awatkins: get base residue type
		chemical::ResidueType const & post_base_type = pose.residue(hbs_post_pos_).type();
		TR<< post_base_type.name() << std::endl;

		std::string const base_name( core::chemical::residue_type_base_name( post_base_type ) );

		// the desired set of variant types:
		utility::vector1< std::string > target_variants( post_base_type.properties().get_list_of_variants() );
		if ( !post_base_type.has_variant_type( chemical::HBS_POST ) ) {
			target_variants.push_back( "HBS_POST" );
		}

		ResidueTypeCOP rsd = ResidueTypeFinder( *restype_set ).residue_base_name( base_name ).variants( target_variants ).get_representative_type();
		//restype_set->make_sure_instantiated( rsd.get_self_ptr() );
		//*new_type = rsd;
		conformation::Residue replace_res_post( *rsd, true );
		replace_res_post.set_all_chi(pose.residue(hbs_post_pos_).chi());
		replace_res_post.residue_connection_partner( 3, hbs_pre_pos_, 2 );
		//replace_res_pre.update_residue_connection_mapping();
		//replace_res_pre.mainchain_torsions(pose.residue(hbs_pre_pos_).mainchain_torsions());

		pose.replace_residue( hbs_post_pos_, replace_res_post, true );
		conformation::idealize_position( hbs_post_pos_, pose.conformation() );
		TR<< replace_res_post.name() << std::endl;

	}// if post

	core::pose::ncbb::add_hbs_constraint( pose, hbs_pre_pos_ );

	//awatkins: need to do all at once at the end because occasionally will get error:
	//awatkins: Unable to handle change in the number of residue connections in the presence of pseudobonds!
	//pose.conformation().detect_bonds();
	//pose.conformation().detect_pseudobonds();
	pose.conformation().update_polymeric_connection( hbs_pre_pos_ );
	pose.conformation().update_polymeric_connection( hbs_post_pos_ );

}


std::string
HbsPatcher::get_name() const {
	return "HbsPatcher";
}

/// @brief
HbsPatcher::HbsPatcher(
	core::Size hbs_seq_position
): Mover(), hbs_pre_pos_(hbs_seq_position), hbs_post_pos_(hbs_seq_position+2)
{
	Mover::type( "HbsPatcher" );

}

HbsPatcher::~HbsPatcher(){}

}//hbs
}//simple_moves
}//protocols
