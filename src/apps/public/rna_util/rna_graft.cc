// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

using namespace core;
using namespace basic;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;



typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Quick superimposer for RNA puzzle stuff. Assumes that poses have one chain each!
// Could easily generalize in the future. Could also add grafting functionality [which I
//  current do with a python script].
//   -- Rhiju, Oct 2012
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( IntegerVector, superimpose_res )
OPT_KEY( IntegerVector, graft_res )
OPT_KEY( IntegerVector, unvirtualize_phosphate_res )
OPT_KEY( Boolean, graft_backbone_only )

////////////////////////////////////////////////////////////////////////////
// by default all 5' phosphates are virtualized -- keep some of them in.
void
unvirtualize_phosphates( pose::Pose & pose, utility::vector1< Size > const & unvirtualize_phosphate_residues  ){
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !unvirtualize_phosphate_residues.has_value( pose.pdb_info()->number( n ) ) ) continue;
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, n );
	}
}


////////////////////////////////////////////////////////////////////////////
void
get_pose_and_numbering( std::string const & pdb_file, pose::Pose & pose, utility::vector1< Size > & resnum )
{

	using namespace core::chemical;
	using namespace core::pose;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	import_pose::pose_from_file( pose, *rsd_set,  pdb_file , core::import_pose::PDB_file);

	core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
	core::pose::rna::virtualize_5prime_phosphates( pose );

	resnum.clear();
	PDBInfoOP pdb_info = pose.pdb_info();
	for ( Size n = 1; n <= pose.size(); n++ )  resnum.push_back( pdb_info->number(n) );

}


////////////////////////////////////////////////////////////////////////////
Size
find_index( Size const val, utility::vector1< Size > const & vec ){
	for ( Size n = 1; n <= vec.size(); n++ ) {
		if ( val == vec[n] ) return n;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
bool
is_subset(  utility::vector1< Size > const & vec_sub,
	utility::vector1< Size  >const & vec_big ){
	for ( Size n = 1; n <= vec_sub.size(); n++ ) {
		if ( !find_index(vec_sub[n], vec_big) ) return false;
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////
std::map< Size, Size >
calculate_res_map( utility::vector1< Size > const & superimpose_resnum,
	utility::vector1< Size > const & resnum1,
	utility::vector1< Size > const & resnum2 ){

	runtime_assert( is_subset( superimpose_resnum, resnum1 ) );
	runtime_assert( is_subset( superimpose_resnum, resnum2 ) );

	std::map< Size, Size > res_map;

	for ( Size n = 1; n <= superimpose_resnum.size(); n++ ) {
		Size const i = find_index( superimpose_resnum[n],  resnum1 );
		runtime_assert( i > 0 );
		Size const j = find_index( superimpose_resnum[n],  resnum2 );
		runtime_assert( j > 0 );
		res_map[i] = j;
	}

	return res_map;
}


////////////////////////////////////////////////////////////////////////////
void
superimpose_pdb( pose::Pose & pose1 /* the 'parent pose'*/,
	pose::Pose & pose2 /*the one that got superimposed*/,
	utility::vector1< Size > const & resnum1,
	utility::vector1< Size > const & resnum2){

	using namespace core::pose;
	using namespace basic::options;

	utility::vector1< Size > superimpose_resnum;
	if ( option[superimpose_res].user() ) {
		// check if all specified residues are actually in the file.
		superimpose_resnum = option[ superimpose_res ]();
	} else {
		// find all shared residues
		for ( Size n = 1; n <= resnum1.size(); n++ ) {
			Size const n2 =  find_index( resnum1[n], resnum2 );
			if ( n2 == 0 ) continue;
			superimpose_resnum.push_back( resnum1[n] );

			if ( pose1.aa(n) != pose2.aa(n2) ) {
				std::cout << "Mismatch in sequence? " << std::endl;
				std::cout << " residue " << resnum1[n]  << " " << pose1.sequence()[n -1] <<  "  at pose position " << n << std::endl;
				std::cout << " residue " << resnum2[n2] << " " << pose2.sequence()[n2-1] <<  "  at pose position " << n2 << std::endl;
				runtime_assert( pose1.aa(n) ==  pose2.aa( n2 ) );
			}
		}
	}

	if ( superimpose_resnum.size() > 0 ) {

		std::map< Size, Size > superimpose_res_map = calculate_res_map( superimpose_resnum, resnum2, resnum1 );

		// now superimpose based on this subset. Copied from RNA_DeNovoProtocol.cc
		id::AtomID_Map< id::AtomID > const & alignment_atom_id_map_native =
			protocols::stepwise::modeler::align::create_alignment_id_map_legacy( pose2, pose1, superimpose_res_map ); // perhaps this should move to toolbox.
		core::scoring::superimpose_pose( pose2, pose1, alignment_atom_id_map_native );
	} else {
		std::cout << "WARNING: no overlap residues found or specified (by -superimpose_res)!!!! " << std::endl;
	}

	unvirtualize_phosphates( pose1, option[ unvirtualize_phosphate_res ]() );
	unvirtualize_phosphates( pose2, option[ unvirtualize_phosphate_res ]() );

}


////////////////////////////////////////////////////////////////////////////
void
output_superimposed_pdb( pose::Pose const & pose2, std::string const & pdb_file2 ){

	using namespace basic::options;

	std::string outfile = pdb_file2;

	std::string::size_type pos = pdb_file2.find( ".pdb", 0 );
	std::string const new_prefix = ".sup.pdb";
	if ( pos == std::string::npos ) {
		utility_exit_with_message(  "If you want to output a lores silent file, better use .out suffix ==> " + pdb_file2 );
	}
	outfile.replace( pos, new_prefix.length(), new_prefix );

	std::cout << "Outputting superimposed pdb to: " << outfile << std::endl;
	pose2.dump_pdb( outfile );

}


////////////////////////////////////////////////////////////////////////////
void
graft_in_positions( pose::Pose const & pose1,
	pose::Pose & pose_target,
	utility::vector1< Size > const & resnum1,
	utility::vector1< Size > const & resnum_target,
	utility::vector1< Size > const & resnum_graft,
	bool const graft_backbone_only){

	using namespace core::conformation;
	using namespace core::id;

	for ( Size i = 1; i <= pose1.size(); i++ ) {

		Size q = find_index( resnum1[i], resnum_target );
		if ( q == 0 ) continue;

		if ( resnum_graft.size() > 0 ) { //only copy over 'graft res'
			Size r = find_index( resnum1[i], resnum_graft );
			if ( r == 0 ) continue;
		}

		Residue const & rsd_i = pose1.residue(i);

		for ( Size ii = 1; ii <= rsd_i.natoms(); ii++ ) {

			if ( rsd_i.is_virtual( ii ) ) continue;
			if ( graft_backbone_only &&
					( ( ii >= rsd_i.first_sidechain_atom() && ii <= rsd_i.nheavyatoms() ) ||
					( ii >= rsd_i.first_sidechain_hydrogen() && ii <= rsd_i.natoms() ) ) ) continue;
			std::string const & atom_name = rsd_i.atom_name( ii ) ;

			if ( !pose_target.residue( q ).has( atom_name )  ) continue;
			Size const qq = pose_target.residue( q ).atom_index( atom_name );
			pose_target.set_xyz(  AtomID(qq, q), rsd_i.xyz( ii ) );

		}

	}

}

////////////////////////////////////////////////////////////////////////////
void
graft_in_positions( pose::Pose const & pose1,
	pose::Pose & pose_target,
	utility::vector1< Size > const & resnum1,
	utility::vector1< Size > const & resnum_target) {
	utility::vector1< Size > resnum_dummy;
	graft_in_positions( pose1, pose_target, resnum1, resnum_target, resnum_dummy, false );

}

////////////////////////////////////////////////////////////////////////////
void
close_loops( core::pose::Pose & pose,
	utility::vector1< Size > const & resnum ) {

	using namespace core::kinematics;

	utility::vector1< Size > cutpoints;
	for ( Size n = 2; n <= resnum.size(); n++ ) {
		if ( resnum[n] == resnum[n-1] + 1 ) {
			if ( ( pose.residue( n-1 ).xyz( " O3'" ) - pose.residue( n ).xyz( " P  " ) ).length() > 2.0 ) {
				FoldTree f = pose.fold_tree();
				f.new_jump( n-1, n, n-1 );
				f.set_jump_atoms( f.num_jump(),
					core::chemical::rna::chi1_torsion_atom( pose.residue_type( n-1 ) ),
					core::chemical::rna::chi1_torsion_atom( pose.residue_type( n   ) )   );
				pose.fold_tree( f );
				correctly_add_cutpoint_variants( pose, n-1 );
				cutpoints.push_back( n-1 );
			}
		}
	}

	protocols::farna::movers::RNA_LoopCloser rna_loop_closer;
	rna_loop_closer.apply( pose, cutpoints );
}

////////////////////////////////////////////////////////////////////////////
void
graft_pdb( pose::Pose const & pose1, pose::Pose const & pose2,
	utility::vector1< Size > const & resnum1,
	utility::vector1< Size > const & resnum2,
	pose::Pose & pose_target,
	utility::vector1< Size > & resnum_target ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace protocols::farna;

	pose_target.clear();
	resnum_target.clear();

	// make a generic RNA pose with the desired sequence. God I hope this works. Need to figure out PDBInfo crap too.
	// what's the sequence?

	utility::vector1< Size > graft_resnum;
	if ( option[ graft_res ].user() ) {
		graft_resnum = option[ graft_res ]();
	} else {
		graft_resnum = resnum2; // copy over everything in pose2.
	}

	std::list< std::pair< Size, char > > resnum_seq_list; // need to use a list, since we can sort it.

	for ( Size n = 1; n <= pose1.size(); n++ ) {
		if ( !find_index( resnum1[n], graft_resnum ) ) {
			resnum_seq_list.push_back( std::make_pair( resnum1[n], pose1.sequence()[ n-1 ] ) );
		}
	}

	for ( Size n = 1; n <= pose2.size(); n++ ) {
		if ( find_index( resnum2[n], graft_resnum ) ) {
			resnum_seq_list.push_back( std::make_pair( resnum2[n], pose2.sequence()[ n-1 ] ) );
		}
	}

	resnum_seq_list.sort();

	std::string sequence_target;
	for ( std::list< std::pair< Size, char > >::const_iterator iter = resnum_seq_list.begin(),
			end = resnum_seq_list.end(); iter != end; ++iter ) {
		resnum_target.push_back( iter->first );
		sequence_target += iter->second;
	}

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	make_pose_from_sequence( pose_target, sequence_target, *rsd_set );

	// Copy in non-virtual atoms from pose1;
	graft_in_positions( pose1, pose_target, resnum1, resnum_target );
	// Copy in non-virtual atoms from pose2;
	graft_in_positions( pose2, pose_target, resnum2, resnum_target, graft_resnum, option[ graft_backbone_only ]() );

	// close loops
	close_loops( pose_target, resnum_target );

	// kind of ad hoc. let's see if it works.
	core::pose::rna::figure_out_reasonable_rna_fold_tree( pose_target );
	core::pose::rna::virtualize_5prime_phosphates( pose_target );

	// chain boundaries -- make sure they are virtualized [this could be semi-dangerous for group I ribozymes with
	// discontinuous numbering -- so don't do it?
	for ( Size n = 1; n <= resnum_target.size(); n++ ) {
		if ( n == 1 || ( resnum_target[n] > resnum_target[n-1] + 1 ) ) {
			add_variant_type_to_pose_residue( pose_target, core::chemical::VIRTUAL_PHOSPHATE, n );
		}
	}

	// Need to fill in PDBInfo!
	PDBInfoOP pdb_info( new PDBInfo( pose_target ) );
	pdb_info->set_numbering( resnum_target );
	pose_target.pdb_info(  pdb_info );
	std::cout << pose_target.annotated_sequence() << std::endl;
	unvirtualize_phosphates( pose_target, option[ unvirtualize_phosphate_res ]() );
	std::cout << pose_target.annotated_sequence() << std::endl;

	std::string outfile = "graft.pdb";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();

	std::cout << "Outputting: " << outfile << std::endl;
	pose_target.dump_pdb( outfile );

}

////////////////////////////////////////////////////////////////////////////
void
rna_superimpose_and_graft_test(){

	using namespace basic::options;

	utility::vector1< std::string > pdb_files = option[ in::file::s ]();
	runtime_assert( pdb_files.size() >= 2 );

	pose::Pose pose1, pose2, pose_target;
	utility::vector1< Size > resnum1, resnum2, resnum_target;

	// note that following will also try to estimate a fold_tree & virtualize 5' phosphates
	// (which will be ignored in superposition).
	get_pose_and_numbering( pdb_files[ 1 ], pose1, resnum1 );

	for ( Size i = 2; i <= pdb_files.size(); i++ ) {
		get_pose_and_numbering( pdb_files[ i ], pose2, resnum2 );
		superimpose_pdb( pose1 /* the 'parent pose'*/,
			pose2 /*the one that got superimposed*/,
			resnum1,
			resnum2); // in principle could store resnum, pdb_file1 inside PDBInfo object!
		graft_pdb( pose1, pose2, resnum1, resnum2, pose_target, resnum_target );

		resnum1 = resnum_target;
		pose1   = pose_target;
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	rna_superimpose_and_graft_test();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;

		utility::vector1< Size > blank_size_vector;

		//Uh, options?
		NEW_OPT( superimpose_res, "residues over which to superimpose second pose onto first pose", blank_size_vector );
		NEW_OPT( graft_res, "residues to graft from second pose", blank_size_vector );
		NEW_OPT( unvirtualize_phosphate_res, "do not virtualize these phosphate (useful for cdiAMP)", blank_size_vector );
		NEW_OPT( graft_backbone_only, "only graft backbone from 2nd PDB", false );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
