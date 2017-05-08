// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseBetaAntiParallelJumpSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @details
/// @author Rhiju Das


//////////////////////////////////
// Unit headers
#include <protocols/stepwise/sampler/protein/util.hh>

// Package headers
#include <protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.hh>
#include <protocols/stepwise/sampler/protein/ProteinBetaAntiParallelStepWiseSampler.hh>
#include <protocols/stepwise/sampler/protein/ProteinFragmentStepWiseSampler.hh>
#include <protocols/stepwise/sampler/protein/ProteinMainChainStepWiseSampler.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClustererSilentBased.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/protein/StepWiseProteinBackboneSampler.hh>
#include <protocols/stepwise/setup/util.hh>

// Protocols headers
#include <protocols/jd2/util.hh>

// Core headers
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// numeric headers
#include <numeric/xyz.functions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.protein.util" );

using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {


StepWiseSamplerSizedOP
get_basic_protein_sampler(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_res_list,
	protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
	protocols::stepwise::modeler::options::StepWiseModelerOptionsCOP options,
	utility::vector1< protocols::stepwise::modeler::protein::InputStreamWithResidueInfoOP > & input_streams ) {

	using namespace protocols::stepwise::modeler::protein;
	using namespace protocols::stepwise::sampler::input_streams;

	StepWiseSamplerSizedOP sampler;
	if ( options->frag_files().size() > 0 ) {
		std::string const frag_file  = options->frag_files()[ 1 ];
		utility::vector1< Size > const & slice_res = working_parameters->working_res_list();
		sampler = StepWiseSamplerSizedOP( new ProteinFragmentStepWiseSampler( frag_file, slice_res, moving_res_list ) );
		if ( input_streams.size() == 1 ) {
			StepWiseSamplerSizedOP sampler_identity( new InputStreamStepWiseSampler( input_streams[1] ) );
			sampler = StepWiseSamplerSizedOP( new StepWiseSamplerSizedComb( sampler_identity /*outer*/, sampler /*inner, fragment generator above*/) );
		}
	} else if ( input_streams.size() == 2 ) {
		// assume that we want to "combine" two streams of poses...
		// This would be the mode if we have a bunch of templates from which we will graft chunks.
		// Or if we have SWA-based little fragments that we want to paste in.
		//   runtime_assert( stepwise_pose_setup_ != 0 );
		StepWiseSamplerSizedOP input_stream_sampler1( new InputStreamStepWiseSampler( input_streams[1] ) );
		StepWiseSamplerSizedOP input_stream_sampler2( new InputStreamStepWiseSampler( input_streams[2] ) );
		sampler = StepWiseSamplerSizedOP( new StepWiseSamplerSizedComb( input_stream_sampler1 /*outer*/, input_stream_sampler2 /*inner*/) );
	} else if ( options->sample_beta() ) {
		if ( moving_res_list.size() !=  1 ) utility_exit_with_message( "Sample beta only works for adding one residue to a beta sheet...");
		sampler = StepWiseSamplerSizedOP( new ProteinBetaAntiParallelStepWiseSampler( pose, moving_res_list[1] ) );
	} else if ( moving_res_list.size() > 0 ) {

		//////////////////////////////////////////////////////////////////////
		//  DEFAULT -- enumeration of conformations for moving residues.
		// Screen that predefines (phi, psi, omega)  for moving residues
		// --> input for later mover carries out all the green packer moves.
		//////////////////////////////////////////////////////////////////////
		StepWiseProteinBackboneSampler backbone_sampler( working_parameters );
		backbone_sampler.set_n_sample( options->n_sample() );
		backbone_sampler.set_native_pose( working_parameters->working_native_pose() );

		//when called in StepWiseMonteCarlo -- had a different convention in original StepWiseAssembly,
		// perhaps should be deprecated. If moving_res_list is 3-4, following samples psi/omega at 2 and phi at 5:
		backbone_sampler.set_expand_loop_takeoff( options->expand_loop_takeoff() );

		///////////////////////////////////////////////////////////////
		// Following has not been tested in a while, may not work:
		backbone_sampler.set_filter_native_big_bins( options->filter_native_big_bins()  );
		std::string const silent_file_centroid = core::io::silent::get_outfile_name_with_tag( options->silent_file(), "_centroid" );
		if ( options->centroid_output() ) backbone_sampler.set_silent_file( silent_file_centroid );
		if ( options->centroid_screen() ) {
			backbone_sampler.setup_centroid_screen( options->centroid_score_diff_cut(), options->centroid_weights(),
				options->nstruct_centroid(), options->ghost_loops() );
		}

		// Could also put loose chainbreak closure check here.
		Pose sampler_pose = pose;
		backbone_sampler.apply( sampler_pose );
		sampler = StepWiseSamplerSizedOP( new ProteinMainChainStepWiseSampler( backbone_sampler.which_torsions(),
			backbone_sampler.main_chain_torsion_set_lists_real(),
			options->choose_random() ) );
		TR << "Using ProteinMainChainStepWiseSampler. Num poses: " << backbone_sampler.main_chain_torsion_set_lists_real().size() << std::endl;
	} else {
		sampler = 0; // no op.
	}

	if ( sampler ) sampler->init();
	return sampler;
}

/////////////////////////////////////////////////////////////////////////
void
do_set_xyz( pose::Pose const & pose, Size const i, pose::Pose & scratch_pose, Size const i_scratch, core::kinematics::Stub const & stub ) {

	using namespace core::id;
	using namespace core::chemical;

	scratch_pose.set_xyz( NamedAtomID( " CA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " CA ", i) ) ) );
	scratch_pose.set_xyz( NamedAtomID( " N  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " N  ", i) ) ) );
	scratch_pose.set_xyz( NamedAtomID( " C  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " C  ", i) ) ) );
	scratch_pose.set_xyz( NamedAtomID( " O  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " O  ", i) ) ) );

	if ( pose.residue_type(i).has( " H  " ) ) {
		scratch_pose.set_xyz( NamedAtomID( " H  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " H  ", i) ) ) );
	} else if ( pose.residue_type(i).has( " H1 " ) ) {
		scratch_pose.set_xyz( NamedAtomID( " H  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " H1 ", i) ) ) );
	}

	if ( pose.residue_type(i).has( "1HA " ) ) {
		scratch_pose.set_xyz( NamedAtomID( "1HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( "1HA ", i) ) ) );
	} else if ( pose.residue_type(i).has( " HA " ) ) {
		scratch_pose.set_xyz( NamedAtomID( "1HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " HA ", i) ) ) );
	}

	if ( pose.residue_type(i).has( "2HA " ) ) {
		scratch_pose.set_xyz( NamedAtomID( "2HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( "2HA ", i) ) ) );
	} else if ( pose.residue_type(i).has( " CB " ) ) {
		scratch_pose.set_xyz( NamedAtomID( "2HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " CB ", i) ) ) );
	}
}


////////////////////////////////////////////////////////
// Move ths out of stepwise_protein_test.cc!!
////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void
generate_beta_database_test() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::pose;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io::silent;

	static ScoreFunctionOP scorefxn = get_score_function();
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Pose pose, scratch_pose;

	// little pose will be used to output beta geometries.
	std::string sequence = "GG";
	core::pose::make_pose_from_sequence( scratch_pose, sequence, *rsd_set );
	FoldTree f( scratch_pose.fold_tree() );

	f.new_jump(1,2,1);
	f.set_jump_atoms( 1, " CA ", " CA ", true /*bKeepStubInResidue. I think this is still protein-centric*/ );
	f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.
	scratch_pose.fold_tree( f );

	pose::remove_variant_type_from_pose_residue( scratch_pose, LOWER_TERMINUS_VARIANT, 1 );
	pose::remove_variant_type_from_pose_residue( scratch_pose, UPPER_TERMINUS_VARIANT, 2 );

	// main loop
	utility::vector1< utility::file::FileName > in_files = jd2::input_pdb_files_from_command_line();

	SilentFileOptions opts;
	SilentFileDataOP silent_file_data( new SilentFileData(opts) );
	std::string const silent_file( option[ out::file::silent ]() );

	Size count( 0 );
	for ( Size n = 1; n <= in_files.size(); n++ ) {

		std::cout << "-------- " << in_files[n] << "----------" << std::endl;
		import_pose::pose_from_file( pose, *rsd_set, in_files[ n ] , core::import_pose::PDB_file);

		// Look for beta pairings (antiparallel for now)
		(*scorefxn)(pose);
		HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
		hbond_options->use_hb_env_dep( false );
		HBondSetOP hbond_set( new hbonds::HBondSet( *hbond_options ) );

		fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

		std::map< std::pair< Size, Size >, bool > donor_to_acceptor_bb_hb;

		// In antiparallel HB, there are two hbonds between the residue backbones:
		//   C=O ... HN  and a NH ... O=C

		for ( Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
			hbonds::HBond const & hbond( hbond_set->hbond( i ) );

			Size const don_res_num = hbond.don_res();
			Size const acc_res_num = hbond.acc_res();

			if ( hbond.don_hatm_is_protein_backbone() && hbond.acc_atm_is_protein_backbone() ) {
				donor_to_acceptor_bb_hb[ std::make_pair( don_res_num, acc_res_num ) ] = true;
			}
		}

		for ( Size i = 1; i <= pose.size(); i++ ) {
			for ( Size j = 1; j <= pose.size(); j++ ) {
				if ( donor_to_acceptor_bb_hb[ std::make_pair( i, j ) ] &&
						donor_to_acceptor_bb_hb[ std::make_pair( j, i ) ] ) {
					//Create a coordinate system at residue i,
					Stub stub( pose.xyz( NamedAtomID( " CA ", i ) ),
						pose.xyz( NamedAtomID( " CA ", i ) ),
						pose.xyz( NamedAtomID( " N  ", i ) ),
						pose.xyz( NamedAtomID( " C  ", i ) ) );

					// spit out coordinates.
					do_set_xyz( pose, i, scratch_pose, 1, stub );
					do_set_xyz( pose, j, scratch_pose, 2, stub );

					count++;

					std::string const tag = "S_"+ObjexxFCL::string_of( count );
					(*scorefxn)( scratch_pose );
					BinarySilentStruct s( opts, scratch_pose, tag );

					silent_file_data->write_silent_struct( s, silent_file, false /*write score only*/ );
					silent_file_data->add_structure( s );

				}
			}
		}

	}

	/////////////////////////////////////////////////////////////////////
	// CLUSTER!
	/////////////////////////////////////////////////////////////////////
	// probably could convert following to StepWiseLegacyClusterer [pose based]
	protocols::stepwise::modeler::align::StepWiseLegacyClustererSilentBased stepwise_clusterer( silent_file_data );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )  max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	//  stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] ); // false by default
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
	stepwise_clusterer.set_score_diff_cut( 10.0 );

	// Do not superimpose inside the clusterer -- trust our superposition over residue 1.
	utility::vector1< Size > calc_rms_residues;
	calc_rms_residues.push_back( 1 );
	calc_rms_residues.push_back( 2 );
	stepwise_clusterer.set_calc_rms_res( calc_rms_residues );

	Real cluster_radius( 0.1 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius );

	stepwise_clusterer.cluster();

	std::string const silent_file_cluster   = core::io::silent::get_outfile_name_with_tag( silent_file, "_cluster" );
	stepwise_clusterer.output_silent_file( silent_file_cluster );

	/////////////////////////////////////////////////////////////////////
	// Output jump information
	/////////////////////////////////////////////////////////////////////
	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	using namespace protocols::stepwise::modeler;
	PoseList const & pose_list = stepwise_clusterer.clustered_pose_list();

	for ( auto const & elem : pose_list ) {
		pose = *(elem.second);
		Stub stub1_forward( pose.xyz( NamedAtomID( " N  ", 1) ),pose.xyz( NamedAtomID( " CA ", 1) ),pose.xyz( NamedAtomID( " C  ", 1) ));
		Stub stub1_reverse( pose.xyz( NamedAtomID( " C  ", 1) ),pose.xyz( NamedAtomID( " CA ", 1) ),pose.xyz( NamedAtomID( " N  ", 1) ));

		Stub stub2_forward( pose.xyz( NamedAtomID( " N  ", 2) ),pose.xyz( NamedAtomID( " CA ", 2) ),pose.xyz( NamedAtomID( " C  ", 2) ));
		Stub stub2_reverse( pose.xyz( NamedAtomID( " C  ", 2) ),pose.xyz( NamedAtomID( " CA ", 2) ),pose.xyz( NamedAtomID( " N  ", 2) ));

		out << "A N N " <<  Jump( stub1_forward, stub2_forward ) << std::endl;
		out << "A N C " <<  Jump( stub1_forward, stub2_reverse ) << std::endl;
		out << "A C N " <<  Jump( stub1_reverse, stub2_forward ) << std::endl;
		out << "A C C " <<  Jump( stub1_reverse, stub2_reverse ) << std::endl;
	}

	out.close();

	std::cout << "Put JUMP transforms in " << outfile << std::endl;

}

} //protein
} //sampler
} //stepwise
} //protocols
