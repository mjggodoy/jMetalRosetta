// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file per_residue_energies.cc
/// @brief simple application for printing out energies from a Pose.
/// @author James Thompson
/// @author Promoted to public by Rocco Moretti (rmoretti@u.washington.edu)

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

#include <string>

void all_pair_energies(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	utility::vector1< utility::vector1< core::Real > > & pairwise_energies
) {
	runtime_assert( pairwise_energies.size() == pose.size() );
	runtime_assert( pairwise_energies.front().size() == pose.size() );

	scorefxn->score(pose);
	utility::vector1< bool > exclude_mask( pose.size(), true );
	for ( unsigned ii = 1; ii <= pose.size(); ++ii ) {
		for ( unsigned jj = ii+1; jj <= pose.size(); ++jj ) {
			if ( scorefxn->are_they_neighbors( pose, ii, jj ) ) {
				utility::vector1< bool > mask = exclude_mask;
				//std::cout << "calculating pairwise energy for " << ii << "," << jj << std::endl;
				mask[ii] = false;
				mask[jj] = false;
				core::Real const total = scorefxn->get_sub_score( pose, mask );
				pairwise_energies[ii][jj] = total;
				pairwise_energies[jj][ii] = total;
			}
		}
	}
}


int
main( int argc, char* argv [] ) {
	try {

		devel::init( argc, argv );

		using core::Size;
		using core::Real;
		using core::pose::Pose;
		using utility::vector1;
		using namespace basic;
		using namespace core::io::silent;
		using namespace core::import_pose::pose_stream;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace ObjexxFCL;

		MetaPoseInputStream input = streams_from_cmd_line();
		core::chemical::ResidueTypeSetCOP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			option[ in::file::residue_type_set ]()
		);

		// set up ScoreFunction, make certain that hydrogen bonding energies
		// are kept in the EnergyGraph.
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		scorefxn->set_energy_method_options( *emopts );

		SilentFileOptions opts; // initialized from the command line
		SilentFileData sfd(opts);
		core::pose::Pose current_pose;
		while ( input.has_another_pose() ) {
			input.fill_pose( current_pose, *rsd_set );

			if ( ! current_pose.pdb_info() ) {
				// Make a default PDB info if one doesn't already exist - this simplifies logic later.
				core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(current_pose) );
				current_pose.pdb_info( new_pdb_info );
			}

			(*scorefxn)(current_pose);
			EnergyMap weights( current_pose.energies().weights() );

			if ( option[ james::debug ]() ) {
				using utility::vector1;
				vector1< vector1< core::Real > > pair_energies( current_pose.size(),
					vector1< Real > (current_pose.size(), 0.0
					) );
				all_pair_energies(current_pose,scorefxn,pair_energies);

				for ( unsigned ii = 1; ii <= current_pose.size(); ++ii ) {
					for ( unsigned jj = ii+1; jj <= current_pose.size(); ++jj ) {
						if ( pair_energies[ii][jj] > 0 ) {
							SilentStructOP ss( new ScoreFileSilentStruct(opts) );
							ss->decoy_tag( "residue_" + string_of(ii) + "_" + string_of(jj) );
							ss->add_string_value( "pose_id", core::pose::tag_from_pose(current_pose) );
							ss->add_energy( "score", pair_energies[ii][jj] );
							sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
						}
					}
				}
			} else {
				for ( Size jj = 1; jj <= current_pose.size(); ++jj ) {
					EnergyMap rsd_energies(
						weights * current_pose.energies().residue_total_energies(jj)
					);

					SilentStructOP ss( new ScoreFileSilentStruct(opts) );
					ss->decoy_tag( "residue_" + string_of(jj) );
					ss->add_string_value( "pose_id", core::pose::tag_from_pose(current_pose) );
					ss->add_string_value( "pdb_id", string_of( current_pose.pdb_info()->number(jj) ) + current_pose.pdb_info()->chain(jj) );
					Real total(0);
					for ( int ii = 1; ii <= n_score_types; ++ii ) {
						if ( weights[ ScoreType(ii) ] != 0.0 ) {
							Real const value( rsd_energies[ ScoreType(ii) ] );
							std::string const scorename( name_from_score_type( ScoreType(ii) ) );
							total += value;
							ss->add_energy( scorename, value );
						}
					} // for n_score_types
					ss->add_energy( "score", total );
					sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
				} // for current_pose.size()
			}
		} // while ( input.has_another_pose() )

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
