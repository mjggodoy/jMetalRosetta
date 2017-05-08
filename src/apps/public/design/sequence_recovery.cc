// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file sequence_recovery.cc
/// @brief A protocol which outputs sequence recovery statistics ala the table in the "Native sequences are close to optimal" paper.
/// @author Ron Jacak
/// @author P. Douglas Renfrew (renfrew@nyu.edu) ( added rotamer recovery, cleanup )

// Unit headers
#include <devel/init.hh>

//project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


// basic headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <basic/options/util.hh>
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Option keys
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <sstream>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/vector0.hh>


static THREAD_LOCAL basic::Tracer TR( "sequence_recovery" );

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;

namespace sequence_recovery {
FileOptionKey const native_pdb_list( "sequence_recovery::native_pdb_list" );
FileOptionKey const redesign_pdb_list( "sequence_recovery::redesign_pdb_list" );
FileOptionKey const parse_taskops_file( "sequence_recovery::parse_taskops_file" );
BooleanOptionKey const rotamer_recovery( "sequence_recovery::rotamer_recovery" );
StringOptionKey const seq_recov_filename( "sequence_recovery::seq_recov_filename" );
StringOptionKey const sub_matrix_filename( "sequence_recovery::sub_matrix_filename" );
IntegerOptionKey const se_cutoff( "sequence_recovery::se_cutoff" );
}

std::string usage_string;

void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream << "No files given: Use either -file:s or -file:l to designate a single pdb or a list of pdbs.\n\n"
		<< "Usage: " << exe
		<< "\n\t-native_pdb_list <list file>"
		<< "\n\t-redesign_pdb_list <list file>"

		<< "\n\t[-seq_recov_filename <file>] file to output sequence recoveries to (default: sequencerecovery.txt)"
		<< "\n\t[-sub_matrix_filename <file>] file to output substitution matrix to (default: submatrix.txt)"
		<< "\n\t[-parse_taskops_file <file>] tagfile which contains task operations to apply before measuring recovery (optional)"
		<< "\n\t[-ignore_unrecognized_res]"

		<< "\n\n";

	usage_string = usage_stream.str();

}


/// @brief load custom TaskOperations according to an xml-like utility::tag file
/// @details The sequence recovery app can only handle taskops that do not use
/// ResidueSelectors, unless they are anonymous (i.e. unnamed) ResidueSelectors
/// that are declared as subtags of TaskOperations.
core::pack::task::TaskFactoryOP setup_tf( core::pack::task::TaskFactoryOP task_factory_ ) {

	using namespace core::pack::task::operation;

	if ( option[ sequence_recovery::parse_taskops_file ].user() ) {
		std::string tagfile_name( option[ sequence_recovery::parse_taskops_file ]() );
		TaskOperationFactory::TaskOperationOPs tops;
		basic::datacache::DataMap dm; // empty data map! any TaskOperation that tries to retrieve data out of this datamap will fail to find it.
		TaskOperationFactory::get_instance()->newTaskOperations( tops, dm, tagfile_name );
		for ( TaskOperationFactory::TaskOperationOPs::iterator it( tops.begin() ), itend( tops.end() ); it != itend; ++it ) {
			task_factory_->push_back( *it );
		}
	} else {
		task_factory_->push_back( TaskOperationCOP( new pack::task::operation::InitializeFromCommandline ) );
	}

	return task_factory_;

}


/// @brief helper method which uses the tenA nb graph in the pose object to fill a vector with nb counts
void fill_num_neighbors( pose::Pose & pose, utility::vector1< core::Size > & num_nbs ) {

	using core::conformation::PointGraph;
	using core::conformation::PointGraphOP;

	PointGraphOP pg( new PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, 10.0 /* ten angstrom distance */ ); // create edges

	num_nbs.resize( pose.size(), 0 );
	for ( core::Size ii=1; ii <= pose.size(); ++ii ) {

		// a PointGraph is a typedef of UpperEdgeGraph< PointGraphVertexData, PointGraphEdgeData >
		// so any of the method in UpperEdgeGraph should be avail. here. The UpperEdgeGraph provides access to nodes
		// via a get_vertex() method, and each vertex can report back how many nbs it has.
		// So something that before was really complicated (nb count calculation) is done in <10 lines of code.
		// the assumption we're making here is that a pose residue position ii is the same index as the point graph vertex
		// that is indeed the case if you look at what the function residue_point_graph_from_pose().
		num_nbs[ ii ] = pg->get_vertex(ii).num_neighbors_counting_self();
	}

	return;
}

/// @brief return the set of residues that are designable based given pose
std::set< Size > fill_designable_set( pose::Pose & pose, pack::task::TaskFactoryOP & tf ) {

	//we need to score the pose for many of the task operations passed from cmd line
	if ( option[ sequence_recovery::parse_taskops_file ].user() ) {
		scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
		(*scorefxn)( pose );
	}
	std::set< Size > designable_set;
	core::pack::task::PackerTaskOP design_task( tf->create_task_and_apply_taskoperations( pose ) );

#ifndef NDEBUG
	TR<< "Task for " << pose.pdb_info()->name() << " is: \n" << *(design_task)  << std::endl;
#endif

	// iterate over all residues
	for ( Size ii = 1; ii<= design_task->total_residue(); ++ii ) {
		if ( design_task->being_designed( ii ) ) {
			designable_set.insert( ii );
		}
	}

	return designable_set;

}


/// @brief iterates over all designed positions and determines identity to native. outputs recoveries to file.
void measure_sequence_recovery( utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses ) {

	// setup main arrays used for calculation
	utility::vector1< core::Size > n_correct( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_core( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_surface( chemical::num_canonical_aas, 0 );

	ObjexxFCL::FArray2D_int sub_matrix( chemical::num_canonical_aas, chemical::num_canonical_aas, 0 );

	Size n_correct_total(0); Size n_total(0);
	Size n_correct_total_core(0); Size n_total_core(0);
	Size n_correct_total_surface(0); Size n_total_surface(0);

	Size surface_exposed_cutoff = option[ sequence_recovery::se_cutoff ];
	Size core_cutoff = 24;

	// iterate through all the structures
	utility::vector1< core::pose::Pose >::iterator native_itr( native_poses.begin() ), native_last( native_poses.end() );
	utility::vector1< core::pose::Pose >::iterator redesign_itr( redesign_poses.begin() ), redesign_last( redesign_poses.end() );

	while ( ( native_itr != native_last ) && (redesign_itr != redesign_last ) ) {

		// get local copies of the poses
		core::pose::Pose native_pose( *native_itr );
		core::pose::Pose redesign_pose( *redesign_itr );

		// figure out the task & neighbor info
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		std::set< Size > design_set;
		utility::vector1< core::Size > num_neighbors;

		// setup what residues we are going to look at...
		setup_tf( task_factory );
		design_set = fill_designable_set( native_pose, task_factory );
		fill_num_neighbors( native_pose, num_neighbors );

		// record native sequence
		// native_sequence vector is sized for the WHOLE pose not just those being designed
		// it doesn't matter because we only iterate over the number of designed positions
		Size const nres( native_pose.size() );
		utility::vector1< chemical::AA > native_sequence( nres );

		// iterate over designable positions
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			if ( ! native_pose.residue(*it).is_protein() ) {
				native_sequence[ *it ] = chemical::aa_unk;
				continue;
			}
			//figure out info about the native pose
			native_sequence[ *it ] = native_pose.residue( *it ).aa();
			n_native[ native_pose.residue(*it).aa() ]++;

			//determine core/surface
			if ( num_neighbors[*it] >= core_cutoff ) {
				n_native_core[ native_pose.residue(*it).aa() ]++;
				n_total_core++;
			}

			if ( num_neighbors[*it] < surface_exposed_cutoff ) {
				n_native_surface[ native_pose.residue(*it).aa() ]++;
				n_total_surface++;
			}

		} // end finding native seq

		/// measure seq recov
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			// don't worry about recovery of non-protein residues
			if ( redesign_pose.residue( *it ).is_protein() ) {
				n_total++;

				// increment the designed count
				n_designed[ redesign_pose.residue(*it).aa() ]++;

				if ( num_neighbors[*it] >= core_cutoff ) { n_designed_core[ redesign_pose.residue(*it).aa() ]++; }
				if ( num_neighbors[*it] < surface_exposed_cutoff ) { n_designed_surface[ redesign_pose.residue(*it).aa() ]++; }

				// then check if it's the same
				if ( native_sequence[ *it ] == redesign_pose.residue(*it).aa() ) {
					n_correct[ redesign_pose.residue(*it).aa() ]++;

					if ( num_neighbors[*it] >= core_cutoff ) {
						n_correct_core[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_core++;
					}
					if ( num_neighbors[*it] < surface_exposed_cutoff ) {
						n_correct_surface[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_surface++;
					}
					n_correct_total++;
				}

				// set the substitution matrix for this go round
				sub_matrix( native_pose.residue(*it).aa(), redesign_pose.residue(*it).aa() )++;
			}

		} // end measure seq reovery

		// increment iterators
		native_itr++; redesign_itr++;
	}

	// open sequence recovery file stream
	utility::io::ozstream outputFile( option[ sequence_recovery::seq_recov_filename ].value() ) ;

	// write header
	outputFile << "Residue\tNo.correct core\tNo.native core\tNo.designed core\tNo.correct/ No.native core\tNo.correct/ No.designed core\t"
		<< "No.correct\tNo.native\tNo.designed\tNo.correct/ No.native\tNo.correct/ No.designed\t"
		<< "Residue\tNo.correct surface\tNo.native surface\tNo.designed surface\tNo.correct/ No.native\tNo.correct/ No.designed" << std::endl;

	// write AA data
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
			<< n_correct_core[ ii ] << "\t" << n_native_core[ ii ] << "\t" << n_designed_core[ ii ] << "\t";

		if ( n_native_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		//if ( n_designed_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";

		outputFile << n_correct[ ii ] << "\t" << n_native[ ii ] << "\t" << n_designed[ ii ] << "\t";
		if ( n_native[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		//if ( n_designed[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
			<< n_correct_surface[ ii ] << "\t" << n_native_surface[ ii ] << "\t" << n_designed_surface[ ii ] << "\t";

		if ( n_native_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		//if ( n_designed_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";

		outputFile << std::endl;
	}

	// write totals
	outputFile << "Total\t"
		<< n_correct_total_core << "\t" << n_total_core << "\t\t" << F(5,3, (float)n_correct_total_core/n_total_core ) << "\t\t"
		<< n_correct_total << "\t" << n_total << "\t\t" << F(5,3, (float)n_correct_total/n_total ) << "\t\tTotal\t"
		<< n_correct_total_surface << "\t" << n_total_surface << "\t\t" << F(5,3, (float)n_correct_total_surface/n_total_surface )
		<< std::endl;


	// output the sequence substitution file
	utility::io::ozstream matrixFile( option[ sequence_recovery::sub_matrix_filename ].value() ) ; //defaults to submatrix.txt

	// write the header
	matrixFile << "AA_TYPE" << "\t" ;
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		matrixFile << "nat_"<<chemical::name_from_aa( chemical::AA(ii) ) << "\t";
	}
	matrixFile<<std::endl;

	// now write the numbers
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) { //redesigns
		matrixFile << "sub_" << chemical::name_from_aa( chemical::AA(ii) );
		for ( Size jj = 1; jj <= chemical::num_canonical_aas; ++jj ) { //natives
			//std::cout<<"Native: "<< jj << " Sub: " << ii << "  Value: "<<sub_matrix( jj, ii ) << std::endl;
			matrixFile<< "\t" << sub_matrix( jj, ii );
		}
		matrixFile << std::endl;
	}

	//  ///output the sequence substitution file with percent of native recovered
	//  utility::io::ozstream matrixFileratio( "submatrix.ratio.txt" ) ; //allow naming later
	//  //write the header
	//  matrixFileratio << "AA_TYPE" << "\t" ;
	//  for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
	//    matrixFileratio << "nat_"<<chemical::name_from_aa( chemical::AA(ii) ) << "\t";
	//  }
	//  matrixFileratio<<std::endl;

	//  //now write the numbers
	//  for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) { //redesigns
	//    matrixFileratio << "sub_" << chemical::name_from_aa( chemical::AA(ii) );
	//    for ( Size jj = 1; jj <= chemical::num_canonical_aas; ++jj ) { //natives
	//      //std::cout<<"Native: "<< jj << " Sub: " << ii << "  Value: "<<sub_matrix( jj, ii ) << std::endl;
	//      matrixFileratio<< "\t"<< F(4,2, (float)sub_matrix( jj, ii )/sub_matrix( jj, jj ) );
	//    }
	//    matrixFileratio << std::endl;
	//  }

}


//@brief method which contains logic for calculating rotamer recovery. not implemented.
void measure_rotamer_recovery( utility::vector1<core::pose::Pose> & /*native_poses*/, utility::vector1<core::pose::Pose> & /*redesign_poses*/ ) {}


//@brief main method for the sequence recovery protocol
int main( int argc, char* argv[] ) {

	try {

		using utility::file::file_exists;
		using utility::file::FileName;

		option.add( sequence_recovery::native_pdb_list, "List of pdb files of the native structures." );
		option.add( sequence_recovery::redesign_pdb_list, "List of pdb files of the redesigned structures." );
		option.add( sequence_recovery::parse_taskops_file, "XML file which contains task operations to apply before measuring recovery (optional)" );
		option.add( sequence_recovery::rotamer_recovery, "Compare the rotamer recovery instead of sequence recovery." ).def( false );
		option.add( sequence_recovery::seq_recov_filename, "Name of file for sequence recovery output." ).def("sequencerecovery.txt");
		option.add( sequence_recovery::sub_matrix_filename, "Name of file substitution matrix output." ).def("submatrix.txt");
		option.add( sequence_recovery::se_cutoff, "Integer for how many nbs a residue must have less than or equal to to be considered surface exposed." ).def( 16 );


		devel::init( argc, argv );

		// changing this so that native_pdb_list and redesign_pdb_list do not have default values. giving these options can lead
		// to users measuring recovery against the wrong set of PDBs.
		if ( argc == 1 || !option[ sequence_recovery::native_pdb_list ].user() || !option[ sequence_recovery::redesign_pdb_list ].user() ) {
			init_usage_prompt( argv[0] );
			utility_exit_with_message_status( usage_string, 1 );
		}

		// read list file. open the file specified by the flag 'native_pdb_list' and read in all the lines in it
		std::vector< FileName > native_pdb_file_names;
		std::string native_pdb_list_file_name( option[ sequence_recovery::native_pdb_list ].value() );
		std::ifstream native_data( native_pdb_list_file_name.c_str() );
		std::string native_line;
		if ( !native_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + native_pdb_list_file_name + '\n' );
		}
		while ( getline( native_data, native_line ) ) {
			native_pdb_file_names.push_back( FileName( native_line ) );
		}

		native_data.close();

		// read list file. open the file specified by the flag 'redesign_pdb_list' and read in all the lines in it
		std::vector< FileName > redesign_pdb_file_names;
		std::string redesign_pdb_list_file_name( option[ sequence_recovery::redesign_pdb_list ].value() );
		std::ifstream redesign_data( redesign_pdb_list_file_name.c_str() );
		std::string redesign_line;
		if ( !redesign_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + redesign_pdb_list_file_name + '\n' );
		}
		while ( getline( redesign_data, redesign_line ) ) {
			redesign_pdb_file_names.push_back( FileName( redesign_line ) );
		}
		redesign_data.close();

		// check that the vectors are the same size. if not error out immediately.
		if ( native_pdb_file_names.size() != redesign_pdb_file_names.size() ) {
			utility_exit_with_message( "Size of native pdb list file: " + native_pdb_list_file_name + " does not equal size of redesign pdb list: " + redesign_pdb_list_file_name + "!\n" );
		}

		// iterate over both FileName vector and read in the PDB files
		utility::vector1< pose::Pose > native_poses;
		utility::vector1< pose::Pose > redesign_poses;

		std::vector< FileName >::iterator native_pdb( native_pdb_file_names.begin() ), native_last_pdb(native_pdb_file_names.end());
		std::vector< FileName >::iterator redesign_pdb( redesign_pdb_file_names.begin() ), redesign_last_pdb(redesign_pdb_file_names.end());

		while ( ( native_pdb != native_last_pdb ) && ( redesign_pdb != redesign_last_pdb ) ) {

			// check to make sure the file exists
			if ( !file_exists( *native_pdb ) ) {
				utility_exit_with_message( "Native pdb " + std::string(*native_pdb) + " not found! skipping" );
			}
			if ( !file_exists( *redesign_pdb ) ) {
				utility_exit_with_message( "Redesign pdb " + std::string(*redesign_pdb) + " not found! skipping" );
			}

			TR << "Reading in poses " << *native_pdb << " and " << *redesign_pdb << std::endl;
			core::pose::Pose native_pose, redesign_pose;
			core::import_pose::pose_from_file( native_pose, *native_pdb , core::import_pose::PDB_file);
			core::import_pose::pose_from_file( redesign_pose, *redesign_pdb , core::import_pose::PDB_file);

			native_poses.push_back( native_pose ); redesign_poses.push_back( redesign_pose );
			native_pdb++; redesign_pdb++;
		}


		if ( option[ sequence_recovery::rotamer_recovery ].value() ) {
			TR << "Measuring rotamer recovery"  << std::endl;
			measure_rotamer_recovery( native_poses, redesign_poses );
		} else {
			TR << "Measuring sequence recovery" << std::endl;
			measure_sequence_recovery( native_poses, redesign_poses );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

