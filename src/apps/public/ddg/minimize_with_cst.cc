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
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <devel/init.hh>

#include <numeric/xyzVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <cstdlib>
#include <string>


#include <core/chemical/ResidueTypeSet.fwd.hh>

//#include <core/scoring/ScoringManager.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>
//#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/Jump.hh>

#include <core/pose/Pose.fwd.hh>


#include <basic/options/option.hh>

//#include <core/util/basic.hh>
//#include <core/io/database/open.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>


//protocols
//#include <protocols/looprelax/looprelax_main.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>

//new options stuff?
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
// C++ headers

//C++ filechek
#include <sys/stat.h>

#include <utility/vector0.hh>


using namespace core;
using namespace std;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace ddg;
using namespace id;

ScoreFunction&
reduce_fa_rep(float fraction_fa_rep, ScoreFunction & s){
	s.set_weight( core::scoring::score_type_from_name("fa_rep"),
		s.get_weight(core::scoring::score_type_from_name("fa_rep"))*fraction_fa_rep);
	return s;
}

void
setup_ca_constraints(pose::Pose & pose, ScoreFunction & s, float const CA_cutoff, float const cst_tol){
	//create constraints for all residues
	//type: HARMONIC
	//static float const CA_cutoff(9.0);
	if ( !basic::options::option[OptionKeys::ddg::sc_min_only]() ) {
		int nres = pose.size();
		if ( basic::options::option[basic::options::OptionKeys::constraints::cst_file].user() ) {
			core::scoring::constraints::ConstraintSetOP cstset( new core::scoring::constraints::ConstraintSet() );
			cstset = core::scoring::constraints::ConstraintIO::read_constraints(option[basic::options::OptionKeys::constraints::cst_file][1],cstset,pose);
			pose.constraint_set(cstset);
		} else {
			for ( int i = 1; i <= nres; i++ ) {
				core::Size atom_i;
				if ( pose.residue(i).has( " CA " ) ) {
					atom_i = pose.residue(i).atom_index(" CA ");
				} else {
					atom_i = pose.residue(i).nbr_atom();
				}
				Vector const CA_i( pose.residue(i).xyz(atom_i) );
				for ( int j = 1; j < i; j++ ) {
					core::Size atom_j;
					if ( pose.residue(j).has( " CA " ) ) {
						atom_j = pose.residue(j).atom_index(" CA ");
					} else {
						atom_j = pose.residue(j).nbr_atom();
					}
					Vector const CA_j( pose.residue(j).xyz(atom_j) );
					Real const CA_dist = (CA_i - CA_j).length();
					if ( CA_dist < CA_cutoff ) {
						std::cout << "c-alpha constraints added to residues " << i << " and " << j << " dist " << CA_dist << " and tol " << cst_tol << std::endl;
						ConstraintCOP cst( ConstraintOP( new AtomPairConstraint( AtomID(atom_i,i),AtomID(atom_j,j),core::scoring::func::FuncOP(new core::scoring::func::HarmonicFunc(CA_dist, cst_tol))) ) );
						pose.add_constraint(cst);
					}
				}
			}
		}

		s.set_weight(atom_pair_constraint, basic::options::option[OptionKeys::ddg::constraint_weight]());
	}
}


void
minimize_with_constraints(pose::Pose & p, ScoreFunction & s,std::string output_tag){
	using namespace basic::options;

	std::string out_pdb_prefix = basic::options::option[OptionKeys::ddg::out_pdb_prefix ]();
	std::string output_directory = basic::options::option[OptionKeys::ddg::output_dir]();

	//silent file capabilities
	bool write_silent_file = false;
	core::io::silent::SilentFileOptions opts; // initialized from the command line
	core::io::silent::SilentFileData sfd(opts);
	std::string silentfilename="";
	if ( basic::options::option[out::file::silent].user() ) {
		sfd.set_filename(basic::options::option[out::file::silent]());
		silentfilename=basic::options::option[out::file::silent]();
		write_silent_file = true;
	}

	struct stat stFileInfo;
	int file_stat = stat((output_directory + "/" + out_pdb_prefix+"."+output_tag+"_0001.pdb").c_str(),
		&stFileInfo);
	if ( file_stat != 0 ) {
		core::optimization::AtomTreeMinimizer min_struc;
		float minimizer_tol = 0.000001;

		core::kinematics::MoveMap mm;
		mm.set_bb(true);
		mm.set_chi(true);
		if ( basic::options::option[OptionKeys::ddg::sc_min_only]() ) {
			minimizer_tol = 0.0001;
			mm.set_bb(false);
		}
		s.show(std::cout, p);
		// This used to be higher, but then it couldn't respond to the adjustment
		// for in case we are only doing sidechain minimization!
		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", minimizer_tol, true /*use_nb_list*/,
			false /*deriv_check_in*/, true /*deriv_check_verbose_in*/);

		options.nblist_auto_update( true );
		options.max_iter(5000); //otherwise, they don't seem to converge


		if ( basic::options::option[OptionKeys::ddg::ramp_repulsive]() ) {
			//set scorefxn fa_rep to 1/10 of original weight and then minimize
			ScoreFunctionOP one_tenth_orig(s.clone());
			reduce_fa_rep(0.1,*one_tenth_orig);
			//min_struc.run(p,mm,s,options);
			min_struc.run(p,mm,*one_tenth_orig,options);
			std::cout << "one tenth repulsive fa_rep score-function" << std::endl;
			one_tenth_orig->show(std::cout, p);

			//then set scorefxn fa_rep to 1/3 of original weight and then minimize
			ScoreFunctionOP one_third_orig(s.clone());
			reduce_fa_rep(0.33,*one_third_orig);
			min_struc.run(p,mm,*one_third_orig,options);
			std::cout << "one third repulsive fa_rep score-function" << std::endl;
			one_third_orig->show(std::cout, p);
			//then set scorefxn fa_rep to original weight and then minimize
		}
		pose::Pose before(p);

		min_struc.run(p,mm,s,options);

		if ( basic::options::option[OptionKeys::ddg::initial_repack]() ) {
			//repack
			std::cout << "repacking" << std::endl;
			ScoreFunctionOP spack = ScoreFunctionFactory::create_score_function("soft_rep_design.wts");
			pack::task::PackerTaskOP repack(pack::task::TaskFactory::create_packer_task(p));
			repack->restrict_to_repacking();
			for ( unsigned int i = 1; i <= p.size(); i++ ) {
				repack->nonconst_residue_task(i).or_include_current(true);
				repack->nonconst_residue_task(i).or_ex1(true);
				repack->nonconst_residue_task(i).or_ex2(true);
			}
			pack::pack_rotamers(p,(*spack),repack);
		}

		s.show(std::cout, p);
		while ( std::abs(s(p)-s(before)) > 1 ) { //make sure furhter minimizations lead to same answer
			std::cout << "running another iteration of minimization. difference is: " << (s(p)-s(before)) << std::endl;
			before = p;
			min_struc.run(p,mm,s,options);
		}

		if ( write_silent_file ) {
			core::io::silent::SilentFileOptions opts; // initialized from the command line
			core::io::silent::BinarySilentStruct ss(opts,p,(out_pdb_prefix+"."+output_tag+"_0001.pdb"));
			sfd.write_silent_struct(ss,silentfilename,false);
		} else {
			p.dump_pdb(output_directory + "/" + out_pdb_prefix+"."+output_tag+"_0001.pdb");
		}
	} else {
		std::cout << "file " << output_directory + "/" + out_pdb_prefix+"."+output_tag+"_0001.pdb" << " already exists! skipping" << std::endl;
	}
}

void
optimize_pose(pose::Pose & p, ScoreFunctionOP scorefxn,std::string output_tag){
	Real cst_tol;
	if ( basic::options::option[OptionKeys::ddg::harmonic_ca_tether].user() ) {
		cst_tol = basic::options::option[ OptionKeys::ddg::harmonic_ca_tether ]();
	} else {
		cst_tol = 0.5;
	}
	if ( !basic::options::option[OptionKeys::ddg::no_constraints].user() || //flag not present
			!basic::options::option[OptionKeys::ddg::no_constraints]() ) { //flag set to false
		setup_ca_constraints(p,(*scorefxn),9.0,cst_tol);
	}
	minimize_with_constraints(p, (*scorefxn),output_tag);
}

bool
already_minimized(std::string query,utility::vector1<std::string> check){
	//std::cout << "query is " << query << std::endl;
	for ( unsigned int i=1; i <= check.size(); i++ ) {
		//std::cout << "DEBUG: tag" << check[i] << std::endl;
		if ( check[i].compare(query) == 0 ) {
			//std::cout << "query " << query << " totall matches " << check[i] << std::endl;
			return true;
		}
	}
	return false;
}

int
main( int argc, char* argv [] )
{
	try {

		using namespace core;
		using namespace core::pose;
		using namespace utility;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::io::silent;
		// options, random initialization. MAKE SURE THIS COMES BEFORE OPTIONS
		devel::init( argc, argv );

		core::chemical::ResidueTypeSetCOP rsd_set;
		if ( option[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}

		// io::pdb::pose_from_file( pose, options::start_file() , core::import_pose::PDB_file); // gets filename from -s option

		ScoreFunctionOP scorefxn = get_score_function();

		/* if(basic::options::option[score::patch].user()){
		scorefxn->apply_patch_from_file(basic::options::option[score::patch]);
		}*/
		bool read_from_silent = false;
		vector1<std::string> files;
		// utility::vector1<std::string> tags;
		utility::vector1<file::FileName> list;
		if ( basic::options::option[in::file::s].user() ) {
			std::cout << "using -s option" << std::endl;
			list = option[in::file::s]();
		} else if ( option[in::file::l].user() ) {
			std::cout << "using -l option " << std::endl;
			list = basic::options::option[ in::file::l ]();
		} else if ( option[in::file::silent].user() ) {
			read_from_silent = true;
			//list = basic::options::option[in::file::silent]();
		}
		for ( unsigned int h=1; h<=list.size(); h++ ) {
			utility::io::izstream pdbs(list[h]);
			std::string fname;
			while ( pdbs >> fname ) {
				files.push_back(fname);
			}
		}
		if ( read_from_silent ) {
			files.push_back(*(option[in::file::silent]().begin()));
		}

		utility::vector1<std::string> tags_done;
		//if outputting silent files as well:
		if ( basic::options::option[out::file::silent].user() ) {
			std::string out = basic::options::option[out::file::silent]();
			SilentFileOptions opts; // initialized from the command line
			SilentFileData out_sfd(out,false,false,"binary",opts);
			if ( utility::file::file_exists(out) ) {
				out_sfd.read_file(out);
				tags_done = out_sfd.tags();
			}
		}
		std::string prefix = basic::options::option[OptionKeys::ddg::out_pdb_prefix ]();
		for ( unsigned int f=1; f<=files.size(); f++ ) {
			pose::Pose pose;
			std::cout << "examining file: " << files[f] << std::endl;
			if ( read_from_silent ) {
				SilentFileOptions opts; // initialized from the command line
				SilentFileData in_sfd(files[f],false,false,"binary",opts);
				in_sfd.read_file(files[f]);
				utility::vector1<std::string> tags = in_sfd.tags();
				for ( unsigned int t=1; t <= tags.size(); t++ ) {
					if ( !already_minimized(prefix+"."+tags[t]+"_0001",tags_done) ) {
						SilentStructOP ss = in_sfd[tags[t]];
						ss->fill_pose(pose,*rsd_set);
						scorefxn->show(std::cout, pose);
						std::cout << "tag assigned to pose: " << ss->decoy_tag() << std::endl;
						optimize_pose(pose,scorefxn,ss->decoy_tag());
					} else {
						std::cout << "tag " << tags[t] << " exists in outfile. skipping!" << std::endl;
					}
				}
			} else {
				core::import_pose::pose_from_file(pose, files[f], core::import_pose::PDB_file);
				std::string output = pose.pdb_info()->name();
				std::string pdb_prefix( utility::string_split( utility::string_split( output, '/' ).back(), '.' ).front() );
				optimize_pose(pose, scorefxn, pdb_prefix);
				//create constraints for all residues
				//type: HARMONIC
				//then minimize
			}
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
