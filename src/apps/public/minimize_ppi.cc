// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk, ragul, skelow, bazzoli (bazzoli@ku.edu)

// Project Headers
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

//Auto Headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

using namespace core;
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace protocols::simple_moves;
using namespace protocols::rigid;
using namespace core::chemical;
using namespace core::conformation;

OPT_KEY( Real, cst_force_constant )
//OPT_KEY( Real, coord_cst_weight )
OPT_KEY( Boolean, print_init )
OPT_KEY( Boolean, print_unbound )
OPT_KEY( Boolean, print_complex )
OPT_KEY( Boolean, iface_rmsd )
OPT_KEY( String, ref_decoy )
OPT_KEY( Boolean, score_only )

static basic::Tracer TR( "apps.pilot.minimize_ppi" );

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue
core::Size lig_res_num;

void define_interface( core::pose::Pose & ref_pose ) {
	lig_res_num =0;
	for ( int j = 1, resnum = ref_pose.size(); j <= resnum; ++j ) {
		if ( !ref_pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		TR << "No ligand given in reference PDB structure.  Cannot identify interface."<<std::endl;
		exit (1);
	}

	TR <<"sele ";

	EnergyGraph & energy_graph(ref_pose.energies().energy_graph());
	for ( utility::graph::Graph::EdgeListIter
			iru  = energy_graph.get_node( lig_res_num )->lower_edge_list_begin(),
			irue = energy_graph.get_node( lig_res_num )->lower_edge_list_end();
			iru != irue; ++iru ) {
		EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
		Size const j( edge->get_first_node_ind() );

		// the pair energies cached in the link
		EnergyMap const & emap( edge->fill_energy_map());
		Real const attr( emap[ fa_atr ] );
		//TR<<"\n"<<j<<": "<<attr<<"\n";
		if ( attr < -.2 ) {
			// create string id to store in set
			std::ostringstream residuestream;
			TR << "resi "<< ref_pose.pdb_info()->number(j)<<" or ";
			residuestream << ref_pose.pdb_info()->chain(j) << ref_pose.pdb_info()->number(j);
			std::string res_id = residuestream.str();
			interface.insert(res_id);
		}
	}

	TR << std::endl;
	TR << lig_res_num<< std::endl;

	//      ref_pose.delete_polymer_residue(lig_res_num);
}

bool
is_interface_heavyatom(
	core::pose::Pose const & pose,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	// ws get residue "key" for set
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();
	core::conformation::Residue const & rsd = pose.residue(resno);
	if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);
	return false;
}

Real
calpha_pdb_superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
	for ( Size ii = 1; ii <= mod_pose.size(); ++ii ) {
		if ( ! mod_pose.residue(ii).has("CA") ) continue;
		if ( ! mod_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= ref_pose.size(); ++jj ) {
			if ( ! ref_pose.residue(jj).has("CA") ) continue;
			if ( ! ref_pose.residue(jj).is_protein() ) continue;
			if ( mod_pose.pdb_info()->chain(ii) != ref_pose.pdb_info()->chain(jj) ) continue;
			if ( mod_pose.pdb_info()->number(ii) != ref_pose.pdb_info()->number(jj) ) continue;
			id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
			id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
			atom_map.set( id1, id2 );
			break;
		}
	}
	return superimpose_pose( mod_pose, ref_pose, atom_map );
}

Real
interface_rmsd(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( Size ii = 1; ii <= ref_pose.size(); ++ii ) {
		if ( ! ref_pose.residue(ii).has("CA") ) continue;
		if ( ! ref_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= mod_pose.size(); ++jj ) {
			if ( ! ref_pose.residue(ii).has("CA") ) continue;
			if ( ! ref_pose.residue(ii).is_protein() ) continue;
			if ( mod_pose.pdb_info()->chain(jj) != ref_pose.pdb_info()->chain(ii) ) continue;
			if ( mod_pose.pdb_info()->number(jj) != ref_pose.pdb_info()->number(ii) ) continue;
			Size num_atoms ( ref_pose.residue(ii).natoms() );
			for ( core::Size i = 1; i <= num_atoms; ++i ) {
				if ( is_interface_heavyatom ( ref_pose, mod_pose, ii, i) ) {
					Size num_atoms2 ( mod_pose.residue(jj).natoms() );
					for ( core::Size j = 1; j <= num_atoms2; ++j ) {
						if ( !ref_pose.residue(ii).atom_name(i).compare(mod_pose.residue(jj).atom_name(j)) ) {
							p1_coords.push_back(ref_pose.residue(ii).xyz(i));
							p2_coords.push_back(mod_pose.residue(jj).xyz(j));
						}
					}
				}
			}
		}
	}
	assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );

}


Real ligand_rmsd (pose::Pose const & pose1, pose::Pose const & pose2) {

	core::Size ref_res_num = 0;
	for ( int j = 1, resnum = pose1.size(); j <= resnum; ++j ) {
		if ( !pose1.residue(j).is_protein() ) {
			ref_res_num = j;
			break;
		}
	}
	if ( ref_res_num == 0 ) {
		std::cout<<"Error, no reference ligand for RMSD calculation" << std::endl;
		exit(1);
	}
	core::conformation::Residue const & pose1_rsd = pose1.conformation().residue(ref_res_num);

	core::Size inp_res_num = 0;
	for ( int j = 1, resnum = pose2.size(); j <= resnum; ++j ) {
		if ( !pose2.residue(j).is_protein() ) {
			inp_res_num = j;
			break;
		}
	}
	if ( inp_res_num == 0 ) {
		std::cout<<"Error, no input ligand for RMSD calculation" << std::endl;
		exit(1);
	}
	core::conformation::Residue const & pose2_rsd = pose2.conformation().residue(inp_res_num);


	//CALCULATE RMSD
	//IMPORTANT:this rmsd calculation does not consider symmetry
	//assert input & reference ligand should have same num of atoms
	if ( !(pose2_rsd.nheavyatoms() == pose1_rsd.nheavyatoms()) ) {
		std::cout<<"Error, input & reference ligand should have same no. of atoms in the same order" << std::endl;
		exit(1);
	} else {
		core::Real rmsd, dist_sum(0.);
		Size j = 1;
		for ( Size i = 1; i <= pose1_rsd.nheavyatoms(); ++i ) {
			if ( j <= pose2_rsd.nheavyatoms() ) {
				assert ( pose1_rsd.atom_name(i) == pose2_rsd.atom_name(j) );
				core::Real x_dist =  ( (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) * (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) );
				core::Real y_dist =  ( (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) * (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) );
				core::Real z_dist =  ( (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) * (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) );
				dist_sum += x_dist + y_dist + z_dist;
			}
			++j;
		}
		rmsd = sqrt(dist_sum/pose1_rsd.nheavyatoms());
		return rmsd;
	}

}


int main( int argc, char * argv [] ){
	try{
		NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.0000001 );
		//NEW_OPT( coord_cst_weight, "coordinate constraint weight", 1 );
		NEW_OPT( print_init, "print the initial complex for debugging", false );
		NEW_OPT( print_unbound, "print the mimized protein for debugging", false );
		NEW_OPT( print_complex, "print the minimized complex", true );
		NEW_OPT( iface_rmsd, "calculate the interface rmsd", false );
		NEW_OPT( ref_decoy, "the structure to compute RMSD and relative score to", "" );
		NEW_OPT( score_only, "compute all scores for the iput complex without minimization", false );

		devel::init(argc, argv);

		//setup scorefxn
		scoring::ScoreFunctionOP scorefxn(get_score_function());
		// scoring::ScoreFunctionOP repack_scorefxn = getScoreFunction();

		std::string const ref_decoy_fname = option[ ref_decoy ];
		// create pose from pdb
		pose::Pose ref_pose;
		std::string outfname;
		utility::io::ozstream outstream;
		if ( option[ iface_rmsd ] ) {
			core::import_pose::pose_from_file( ref_pose, ref_decoy_fname , core::import_pose::PDB_file);
			define_interface( ref_pose );
			TR << "Defined interface" << std::endl;
			if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
				outfname = "minrmsd." + option[ OptionKeys::out::output_tag ]() + ".out";
			} else {
				outfname = "minrmsd.out";
			}
			//std::cout<<outfname<<" output_tag: "<<option[ OptionKeys::out::output_tag ]()<<std::endl;
			outstream.open(outfname, std::ios::out);
		}

		//Register calculators
		std::string sasa_calc_name = "sasa";
		std::string hbond_calc_name = "hbond";
		std::string packstat_calc_name = "packstat";
		std::string burunsat_calc_name = "burunsat";
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator(new core::pose::metrics::simple_calculators::SasaCalculatorLegacy);
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );
		core::pose::metrics::PoseMetricCalculatorOP hb_calc(new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator());
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );
		core::pose::metrics::PoseMetricCalculatorOP packstat_calc(new protocols::toolbox::pose_metric_calculators::PackstatCalculator());
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calc );
		core::pose::metrics::PoseMetricCalculatorOP burunsat_calc(new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name));
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );

		//Description of output scores
		std::cout << "Description of output scores"<< std::endl;
		std::cout << "Interface_Scores:TAG" <<"  "<< "input_pdb_name" <<"  " << "bound_energy" <<" " << "Interface_Energy" <<" "<< "Total_BSA" <<" "<< "Interface_HB" <<"  "<< "Total_packstats" <<" "<< "Interface_unsat" << std::endl;

		for ( core::Size f=1; f <= basic::options::start_files().size(); f++ ) {
			//setup the bound pose
			pose::Pose bound_pose;
			std::string const input_pdb_name = basic::options::start_files().at(f);
			core::import_pose::pose_from_file( bound_pose, input_pdb_name , core::import_pose::PDB_file);
			pose::Pose pre_min_darc_pose = bound_pose;
			//create tag for output filename
			int pfounddir = input_pdb_name.find_last_of("/\\");
			int pfounddot = input_pdb_name.find_last_of(".");
			std::string tag = input_pdb_name.substr((pfounddir+1),(pfounddot-(pfounddir+1)));
			std::string init_pdb = "init_" + tag + ".pdb";
			std::string mini_pdb = "mini_" + tag + ".pdb";
			std::string unbo_pdb = "unbo_" + tag + ".pdb";

			int nres = bound_pose.size();
			Real coord_sdev( option[ OptionKeys::cst_force_constant ] );
			//take reciprocal and sqrt to pass as force constant
			coord_sdev = sqrt(1/coord_sdev);
			//std::cout<<" coord sdev "<< coord_sdev <<std::endl;
			Real cst_weight = 1;


			if ( !option[ score_only ] ) {

				//Initial score
				(*scorefxn)(bound_pose);
				TR << "Initial score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;
				if ( option [ print_init ] ) {
					bound_pose.dump_scored_pdb( init_pdb, *scorefxn );
				}

				ConstraintSetOP cst_set( new ConstraintSet() );
				core::scoring::func::HarmonicFuncOP spring(
					new core::scoring::func::HarmonicFunc( 0 /*mean*/, coord_sdev /*std-dev*/));
				conformation::Conformation const & conformation( bound_pose.conformation() );

				// jk we need an anchor in order to use CoordinateConstraint !!!
				Size const my_anchor = 1;
				core::kinematics::FoldTree fold_tree=bound_pose.fold_tree();
				core::kinematics::FoldTree rerooted_fold_tree = fold_tree;
				rerooted_fold_tree.reorder( my_anchor );
				bound_pose.fold_tree( rerooted_fold_tree);

				for ( int i=1; i <= nres; i++ ) {
					//Residue const  & reside = pose.residue( i );
					Residue const & nat_i_rsd( bound_pose.residue(i) );
					for ( Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {

						AtomID CAi ( ii, i );
						ConstraintOP cst(new CoordinateConstraint(
							CAi, AtomID(1,my_anchor), conformation.xyz( CAi ), spring ));

						cst_set->add_constraint(cst);
					}
				}
				bound_pose.constraint_set( cst_set );
				scorefxn->set_weight( coordinate_constraint, cst_weight );

				TR << "Starting minimization...." << std::endl;
				//AtomTreeMinimizer minimizer;
				AtomTreeMinimizer minimizer;
				MinimizerOptions min_options( "dfpmin", 0.00001, true, false );
				kinematics::MoveMap mm_all;
				mm_all.set_chi( true );
				mm_all.set_bb( true );
				mm_all.set_jump( true );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				(*scorefxn)(bound_pose);
				TR << "Post minimization 1 constrained score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;
				//bound_pose.dump_scored_pdb( "1.pdb", *scorefxn );

				/* // Setup packer task for repacking
				pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( bound_pose ));
				base_packer_task->set_bump_check( false );
				base_packer_task->initialize_from_command_line();
				base_packer_task->or_include_current( true );
				for ( Size ii = 1; ii <= bound_pose.size(); ++ii ) {
				base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
				}
				// First repack
				pack::pack_rotamers( bound_pose, *scorefxn, base_packer_task );
				// Report Scores
				(*scorefxn)(bound_pose);
				//bound_pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
				TR << "Score after repacking once: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

				// iterate over minimizing and repacking
				for ( Size iter = 1; iter <= 5; ++iter ) {
				cst_weight = cst_weight/2;
				scorefxn->set_weight( coordinate_constraint, cst_weight );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				(*scorefxn)(bound_pose);
				TR << "Current score after minimizing: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
				pack::pack_rotamers( bound_pose, *scorefxn, base_packer_task );
				(*scorefxn)(bound_pose);
				TR << "Current score after repacking: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
				}
				//bound_pose.dump_scored_pdb( "2.pdb", *scorefxn ); */

				// final minimization
				bound_pose.remove_constraints();
				(*scorefxn)(bound_pose);
				if ( option[ print_complex ] ) {
					//align minimized pose to the original docked pose and dump pdb complex and ligand
					protocols::simple_moves::SuperimposeMoverOP sp_mover(new protocols::simple_moves::SuperimposeMover());
					sp_mover->set_reference_pose( pre_min_darc_pose, 1, (pre_min_darc_pose.size()-1) );
					sp_mover->set_target_range( 1, (bound_pose.size()-1) );
					sp_mover->apply( bound_pose );
					bound_pose.dump_scored_pdb( mini_pdb, *scorefxn );
				}
				(*scorefxn)(bound_pose);
				TR << "Final score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
				TR << "Successfully finished minimizing input." << std::endl;
			}
			(*scorefxn)(bound_pose);

			//setup the unbound pose
			core::pose::Pose unbound_pose = bound_pose;
			core::Real const unbound_dist = 80.;
			//Size const rb_jump = 1; // use the first jump as the one between partners
			Size const rb_jump = bound_pose.num_jump(); // use the LAST jump as the one between partners
			protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose, rb_jump );
			trans_mover.trans_axis( trans_mover.trans_axis() );
			trans_mover.step_size(unbound_dist);
			trans_mover.apply( unbound_pose );
			(*scorefxn)(unbound_pose);
			if ( option[ print_unbound ] ) {
				unbound_pose.dump_pdb( unbo_pdb );
			}

			// define containers for metrics for total complex
			basic::MetricValue<Real> tot_sasa_mval;
			basic::MetricValue<Size> tot_hb_mval;
			basic::MetricValue<Real> tot_packstat_mval;
			basic::MetricValue<Size> tot_unsat_mval;

			// calculate and store total metrics for bound and unbound poses
			//core::Real bound_energy = 0.0, unbound_energy = 0.0, Interface_Energy = 0.0;
			//core::Real bound_sasa = 0.0, unbound_sasa = 0.0, Total_BSA = 0.0;
			//core::Size  bound_hb = 0,   unbound_hb = 0, Interface_HB = 0;
			core::Real /*bound_packstat = 0.0, */unbound_packstat = 0.0/*, Total_packstats = 0.0*/;
			//core::Size  bound_unsat = 0, unbound_unsat = 0, Interface_unsat = 0;

			//calculate interface Energy
			core::Real bound_energy = bound_pose.energies().total_energies()[ total_score ];
			core::Real unbound_energy = unbound_pose.energies().total_energy();
			core::Real Interface_Energy = bound_energy - unbound_energy;

			//delta sasa calculation
			bound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
			core::Real bound_sasa = tot_sasa_mval.value();
			unbound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
			core::Real unbound_sasa = tot_sasa_mval.value();
			core::Real Total_BSA = unbound_sasa - bound_sasa;

			//interface hb calculation
			bound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
			core::Size bound_hb = tot_hb_mval.value();
			unbound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
			core::Size unbound_hb = tot_hb_mval.value();
			core::Size Interface_HB = bound_hb - unbound_hb;

			//packstat calculation
			bound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
			core::Real bound_packstat = tot_packstat_mval.value();
			//unbound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
			//unbound_packstat = tot_packstat_mval.value();
			core::Real Total_packstats = bound_packstat - unbound_packstat;

			//unsat polar calculation
			bound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
			core::Size bound_unsat = tot_unsat_mval.value();
			unbound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
			core::Size unbound_unsat = tot_unsat_mval.value();
			core::Size Interface_unsat = bound_unsat - unbound_unsat;

			//print RMSD between input and minimized ligands (only after superposition)
			if ( option[print_complex] ) {
				std::cout << "computing post minization ligand RMSD for " << input_pdb_name << ":"
					<< ligand_rmsd(pre_min_darc_pose, bound_pose) << std::endl;
			}

			std::cout << "Interface_Scores:"<< tag <<"\t"<< input_pdb_name <<"\t" << bound_energy <<"\t" << Interface_Energy <<"\t"<< Total_BSA <<"\t"<< Interface_HB <<"\t"<< Total_packstats <<"\t"<< Interface_unsat << std::endl;
			if ( option[ iface_rmsd ] ) {
				/*core::Real CA_rms = */calpha_pdb_superimpose_pose( unbound_pose, ref_pose);
				core::Real CA_rms = core::scoring::CA_rmsd( unbound_pose, ref_pose );
				std::cout << "superimpose to native. Rms to native: " << CA_rms << std::endl;
				core::Real heavyatom_rms = interface_rmsd( ref_pose, unbound_pose );
				std::cout << "Interface rmsd: " << heavyatom_rms << std::endl;
				outstream << mini_pdb << ' ' << CA_rms << ' ' << heavyatom_rms <<std::endl;
			}

		}//end for loop of all decoys

		outstream.close();
		outstream.clear();

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
