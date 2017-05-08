// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/oop/OopRandomPuckMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ncbb/oop/OopDockDesignProtocol.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/ncbb/oop/OopDockDesignProtocolCreator.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <string>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::oop;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR( "ODDP" );

namespace protocols {
namespace ncbb {
namespace oop {

OopDockDesignProtocol::OopDockDesignProtocol():
	mc_temp_ ( 1.0),
	pert_mc_temp_ (0.8),
	pert_dock_rot_mag_ (1.0),
	pert_dock_trans_mag_ (0.5),
	pert_pep_small_temp_ (0.8),
	pert_pep_small_H_ (2.0),
	pert_pep_small_L_ (2.0),
	pert_pep_small_E_ (2.0),
	pert_pep_shear_temp_ (0.8),
	pert_pep_shear_H_ (2.0),
	pert_pep_shear_L_ (2.0),
	pert_pep_shear_E_ (2.0),

	pert_pep_num_rep_ (100),
	pert_num_ (10),
	dock_design_loop_num_ (10),

	no_design_(false),
	final_design_min_ (true),
	use_soft_rep_ (false),
	mc_initial_pose_ (false),
	oop_design_first_ (false),

	pymol_ (false),
	keep_history_ (false)
{
	Mover::type("OopDockDesignProtocol");

	score_fxn_ = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

OopDockDesignProtocol::OopDockDesignProtocol(
	scoring::ScoreFunctionOP score_function,
	core::Real const mc_temp,
	core::Real const pert_mc_temp,
	core::Real const pert_dock_rot_mag,
	core::Real const pert_dock_trans_mag,
	core::Real const pert_pep_small_temp,
	core::Real const pert_pep_small_H,
	core::Real const pert_pep_small_L,
	core::Real const pert_pep_small_E,
	core::Real const pert_pep_shear_temp,
	core::Real const pert_pep_shear_H,
	core::Real const pert_pep_shear_L,
	core::Real const pert_pep_shear_E,

	core::Size const pert_pep_num_rep,
	core::Size const pert_num,
	core::Size const dock_design_loop_num,

	bool const no_design,
	bool const final_design_min,
	bool const use_soft_rep,
	bool const mc_initial_pose,
	bool const oop_design_first,

	bool const pymol,
	bool const keep_history

): Mover("OopDockDesignProtocol"),
	score_fxn_(score_function),
	mc_temp_ (mc_temp),
	pert_mc_temp_(pert_mc_temp),
	pert_dock_rot_mag_(pert_dock_rot_mag),
	pert_dock_trans_mag_(pert_dock_trans_mag),
	pert_pep_small_temp_(pert_pep_small_temp),
	pert_pep_small_H_(pert_pep_small_H),
	pert_pep_small_L_(pert_pep_small_L),
	pert_pep_small_E_(pert_pep_small_E),
	pert_pep_shear_temp_(pert_pep_shear_temp),
	pert_pep_shear_H_(pert_pep_shear_H),
	pert_pep_shear_L_(pert_pep_shear_L),
	pert_pep_shear_E_(pert_pep_shear_E),
	pert_pep_num_rep_(pert_pep_num_rep),
	pert_num_(pert_num),
	dock_design_loop_num_(dock_design_loop_num),
	no_design_(no_design),
	final_design_min_(final_design_min),
	use_soft_rep_(use_soft_rep),
	mc_initial_pose_(mc_initial_pose),
	oop_design_first_(oop_design_first),
	pymol_(pymol),
	keep_history_(keep_history)
{
}

OopDockDesignProtocol::OopDockDesignProtocol(
	scoring::ScoreFunctionOP score_function,
	core::Real const mc_temp,
	core::Real const pert_dock_rot_mag,
	core::Real const pert_dock_trans_mag,
	core::Size const dock_design_loop_num,
	bool const no_design,
	bool const final_design_min,
	bool const pymol,
	bool const keep_history
): Mover("OopDockDesignProtocol"),
	score_fxn_(score_function),
	mc_temp_ (mc_temp),
	pert_mc_temp_ (0.8),
	pert_dock_rot_mag_ (pert_dock_rot_mag),
	pert_dock_trans_mag_ (pert_dock_trans_mag),
	pert_pep_small_temp_ (0.8),
	pert_pep_small_H_ (2.0),
	pert_pep_small_L_ (2.0),
	pert_pep_small_E_ (2.0),
	pert_pep_shear_temp_ (0.8),
	pert_pep_shear_H_ (2.0),
	pert_pep_shear_L_ (2.0),
	pert_pep_shear_E_ (2.0),

	pert_pep_num_rep_ (100),
	pert_num_ (10),
	dock_design_loop_num_(dock_design_loop_num),

	no_design_(no_design),
	final_design_min_(final_design_min),
	use_soft_rep_ (false),
	mc_initial_pose_ (false),
	oop_design_first_ (false),

	pymol_(pymol),
	keep_history_(keep_history)
{}

void
OopDockDesignProtocol::apply(
	core::pose::Pose & pose
) {
	// create score function

	scoring::ScoreFunctionOP soft_score_fxn  = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*soft_score_fxn);
	soft_score_fxn->set_etable( FA_STANDARD_SOFT );

	scoring::ScoreFunctionOP pert_score_fxn;
	if ( use_soft_rep_  ) {
		pert_score_fxn = soft_score_fxn;
	} else {
		pert_score_fxn = score_fxn_;
	}

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn_, mc_temp_ ) );

	/*********************************************************************************************************************
	Pertubation Phase
	**********************************************************************************************************************/

	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *pert_score_fxn, pert_mc_temp_ ) );

	/*********************************************************
	Docking Setup
	**********************************************************/

	// create a rigid body mover to move the peptide around in the pocket
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, pert_dock_rot_mag_,  pert_dock_trans_mag_) );

	/*********************************************************
	Peptide Setup
	**********************************************************/

	// get peptide start and end positions
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.size() );
	TR << "pep_start: " << pep_start << " pep_end: " << pep_end << std::endl;

	// create movemap for peptide
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep_start, pep_end);

	////kdrew: automatically find oop positions
	utility::vector1< core::Size > oop_seq_positions = core::pose::ncbb::initialize_oops(pose);

	for ( Size i = 1; i <= oop_seq_positions.size(); i++ ) {
		pert_pep_mm->set_bb( oop_seq_positions[i], false );

		if ( score_fxn_->has_zero_weight( core::scoring::atom_pair_constraint ) ) {
			score_fxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		}
	}


	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, pert_pep_small_temp_, 1 ) );
	pert_pep_small->angle_max( 'H', pert_pep_small_H_ );
	pert_pep_small->angle_max( 'L', pert_pep_small_L_ );
	pert_pep_small->angle_max( 'E', pert_pep_small_E_ );

	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, pert_pep_shear_temp_ , 1 ) );
	pert_pep_shear->angle_max( 'H', pert_pep_shear_H_ );
	pert_pep_shear->angle_max( 'L', pert_pep_shear_L_ );
	pert_pep_shear->angle_max( 'E', pert_pep_shear_E_ );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	//pert_pep_random->add_mover( pert_pep_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, pert_pep_num_rep_ ) );

	/*********************************************************
	OopPuck Setup
	**********************************************************/

	//oop::OopRandomSmallMoverOP opm_small( new oop::OopRandomSmallMover( option[ oddm::oop_positions ].value(), 2.0 ) );
	//oop::OopRandomPuckMoverOP opm_puck( new oop::OopRandomPuckMover( option[ oddm::oop_positions ].value() ) );
	simple_moves::oop::OopRandomSmallMoverOP opm_small( new simple_moves::oop::OopRandomSmallMover( oop_seq_positions, 2.0 ) );
	simple_moves::oop::OopRandomPuckMoverOP opm_puck( new simple_moves::oop::OopRandomPuckMover( oop_seq_positions ) );

	/******************************************************************************
	Rotamer Trials Setup
	*******************************************************************************/

	// create a task factory and task operations
	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP pert_tf( new TaskFactory() );
	pert_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	//operation::RestrictToInterfaceOP pert_rtio( new operation::RestrictToInterface(1, 2) ); //magic numbers: assume chains 1 and 2
	//pert_tf->push_back( pert_rtio );

	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt( new simple_moves::EnergyCutRotamerTrialsMover( pert_score_fxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 0.5 );
	pert_random->add_mover( opm_small, 0.5 );
	pert_random->add_mover( opm_puck, 0.1 );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	/*********************************************************************************************************************
	Design Min Phase
	**********************************************************************************************************************/
	/*********************************************************
	Design Setup
	**********************************************************/

	// create a task factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );


	if ( no_design_ ) {
		core::pack::task::operation::RestrictToRepackingOP rtrp( new core::pack::task::operation::RestrictToRepacking() );
		desn_tf->push_back( rtrp );
	}

	/*
	desn_tf->push_back( pert_rtio ); //not bothering to construct it a second time
	*/

	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( pert_score_fxn );

	/*********************************************************
	Minimize Setup
	**********************************************************/

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
	//kdrew: set backbone of target false and backbone of oop true, decide whether to do this or not
	desn_mm->set_bb( false );
	desn_mm->set_bb_true_range( pep_start, pep_end );
	//desn_mm->set_bb( true );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );

	//kdrew: loop through oop positions and set VZP and VYP torsions to minimize
	/*
	for( Size i = 1; i <= oop_seq_positions.size(); i++  )
	{
	AtomID aidVZP( pose.residue(oop_seq_positions[i]).atom_index("VZP"), oop_seq_positions[i] );
	AtomID aidCYP( pose.residue(oop_seq_positions[i]).atom_index("CYP"), oop_seq_positions[i] );
	AtomID aidN( pose.residue(oop_seq_positions[i]).atom_index("N"), oop_seq_positions[i] );
	AtomID aidCA( pose.residue(oop_seq_positions[i]).atom_index("CA"), oop_seq_positions[i] );
	DOF_ID const vzp_tor_dof_id( pose.atom_tree().torsion_angle_dof_id( aidVZP, aidCYP, aidN, aidCA ) );
	desn_mm->set( vzp_tor_dof_id, true );

	AtomID aidVYP( pose.residue(oop_seq_positions[i]+1).atom_index("VYP"), oop_seq_positions[i]+1 );
	AtomID aidCZP( pose.residue(oop_seq_positions[i]+1).atom_index("CZP"), oop_seq_positions[i]+1 );
	AtomID aidNpost( pose.residue(oop_seq_positions[i]+1).atom_index("N"), oop_seq_positions[i]+1 );
	AtomID aidCApost( pose.residue(oop_seq_positions[i]+1).atom_index("CA"), oop_seq_positions[i]+1 );
	DOF_ID const vyp_tor_dof_id( pose.atom_tree().torsion_angle_dof_id( aidVYP, aidCZP, aidNpost, aidCApost ) );
	desn_mm->set( vyp_tor_dof_id, true );
	}
	*/

	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min( new TaskAwareMinMover( desn_min, desn_tf ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a sequence mover to hold pack rotamers and minimization movers
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pr );
	desn_sequence->add_mover( desn_ta_min );

	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	if ( pymol_ ) {
		protocols::moves::PyMOLObserverOP pymover = protocols::moves::AddPyMOLObserver(pose, keep_history_ );
	}

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= Size( dock_design_loop_num_ ); ++k ) {

		pert_mc->reset(pose);
		//kdrew: a quick design/repack prior to pertubation, often the initial structure given is aligned to hotspot Ca Cb vector
		//kdrew: and do not want to perturb away until designed in hotspot residue
		if ( k == 1 && oop_design_first_ ) {
			desn_sequence->apply( pose );
		}

		// pert loop
		for ( Size j = 1; j <= Size( pert_num_ ); ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
			curr_job->add_string_real_pair( "ENERGY_PERT (pert score)", (*pert_score_fxn)(pose) );
		}
		pert_mc->recover_low( pose );
		curr_job->add_string_real_pair( "ENERGY_PERT (pert score) recovered low", (*pert_score_fxn)(pose) );

		// design
		TR << "DESIGN: " << k << std::endl;
		desn_sequence->apply( pose );
		curr_job->add_string_real_pair( "ENERGY_DESN (hard score)", (*score_fxn_)(pose) );

		//kdrew: reset mc after first cycle if not considering initial pose
		if ( !mc_initial_pose_  && k == 1 ) {
			mc->reset(pose);
			TR<< "after mc->reset" << std::endl;
			mc->show_state();
		}

		TR<< "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR<< "post mc->boltzmann" << std::endl;
		mc->show_state();

	}//dock_design for loop

	mc->recover_low( pose );

	curr_job->add_string_real_pair( "ENERGY_FINAL (pert score) ", (*pert_score_fxn)(pose) );
	curr_job->add_string_real_pair( "ENERGY_FINAL (hard score) ", (*score_fxn_)(pose) );

	TR << "Ending main loop..." << std::endl;
	TR << "Checking pose energy..." << std::endl;

	// calc energy
	//energy_complex = (*score_fxn_)(pose);
	// we don't actually apply a cutoff here, it seems.
	TR << "Energy less than cutoff, doing final design and running filters..." << std::endl;

	if ( final_design_min_ ) {
		final_design_min( pose, score_fxn_, desn_tf );
	}

	calculate_statistics( curr_job, pose, score_fxn_ );

}

protocols::moves::MoverOP
OopDockDesignProtocol::clone() const {
	return protocols::moves::MoverOP( new OopDockDesignProtocol (
		score_fxn_,
		mc_temp_,
		pert_mc_temp_,
		pert_dock_rot_mag_,
		pert_dock_trans_mag_,
		pert_pep_small_temp_,
		pert_pep_small_H_,
		pert_pep_small_L_,
		pert_pep_small_E_,
		pert_pep_shear_temp_,
		pert_pep_shear_H_,
		pert_pep_shear_L_,
		pert_pep_shear_E_,
		pert_pep_num_rep_,
		pert_num_,
		dock_design_loop_num_,
		no_design_,
		final_design_min_,
		use_soft_rep_,
		mc_initial_pose_,
		oop_design_first_,
		pymol_,
		keep_history_
		) );
}

void
OopDockDesignProtocol::parse_my_tag
(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	init_common_options( tag, data, score_fxn_, mc_temp_, pert_mc_temp_, pert_dock_rot_mag_, pert_dock_trans_mag_, pert_pep_small_temp_, pert_pep_small_H_, pert_pep_small_L_, pert_pep_small_E_, pert_pep_shear_temp_,pert_pep_shear_H_, pert_pep_shear_L_, pert_pep_shear_E_,pert_pep_num_rep_, pert_num_, dock_design_loop_num_, no_design_, final_design_min_, use_soft_rep_, mc_initial_pose_, pymol_, keep_history_ );

	if ( tag->hasOption( "oop_design_first") ) {
		oop_design_first_ = tag->getOption<bool>("oop_design_first", oop_design_first_);
	} else {
		oop_design_first_ = false;
	}
}

// MoverCreator
// XRW TEMP moves::MoverOP OopDockDesignProtocolCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new OopDockDesignProtocol() );
// XRW TEMP }

// XRW TEMP std::string OopDockDesignProtocolCreator::keyname() const {
// XRW TEMP  return OopDockDesignProtocol::mover_name();
// XRW TEMP }

// XRW TEMP std::string OopDockDesignProtocol::mover_name(){
// XRW TEMP  return "OopDockDesign";
// XRW TEMP }

std::string OopDockDesignProtocol::get_name() const {
	return mover_name();
}

std::string OopDockDesignProtocol::mover_name() {
	return "OopDockDesign";
}

void OopDockDesignProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	add_attributes_for_common_options( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "oop_design_first", xsct_rosetta_bool, "First optimize sequence of the oop ligand", "false" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string OopDockDesignProtocolCreator::keyname() const {
	return OopDockDesignProtocol::mover_name();
}

protocols::moves::MoverOP
OopDockDesignProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP( new OopDockDesignProtocol );
}

void OopDockDesignProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	OopDockDesignProtocol::provide_xml_schema( xsd );
}


}
}
}

