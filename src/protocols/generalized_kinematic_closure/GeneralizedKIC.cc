// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/GeneralizedKIC.cc
/// @brief  Kinematic closure of arbitrary segments that could go through side-chains (e.g. disulfides).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// BOINC includes -- keep these first:
#ifdef BOINC
#include <utility/io/izstream.hh>
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC

// Unit Headers
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKICCreator.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>

#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::T;
using basic::Error;
using basic::Warning;


namespace protocols {
namespace generalized_kinematic_closure {

static THREAD_LOCAL basic::Tracer TR( "protocols.generalized_kinematic_closure.GeneralizedKIC" );

// XRW TEMP std::string
// XRW TEMP GeneralizedKICCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return GeneralizedKIC::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP GeneralizedKICCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new GeneralizedKIC );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP GeneralizedKIC::mover_name()
// XRW TEMP {
// XRW TEMP  return "GeneralizedKIC";
// XRW TEMP }

/// @brief Constructor for GeneralizedKIC mover.
///
GeneralizedKIC::GeneralizedKIC():
	Mover("GeneralizedKIC"),
	loopresidues_(),
	tailresidues_(),
	lower_anchor_connID_(0),
	upper_anchor_connID_(0),
	build_ideal_geometry_(false),
	last_run_successful_(false),
	atomlist_(),
	pivot_1_rsd_(0),
	pivot_1_atmname_(""),
	pivot_2_rsd_(0),
	pivot_2_atmname_(""),
	pivot_3_rsd_(0),
	pivot_3_atmname_(""),
	perturberlist_(),
	filterlist_(),
	selector_( selector::GeneralizedKICselectorOP( new selector::GeneralizedKICselector ) ),
	n_closure_attempts_(100),
	min_solution_count_(0),
	rosettascripts_filter_(),
	rosettascripts_filter_exists_(false),
	pre_selection_mover_(),
	pre_selection_mover_exists_(false),
	ntries_before_giving_up_(0),
	solutions_(),
	attach_boinc_ghost_observer_(false),
	low_memory_mode_(false),
	dont_fail_if_no_solution_found_(false),
	bondangle_solutions_(),
	bondlength_solutions_(),
	correct_polymer_dependent_atoms_(false)
	//TODO -- make sure above data are copied properly when duplicating this mover.
{}

/// @brief Copy constructor for GeneralizedKIC mover.
///
GeneralizedKIC::GeneralizedKIC( GeneralizedKIC const &src ):
	Mover("GeneralizedKIC"),
	loopresidues_(src.loopresidues_),
	tailresidues_(src.tailresidues_),
	lower_anchor_connID_(src.lower_anchor_connID_),
	upper_anchor_connID_(src.upper_anchor_connID_),
	build_ideal_geometry_(src.build_ideal_geometry_),
	last_run_successful_(src.last_run_successful_),
	atomlist_(src.atomlist_),
	pivot_1_rsd_(src.pivot_1_rsd_),
	pivot_1_atmname_(src.pivot_1_atmname_),
	pivot_2_rsd_(src.pivot_2_rsd_),
	pivot_2_atmname_(src.pivot_2_atmname_),
	pivot_3_rsd_(src.pivot_3_rsd_),
	pivot_3_atmname_(src.pivot_3_atmname_),
	perturberlist_(), //Will be cloned
	filterlist_(), //Will be cloned
	selector_( src.selector_->clone() ), //Cloned!
	n_closure_attempts_(src.n_closure_attempts_),
	min_solution_count_(src.min_solution_count_),
	rosettascripts_filter_(src.rosettascripts_filter_), //Not cloned!
	rosettascripts_filter_exists_(src.rosettascripts_filter_exists_),
	pre_selection_mover_(src.pre_selection_mover_), //Not cloned!
	pre_selection_mover_exists_(src.pre_selection_mover_exists_),
	ntries_before_giving_up_(src.ntries_before_giving_up_),
	solutions_(), //Copied below.
	attach_boinc_ghost_observer_(src.attach_boinc_ghost_observer_),
	low_memory_mode_(src.low_memory_mode_),
	dont_fail_if_no_solution_found_(src.dont_fail_if_no_solution_found_),
	bondangle_solutions_(src.bondangle_solutions_),
	bondlength_solutions_(src.bondlength_solutions_),
	correct_polymer_dependent_atoms_(src.correct_polymer_dependent_atoms_)
	//TODO -- make sure above data are copied properly when duplicating this mover.
{
	//Clone elements in the perturber list
	perturberlist_.clear();
	for ( core::Size i=1, imax=src.perturberlist_.size(); i<=imax; ++i ) {
		perturberlist_.push_back( src.perturberlist_[i]->clone() );
	}

	//Clone elements in the filter list
	filterlist_.clear();
	for ( core::Size i=1, imax=src.filterlist_.size(); i<=imax; ++i ) {
		filterlist_.push_back( src.filterlist_[i]->clone() );
	}

	//Clone elements (poses) in solutions list.  Note -- potentially memory-expensive!.
	solutions_.clear();
	for ( core::Size i=1; i<=src.solutions_.size(); ++i ) {
		solutions_.push_back( src.solutions_[i]->clone() );
	}

	return;
}

/// @brief Destructor for GeneralizedKIC mover.
///
GeneralizedKIC::~GeneralizedKIC() = default;

/// @brief Clone operator to create a pointer to a fresh GeneralizedKIC object that copies this one.
///
protocols::moves::MoverOP GeneralizedKIC::clone() const {
	return protocols::moves::MoverOP( new GeneralizedKIC( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh GeneralizedKIC object that does NOT copy this one.
///
protocols::moves::MoverOP GeneralizedKIC::fresh_instance() const {
	return protocols::moves::MoverOP( new GeneralizedKIC );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void GeneralizedKIC::apply ( core::pose::Pose & pose )
{
	if ( preselection_mover_exists() && low_memory_mode() && TR.Warning.visible() ) {
		TR.Warning << "*******************************************************************************" << std::endl;
		TR.Warning << "WARNING!  Use of a preselection mover with GeneralizedKIC's low-memory mode is not advisable.  To save memory, in low-memory mode GeneralizedKIC stores only mainchain DoFs of solutions, not full solution poses.  As a result, preselection movers' effects are not stored.  The preselection mover is reapplied after a solution is selected, but if the preselection mover has any stochastic behaviour, it is not guaranteed to produce the same result." << std::endl;
		TR.Warning << "*******************************************************************************" << std::endl;
		TR.Warning.flush();
	}

	clear_stored_solutions(); //Wipe any solution poses that have been stored.

	if ( loopresidues_.size()==0 ) {
		utility_exit_with_message("Error!  GeneralizedKIC apply() function called without setting a loop to be closed.\n");
	}

	check_loop_residues_sensible( pose ); //Check that loop residues haven't been specified multiple times.
	check_tail_residues_sensible( pose ); //Check that tail residues haven't been specified multiple times, and that the tail residues don't overlap with loop residues.

	infer_anchor_connIDs(pose);

	core::pose::Pose perturbedloop_pose; //A pose in which the perturbed loop will be stored.
	utility::vector1 < std::pair < core::Size, core::Size > > residue_map; //Map of < residue in perturbedloop_pose, residue in pose >
	utility::vector1 < std::pair < core::Size, core::Size > > tail_residue_map; //Map of < tail residue in perturbedloop_pose, residue in pose >.
	//Tail residues are residues that are not in the loop to be closed, but which "come along for the ride".

	addloweranchor(perturbedloop_pose, pose);
	addloopgeometry(perturbedloop_pose, pose, build_ideal_geometry_, residue_map);
	addupperanchor(perturbedloop_pose, pose);
	addtailgeometry(perturbedloop_pose, pose, build_ideal_geometry_, residue_map, tail_residue_map);

	core::Size solution_index(0); //Storage for the index of the solution chosen by the selector.

	//Storage for the torsions, bondangles, and bondlengths vectors for the selected solution in low-memory mode:
	utility::vector1 < core::Real > selected_torsions;
	utility::vector1 < core::Real > selected_bondangles;
	utility::vector1 < core::Real > selected_bondlengths;

	last_run_successful_ = doKIC(perturbedloop_pose, pose, residue_map, tail_residue_map, solution_index, selected_torsions, selected_bondangles, selected_bondlengths );
	if ( rosettascripts_filter_exists_ ) rosettascripts_filter_->set_value(last_run_successful_);

	if ( last_run_successful_ && solution_index!=0 ) {
		if ( TR.visible() ) { TR << "Closure successful." << std::endl; TR.flush(); }
		if ( low_memory_mode() ) {
			set_last_move_status( set_final_solution( pose, perturbedloop_pose, residue_map, tail_residue_map, selected_torsions, selected_bondangles, selected_bondlengths ) );
			if ( get_last_move_status() != protocols::moves::MS_SUCCESS ) { set_last_move_status( protocols::moves::FAIL_RETRY ); }
		} else {
			set_final_solution( pose, solution_index ); //Set the pose to the solution.
			set_last_move_status( protocols::moves::MS_SUCCESS );
		}
	} else {
		if ( TR.visible() ) { TR << "Closure unsuccessful." << std::endl; TR.flush();}
		if ( !dont_fail_if_no_solution_found_ ) {
			set_last_move_status( protocols::moves::FAIL_RETRY );
		}
	}

	clear_stored_solutions(); //Wipe any solution poses that have been stored.

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("GeneralizedKIC").
// XRW TEMP std::string GeneralizedKIC::get_name() const{
// XRW TEMP  return "GeneralizedKIC";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
GeneralizedKIC::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & /*pose*/
) {
	using namespace core::id;
	using namespace protocols::filters;

	TR << "Parsing options for GeneralizedKIC (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	if ( tag->hasOption("low_memory_mode") ) {
		set_low_memory_mode( tag->getOption<bool>( "low_memory_mode" , false ) );
	}

	set_correct_polymer_dependent_atoms( tag->getOption<bool>( "correct_polymer_dependent_atoms", correct_polymer_dependent_atoms() ) );

	set_dont_fail_if_no_solution_found( tag->getOption<bool>( "dont_fail_if_no_solution_found", dont_fail_if_no_solution_found_ ) );

	runtime_assert_string_msg( tag->hasOption("selector"), "RosettaScript parsing error: the <GeneralizedKIC> mover must have a selector specfied with the selector=<selector_name> statement." );
	set_selector_type( tag->getOption<std::string>("selector", "") );
	if ( tag->hasOption("selector_scorefunction") ) {
		set_selector_scorefunction( protocols::rosetta_scripts::parse_score_function( tag, "selector_scorefunction", data_map )->clone() );
	}
	if ( tag->hasOption("selector_kbt") ) {
		set_selector_kbt( tag->getOption<core::Real>("selector_kbt", 1.0) );
	}
	set_ntries_before_giving_up(tag->getOption<core::Size>("stop_if_no_solution", 0)); //Number of tries to make before stopping if no solution has been found yet.
	if ( tag->hasOption("closure_attempts") ) set_closure_attempts( tag->getOption<core::Size>("closure_attempts", 100) );

	//Check for depreciated option:
	runtime_assert_string_msg ( !tag->hasOption("stop_when_solution_found"),
		"Error encountered when parsing options for GeneralizedKIC: The \"stop_when_solution_found\" option has been removed.  Instead, you can use \"stop_when_n_solutions_found\", which takes an integer rather than a boolean and stops after at least N solutions have been found." );
	if ( tag->hasOption("stop_when_n_solutions_found") ) set_min_solution_count( tag->getOption<core::Size>("stop_when_n_solutions_found", 0) );

	if ( tag->hasOption("pre_selection_mover") ) {
		protocols::moves::MoverOP curmover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "pre_selection_mover" ), movers );
		set_preselection_mover(curmover);
		TR << "GeneralizedKIC mover \"" << tag->getOption< std::string >("name", "") << "\" has been assigned mover \"" << tag->getOption< std::string >("pre_selection_mover") << "\" as a pre-selection mover that will be applied to all solutions prior to applying the selector." << std::endl;
	}
	if ( tag->hasOption( "contingent_filter" ) ) {
		FilterOP curfilter = protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "contingent_filter" ), filters  );
		runtime_assert_string_msg( curfilter != nullptr, "Invalid filter specified with contingent_filter tag in GeneralizedKIC." );
		rosettascripts_filter_exists_=true;
		rosettascripts_filter_ = utility::pointer::dynamic_pointer_cast< ContingentFilter >(curfilter);
		runtime_assert_string_msg( rosettascripts_filter_ != nullptr, "Only a ContingentFilter can be passed to GeneralizedKIC with the contingent_filter tag." );
		rosettascripts_filter_->set_value(true); //By default, the contingent filter is set to "True".
		TR << "GeneralizedKIC mover \"" << tag->getOption<std::string>("name", "") << "\" linked to ContingentFilter filter \"" << tag->getOption< std::string >("contingent_filter") << "\".  The filter's value will be set by the success of the kinematic closure." << std::endl;
	}

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( auto const & branch_tag : branch_tags ) {
		if ( branch_tag->getName() == "AddResidue" ) {
			runtime_assert_string_msg( branch_tag->hasOption("res_index"), "RosettaScript parsing error: the <AddResidue> group within a <GeneralizedKIC> block must include a \"res_index=<index>\" statement.");
			add_loop_residue( branch_tag->getOption<core::Size>("res_index", 0) );
		} else if ( branch_tag->getName() == "AddTailResidue" ) {
			runtime_assert_string_msg( branch_tag->hasOption("res_index"), "RosettaScript parsing error: the <AddTailResidue> group within a <GeneralizedKIC> block must include a \"res_index=<index>\" statement.");
			add_tail_residue( branch_tag->getOption<core::Size>("res_index", 0) );
		} else if ( branch_tag->getName() == "SetPivots" ) {
			runtime_assert_string_msg( branch_tag->hasOption("res1"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"res1=<index>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("res2"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"res2=<index>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("res3"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"res3=<index>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("atom1"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"atom1=<atom_name>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("atom2"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"atom2=<atom_name>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("atom3"), "RosettaScript parsing error: the <SetPivots> group within a <GeneralizedKIC> block must include a \"atom3=<atom_name>\" statement.");
			set_pivot_atoms(
				branch_tag->getOption<core::Size>("res1", 0),
				branch_tag->getOption<std::string>("atom1", ""),
				branch_tag->getOption<core::Size>("res2", 0),
				branch_tag->getOption<std::string>("atom2", ""),
				branch_tag->getOption<core::Size>("res3", 0),
				branch_tag->getOption<std::string>("atom3", "")
			);
		} else if  ( branch_tag->getName() == "SampleCisPeptideBond" ) {
			add_perturber("sample_cis_peptide_bond");
			add_value_to_perturber_value_list( branch_tag->getOption("cis_prob", 0.1) );
			//Loop through the sub-tags to find out what information we're adding to this perturber:
			utility::vector1< utility::tag::TagCOP > const subbranch_tags( branch_tag->getTags() );
			for ( auto const & subbranch_tag : subbranch_tags ) {
				if ( subbranch_tag->getName() == "AddResidue" ) { //Adding a residue to the PERTURBER'S list of residues.
					runtime_assert_string_msg( subbranch_tag->hasOption("index"), "RosettaScript parsing error: the <AddResidue> subgroup within the <SampleCisProline> block in a <GeneralizedKIC> block must include a \"index=<residue_index>\" statement.");
					add_residue_to_perturber_residue_list ( subbranch_tag->getOption<core::Size>("index", 0) );
				}
			}
		} else if ( branch_tag->getName() == "CloseBond" ) {
			runtime_assert_string_msg( branch_tag->hasOption("res1"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"res1=<index>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("res2"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"res2=<index>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("atom1"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"atom1=<atom_name>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("atom2"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"atom2=<atom_name>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("bondlength"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"bondlength=<value>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("angle1"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"angle1=<value>\" statement.");
			runtime_assert_string_msg( branch_tag->hasOption("angle2"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"angle2=<value>\" statement.");
			//runtime_assert_string_msg( (*tag_it)->hasOption("torsion"), "RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must include a \"torsion=<value>\" statement.");
			if ( branch_tag->hasOption("randomize_flanking_torsions") ) {
				runtime_assert_string_msg( branch_tag->hasOption("prioratom_res") && branch_tag->hasOption("prioratom"),
					"RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must define the atom before the bond with the \"prioratom_res=<value>\" and \"prioratom=<atom_name>\" statements if the \"randomize_flanking_torsions\" statement is used.");
				runtime_assert_string_msg( branch_tag->hasOption("followingatom_res") && branch_tag->hasOption("followingatom"),
					"RosettaScript parsing error: the <CloseBond> group within a <GeneralizedKIC> block must define the atom after the bond with the \"followingatom_res=<value>\" and \"followingatom=<atom_name>\" statements if the \"randomize_flanking_torsions\" statement is used.");
			}
			close_bond (
				branch_tag->getOption<core::Size>("res1", 0),
				branch_tag->getOption<std::string>("atom1", ""),
				branch_tag->getOption<core::Size>("res2", 0),
				branch_tag->getOption<std::string>("atom2", ""),
				branch_tag->getOption<core::Size>("prioratom_res", 0),
				branch_tag->getOption<std::string>("prioratom", ""),
				branch_tag->getOption<core::Size>("followingatom_res", 0),
				branch_tag->getOption<std::string>("followingatom", ""),
				branch_tag->getOption<core::Real>("bondlength", 0.0),
				branch_tag->getOption<core::Real>("angle1", 0.0),
				branch_tag->getOption<core::Real>("angle2", 0.0),
				branch_tag->getOption<core::Real>("torsion", 0.0),
				!branch_tag->hasOption("torsion"), //If the torsion option is not specified, randomize this torsion
				branch_tag->getOption<bool>("randomize_flanking_torsions", false)
			);
		} else if ( branch_tag->getName() == "AddPerturber" ) { //If we're adding a perturber, need to loop through sub-tags.
			runtime_assert_string_msg( branch_tag->hasOption("effect"), "RosettaScript parsing error: the <AddPerturber> group within a <GeneralizedKIC> block must include a \"effect=<perturber_effect_type>\" statement." );

			std::string effect( branch_tag->getOption<std::string>("effect", "") );
			add_perturber( effect );

			//If this is a perturber that uses torsion bin transition probabilities, we need to initialize the
			//BinTransitionCalculator object and load a bin_params file:
			if ( effect=="randomize_backbone_by_bins" || effect=="perturb_backbone_by_bins" || effect=="set_backbone_bin" ) {
				load_perturber_bin_params( branch_tag->getOption<std::string>("bin_params_file", "ABBA") );
			}
			if ( effect=="perturb_backbone_by_bins" ) {
				core::Size const iterationcount( branch_tag->getOption<core::Size>("iterations", 1) );
				set_perturber_iterations( iterationcount );
				if ( TR.visible() ) TR << "Set iterations for perturb_backbone_by_bins GeneralizedKICperturber to " << iterationcount << "." << std::endl;
				bool const mustswitch( branch_tag->getOption<bool>("must_switch_bins", false) );
				set_perturber_must_switch_bins( mustswitch );
				if ( TR.visible() ) TR << "The perturb_backbone_by_bins GeneralizedKICperturber " << (mustswitch ? "must always " : "need not necessarily ") << "switch torsion bins." << std::endl;
			}
			if ( effect=="set_backbone_bin" ) { //Get the bin that we're setting these residues to lie in.
				if ( !branch_tag->hasOption( "bin" ) ) utility_exit_with_message( "RosettaScript parsing error: the set_backbone_bin GeneralizedKICperturber requires a \"bin=\" option to be specified." );
				std::string const bin( branch_tag->getOption<std::string>("bin", "") );
				set_perturber_bin( bin );
				if ( TR.visible() ) TR << "The set_backbone_bin GeneralizedKICperturber was set to use bin " << bin << "." << std::endl;
			}

			//If this is the randomize_alpha_backbone_by_rama perturber, the user may have specified a custom Ramachandran map:
			if ( effect == "randomize_alpha_backbone_by_rama" ) {
				if ( branch_tag->hasOption("custom_rama_table") ) {
					std::string const custom_table( branch_tag->getOption<std::string>("custom_rama_table") );
					set_perturber_custom_rama_table( custom_table );
					if ( TR.visible() ) TR << "The perturb_alpha_backbone_by_rama perturber was set to use the \"" << custom_table << "\" custom Ramachandran table for sampling." << std::endl;
				}
			}

			//Loop through the sub-tags to find out what information we're adding to this perturber:
			utility::vector1< utility::tag::TagCOP > const subbranch_tags( branch_tag->getTags() );
			for ( auto const & subbranch_tag : subbranch_tags ) {
				if ( subbranch_tag->getName() == "AddResidue" ) { //Adding a residue to the PERTURBER'S list of residues.
					runtime_assert_string_msg( subbranch_tag->hasOption("index"), "RosettaScript parsing error: the <AddResidue> subgroup within the <AddPerturber> block in a <GeneralizedKIC> block must include a \"index=<residue_index>\" statement.");
					add_residue_to_perturber_residue_list ( subbranch_tag->getOption<core::Size>("index", 0) );
				} else if ( subbranch_tag->getName() == "AddAtoms" ) { //Adding a set of one, two, three, or four atoms to the PERTURBER'S list of atoms.
					runtime_assert_string_msg( subbranch_tag->hasOption("atom1") && subbranch_tag->hasOption("res1"), "RosettaScript parsing error: the <AddAtoms> subgroup within the <AddPerturber> block in a <GeneralizedKIC> block must define at least one atom.");
					utility::vector1 < core::id::NamedAtomID > atomset;
					if ( subbranch_tag->hasOption("atom1") && subbranch_tag->hasOption("res1") ) {
						atomset.push_back( NamedAtomID( subbranch_tag->getOption<std::string>("atom1"), subbranch_tag->getOption<core::Size>("res1") ) );
					}
					if ( subbranch_tag->hasOption("atom2") && subbranch_tag->hasOption("res2") ) {
						atomset.push_back( NamedAtomID( subbranch_tag->getOption<std::string>("atom2"), subbranch_tag->getOption<core::Size>("res2") ) );
					}
					if ( subbranch_tag->hasOption("atom3") && subbranch_tag->hasOption("res3") ) {
						atomset.push_back( NamedAtomID( subbranch_tag->getOption<std::string>("atom3"), subbranch_tag->getOption<core::Size>("res3") ) );
					}
					if ( subbranch_tag->hasOption("atom4") && subbranch_tag->hasOption("res4") ) {
						atomset.push_back( NamedAtomID( subbranch_tag->getOption<std::string>("atom4"), subbranch_tag->getOption<core::Size>("res4") ) );
					}
					add_atomset_to_perturber_atomset_list(atomset);
				} else if ( subbranch_tag->getName() == "AddValue" ) { //Adding a residue to the PERTURBER'S list of residues.
					runtime_assert_string_msg( subbranch_tag->hasOption("value"), "RosettaScript parsing error: the <AddValue> subgroup within the <AddPerturber> block in a <GeneralizedKIC> block must include a \"value=<real_value>\" statement.");
					add_value_to_perturber_value_list ( subbranch_tag->getOption<core::Real>("value", 0.0) );
				}
			}

		} else if ( branch_tag->getName() == "AddFilter" ) { //If we're adding a filter, need to loop through sub-tags.
			runtime_assert_string_msg( branch_tag->hasOption("type"), "RosettaScript parsing error: the <AddFilter> group within a <GeneralizedKIC> block must include a \"type=<filter_type>\" statement." );
			std::string const filtertype( branch_tag->getOption<std::string>("type", "") );
			add_filter( filtertype );

			//Parse filter-specific options:
			if ( filtertype=="backbone_bin" ) {
				//Residue number:
				runtime_assert_string_msg( branch_tag->hasOption("residue"), "RosettaScript parsing error: when adding a backbone_bin filter, the <AddFilter> group within a <GeneralizedKIC> block must have a \"residue=(&int)\" statement." );
				core::Size const resnum( branch_tag->getOption<core::Size>("residue",0) );
				set_filter_resnum(resnum);
				if ( TR.visible() ) TR << "Set the residue number for backbone_bin filter to " << resnum << "." << std::endl;

				//Bin transition probabilities params file:
				load_filter_bin_params( branch_tag->getOption<std::string>("bin_params_file", "ABBA") );

				//Bin name:
				runtime_assert_string_msg( branch_tag->hasOption("bin"), "RosettaScript parsing error: when adding a backbone_bin filter, the <AddFilter> group within a <GeneralizedKIC> block must have a \"bin=\" statement." );
				std::string const binname( branch_tag->getOption<std::string>("bin", "")  );
				set_filter_bin(binname);
				if ( TR.visible() ) TR << "Set the bin name for the backbone_bin filter to " << binname << "." << std::endl;
			} else if ( filtertype == "alpha_aa_rama_check" || filtertype == "rama_prepro_check" ) {
				//Residue number:
				runtime_assert_string_msg( branch_tag->hasOption("residue"), "RosettaScript parsing error: when adding an alpha_aa_rama_check or rama_prepro_check filter, the <AddFilter> group within a <GeneralizedKIC> block must have a \"residue=(&int)\" statement." );
				core::Size const resnum( branch_tag->getOption<core::Size>("residue",0) );
				set_filter_resnum(resnum);
				if ( TR.visible() ) TR << "Set the residue number for " << filtertype << " filter to " << resnum << "." << std::endl;

				core::Real const cutoff_energy( branch_tag->getOption<core::Real>("rama_cutoff_energy", 0.3) );
				set_filter_rama_cutoff_energy( cutoff_energy );
				if ( TR.visible() ) TR << "Set the " << (filtertype == "rama_prepro_check" ? "rama_prepro" : "rama" ) << " term cutoff energy for the " << filtertype << " filter to " << cutoff_energy << std::endl;
			}

			//Loop through the sub-tags to find out what information we're adding to this filter, if any:
			utility::vector1< utility::tag::TagCOP > const subbranch_tags( branch_tag->getTags() );
			for ( auto const & subbranch_tag : subbranch_tags ) {
				if ( subbranch_tag->getName() == "AddFilterParameterReal" ) { //Adding a real-valued filter parameter.
					runtime_assert_string_msg( subbranch_tag->hasOption("value"), "RosettaScript parsing error: the <AddFilterParameterReal> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"value=<real value>\" statement.");
					runtime_assert_string_msg( subbranch_tag->hasOption("name"), "RosettaScript parsing error: the <AddFilterParameterReal> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"name=<string>\" statement.");
					add_filter_parameter(subbranch_tag->getOption<std::string>("name"), subbranch_tag->getOption<core::Real>("value") );
				} else if ( subbranch_tag->getName() == "AddFilterParameterInteger" ) { //Adding an integer-valued filter parameter.
					runtime_assert_string_msg( subbranch_tag->hasOption("value"), "RosettaScript parsing error: the <AddFilterParameterInteger> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"value=<integer value>\" statement.");
					runtime_assert_string_msg( subbranch_tag->hasOption("name"), "RosettaScript parsing error: the <AddFilterParameterInteger> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"name=<string>\" statement.");
					add_filter_parameter(subbranch_tag->getOption<std::string>("name"), subbranch_tag->getOption<core::Size>("value") );
				} else if ( subbranch_tag->getName() == "AddFilterParameterBoolean" ) { //Adding a Boolean-valued filter parameter.
					runtime_assert_string_msg( subbranch_tag->hasOption("value"), "RosettaScript parsing error: the <AddFilterParameterBoolean> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"value=<Boolean value>\" statement.");
					runtime_assert_string_msg( subbranch_tag->hasOption("name"), "RosettaScript parsing error: the <AddFilterParameterBoolean> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"name=<string>\" statement.");
					add_filter_parameter(subbranch_tag->getOption<std::string>("name"), subbranch_tag->getOption<bool>("value") );
				} else if ( subbranch_tag->getName() == "AddFilterParameterString" ) { //Adding a string-valued filter parameter.
					runtime_assert_string_msg( subbranch_tag->hasOption("value"), "RosettaScript parsing error: the <AddFilterParameterString> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"value=<string>\" statement.");
					runtime_assert_string_msg( subbranch_tag->hasOption("name"), "RosettaScript parsing error: the <AddFilterParameterString> subgroup within the <AddFilter> block in a <GeneralizedKIC> block must include a \"name=<string>\" statement.");
					add_filter_parameter(subbranch_tag->getOption<std::string>("name"), subbranch_tag->getOption<std::string>("value") );
				}
			}

		} else if ( branch_tag->getName() == "AddAtomPairDistanceFilter" ) { //Shorthand for adding an atom_pair_distance filter
			runtime_assert_string_msg(branch_tag->hasOption("atom1"), "RosettaScript parsing error: the <AddAtomPairDistanceFilter> block in a <GeneralizedKIC> block must include an \"atom1=<string>\" statement." );
			runtime_assert_string_msg(branch_tag->hasOption("atom2"), "RosettaScript parsing error: the <AddAtomPairDistanceFilter> block in a <GeneralizedKIC> block must include an \"atom2=<string>\" statement." );
			runtime_assert_string_msg(branch_tag->hasOption("res1"), "RosettaScript parsing error: the <AddAtomPairDistanceFilter> block in a <GeneralizedKIC> block must include an \"res1=<integer value>\" statement." );
			runtime_assert_string_msg(branch_tag->hasOption("res2"), "RosettaScript parsing error: the <AddAtomPairDistanceFilter> block in a <GeneralizedKIC> block must include an \"res2=<integer value>\" statement." );
			runtime_assert_string_msg(branch_tag->hasOption("distance"), "RosettaScript parsing error: the <AddAtomPairDistanceFilter> block in a <GeneralizedKIC> block must include an \"distance=<real value>\" statement." );
			add_filter( "atom_pair_distance" );
			add_filter_parameter( "atom1", branch_tag->getOption<std::string>("atom1") );
			add_filter_parameter( "atom2", branch_tag->getOption<std::string>("atom2") );
			add_filter_parameter( "res1", branch_tag->getOption<core::Size>("res1") );
			add_filter_parameter( "res2", branch_tag->getOption<core::Size>("res2") );
			add_filter_parameter( "distance", branch_tag->getOption<core::Real>("distance") );
			if ( branch_tag->hasOption("greater_than") ) add_filter_parameter( "greater_than", branch_tag->getOption<bool>("greater_than") );
		}

	}

	if ( TR.visible() ) TR.flush();

	return;
}


/// @brief Add a residue (by index in the pose) to the list of residues making up the loop to be closed.
void GeneralizedKIC::add_loop_residue( core::Size const residue_index )
{
	runtime_assert_string_msg(residue_index!=0, "Cannot add residue index 0 to the list of loop residues!");
	loopresidues_.push_back(residue_index);
	TR.Debug << "Added residue " << residue_index << " to the list of residues to close by GeneralizedKIC.  List now has " << loopresidues_.size() << " residues in it." << std::endl;
	return;
}

/// @brief Add a residue (by index in the pose) to the list of residues making up the tails attached to
/// the loop to be closed.
/// @details  "Tails" are residues that are not part of the loop to be closed, but which are attached to
/// the loop and which "come along for the ride" as the loop moves.
void GeneralizedKIC::add_tail_residue( core::Size const residue_index )
{
	runtime_assert_string_msg(residue_index!=0, "Cannot add residue index 0 to the list of tail residues!");
	tailresidues_.push_back(residue_index);
	TR.Debug << "Added residue " << residue_index << " to the list of tail residues attached to the loop to close by GeneralizedKIC.  List now has " << tailresidues_.size() << " residues in it." << std::endl;
	return;
}

/// @brief Check that loop residues haven't been specified multiple times.
/// @details Also checks that the loop residues are all within the pose.
void GeneralizedKIC::check_loop_residues_sensible( core::pose::Pose const &pose)
{
	if ( loopresidues_.size()==0 ) return;
	core::Size const nres = pose.size();

	for ( core::Size ir=1, irmax=loopresidues_.size(); ir<=irmax; ++ir ) {
		runtime_assert_string_msg( loopresidues_[ir]>0 && loopresidues_[ir]<=nres, "One or more of the loop residue indices is not within the range of pose residue indices.  (This doesn't make sense -- check your input.)");
		if ( ir>1 ) {
			for ( core::Size jr=1; jr<ir; ++jr ) {
				runtime_assert_string_msg( loopresidues_[ir]!=loopresidues_[jr], "One or more of the loop residues was specified more than once.  (This doesn't make sense -- check your input.)" );
			}
		}
	}

	return;
}

/// @brief Check that tail residues haven't been specified multiple times,
/// and that the tail residues don't overlap with loop residues.
/// @details Also checks that the tail residues are all within the pose.
void GeneralizedKIC::check_tail_residues_sensible( core::pose::Pose const &pose )
{
	if ( tailresidues_.size()==0 ) return;
	core::Size const nres = pose.size();
	core::Size const nloopresidues=loopresidues_.size();

	for ( core::Size ir=1, irmax=tailresidues_.size(); ir<=irmax; ++ir ) {
		runtime_assert_string_msg( tailresidues_[ir]>0 && tailresidues_[ir]<=nres, "One or more of the tail residue indices is not within the range of pose residue indices.  (This doesn't make sense -- check your input.)");
		if ( ir>1 ) {
			for ( core::Size jr=1; jr<ir; ++jr ) {
				runtime_assert_string_msg( tailresidues_[ir]!=tailresidues_[jr], "One or more of the tail residues was specified more than once.  (This doesn't make sense -- check your input.)" );
			}
		}
		if ( nloopresidues>0 ) {
			for ( core::Size jr=1; jr<=nloopresidues; ++jr ) {
				runtime_assert_string_msg( tailresidues_[ir]!=loopresidues_[jr], "One or more tail residues are ALSO part of the loop to be closed.  (This doesn't make sense -- check your input.)" );
			}
		}
	}

	return;
}

/// @brief Clears the list of loop residues
void GeneralizedKIC::clear_loop_residues()
{
	loopresidues_.clear();
}

/// @brief Clears the list of perturber residues for perturber perturber_idx
void
GeneralizedKIC::clear_perturber_residue_list( core::Size const perturber_idx )
{
	perturberlist_[perturber_idx]->clear_residues();
}

/// @brief Clears the list of perturber residues
void
GeneralizedKIC::clear_perturber_residue_lists()
{
	for ( core::Size i=1; i<=perturberlist_.size(); ++i ) {
		clear_perturber_residue_list( i );
	}
}

// @brief Function to set the effect of this mover on parts of the pose that are covalently attached to
// the loop to be closed, but which aren't part of it.  Settings are:
// 0 -- Moves the loop only; can pull apart covalent bonds to anything outside of the loop that isn't
//      an anchor point.
// 1 -- Moves the loop and anything downstream of the loop in the foldtree.  Can still pull apart
//    connections to non-child geometry.
// CURRENTLY DEPRECATED
/*void GeneralizedKIC::set_mover_effect_on_bonded_geometry( core::Size const effect )
{
if(effect > 1) {
utility_exit_with_message("Error!  GeneralizedKIC mover effect mode can only be set to 0 or 1.\n");
}
effect_on_bonded_geometry_ = effect;

return;
}*/


/// @brief Function to set the pivot atoms for kinematic closure:
void GeneralizedKIC::set_pivot_atoms(
	core::Size const rsd1,
	std::string const &at1,
	core::Size const rsd2,
	std::string const &at2,
	core::Size const rsd3,
	std::string const &at3
) {
	pivot_1_rsd_ = rsd1;
	pivot_2_rsd_ = rsd2;
	pivot_3_rsd_ = rsd3;

	pivot_1_atmname_ = at1;
	pivot_2_atmname_ = at2;
	pivot_3_atmname_ = at3;

	return;
}

/// @brief Tells GeneralizedKIC to close a bond, setting bond length, bond angle, and bond torsion values.  This
/// actually just adds appropriate set_dihedral, set_bondangle, and set_bondlength perturbers to the perturber
/// list.  Note that subsequent perturbers OR the closure itself can overwrite the bond length, bond angle, or
/// torsion angles set here.
/// @details
/// @param[in] rsd1 -- The index of the first atom's residue (indexed based on residue indices in the original pose).
/// @param[in] at1 -- The name of the first atom defining the bond to be closed.
/// @param[in] rsd2 -- The index of the second atom's residue (indexed based on residue indices in the original pose).
/// @param[in] at1 -- The name of the second atom defining the bond to be closed.
/// @param[in] bondlength -- The length of the bond between the two atoms.
/// @param[in] bondangle1 -- The bond angle defined by (atom preceding at1 in the chain to be closed), (atm1), (atm2).
/// @param[in] bondangle1 -- The bond angle defined by (atm1), (atm2), (atom following at2 in the chain to be closed).
/// @param[in] torsion -- The torsion angle defined by (atom preceding at1 in the chain to be closed), (atm1), (atm2), (atom following at2 in the chain to be closed).
void GeneralizedKIC::close_bond (
	core::Size const rsd1,
	std::string const &at1,
	core::Size const rsd2,
	std::string const &at2,
	core::Size const rsd1_before,
	std::string const &at1_before,
	core::Size const rsd2_after,
	std::string const &at2_after,
	core::Real const &bondlength,
	core::Real const &bondangle1,
	core::Real const &bondangle2,
	core::Real const &torsion,
	bool const randomize_this_torsion,
	bool const randomize_flanking_torsions
) {
	using namespace core::id;
	utility::vector1 < NamedAtomID > atomset;

	add_perturber("set_bondlength");
	atomset.clear();
	atomset.push_back( NamedAtomID(at1, rsd1) );
	atomset.push_back( NamedAtomID(at2, rsd2) );
	add_atomset_to_perturber_atomset_list(atomset);
	add_value_to_perturber_value_list(bondlength);

	add_perturber("set_dihedral");
	add_atomset_to_perturber_atomset_list(atomset); //If only two atoms are given, the perturber infers the other two needed to define the dihedral as the upstream and downstream atoms in the atom list.
	add_value_to_perturber_value_list( ( randomize_this_torsion ? numeric::random::rg().uniform()*360.0 : torsion) );
	if ( randomize_flanking_torsions ) {
		utility::vector1 < NamedAtomID > atomset2;
		atomset2.push_back( NamedAtomID(at1_before, rsd1_before) );
		atomset2.push_back( NamedAtomID(at1, rsd1) );
		add_atomset_to_perturber_atomset_list(atomset2);
		add_value_to_perturber_value_list(numeric::random::rg().uniform()*360.0);
		atomset2.clear();
		atomset2.push_back( NamedAtomID(at2, rsd2) );
		atomset2.push_back( NamedAtomID(at2_after, rsd2_after) );
		add_atomset_to_perturber_atomset_list(atomset2);
		add_value_to_perturber_value_list(numeric::random::rg().uniform()*360.0);
	}

	add_perturber("set_bondangle");
	atomset.clear();
	atomset.push_back( NamedAtomID(at1, rsd1) ); //If only one atom is given, the perturber infers the other two needed to define the bond angle as the upstream and downstream atoms in the atom list.
	add_atomset_to_perturber_atomset_list(atomset);
	atomset.clear();
	atomset.push_back( NamedAtomID(at2, rsd2) );
	add_atomset_to_perturber_atomset_list(atomset);
	add_value_to_perturber_value_list(bondangle1);
	add_value_to_perturber_value_list(bondangle2);

	return;
}

/// @brief Set the selector (the algorithm controlling how a solution will be chosen
/// from among the solutions passing filters).
void GeneralizedKIC::set_selector_type ( selector::selector_type const &stype) {
	selector_->set_selector_type( stype );
	return;
}

/// @brief Set the selector (the algorithm controlling how a solution will be chosen
/// from among the solutions passing filters).  This sets the selector by name.
void GeneralizedKIC::set_selector_type ( std::string const &stypename) {
	selector_->set_selector_type( stypename );
	return;
}

/// @brief Set the selector's scorefunction.
///
void GeneralizedKIC::set_selector_scorefunction ( core::scoring::ScoreFunctionOP sfxn ) {
	selector_->set_scorefunction( sfxn );
	return;
}

/// @brief Set the selector's Boltzmann temperature value.
///
void GeneralizedKIC::set_selector_kbt ( core::Real const &kbt ) {
	selector_->set_boltzmann_temp( kbt );
	return;
}


/// @brief Add a new perturber to the list of perturbers.
void GeneralizedKIC::add_perturber () {
	perturber::GeneralizedKICperturberOP newperturber( new perturber::GeneralizedKICperturber );
	perturberlist_.push_back( newperturber );
	return;
}


/// @brief Add a new perturber to the list of perturbers, setting the effect.
void GeneralizedKIC::add_perturber ( perturber::perturber_effect const &effect ) {
	add_perturber();
	perturberlist_[perturberlist_.size()]->set_perturber_effect(effect);
	return;
}


/// @brief Add a new perturber to the list of perturbers, setting the effect by effect name string.
void GeneralizedKIC::add_perturber ( std::string const &effectname ) {
	add_perturber();
	perturberlist_[perturberlist_.size()]->set_perturber_effect(effectname);
	return;
}

/// @brief Set a perturber's effect.
/// @details
///
/// @param[in] perturber_index -- The index in the list of perturbers already added.
/// @param[in] effect -- The perturber effect type, based on the perturber::perturber_effect enum (e.g. set_dihedral, randomize_backbone, etc.).
void GeneralizedKIC::set_perturber_effect ( core::Size const perturber_index, perturber::perturber_effect const &effect )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_effect() is out of range." );
	perturberlist_[perturber_index]->set_perturber_effect(effect);
	return;
}

/// @brief Initialize a perturber's BinTransitionCalculator object, and load a bin_params file.
///
void GeneralizedKIC::load_perturber_bin_params( core::Size const perturber_index, std::string const &bin_params_file )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::load_perturber_bin_params() is out of range." );
	perturberlist_[perturber_index]->load_bin_params( bin_params_file );
	return;
} //load_perturber_bin_params

/// @brief Initialize a perturber's BinTransitionCalculator object, and load a bin_params file.
/// @details This acts on the last perturber in the perturber list.
void GeneralizedKIC::load_perturber_bin_params( std::string const &bin_params_file )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::load_perturber_bin_params().");
	load_perturber_bin_params( perturberlist_.size(), bin_params_file );
	return;
} //load_perturber_bin_params

/// @brief Set the number of iterations for a perturber.
///
void GeneralizedKIC::set_perturber_iterations( core::Size const perturber_index, core::Size const val )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_iterations() is out of range." );
	perturberlist_[perturber_index]->set_iterations( val );
	return;
} //set_perturber_iterations

/// @brief Set the number of iterations for a perturber.
/// @details This acts on the last perturber in the perturber list.
void GeneralizedKIC::set_perturber_iterations( core::Size const val )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_perturber_iterations().");
	set_perturber_iterations( perturberlist_.size(), val );
	return;
} //set_perturber_iterations

/// @brief Set whether the perturb_backbone_by_bins perturber requires residues to change their torsion
/// bins every move, or whether they can stay within the same bin.
void GeneralizedKIC::set_perturber_must_switch_bins( core::Size const perturber_index, bool const val )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_must_switch_bins() is out of range." );
	perturberlist_[perturber_index]->set_must_switch_bins( val );
	return;
} //set_perturber_must_switch_bins

/// @brief Set whether the perturb_backbone_by_bins perturber requires residues to change their torsion
/// bins every move, or whether they can stay within the same bin.
/// @details This acts on the last perturber in the perturber list.
void GeneralizedKIC::set_perturber_must_switch_bins( bool const val )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_perturber_must_switch_bins().");
	set_perturber_must_switch_bins( perturberlist_.size(), val );
	return;
} //set_perturber_must_switch_bins

/// @brief Set the bin for the set_backbone_bin perturber.
///
void GeneralizedKIC::set_perturber_bin( core::Size const perturber_index, std::string const &bin )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_bin() is out of range." );
	perturberlist_[perturber_index]->set_bin( bin );
	return;
} //set_perturber_bin

/// @brief Set the bin for the set_backbone_bin perturber.
/// @details This acts on the last perturber in the perturber list.
void GeneralizedKIC::set_perturber_bin( std::string const &bin )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_bin().");
	set_perturber_bin(perturberlist_.size(), bin);
	return;
} //set_perturber_bin

/// @brief Set the custom Ramachandran table for the perturb_alpha_backbone_by_rama perturber.
///
void GeneralizedKIC::set_perturber_custom_rama_table( core::Size const perturber_index, std::string const &table_name ) {
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_custom_rama_table() is out of range." );
	perturberlist_[perturber_index]->set_custom_rama_table( table_name );
	return;
}

/// @brief Set the curstom Ramachandran table for the perturb_alpha_backbone_by_rama perturber.
/// @details This acts on the last perturber in the perturber list.
void GeneralizedKIC::set_perturber_custom_rama_table( std::string const &table_name ) {
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_perturber_custom_rama_table().");
	set_perturber_custom_rama_table(perturberlist_.size(), table_name);
	return;
}

/// @brief Set the custom Ramachandran table for the perturb_alpha_backbone_by_rama perturber.
/// @details This version works by Rama_Table_Type.
void GeneralizedKIC::set_perturber_custom_rama_table( core::Size const perturber_index, core::scoring::Rama_Table_Type const table_type ) {
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_custom_rama_table() is out of range." );
	perturberlist_[perturber_index]->set_custom_rama_table( table_type );
	return;
}

/// @brief Set the curstom Ramachandran table for the perturb_alpha_backbone_by_rama perturber.
/// @details This acts on the last perturber in the perturber list, and works by Rama_Table_Type.
void GeneralizedKIC::set_perturber_custom_rama_table( core::scoring::Rama_Table_Type const table_type ) {
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_perturber_custom_rama_table().");
	set_perturber_custom_rama_table(perturberlist_.size(), table_type);
	return;
}

/// @brief Set whether the perturber's generated poses should be used for BOINC graphics.
/// @details Does nothing outside of the BOINC build.
void GeneralizedKIC::set_perturber_attach_boinc_ghost_observer( core::Size const perturber_index, bool const setting )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::set_perturber_attach_boinc_ghost_observer() is out of range." );
	perturberlist_[perturber_index]->set_attach_boinc_ghost_observer( setting );
	return;
}

/// @brief Set whether the perturber's generated poses should be used for BOINC graphics.
/// @details Does nothing outside of the BOINC build.  This version acts on the last perturber in the perturber list.
void GeneralizedKIC::set_perturber_attach_boinc_ghost_observer( bool const setting )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::set_perturber_attach_boinc_ghost_observer().");
	set_perturber_attach_boinc_ghost_observer(perturberlist_.size(), setting);
	return;
}


/// @brief Add a value to the list of values that a perturber takes.
void GeneralizedKIC::add_value_to_perturber_value_list ( core::Size const perturber_index, core::Real const &val )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::add_value_to_perturber_value_list() is out of range." );
	perturberlist_[perturber_index]->add_inputvalue(val);
	return;
}


/// @brief Add a value to the list of values that a perturber takes.  This operates on the last perturber in the perturber list.
void GeneralizedKIC::add_value_to_perturber_value_list ( core::Real const &val )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::add_value_to_perturber_value_list().");
	add_value_to_perturber_value_list( perturberlist_.size(), val);
	return;
}

/// @brief Add a residue to the list of residues that a perturber takes.  Note that residue_index is based on indices of the ORIGINAL POSE,
/// not the loop in isolation.
void GeneralizedKIC::add_residue_to_perturber_residue_list ( core::Size const perturber_index, core::Size const residue_index )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::add_residue_to_perturber_residue_list() is out of range." );
	perturberlist_[perturber_index]->add_residue(residue_index);
	return;
}

/// @brief Add a residue to the list of residues that a perturber takes.  Note that residue_index is based on indices of the ORIGINAL POSE,
/// not the loop in isolation.  This version acts on the last perturber added.
void GeneralizedKIC::add_residue_to_perturber_residue_list ( core::Size const residue_index )
{
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::add_residue_to_perturber_residue_list().");
	add_residue_to_perturber_residue_list( perturberlist_.size(), residue_index );
	return;
}


/// @brief Add a set of AtomIDs to the list of sets of AtomIDs that a perturber takes.
void GeneralizedKIC::add_atomset_to_perturber_atomset_list ( core::Size const perturber_index, utility::vector1 < core::id::NamedAtomID > const &atomset )
{
	runtime_assert_string_msg( perturber_index <= perturberlist_.size() && perturber_index > 0, "The perturber index provided to GeneralizedKIC::add_atomset_to_perturber_atomset_list() is out of range." );
	perturberlist_[perturber_index]->add_atom_set(atomset);
	return;
}


/// @brief Add a set of AtomIDs to the list of sets of AtomIDs that a perturber takes.  This operates on the last perturber in the perturber list.
void GeneralizedKIC::add_atomset_to_perturber_atomset_list ( utility::vector1 < core::id::NamedAtomID > const &atomset ) {
	runtime_assert_string_msg(perturberlist_.size()>0, "No perturbers specified.  Aborting from GeneralizedKIC::add_atomset_to_perturber_atomset_list().");
	add_atomset_to_perturber_atomset_list( perturberlist_.size(), atomset );
	return;
}


/// @brief Add a new filter to the list of filters.
void GeneralizedKIC::add_filter () {
	filter::GeneralizedKICfilterOP newfilter( new filter::GeneralizedKICfilter );
	filterlist_.push_back( newfilter );
	return;
}


/// @brief Add a new filter to the list of filters, setting the filter type.
void GeneralizedKIC::add_filter ( filter::filter_type const &filtertype ) {
	add_filter();
	filterlist_[filterlist_.size()]->set_filter_type(filtertype);
	return;
}


/// @brief Add a new filter to the list of filters, setting the filter type by name.
/// @details See src/protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.cc for
/// the list of filter type names.
void GeneralizedKIC::add_filter ( std::string const &filtertypename ) {
	add_filter();
	filterlist_[filterlist_.size()]->set_filter_type(filtertypename);
	return;
}


/// @brief Add a real-valued parameter to a filter's parameter list.
void GeneralizedKIC::add_filter_parameter (core::Size const filter_index, std::string const &param_name, core::Real const &value)
{
	filterlist_[filter_index]->add_filter_param(param_name, value);
	return;
}


/// @brief Add an integer-valued parameter to a filter's parameter list.
void GeneralizedKIC::add_filter_parameter (core::Size const filter_index, std::string const &param_name, core::Size const value)
{
	filterlist_[filter_index]->add_filter_param(param_name, value);
	return;
}


/// @brief Add a Boolean-valued parameter to a filter's parameter list.
void GeneralizedKIC::add_filter_parameter (core::Size const filter_index, std::string const &param_name, bool const value)
{
	filterlist_[filter_index]->add_filter_param(param_name, value);
	return;
}


/// @brief Add a string-valued parameter to a filter's parameter list.
void GeneralizedKIC::add_filter_parameter (core::Size const filter_index, std::string const &param_name, std::string const &value)
{
	filterlist_[filter_index]->add_filter_param(param_name, value);
	return;
}


/// @brief Add a real-valued parameter to the last filter's parameter list.
void GeneralizedKIC::add_filter_parameter (std::string const &param_name, core::Real const &value)
{
	add_filter_parameter(filterlist_.size(), param_name, value);
	return;
}


/// @brief Add an integer-valued parameter to the last filter's parameter list.
void GeneralizedKIC::add_filter_parameter (std::string const &param_name, core::Size const value)
{
	add_filter_parameter(filterlist_.size(), param_name, value);
	return;
}


/// @brief Add a Boolean-valued parameter to the last filter's parameter list.
void GeneralizedKIC::add_filter_parameter (std::string const &param_name, bool const value)
{
	add_filter_parameter(filterlist_.size(), param_name, value);
	return;
}


/// @brief Add a string-valued parameter to the last filter's parameter list.
void GeneralizedKIC::add_filter_parameter (std::string const &param_name, std::string const &value)
{
	add_filter_parameter(filterlist_.size(), param_name, value);
	return;
}

/// @brief Set the residue number that a backbone_bin filter is acting on.
///
void GeneralizedKIC::set_filter_resnum( core::Size const filter_index, core::Size const value )
{
	if ( filter_index < 1 || filter_index >filterlist_.size() ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_resnum(): Filter index is out of range.\n" );
	}
	filterlist_[filter_index]->set_resnum(value);
	return;
}

/// @brief Set the residue number that a backbone_bin filter is acting on.
/// @details This version acts on the last filter in the filter list.
void GeneralizedKIC::set_filter_resnum( core::Size const value )
{
	if ( filterlist_.size() < 1 ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_resnum(): No filters have been defined!\n" );
	}
	set_filter_resnum(filterlist_.size(), value);
	return;
}

/// @brief Set the bin name that a backbone_bin filter is looking for.
///
void GeneralizedKIC::set_filter_bin (
	core::Size const filter_index,
	std::string const &name_in
) {
	if ( filter_index < 1 || filter_index >filterlist_.size() ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_bin(): Filter index is out of range.\n" );
	}
	filterlist_[filter_index]->set_binname(name_in);
	return;
}

/// @brief Set the bin name that a backbone_bin filter is looking for.
/// @details This version acts on the last filter in the filter list.
void GeneralizedKIC::set_filter_bin( std::string const &name_in )
{
	if ( filterlist_.size() < 1 ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_bin(): No filters have been defined!\n" );
	}
	set_filter_bin(filterlist_.size(), name_in);
	return;
}

/// @brief Set the rama term cutoff energy for the alpha_aa_rama_check and rama_prepro_check filters.
///
void GeneralizedKIC::set_filter_rama_cutoff_energy(
	core::Size const filter_index,
	core::Real const &cutoff_energy
) {
	if ( filter_index < 1 || filter_index >filterlist_.size() ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_rama_cutoff_energy(): Filter index is out of range.\n" );
	}
	filterlist_[filter_index]->set_rama_cutoff_energy(cutoff_energy);
	return;
}

/// @brief Set the rama term cutoff energy for the alpha_aa_rama_check and rama_prepro_check filters.
/// @details This version acts on the last filter in the filter list.
void GeneralizedKIC::set_filter_rama_cutoff_energy( core::Real const &cutoff_energy )
{
	if ( filterlist_.size() < 1 ) {
		utility_exit_with_message( "In GeneralizedKIC::set_filter_rama_cutoff_energy(): No filters have been defined!\n" );
	}
	set_filter_rama_cutoff_energy( filterlist_.size(), cutoff_energy );
	return;
}

/// @brief Set whether the filter's generated poses should be used for BOINC graphics.
/// @details Does nothing outside of the BOINC build.
void GeneralizedKIC::set_filter_attach_boinc_ghost_observer( core::Size const filter_index, bool const setting )
{
	runtime_assert_string_msg( filter_index <= filterlist_.size() && filter_index > 0, "The filter index provided to GeneralizedKIC::set_filter_attach_boinc_ghost_observer() is out of range." );
	filterlist_[filter_index]->set_attach_boinc_ghost_observer( setting );
	return;
}

/// @brief Set whether the filter's generated poses should be used for BOINC graphics.
/// @details Does nothing outside of the BOINC build.  This version acts on the last filter in the filter list.
void GeneralizedKIC::set_filter_attach_boinc_ghost_observer( bool const setting )
{
	runtime_assert_string_msg(filterlist_.size()>0, "No filters specified.  Aborting from GeneralizedKIC::set_filter_attach_boinc_ghost_observer().");
	set_filter_attach_boinc_ghost_observer(filterlist_.size(), setting);
	return;
}

/// @brief Initialize a filter's BinTransitionCalculator object, and load a bin_params file.
///
void GeneralizedKIC::load_filter_bin_params( core::Size const filter_index, std::string const &bin_params_file )
{
	runtime_assert_string_msg( filter_index <= filterlist_.size() && filter_index > 0, "The filter index provided to GeneralizedKIC::load_filter_bin_params() is out of range." );
	filterlist_[filter_index]->load_bin_params( bin_params_file );
	return;
} //load_filter_bin_params

/// @brief Initialize a filter's BinTransitionCalculator object, and load a bin_params file.
/// @details This acts on the last filter in the filter list.
void GeneralizedKIC::load_filter_bin_params( std::string const &bin_params_file )
{
	runtime_assert_string_msg(filterlist_.size()>0, "No filters specified.  Aborting from GeneralizedKIC::load_filter_bin_params().");
	load_filter_bin_params( filterlist_.size(), bin_params_file );
	return;
} //load_filter_bin_params

/// @brief Set the number of closure attempts.
/// @details Perturbation, closure, and filtering is carried out for every closure
/// attempt.  Successful closures from ALL attempts are then selected from by
/// selectors.
void GeneralizedKIC::set_closure_attempts( core::Size const attempts) {
	TR.Debug << "Closure attempts set to " << attempts << "." << std::endl;
	n_closure_attempts_=attempts;
	return;
}

/// @brief Sets number of tries before giving up.
/// @details If this is set to 0, then no such check is made.
/// The algorithm tries n_closure_attempts_ times if and only if at least one solution is found in the first
/// ntries_before_giving_up_ attempts.
void GeneralizedKIC::set_ntries_before_giving_up ( core::Size const ntries ) {
	if ( ntries!=0 ) TR.Debug << "The algorithm will give up if no closed solution is found in the first " << ntries << " attempts." << std::endl;
	ntries_before_giving_up_=ntries;
	return;
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


/// @brief Gets a path through atoms in the residue from one connection atom to another within a residue.
/// @details  This is not necessarily the shortest path, geometrically, but hopefully it will do for our purposes.  (I
/// didn't want to reimplement Dijkstra's algorithm, since that would tax my programming abilities and would result in
/// unnecessarily slow performance.)  This algorithm works by stepping back (child->parent) from the higher-index atom
/// until (a) it finds the lower-index atom, or (b) it reaches the root.  If it reaches the root, it then starts at the
/// lower-index atom and traces back until it reaches an atom already in the path.  The chosen path is therefore the path
/// from one atom, to the nearest common ancestor, to the other atom.
///
/// @param[in] first_atom -- The index of the first atom in the path.
/// @param[in] second_atom -- The index of the first atom in the path.
/// @param[in] rsd -- The residue object (const instance).
/// @param[out] path_indices -- The list of indices of the atoms found in the path.
///
/// @author Vikram K. Mulligan
void GeneralizedKIC::get_path (
	core::Size const first_atom,
	core::Size const second_atom,
	core::conformation::Residue const &rsd,
	utility::vector1 <core::Size> &path_indices
) {
	runtime_assert_string_msg((first_atom > 0 && first_atom <= rsd.natoms()), "The first atom index passed to GeneralizedKIC::get_path is out of range for the residue.");
	runtime_assert_string_msg((second_atom > 0 && second_atom <= rsd.natoms()), "The second atom index passed to GeneralizedKIC::get_path is out of range for the residue.");
	utility::vector1 <core::Size> temppath; //Temporary storage of path indices.  This vector may have to be reversed.

	core::Size const biggeratom = (first_atom > second_atom ? first_atom : second_atom);
	core::Size const smalleratom = (first_atom > second_atom ? second_atom : first_atom);
	core::Size curatom = biggeratom;

	//TR << "GeneralizedKIC::get_path():  Starting with atom " << rsd.atom_name(curatom) << " (atom " << curatom << ")." << std::endl; //DELETE ME

	bool found=false;

	//Search back from the higher-numbered atom
	while ( true ) {
		temppath.push_back(curatom);
		if ( curatom == smalleratom ) {
			found = true;
			break;
		}
		if ( curatom == 1 ) break;
		curatom = rsd.type().atom_base(curatom); //Move to this atom's parent if we're still searching.
		//TR << "GeneralizedKIC::get_path():  Moving to atom " << rsd.atom_name(curatom) << " (atom " << curatom << ")." << std::endl; //DELETE ME
	}

	//If we didn't find the lower-numbered atom, search back from the lower-numbered atom until we find a common ancestor:
	if ( !found ) {
		//TR << "GeneralizedKIC::get_path():  Path not yet found.  Now tracing back from smaller-numbered atom." << std::endl; //DELETE ME

		curatom = smalleratom;
		utility::vector1 < core::Size > temppath2;

		//TR << "GeneralizedKIC::get_path():  Starting with atom " << rsd.atom_name(curatom) << " (atom " << curatom << ")." << std::endl; //DELETE ME

		while ( true ) {
			temppath2.push_back(curatom);
			for ( core::Size i=1, imax=temppath.size(); i<=imax; ++i ) {
				if ( temppath[i]==curatom ) {
					temppath.resize(i-1); //If we've found a common ancestor, remove anything in the path that traces further back.
					for ( core::Size j=temppath2.size(); j>=1; --j ) { //Append the path from the common ancestor to the lower-numbered atom.
						temppath.push_back(temppath2[j]);
					}
					found=true;
					break;
				}
			}
			if ( found ) break;
			curatom = rsd.type().atom_base(curatom); //Move to this atom's parent if we're still searching.

			//TR << "GeneralizedKIC::get_path():  Path not yet found.  Now tracing back from smaller-numbered atom." << std::endl; //DELETE ME
		}
	}

	//At this point, temppath traces from the higher-numbered to the lower-numbered atom.  Now to copy it over:
	if ( first_atom > second_atom ) path_indices = temppath;
	else { //Copy backwards, if necessary.
		path_indices.clear();
		for ( core::Size i=temppath.size(); i>=1; --i ) {
			path_indices.push_back(temppath[i]);
		}
	}

	return;
}


/// @brief Function to get the FIRST connection in res_with_connection that connects it to other_res:
core::Size GeneralizedKIC::get_connection(
	core::Size const res_with_connection,
	core::Size const other_res,
	core::pose::Pose const &pose
) {
	core::Size nres = pose.size();
	runtime_assert_string_msg(res_with_connection <= nres, "GeneralizedKIC::get_connection got a connection residue that's not in the pose.");
	runtime_assert_string_msg(other_res <= nres, "GeneralizedKIC::get_connection got a residue that's not in the pose.");

	utility::vector1< core::Size > connections = pose.residue(res_with_connection).connections_to_residue(other_res);

	runtime_assert_string_msg(connections.size() > 0 , "The residues passed to GeneralizedKIC::get_connection aren't connected.");

	return connections[1];
}

/// @brief Function that returns true if two residues in a pose have a direct geometric connection,
/// and false otherwise.
bool GeneralizedKIC::has_geometric_connection (
	core::Size const residue1,
	core::Size const residue2,
	core::pose::Pose const &pose
) {
	return pose.residue( residue1 ).is_bonded( residue2 );
}

/// @brief As the list of residues in the loop to be closed is updated, we need to figure out how
/// that loop is connected to the geometry outside of the loop (i.e. what's considred the connection
/// to stationary geometry).  This function loops through all connIDs on the terminal residues of
/// the loop to be closed and picks the first one that links to geometry not in the loop as the
/// anchor connnection.  TODO: Add a manual override to specifiy that a different connection is the
/// anchor.
void GeneralizedKIC::infer_anchor_connIDs(core::pose::Pose const &pose)
{
	//Initialize these to zero:
	lower_anchor_connID_ = 0;
	upper_anchor_connID_ = 0;

	core::Size const loopsize = loopresidues_.size();
	if ( loopsize == 0 ) return;

	core::Size const lower_res = loopresidues_[1];
	core::Size const upper_res = loopresidues_[loopsize];
	core::Size const lower_connID_count = pose.residue(lower_res).n_possible_residue_connections();
	core::Size const upper_connID_count = pose.residue(upper_res).n_possible_residue_connections();

	for ( core::Size i=1; i<=lower_connID_count; ++i ) { //Loop through all of the connections made by the first residue in the loop
		if ( pose.residue(lower_res).connection_incomplete(i) ) continue; //If this connection isn't bonded to anything, go on to the next.
		core::Size const other_res = pose.residue(lower_res).residue_connection_partner(i);
		if ( other_res!=0 && !is_in_list(other_res, loopresidues_) && !is_in_list(other_res, tailresidues_) ) { //If this connection is to something outside the loop (that isn't a tail residue), then we can set the lower_anchor_connID_ to this connection's ID.
			lower_anchor_connID_ = i;
			break;
		}
	}

	for ( core::Size i=1; i<=upper_connID_count; ++i ) { //Loop through all of the connections made by the last residue in the loop
		if ( pose.residue(upper_res).connection_incomplete(i) ) continue; //If this connection isn't bonded to anything, go on to the next.
		core::Size const other_res = pose.residue(upper_res).residue_connection_partner(i);
		if ( other_res!=0 && !is_in_list(other_res, loopresidues_) && !is_in_list(other_res, tailresidues_) ) { //If this connection is to something outside the loop (that isn't a tail residue), then we can set the lower_anchor_connID_ to this connection's ID -- with one caveat:
			if ( loopsize == 1 && lower_anchor_connID_ == i ) continue; //If this is a single-residue loop, then we don't want the same connection ID for the lower anchor and for the upper anchor.
			upper_anchor_connID_ = i;
			break;
		}
	}

	runtime_assert_string_msg( lower_anchor_connID_ !=0, "Unable to auto-assign a lower anchor for the loop to be closed.  Check that the first residue is connected to something that isn't in the loop to be closed and isn't a tail residue." );
	runtime_assert_string_msg( upper_anchor_connID_ !=0, "Unable to auto-assign an upper anchor for the loop to be closed.  Check that the final residue is connected to something that isn't in the loop to be closed and isn't a tail residue." );

	TR.Debug << "Auto-set lower_anchor_connID_ to " << lower_anchor_connID_ << " and upper_anchor_connID_ to " << upper_anchor_connID_ << "." << std::endl;
	TR.Debug.flush();

	return;
}

/// @brief Find the residue that is the anchor of the lower end of the loop that we're about to close and add it to the loop pose by a jump.
void GeneralizedKIC::addloweranchor(
	core::pose::Pose &perturbedloop_pose,
	core::pose::Pose const &pose
) {
	perturbedloop_pose.clear();

	if ( lower_anchor_connID_ == 0 ) {
		//TODO -- if the lower anchor residue isn't connected to anything but the loop.
		//For now, exit with an error.
		utility_exit_with_message("Error!  GeneralizedKIC cannot operate on a loop segment that isn't connected at both its ends to anything else (though this functionality may be added in the future).\n");
	} else {
		core::Size preceding_residue = pose.residue(loopresidues_[1]).residue_connection_partner(lower_anchor_connID_);
		core::conformation::ResidueOP rsd = pose.residue(preceding_residue).clone();
		perturbedloop_pose.append_residue_by_jump((*rsd), 1, "", "", true);
	}

	return;
}


/// @brief Add the loop geometry from the starting pose to the temporary pose used for kinematic closure.  This will build ideal geometry if build_ideal==true (NOT YET TESTED).
void GeneralizedKIC::addloopgeometry(
	core::pose::Pose &perturbedloop_pose,
	core::pose::Pose const &pose,
	bool const build_ideal,
	//core::Size const effect_on_bonded_geom,
	utility::vector1 < std::pair < core::Size, core::Size > > &residue_map
) {
	core::Size const loopsize=loopresidues_.size();

	if ( loopsize==0 ) {
		utility_exit_with_message("Error!  GeneralizedKIC addloopgeometry() function called without setting a loop to be closed.  (This shouldn't actually be possible...)\n");
	}

	core::Size lastres = perturbedloop_pose.size();

	for ( core::Size ir=1; ir<=loopsize; ++ir ) { //Loop through all of the loop residues to add.
		core::conformation::ResidueOP rsd = pose.residue(loopresidues_[ir]).clone();

		core::Size con=0, anchorres=0, anchorres_pose=0, anchorcon=0; //The connection id for this residue, the residue to which we're attaching it, and that residue's connection id.

		//Set anchorres:
		anchorres=lastres;

		//Set anchorres_pose:
		if ( ir>1 ) anchorres_pose=loopresidues_[ir-1]; else anchorres_pose=pose.residue(loopresidues_[1]).residue_connection_partner(lower_anchor_connID_);
		//TODO -- handle case where there is no anchor residue.

		//Set con:
		if ( ir==1 ) con=lower_anchor_connID_;
		else { // Find the first connection ID that links to the previous residue
			for ( core::Size i=1, imax=pose.residue(loopresidues_[ir]).n_possible_residue_connections(); i<=imax; ++i ) {
				if ( pose.residue(loopresidues_[ir]).residue_connection_partner(i)==anchorres_pose ) {
					con=i;
					break;
				}
			}
		}

		//Set anchorcon:
		for ( core::Size i=1, imax=pose.residue(anchorres_pose).n_possible_residue_connections(); i<=imax; ++i ) {
			if ( pose.residue(anchorres_pose).residue_connection_partner(i)==loopresidues_[ir] ) {
				anchorcon=i;
				break;
			}
		}

		//TR << "con=" << con << " anchorres=" << anchorres << " anchorcon=" << anchorcon << " anchorres_pose=" << anchorres_pose << std::endl; //DELETE ME!

		perturbedloop_pose.append_residue_by_bond((*rsd), build_ideal, con, anchorres, anchorcon, false, false);

		//Keep track of which residue in the temporary pose corresponds to which residue in the real pose:
		residue_map.push_back( std::pair< core::Size, core::Size >( perturbedloop_pose.size(), loopresidues_[ir] ) );
		lastres=perturbedloop_pose.size();

	}

	return;
}

/// @brief Add the tail geometry from the starting pose to the temporary pose used for kinematic closure.  This will build ideal geometry if build_ideal==true (NOT YET TESTED).
/// @details The tails are residues that are not part of the loop to be closed, but which are attached to them and which "come along for the ride" as the loop to be closed moves.
/// This function MUST be called AFTER addloopgeometry().
void GeneralizedKIC::addtailgeometry(
	core::pose::Pose &perturbedloop_pose,
	core::pose::Pose const &pose,
	bool const build_ideal,
	//core::Size const effect_on_bonded_geom,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
	utility::vector1 < std::pair < core::Size, core::Size > > &tail_residue_map
) {
	core::Size const tailsize = tailresidues_.size(); //Number of tail residues in the list.
	core::Size const loopsize = loopresidues_.size() + 2; //Number of loop residues + 2 anchor residues.
	if ( tailsize==0 ) return; //If there are no tail residues, do nothing.

	// Tail residues could be specified in any order -- not necessarily starting with what's connected to the loop.
	// For this reason, we'll start by looping through and adding the tail residues that are directly attached to the loop, then the ones that are attached to those, and so on and so forth.
	utility::vector1 <bool> tail_residue_added;
	tail_residue_added.resize( tailsize , false ); //A vector of bools for whether residues have already been added to the pose or not.
	core::Size residues_to_add = tailsize; //The number of residues we still have to add.

	while ( residues_to_add>0 ) {
		bool added_a_residue_this_round=false; //Keeps track of whether residues were added.  If a residue isn't added in a given round, it means that there exist residues in the list that aren't connected to the loop in any way.

		bool breaknow=false;
		for ( core::Size ir=1; ir<=tailsize; ++ir ) { //Loop through the list of tail residues.
			if ( tail_residue_added[ir] ) continue; //Skip ahead if we've already added this residue.

			//If this is a residue that hasn't yet been added...
			for ( core::Size jr=2, jrmax=perturbedloop_pose.size(); jr<=jrmax; ++jr ) { //Loop through the residues already added, checking for a connection to this residue.  Start at 2 to skip the first anchor.
				if ( jr==loopsize ) continue; //Skip the second anchor.
				else { //going through loop residues or tail residues
					if ( pose.residue( tailresidues_[ir] ).is_bonded( get_original_pose_rsd(jr, ( (jr<loopsize)?residue_map:tail_residue_map )  ) ) ) {
						--residues_to_add;
						tail_residue_added[ir]=true;
						//Add the residue:

						core::conformation::ResidueOP rsd = pose.residue(tailresidues_[ir]).clone();
						core::Size con=0, anchorres_pose=get_original_pose_rsd(jr, ( (jr<loopsize)?residue_map:tail_residue_map ) ), anchorcon=0; //The connection ID for this residue, the other residue's index in the original pose, and that residue's connection ID to this residue.

						//Set con:
						for ( core::Size i=1, imax=pose.residue(tailresidues_[ir]).n_possible_residue_connections(); i<=imax; ++i ) { //Loop through this residue's connections
							if ( pose.residue(tailresidues_[ir]).residue_connection_partner(i)==anchorres_pose ) {
								con=i;
								break;
							}
						}
						runtime_assert_string_msg(con!=0, "Internal error (con==0) in GeneralizedKIC::addtailgeometry().  This shouldn't happen.  Consult a developer or an exorcist.");

						//Set anchorcon:
						for ( core::Size i=1, imax=pose.residue(anchorres_pose).n_possible_residue_connections(); i<=imax; ++i ) {
							if ( pose.residue(anchorres_pose).residue_connection_partner(i)==tailresidues_[ir] ) {
								anchorcon=i;
								break;
							}
						}
						runtime_assert_string_msg(anchorcon!=0, "Internal error (anchorcon==0) in GeneralizedKIC::addtailgeometry().  This shouldn't happen.  Consult a developer or an exorcist.");

						perturbedloop_pose.append_residue_by_bond((*rsd), build_ideal, con, jr, anchorcon, false, false);

						tail_residue_map.push_back( std::pair<core::Size,core::Size>( perturbedloop_pose.size(), tailresidues_[ir] )  );

						added_a_residue_this_round=true;
						breaknow=true;
						break;
					}
				}
			} //for loop through residues already added.
			if ( breaknow ) break;
		}

		runtime_assert_string_msg(added_a_residue_this_round, "Tail residues were specified that are neither directly connected to the loop to be closed, nor connected to that loop through other tail residues.  (Check your input!)");

	}

	return;
}


/// @brief Find the residue that is the anchor of the upper end of the loop that we're about to close and add it to the loop pose by a bond.
void GeneralizedKIC::addupperanchor(
	core::pose::Pose &perturbedloop_pose,
	core::pose::Pose const &pose
) {
	core::Size const lastres = perturbedloop_pose.size();
	if ( lastres==0 ) utility_exit_with_message("Error!  Empty pose passed to GeneralizedKIC::addupperanchor.  (This shouldn't be possible.)\n");

	if ( upper_anchor_connID_ == 0 ) {
		//TODO -- if the upper anchor residue isn't connected to anything but the loop.
		//For now, exit with an error.
		utility_exit_with_message("Error!  GeneralizedKIC cannot operate on a loop segment that isn't connected at both its ends to anything else (though this functionality may be added in the future).\n");
	} else {
		core::Size following_residue = pose.residue(loopresidues_[loopresidues_.size()]).residue_connection_partner(upper_anchor_connID_);
		core::conformation::ResidueOP rsd = pose.residue(following_residue).clone();

		//Find the connection point on this residue:
		core::Size con = 0;
		for ( core::Size i=1, imax=pose.residue(following_residue).n_possible_residue_connections(); i<=imax; ++i ) {
			if ( pose.residue(following_residue).residue_connection_partner(i)==loopresidues_[loopresidues_.size()] ) {
				con = i;
				break;
			}
		}

		//Append the residue:
		perturbedloop_pose.append_residue_by_bond((*rsd), false, con, lastres, upper_anchor_connID_, false, false);
	}

	return;
}

/// @brief Do the actual kinematic closure.
/// @details Inputs are pose (the loop to be closed), original_pose (the reference pose, unchanged by operation), residue_map (the mapping of
/// residues from pose to original_pose) and tail_residue_map (the mapping of tail residues from pose to original_pose).  Output is the index
/// of the solution in the solutions_ vector.
/// @notes The selected_torsions, selected_bondangles, and selected_bondlengths vectors are output only in low-memory mode.
bool
GeneralizedKIC::doKIC(
	core::pose::Pose const &pose,
	core::pose::Pose const &original_pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
	utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map,
	core::Size &solution_index,
	utility::vector1 < core::Real > & selected_torsions,
	utility::vector1 < core::Real > & selected_bondangles,
	utility::vector1 < core::Real > & selected_bondlengths
) {
	using namespace core::id;
	using namespace numeric::kinematic_closure;
	using namespace numeric::conversions;

	//Figure out what atoms make up the chain of atoms to be closed.  This will be stored in the class variable atomlist_:
	generate_atomlist(pose, residue_map);

	//Set the pivot atoms:
	utility::vector1< core::Size > pivots; //list of 3 atoms to be used as pivots (indices)
	pick_pivots(original_pose, residue_map, pivots);

	//Vars for closure:
	utility::vector1 < utility::vector1 < utility::vector1 <core::Real> > > t_ang; //Resulting torsion angles.  This is a 3-matrix of [closure attempt #][solution #][torsion #].
	utility::vector1 < utility::vector1 < utility::vector1 <core::Real> > > b_ang;  //Resulting bond angles.  This is a 3-matrix of [closure attempt #][solution #][bond angle #].
	utility::vector1 < utility::vector1 < utility::vector1 <core::Real> > > b_len; //Resulting bond lengths.  This is a 3-matrix of [closure attempt #][solution #][bond length #].
	utility::vector1 <core::Size> nsol_for_attempt; //Number of solutions for each attempt.
	utility::vector1 < core::Real > energies_for_solution; //Energy for each solution, used only in low-memory mode.

	//Total solutions found:
	core::Size total_solution_count=0;

	for ( core::Size iattempt = 1; (n_closure_attempts_>0 ? iattempt<=n_closure_attempts_ : true); ++iattempt ) { //Loop for closure attempts
		TR.Debug << "Generalized kinematic closure attempt " << iattempt << "." << std::endl;

		t_ang.push_back( utility::vector1 < utility::vector1<core::Real> > () );
		b_ang.push_back( utility::vector1 < utility::vector1<core::Real> > () );
		b_len.push_back( utility::vector1 < utility::vector1<core::Real> > () );
		nsol_for_attempt.push_back(0);

		//Translate atomlist_ into the data vectors that the kinematic closure bridgeObjects() function uses:
		utility::vector1< utility::vector1< core::Real > > atoms; //atom xyz
		utility::vector1< core::Real > dt; //desired torsions for each atom
		utility::vector1< core::Real > da; //desired bond angle for each atom
		utility::vector1< core::Real > db; //desired bond length for each atom
		utility::vector1< core::Size > order; //use 1, 2, 3
		generate_bridgeobjects_data_from_atomlist(atoms, dt, da, db, order);

		int nsol=0; //The number of solutions found.  This has to be an int.

		//Perturb the lists of desired torsions, bond angles, and bond lengths using the perturberlist_.
		apply_perturbations(pose, original_pose, residue_map, tail_residue_map, dt, da, db);

		//Do the actual closure:
		bridgeObjects(atoms, dt, da, db, pivots, order, t_ang[iattempt], b_ang[iattempt], b_len[iattempt], nsol);

		//Filter solutions found.  This decrements nsol as solutions are eliminated, and deletes solutions from the t_ang, b_ang, and b_len vectors.
		if ( nsol>0 ) filter_solutions(original_pose, pose, residue_map, tail_residue_map, atomlist_, t_ang[iattempt], b_ang[iattempt], b_len[iattempt], nsol );

		//Flag certain solutions for deletion based on whether preselection movers fail.
		utility::vector1< bool > solutions_to_delete( nsol, false );

		//Only at this point do we actually build poses:
		for ( core::Size j=1, jmax=nsol; j<=jmax; ++j ) {
			core::pose::PoseOP curpose( original_pose.clone() ); //Clone the original pose.
			core::pose::PoseOP looppose( pose.clone() ); //Clone the loop pose.
			set_loop_pose( *looppose, atomlist_, t_ang[iattempt][j], b_ang[iattempt][j], b_len[iattempt][j]);
			copy_loop_pose_to_original( *curpose, *looppose, residue_map, tail_residue_map);
			if ( correct_polymer_dependent_atoms() ) correct_polymer_dependent_atoms_in_pose_segment( *curpose, residue_map );

			//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
			if ( attach_boinc_ghost_observer() ) {
				protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( *curpose );
				protocols::boinc::Boinc::update_graphics_current_ghost( *curpose );
				//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
				//std::cerr.flush();
			}
#endif

			//Apply preselection movers.
			if ( preselection_mover_exists() ) {
				TR << "Applying pre-selection mover to solution " << j << " from attempt " << iattempt << "." << std::endl;
				pre_selection_mover_->apply( *curpose );
				if ( pre_selection_mover_->get_last_move_status()!=protocols::moves::MS_SUCCESS ) {
					TR << "Preselection mover failed!  Continuing to next solution." << std::endl;
					--nsol;
					solutions_to_delete[j] = true;
					continue;
				}
			}
			if ( !low_memory_mode() ) {
				add_solution(curpose);
				TR << "Storing solution " << j << " from attempt " << iattempt << "." << std::endl;
			} else {
				( *( selector_->scorefunction() ) )( *curpose );
				energies_for_solution.push_back( curpose->energies().total_energy() );
			}
		}

		//Delete solutions for which the preselection mover failed:
		core::Size solcount(0);
		for ( auto it=t_ang[iattempt].begin(); it!=t_ang[iattempt].end(); ) {
			++solcount;
			if ( solutions_to_delete[solcount] ) { t_ang[iattempt].erase( it ); }
			else { ++it; }
		}
		solcount=0;
		for ( auto it=b_ang[iattempt].begin(); it!=b_ang[iattempt].end(); ) {
			++solcount;
			if ( solutions_to_delete[solcount] ) { b_ang[iattempt].erase( it ); }
			else { ++it; }
		}
		solcount=0;
		for ( auto it=b_len[iattempt].begin(); it!=b_len[iattempt].end(); ) {
			++solcount;
			if ( solutions_to_delete[solcount] ) { b_len[iattempt].erase( it ); }
			else { ++it; }
		}

		total_solution_count += static_cast<core::Size>(nsol);
		if ( !low_memory_mode() ) {
			debug_assert(total_solution_count == total_stored_solutions());
		}
		nsol_for_attempt[iattempt] = static_cast<core::Size>(nsol);

		if ( min_solution_count()!=0 && total_solution_count >= min_solution_count() ) {
			TR.Debug << "Total solutions=" << total_solution_count << " and min_solution_count=" << min_solution_count() << ".  Breaking from solution-seeking loop." << std::endl;
			break;
		}

		if ( get_ntries_before_giving_up()==iattempt && total_solution_count==0 ) {
			TR.Debug << "GenKIC is set to give up after " << get_ntries_before_giving_up() << " tries without a solution.  This is attempt " << iattempt << " and solution count is " << total_solution_count << ".  Giving up." << std::endl;
			break; //Give up now if we've attempted a certain number of tries and found nothing.
		}

	} //End loop for closure attempts

	// Apply the selector here.  This ultimately picks a single solution and sets pose (the loop pose) to the conformation for that solution, but doesn't touch the original pose.
	bool selector_success(false);
	solution_index=0;
	if ( total_solution_count>0 ) {
		solution_index=select_solution ( pose, original_pose, residue_map, tail_residue_map, atomlist_, t_ang, b_ang, b_len, nsol_for_attempt, total_solution_count, energies_for_solution );
		selector_success=(solution_index!=0);
	}

	// In low-memory mode, populate vectors for the solution's conformation:
	if ( solution_index && low_memory_mode() ) {
		core::Size counter(0);
		bool breaknow(false);
		for ( core::Size i=1, imax=t_ang.size(); i<=imax; ++i ) {
			for ( core::Size j=1, jmax=t_ang[i].size(); j<=jmax; ++j ) {
				++counter;
				if ( counter == solution_index ) {
					selected_torsions=t_ang[i][j];
					selected_bondangles=b_ang[i][j];
					selected_bondlengths=b_len[i][j];
					breaknow=true;
					break;
				}
			}
			if ( breaknow ) break;
		}
	}

	return (total_solution_count>0 && selector_success);
}


/// @brief Generate the list of atomIDs for the chain that will be closed by kinematic closure.
void GeneralizedKIC::generate_atomlist(
	core::pose::Pose const &pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
) {
	using namespace core::id;

	atomlist_.clear();

	runtime_assert_string_msg((pose.size() > 0 && residue_map.size() > 0), "Uninitialized data passed to GeneralizedKIC::generate_atomlist.");

	{ //Scope 1: start with the AtomID of the lower anchor connection point and its 2 parents:
		core::Size firstconindex = get_connection( 1, residue_map[1].first, pose );
		core::Size firstconatomindex = pose.residue(1).residue_connect_atom_index(firstconindex);
		core::Size firstconatomindex_parent = (firstconatomindex==1 ? 2 : pose.residue(1).atom_base(firstconatomindex) );
		core::Size firstconatomindex_grandparent = (firstconatomindex==1 ? 3 : pose.residue(1).atom_base(firstconatomindex_parent) );

		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(firstconatomindex_grandparent, 1), pose.residue(1).xyz(firstconatomindex_grandparent) ) );
		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(firstconatomindex_parent, 1), pose.residue(1).xyz(firstconatomindex_parent) ) );
		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(firstconatomindex, 1), pose.residue(1).xyz(firstconatomindex) ) );
	}

	for ( core::Size ir=1, irmax=residue_map.size(); ir<=irmax; ++ir ) { //Loop through all of the residues in the loop
		core::Size const curres = residue_map[ir].first;
		core::Size const prevres = (ir==1 ? 1 : residue_map[ir-1].first);
		core::Size const nextres = (ir==irmax ? irmax+2 /*The second anchor*/ : residue_map[ir+1].first);
		//TR << "curres=" << curres << " nextres=" << nextres << " prevres=" << prevres << std::endl; //DELETE ME

		core::Size const prevconindex = get_connection(curres, prevres, pose);
		core::Size const nextconindex = get_connection(curres, nextres, pose);
		core::Size const prevconatomindex = pose.residue(curres).residue_connect_atom_index(prevconindex);
		core::Size const nextconatomindex = pose.residue(curres).residue_connect_atom_index(nextconindex);

		//TR << prevconatomindex << " " << nextconatomindex << std::endl ; //DELETE ME

		utility::vector1 < core::Size > pathindices;
		get_path(prevconatomindex, nextconatomindex, pose.residue(curres), pathindices);

		runtime_assert_string_msg(pathindices.size()!=0, "In GeneralizedKIC::generate_atomlist, the atom path could not be generated (pathindices.size()==0).");

		for ( core::Size i=1, imax=pathindices.size(); i<=imax; ++i ) {
			atomlist_.push_back(
				std::pair<AtomID, numeric::xyzVector<core::Real> > (
				AtomID(pathindices[i], curres),
				pose.residue(curres).xyz(pathindices[i])
				)
			);
		}
	} //Looping through the residues of the loop

	{ //Scope 2: end with the AtomID of 3 atoms of the upper anchor connection point:
		core::Size nres = loopresidues_.size()+2;
		core::Size lastconindex = get_connection( nres, residue_map[residue_map.size()].first, pose );
		core::Size lastconatomindex = pose.residue(nres).residue_connect_atom_index(lastconindex);
		core::Size lastconatomindex_parent = (lastconatomindex==1 ? 2 : pose.residue(nres).atom_base(lastconatomindex) );
		core::Size lastconatomindex_grandparent = (lastconatomindex==1 ? 3 : pose.residue(nres).atom_base(lastconatomindex_parent) );

		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(lastconatomindex, nres), pose.residue(nres).xyz(lastconatomindex) ) );
		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(lastconatomindex_parent, nres), pose.residue(nres).xyz(lastconatomindex_parent) ) );
		atomlist_.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> >(AtomID(lastconatomindex_grandparent, nres), pose.residue(nres).xyz(lastconatomindex_grandparent) ) );
	}


	//Comment out the following -- for testing only:
	//TR << "Index\tRes\tAtom\tX\tY\tZ" << std::endl;
	//for(core::Size i=1, imax=atomlist_.size(); i<=imax; ++i) {
	// TR << i << "\t" << atomlist_[i].first.rsd() << "\t" << pose.residue(atomlist_[i].first.rsd()).atom_name(atomlist_[i].first.atomno())
	//   << "\t" << atomlist_[i].second[0] << "\t" << atomlist_[i].second[1] << "\t" << atomlist_[i].second[2] <<  std::endl;
	//}

	return;
}


/// @brief Generate the numeric::kinematic_closure::bridgeObjects data from the atomlist_ object.
void GeneralizedKIC::generate_bridgeobjects_data_from_atomlist(
	utility::vector1< utility::vector1< core::Real > > &atoms, //atom xyz
	utility::vector1< core::Real > &dt, //desired torsions for each atom
	utility::vector1< core::Real > &da, //desired bond angle for each atom
	utility::vector1< core::Real > &db, //desired bond length for each atom
	utility::vector1< core::Size > &order //use 1, 2, 3
) {
	using namespace numeric::kinematic_closure;
	core::Size const atomcount=atomlist_.size();

	//Create the "atoms" object (list of x,y,z coordinates of atoms):
	atoms.clear();
	for ( core::Size i=1; i<=atomcount; ++i ) {
		utility::vector1 <core::Real> xyz;
		xyz.resize(3);
		xyz[1]=atomlist_[i].second[0];
		xyz[2]=atomlist_[i].second[1];
		xyz[3]=atomlist_[i].second[2];
		atoms.push_back(xyz);
	}

	//Calculate dt, da, db:
	utility::vector1<utility::vector1<core::Real> > q0 (3); //Used by numeric::kinematic_closure::chainTORS
	utility::vector1<core::Real> r0 (3); //Used by numeric::kinematic_closure::chainTORS
	chainTORS( atomcount, atoms, dt, da, db, r0, q0 );

	//for(core::Size i=1, imax=atomcount; i<=imax; ++i) TR << "x=" << atoms[i][1] << " y=" << atoms[i][2] << " z=" << atoms[i][3] << " dt=" << dt[i] << " da=" << da[i] << " db=" << db[i] << std::endl; //DELETE ME

	order.resize(3); order[1]=1; order[2]=2; order[3]=3; //Order = 1,2,3

	return;
}


/// @brief Given a residue_map vector of pairs, where each pair is < residue_index_in_perturbedloop_pose, residue_index_in_original_pose >,
/// and a residue index in the original pose, return the corresponding residue index in the perturbed loop pose.
core::Size GeneralizedKIC::get_perturbedloop_rsd ( core::Size const original_pose_rsd, utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map )
{
	for ( core::Size i=1, imax=residue_map.size(); i<=imax; ++i ) {
		if ( residue_map[i].second == original_pose_rsd ) return residue_map[i].first;
	}

	utility_exit_with_message("Error in GeneralizedKIC::get_perturbedloop_rsd.");
	return 0;
}


/// @brief Pick the pivots for kinematic closure.
void GeneralizedKIC::pick_pivots(
	core::pose::Pose const &original_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < core::Size > &pivots
) {
	using namespace protocols::generalized_kinematic_closure;
	//TODO -- I THINK THE FIRST PIVOT HAS TO BE ATOM 5.  CHECK THIS.

	//runtime_assert_string_msg( pivot_1_rsd_ == residue_map[1].second, "The first pivot atom for kinematic closure must be in the first residue in the loop to be closed." );
	//runtime_assert_string_msg( pivot_3_rsd_ == residue_map[residue_map.size()].second, "The third pivot atom for kinematic closure must be in the last residue in the loop to be closed." );
	runtime_assert_string_msg( original_pose.residue(pivot_1_rsd_).has(pivot_1_atmname_), "The first pivot atom was not found in the first residue in the loop to be closed." );
	runtime_assert_string_msg( original_pose.residue(pivot_2_rsd_).has(pivot_2_atmname_), "The second pivot atom was not found in the last residue in the loop to be closed." );
	runtime_assert_string_msg( original_pose.residue(pivot_3_rsd_).has(pivot_3_atmname_), "The third pivot atom was not found in the last residue in the loop to be closed." );

	core::Size const pivot_1_atmindex = original_pose.residue(pivot_1_rsd_).atom_index(pivot_1_atmname_);
	core::Size const pivot_2_atmindex = original_pose.residue(pivot_2_rsd_).atom_index(pivot_2_atmname_);
	core::Size const pivot_3_atmindex = original_pose.residue(pivot_3_rsd_).atom_index(pivot_3_atmname_);

	pivots.clear();
	pivots.resize(3);
	pivots[1]=0; pivots[2]=0; pivots[3]=0;

	core::Size totalatoms = atomlist_.size(); //Total number of atoms in the chain to be closed

	for ( core::Size ia=5; ia<=totalatoms-4; ia++ ) { //Loop through the atom list, ignoring the atoms from the anchor residues
		core::Size curres = get_original_pose_rsd(atomlist_[ia].first.rsd(), residue_map);
		core::Size curat = atomlist_[ia].first.atomno();
		//TR << "curres=" << curres << " curat=" << curat << " original_pose.residue(curres).atom_name(curat)=" << original_pose.residue(curres).atom_name(curat) << std::endl; //DELETE ME
		if ( (curres==pivot_1_rsd_) && (curat==pivot_1_atmindex) ) pivots[1]=ia;
		if ( (curres==pivot_2_rsd_) && (curat==pivot_2_atmindex) ) pivots[2]=ia;
		if ( (curres==pivot_3_rsd_) && (curat==pivot_3_atmindex) ) pivots[3]=ia;
	}

	runtime_assert_string_msg( pivots[1] > 4 && pivots[1] < totalatoms-3, "The first pivot atom was not found in the loop to be closed!" );
	runtime_assert_string_msg( pivots[2] > 4 && pivots[2] < totalatoms-3, "The second pivot atom was not found in the loop to be closed!" );
	runtime_assert_string_msg( pivots[3] > 4 && pivots[3] < totalatoms-3, "The third pivot atom was not found in the loop to be closed!" );

	runtime_assert_string_msg( pivots[2] > pivots[1], "The second pivot cannot be before the first pivot." );
	runtime_assert_string_msg( pivots[2] - pivots[1] > 2, "There must be at least two atoms between the first and second pivot atoms." );
	runtime_assert_string_msg( pivots[3] > pivots[2], "The third pivot cannot be before the second pivot." );
	runtime_assert_string_msg( pivots[3] - pivots[2] > 2, "There must be at least two atoms between the second and third pivot atoms." );

	prune_extra_atoms(pivots);

	//TR << "pivots[1]=" << pivots[1] << " pivots[2]=" << pivots[2] << " pivots[3]=" << pivots[3] << std::endl; //DELETE ME

	return;
}

/// @brief Apply the list of perturbers (everything in perturberlist_) to alter the desired torsion, desired angle, and
///        desired bond length lists prior to calling bridgeObjects.  Note that later perturbers might overwrite earlier
///        perturbers' effects.
/// @details
///
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] original_pose -- A pose consisting of the full, original structure.
/// @param[in] residue_map -- The mapping of (residue in loop_pose, residue in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue in loop_pose, tail residue in original_pose).
/// @param[in,out] torsions -- The desired torsion angles, potentially altered or overwritten by the perturbers.
/// @param[in,out] bondangles -- The desired bond angles, potentially altered or overwritten by the perturbers.
/// @param[in,out] bondlenghts -- The desired bond lengths, potentially altered or overwritten by the perturbers.
void GeneralizedKIC::apply_perturbations(
	core::pose::Pose const &loop_pose,
	core::pose::Pose const &original_pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
	utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map,
	utility::vector1 < core::Real > &torsions,
	utility::vector1 < core::Real > &bondangles,
	utility::vector1 < core::Real > &bondlengths
) const {
	core::Size const nperturbers = perturberlist_.size();
	if ( nperturbers==0 ) {
		TR.Warning << "Warning!  No perturbers specified for GeneralizedKIC mover." << std::endl;
		return;
	}

	for ( core::Size i=1; i<=nperturbers; ++i ) { //Loop through, applying each perturber in turn.
		//TR << "Applying perturber " << i << "." << std::endl;
		if ( perturberlist_[i]->get_perturber_effect() == perturber::perturb_dihedral_bbg ) {
			perturberlist_[i]->init_bbgmover(loop_pose, residue_map);
		}
		perturberlist_[i]->apply(original_pose, loop_pose, residue_map, tail_residue_map, atomlist_, torsions, bondangles, bondlengths);
	}

	return;
}

/// @brief Apply filters to the list of solutions, and proceed to eliminate any and all that fail to pass filters.
/// @details This removes entries from the torsions, bondangles, and bondlengths lists, and decrements nsol (the number of solutions) appropriately.
/// @param[in] original_pose -- A pose consisting of the full, original structure.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue in loop_pose, residue in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue in loop_pose, tail residue in original_pose).
/// @param[in] atomlist -- The list of atomIDs in the chain that was closed, plus their xyz coordinates from the original pose.
/// @param[in,out] torsions -- The torsion angles returned by bridgeObjects, as a matrix of [solution #][torsion index].  Columns can be deleted by filters.
/// @param[in,out] bondangles -- The bond angles returned by bridgeObjects, as a matrix of [solution #][bondangle index].  Columns can be deleted by filters.
/// @param[in,out] bondlenghts -- The bond lengths returned by bridgeObjects, as a matrix of [solution #][bondlength index].  Columns can be deleted by filters.
/// @param[in,out] nsol -- The number of solutions, which can be decremented as filters delete solutions.  Note that this must be an int, not a core::Size.
void GeneralizedKIC::filter_solutions(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 <utility::vector1<core::Real> > &torsions,
	utility::vector1 <utility::vector1<core::Real> > &bondangles,
	utility::vector1 <utility::vector1<core::Real> > &bondlengths,
	int &nsol
) const {
	TR.Debug << "Applying filters to GeneralizedKIC solutions." << std::endl;  TR.Debug.flush();

	//Double-check that the value of nsol matches the size of the torsions, bondangles, and bondlenghts arrays (skipped in release mode, I think):
	assert( torsions.size()==static_cast<core::Size>(nsol));
	assert( bondangles.size()==static_cast<core::Size>(nsol));
	assert( bondlengths.size()==static_cast<core::Size>(nsol));

	if ( nsol==0 ) return;
	core::Size const nfilters = filterlist_.size(); //Number of filters
	if ( nfilters == 0 ) {
		TR.Warning << "Warning!  No filters specified for GeneralizedKIC mover." << std::endl;  TR.Warning.flush();
		return;
	}

	for ( core::Size ifilter=1; ifilter<=nfilters; ++ifilter ) { //Loop through all filters
		if ( nsol==0 ) break;
		TR.Debug << "Applying " << filterlist_[ifilter]->get_this_filter_type_name() << " filter." << std::endl;
		for ( core::Size isolution=1; isolution<=static_cast<core::Size>(nsol); ++isolution ) { //Loop through all of the solutions
			//TR << "Applying filter " << ifilter << " to solution " << isolution << "." << std::endl;  TR.flush(); //DELETE ME
			if ( !filterlist_[ifilter]->apply(original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions[isolution], bondangles[isolution], bondlengths[isolution]) ) { //If this filter returns false
				//TR << "Filter returned false.  Erasing solution " << isolution << "." << std::endl;  TR.flush(); //DELETE ME
				torsions.erase( torsions.begin()+isolution-1 );
				bondangles.erase( bondangles.begin()+isolution-1 );
				bondlengths.erase( bondlengths.begin()+isolution-1 );
				--nsol; //Decrement
				--isolution; //Decrement
				if ( nsol==0 ) break;
			} //else {
			//TR << "Filter returned true." << std::endl;  TR.flush(); //DELETE ME
			//}
		}
	}

	return;

}

/// @brief Applies the selector to choose a solution.
/// @details  If the selector could not select a solution (e.g. if the preselection mover returned failed status for every solution), this function returns 0;
/// otherwise, returns the index of the solution in the solutions_ vector.
/// @param[in,out] pose -- The loop to be closed.
/// @param[in] original_pose -- The original pose.  Can be used for reference by selectors.
/// @param[in] residue_map -- Mapping of (loop residue, original pose residue).
/// @param[in] tail_residue_map -- Mapping of (tail residue index in pose, tail residue index in original_pose).
/// @param[in] atomlist -- The list of (AtomID, original XYZ coordinates of atoms) representing the chain that was closed.
/// @param[in] torsions -- Matrix of [closure attempt #][solution #][torsion #] with torsion values for each torsion angle in the chain.  A selector will pick one solution.
/// @param[in] bondangles -- Matrix of [closure attempt #][solution #][angle #] with bond angle values for each bond angle in the chain.  A selector will pick one solution.
/// @param[in] bondlengths -- Matrix of [closure attempt #][solution #][bondlength #] with bond length for each bond in the chain.  A selector will pick one solution.
/// @param[in] nsol_for_attempt -- List of the number of solutions for each attempt.
/// @param[in] total_solutions -- Total number of solutions found.
/// @param[in] energies_for_solution -- Vectors of the energies from each attempt.  Used only in low-memory mode.
core::Size GeneralizedKIC::select_solution (
	core::pose::Pose const &pose,
	core::pose::Pose const &original_pose, //The original pose
	utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map, //mapping of (loop residue, original pose residue)
	utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map, //mapping of (tail residue index in pose, tail residue index in original_pose)
	utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, //torsions for each atom
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, //bond angle for each atom
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, //bond length for each atom
	utility::vector1 <core::Size> const &nsol_for_attempt,
	core::Size const total_solutions,
	utility::vector1 <core::Real> const &energies_for_solution
) const {
	return selector_->apply(pose, original_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles,
		bondlengths, nsol_for_attempt, total_solutions, solutions_, low_memory_mode(), energies_for_solution );
}

/// @brief Trims extra atoms from the start and end of the atom list, if the first and last pivots are not the fifth and fifth-last atoms, respectively.
///
void GeneralizedKIC::prune_extra_atoms( utility::vector1 <core::Size> &pivots )
{
	core::Size const totalatoms = atomlist_.size();
	runtime_assert_string_msg(pivots.size()==3 && pivots[1]>4 && pivots[3]<(totalatoms-3), "Internal error in GeneralizedKIC::prune_extra_atoms().  Consult a programmer or an exorcist.  This shouldn't happen.");

	core::Size extra_at_start = pivots[1]-5;
	core::Size extra_at_end = totalatoms-4-pivots[3];

	if ( extra_at_start > 0 ) {
		//TR.Debug << "Removing first " << extra_at_start << " atoms from the atom list." << std::endl ; //DELETE ME
		for ( core::Size i=1; i<=extra_at_start; ++i ) atomlist_.erase( atomlist_.begin() ); //Erase the first extra_at_start atoms
		for ( core::Size i=1; i<4; ++i ) pivots[i] -= extra_at_start; //Update the pivot indices
	}

	if ( extra_at_end > 0 ) {
		//TR.Debug << "Removing last " << extra_at_end << " atoms from the atom list." << std::endl ; //DELETE ME
		for ( core::Size i=1; i<=extra_at_end; ++i ) atomlist_.erase(atomlist_.end()-1); //Erase the last extra_at_end atoms
	}

	return;
}

/// @brief Clear the stored solution poses.
///
void
GeneralizedKIC::clear_stored_solutions() {
	solutions_.clear();
}

/// @brief Sets the mover that will be applied to all solutions that pass filters prior to applying the selector.
///
void GeneralizedKIC::set_preselection_mover ( protocols::moves::MoverOP mover )
{
	pre_selection_mover_ = mover;
	pre_selection_mover_exists_ = true;
	return;
}

/// @brief Add a solution to the solutions list.
/// @details Stores an owning pointer to a pose.  Not used in low-memory mode (see function overload).
void
GeneralizedKIC::add_solution(
	core::pose::PoseOP pose_in
) {
	debug_assert( !low_memory_mode() );
	solutions_.push_back(pose_in);
	return;
}

/// @brief Sets the pose to a particular solution from the solutions list.
/// @details This is for regular mode, so it just grabs this from the stored solution pose vector.
void
GeneralizedKIC::set_final_solution(
	core::pose::Pose &output_pose,
	core::Size const solution_index
) const {
	debug_assert( !low_memory_mode() );
	output_pose = *solutions_[solution_index];
}

/// @brief Sets the pose to a particular solution from the solutions list.
/// @details This is for low-memory mode, so it constructs the pose from DoF vectors and reapplies preselection movers.
/// The pose is the original pose (input/output), and the looppose is the loop alone (input/output).
protocols::moves::MoverStatus
GeneralizedKIC::set_final_solution(
	core::pose::Pose &pose,
	core::pose::Pose &looppose,
	utility::vector1 < std::pair < core::Size, core::Size > > const & residue_map,
	utility::vector1 < std::pair < core::Size, core::Size > > const & tail_residue_map,
	utility::vector1< core::Real > const &torsions,
	utility::vector1< core::Real > const &bondangles,
	utility::vector1< core::Real > const &bondlengths
) const {
	debug_assert( low_memory_mode() );
	set_loop_pose( looppose, atomlist_, torsions, bondangles, bondlengths);
	copy_loop_pose_to_original( pose, looppose, residue_map, tail_residue_map);
	if ( correct_polymer_dependent_atoms() ) correct_polymer_dependent_atoms_in_pose_segment( pose, residue_map );
	if ( preselection_mover_exists() ) {
		TR << "Re-applying pre-selection mover to solution." << std::endl;
		pre_selection_mover_->apply( pose );
		if ( pre_selection_mover_->get_last_move_status() != protocols::moves::MS_SUCCESS && TR.Warning.visible() ) {
			TR.Warning << "Warning: in re-applying the preselection mover, the mover failed.  Setting GeneralizedKIC mover status to failure." << std::endl;
			TR.Warning.flush();
		}
		return pre_selection_mover_->get_last_move_status();
	}
	return protocols::moves::MS_SUCCESS;
}

std::string GeneralizedKIC::get_name() const {
	return mover_name();
}

std::string GeneralizedKIC::mover_name() {
	return "GeneralizedKIC";
}

std::string gen_kic_subelement_mangler( std::string const & element_name )
{
	return "gen_kic_subelement_" + element_name + "_type";
}

std::string gen_kic_sample_cis_subelement_mangler( std::string const & element_name )
{
	return "gen_kic_sample_cis_subelement_" + element_name + "_type";
}

std::string gen_kic_add_pert_subelement_mangler( std::string const & element_name )
{
	return "gen_kic_add_pert_subelement_" + element_name + "_type";
}

std::string gen_kic_add_filter_subelement_mangler( std::string const & element_name )
{
	return "gen_kic_add_filter_subelement_" + element_name + "_type";
}

void GeneralizedKIC::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	selector::GeneralizedKICselector::define_valid_selector_name_enumeration( xsd );
	perturber::GeneralizedKICperturber::define_valid_perturber_name_enumeration( xsd );
	filter::GeneralizedKICfilter::define_valid_filter_name_enumeration( xsd );

	typedef utility::tag::XMLSchemaAttribute Attr;
	AttributeList attlist;

	attlist + Attr( "low_memory_mode", xsct_rosetta_bool, "The low_memory_mode option can be used to limit the amount of information about each solution found that is stored, in order to reduce memory consumption. In this mode, only loop degree-of-freedom values are stored prior to selection. The drawback, however, is it can be risk to use low-memory mode with a preselection mover. Any changes to the pose made by the preselection mover will not be stored, so the preselection mover will have to be re-applied to the selected pose after selection. This costs additional time. In addition, if the preselection mover has any stochastic component to its behaviour, then the second application may not produce identical results to the first. This means that what might have been the lowest-energy pose when the preselection mover was first applied could now be a relatively high-energy pose, for example. For this reason, be judicious in the use of low-memory mode. The preferred course is to stop GenKIC when N solutions have been found to limit memory usage." )
		+ Attr( "correct_polymer_dependent_atoms", xsct_rosetta_bool, "If set to true, atoms whose positions are dependent on polymer bonds, such as amino acid \"H\" and \"O\" atoms, will have their positions set to ideal coordinates after closure.  False by default." )
		+ Attr( "dont_fail_if_no_solution_found", xsct_rosetta_bool, "By default, the GeneralizedKIC mover returns failure status if it fails to find a closed solution that passes filters.  Certain usage cases require that it does not fail.  If this option is set to true, it prevents GeneralizedKIC from returning a failure status, even if it finds no solution (in which case the input pose, which may have an open chain of atoms, is returned as the output pose).  This is necessary, for example, when GeneralizedKIC is used in conjunction with a ContingentFilter." )
		+ Attr::required_attribute( "selector", "genkic_selector_name", "The GeneralizedKICselector to use to pick a solution.  Solutions may be picked randomly (\"random_selector\"), by energy (\"lowest_energy_selector\"), or by other criteria.  See the Rosetta help wiki for more information on GeneralizedKICselectors, and for the full range of options." )
		+ Attr( "selector_scorefunction", xs_string, "A scorefunction that the GeneralizedKICselector will use to select the solution to return.")
		+ Attr( "selector_kbt", xsct_real, "If the GeneralizedKICselector is set to \"boltzmann_energy_selector\", this is the Boltzmann temperature used when selecting a solution randomly, weighted by the Boltzmann probability of that solution." )
		+ Attr::attribute_w_default( "stop_if_no_solution", xsct_non_negative_integer, "If this option is used and set to a positive integer N, GeneralizedKIC stops looking for closure solutions if no solution was found in the first N attempts.", "0" )
		+ Attr( "closure_attempts", xsct_non_negative_integer, "The maximum number of times to try to find a solution.  Set this to zero to keep trying indefinitely (in which case the \"stop_when_n_solutions_found\" option should be used to prevent the mover from iterating forever)." )
		+ Attr( "stop_when_n_solutions_found", xsct_non_negative_integer, "If this option is used and set to a positive integer N, GeneralizedKIC will stop looking for closure solutions as soon as at least N solutions have been found that pass filters.  Note that a GeneralizedKICselector must be specified even if this option is set to 1, because a single attempt can yield up to 16 solutions." )
		+ Attr( "pre_selection_mover", xs_string, "As an optional alternative to restricting oneself to backbone-only score terms, GeneralizedKIC permits the user to specify a pre-selection mover defined prior to the GeneralizedKIC block. This permits side-chain moves, such as repacking or side-chain minimization, prior to scoring with the full energy function. Note that this mover will be applied to all candidate solutions passing filters, which makes this option potentially very computationally expensive!\nIf a pre-selection mover is specified, the mover may alter geometry outside of the loop to be closed, or may alter the FoldTree, in which case these alterations will carry through to the final pose returned by GeneralizedKIC. TaskOperations and MoveMaps may be used appropriately in the definition of the mover in question to restrict its effects to the loop on which GeneralizedKIC is operating if the user does not wish to alter geometry outside of this loop. If one wishes to apply more than one pre-selection mover (for example, a PackRotamersMover followed by a sidechain-only MinMover, to repack and minimize side-chains), they may be combined in a ParsedProtocol mover.\nAs a final note, a mover may return failure status, or may encapsulate filters that return failure status, in which case the solution to which the mover has been applied will be discarded." )
		+ Attr( "contingent_filter", xs_string, " A ContingentFilter can also be used to record whether the mover failed. The ContingentFilter is a specialized filter that has its value set by a mover. GeneralizedKIC can set the value of a ContingentFilter, specified using the contingent_filter flag, to true or false depending on whether the closure was successful or unsuccessful. Subsequent application of the filter, possibly at a later point, can then abort trajectories involving unsuccessful loop closure.  Note that this is a largely deprecated feature that has been retained since a few users have found it useful in some cases; for most usage cases, one can simply rely on the fact that the mover returns failure status if no closed solution is found.  It will likely be necessary to use the \"dont_fail_if_no_solution_found\" option if a ContingentFilter is used." );

	// AddResidue sub-element
	AttributeList add_residue_attributes;
	add_residue_attributes + Attr::required_attribute( "res_index", xsct_non_negative_integer, "The \"res_index\" option within an \"AddResidue\" sub-tag specifies the index, in Rosetta numbering, of the residue being added to the chain of residues to be closed." );
	XMLSchemaComplexTypeGenerator add_residue_ct_gen;
	add_residue_ct_gen.element_name( "AddResidue" )
		.complex_type_naming_func( & gen_kic_subelement_mangler )
		.add_attributes( add_residue_attributes )
		.description( "\"AddResidue\" sub-tags are used to specify the residues that make up the chain of atoms to be closed, in order of connectivity." )
		.write_complex_type_to_schema( xsd );

	{
		// AddTailResidue sub-element
		AttributeList add_tail_residue_attributes;
		add_tail_residue_attributes + Attr::required_attribute( "res_index", xsct_non_negative_integer, "The \"res_index\" option within an \"AddTailResidue\" sub-tag specifies the index, in Rosetta numbering, of the residue being added to the chain of residues constituting a \"tail\" that is intended to move with the chain to be closed." );
		XMLSchemaComplexTypeGenerator add_tail_residue_ct_gen;
		add_tail_residue_ct_gen.element_name( "AddTailResidue" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( add_tail_residue_attributes )
			.description( "Loops to be closed might have other chains of residues attached to them.  By default, GeneralizedKIC fails to respect covalent bonds to any but the end residues in the chain to be closed, and can tear these bonds apart as residues move.  Users can optionally specify a list of \"tail\" residues that are covalently attached to, and move with, the chain of residues to be closed.  For example, one might have a series of glycans attached to a serine residue in a protein loop that's being closed with GenKIC.  If these residues are specified as tail residues, then they will remain attached to the serine residue even as the loop that contains the serine residue moves." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// SetPivots sub-element
		AttributeList set_pivots_attributes;
		set_pivots_attributes
			+ Attr::required_attribute( "res1", xsct_non_negative_integer, "The first pivot residue (by Rosetta index). Currently, due to hard-coded assumptions in the kinematic closure numerical library, the first pivot must be the second atom in the chain to be closed." )
			+ Attr::required_attribute( "atom1", xs_string, "The first pivot atom (by atom name). Currently, due to hard-coded assumptions in the kinematic closure numerical library, the first pivot must be the second atom in the chain to be closed." )
			+ Attr::required_attribute( "res2", xsct_non_negative_integer, "The middle pivot residue (by Rosetta index)." )
			+ Attr::required_attribute( "atom2", xs_string, "The middle pivot atom (by atom name)." )
			+ Attr::required_attribute( "res3", xsct_non_negative_integer, "The last pivot residue (by Rosetta index). Currently, due to hard-coded assumptions in the kinematic closure numerical library, the last pivot must be the second-to-last atom in the chain to be closed." )
			+ Attr::required_attribute( "atom3", xs_string, "The last pivot atom (by atom name). Currently, due to hard-coded assumptions in the kinematic closure numerical library, the last pivot must be the second-to-last atom in the chain to be closed." );
		XMLSchemaComplexTypeGenerator set_pivots_ct_gen;
		set_pivots_ct_gen.element_name( "SetPivots" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( set_pivots_attributes )
			.description( "Pivots are atoms in the chain of atoms to be closed that are flanked by bonds whose dihedral values will be solved for analytically by the closure algorithm in order to close the loop. Currently, due to hard-coded assumptions in the kinematic closure numerical library, the first pivot must be the second atom in the chain to be closed, and the last pivot must be the second-to-last atom in the chain to be closed." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// SampleCisPeptideBond sub-element
		// It itself has its own sub-element -- AddResidue
		AttributeList sample_cis_peptide_bond_add_residue_attributes;
		sample_cis_peptide_bond_add_residue_attributes
			+ Attr::required_attribute( "index", xsct_non_negative_integer, "The index, in Rosetta numbering, of a residue whose omega value should sample 0." );
		XMLSchemaSimpleSubelementList samp_cis_add_residue_subelement;
		samp_cis_add_residue_subelement.complex_type_naming_func( & gen_kic_sample_cis_subelement_mangler )
			.add_simple_subelement( "AddResidue", sample_cis_peptide_bond_add_residue_attributes, "\"AddResidue\" sub-sub-tags are used within a \"SampleCisPeptideBond\" sub-tag to list all of the residues at which cis-peptide bonds should be sampled." );

		AttributeList sample_cis_peptide_bond_attributes;
		sample_cis_peptide_bond_attributes
			+ Attr::attribute_w_default( "cis_prob", xsct_real, "The probability of sampling a cis-peptide bond at a given position, where 1.0 means that cis is sampled 100% of the time, and 0.0 means that trans is sampled 100% of the time.", "0.1" );

		XMLSchemaComplexTypeGenerator sample_cis_peptide_bond_ct_gen;
		sample_cis_peptide_bond_ct_gen.element_name( "SampleCisPeptideBond" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( sample_cis_peptide_bond_attributes )
			.set_subelements_repeatable( samp_cis_add_residue_subelement )
			.description( "A \"SampleCisPeptideBond\" sub-tag is a shorthand for adding a sample_cis_peptide_bond GeneralizedKICperturber.  It tells the GeneralizedKIC mover to sample a cis-peptide bond at a particular residue some subset of the time." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// CloseBond subelement
		AttributeList close_bond_attributes;
		close_bond_attributes

			+ Attr::required_attribute( "res1", xsct_non_negative_integer, "The residue (Rosetta numbering) containing the first atom in the bond to be closed." )
			+ Attr::required_attribute( "atom1", xs_string, "The first atom (given as an atom name string) in the bond to be closed." )
			+ Attr::required_attribute( "res2", xsct_non_negative_integer, "The residue (Rosetta numbering) containing the second atom in the bond to be closed." )
			+ Attr::required_attribute( "atom2", xs_string, "The second atom (given as an atom name string) in the bond to be closed." )
			+ Attr::required_attribute( "bondlength", xsct_real, "The length of the bond to be closed." )
			+ Attr::required_attribute( "angle1", xsct_real, "The bond angle defined by B1, A1, A2, where A1 and A2 are the atoms in the bond to be closed and B1 is the last atom in the chain to be closed that's prior to the bond to be closed." )
			+ Attr::required_attribute( "angle2", xsct_real, "The bond angle defined by A1, A2, B2, where A1 and A2 are the atoms in the bond to be closed and B2 is the first atom in the chain to be closed that's past the bond to be closed." )
			+ Attr( "randomize_flanking_torsions", xsct_rosetta_bool, "If true, bonds flanking the bond to be closed will have their torsion values randomized.  If present, requires \"prioratom_res\" and \"followingatom_res\" attributes to be set." )
			+ Attr( "prioratom_res", xsct_non_negative_integer, "The atom preceding the bond to be closed.  This is only required if the \"randomize_flanking_torsions\" option is used." )
			+ Attr( "prioratom", xs_string, "The atom preceding the bond to be closed.  This is only required if the \"randomize_flanking_torsions\" option is used." )
			+ Attr( "followingatom_res", xsct_non_negative_integer, "The atom following the bond to be closed.  This is only required if the \"randomize_flanking_torsions\" option is used." )
			+ Attr( "followingatom", xs_string, "The atom following the bond to be closed.  This is only required if the \"randomize_flanking_torsions\" option is used." )
			+ Attr( "torsion", xsct_real, "An optional torsion value for the bond to be closed.  If not specified, the input value persists." );

		XMLSchemaComplexTypeGenerator close_bond_ct_gen;
		close_bond_ct_gen.element_name( "CloseBond" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( close_bond_attributes )
			.description( "A \"CloseBond\" sub-tag is a shorthand for set_bondangle and set_bondlength GeneralizedKICperturbers, as well as an optional pair of randomize_dihedral perturbers for the flanking dihedral angles, and an optional set_torsion perturber.  It is intended to be used to set the ideal geometry for an open bond in a loop that GeneralizedKIC should close." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// AddPerturber subelement
		AttributeList add_perturber_attributes;
		add_perturber_attributes + Attr::required_attribute( "effect", "genkic_perturber_name", "The name of the GeneralizedKICperturber (e.g. \"perturb_dihedral\", \"set_bondangle\", \"randomize_alpha_backbone_by_rama\", etc.).  See the Rosetta documentation wiki for a full list of available perturbers." )
			+ Attr::attribute_w_default( "bin_params_file", xs_string, "The \"randomize_backbone_by_bins\" and \"perturb_backbone_by_bins\" perturbers require a bin definition file, set using this option.", "ABBA" )
			+ Attr::attribute_w_default( "iterations", xsct_non_negative_integer, "The \"perturb_backbone_by_bins\" perturber requires that a number of iterations be set. If the number of iterations is set to 1, only one residue, randomly chosen from the list provided, will have its mainchain torsions perturbed. If it is set higher, the algorithm iteratively picks a residue at random and perturbs it.", "1" )
			+ Attr::attribute_w_default( "must_switch_bins", xsct_rosetta_bool, "This option is only used by the \"perturb_backbone_by_bins\" perturber.  If \"must_switch_bins\" is set to true, the chosen residue is forced into a different torsion bin; if false, it has some probability of remaining in the same torson bin, in which case its mainchain torsion values will be chosen randomly (in a biased manner, if possible) from within that bin.", "false" )
			+ Attr( "bin", xs_string, "The \"set_backbone_bin\" perturber draws random mainchain torsion values from a mainchain torsion bin defined in a bin params file.  The \"bin\" option names the bin from which the torsion values will be drawn.  Required for this perturber only." ) /*NOT an enum -- a bin name is read from a bin params file*/
			+ Attr( "custom_rama_table", xs_string, "The \"randomize_alpha_backbone_by_rama\" perturber may optionally draw random phi and psi values from one of several custom Ramachandran tables that offer either flattened or more stringent distributions.  These are listed on the Rosetta help wiki." );

		// This one has its own subelements!
		//   1. AddResidue
		AttributeList add_pert_subelement_add_residue_attrs;
		add_pert_subelement_add_residue_attrs
			+ Attr::required_attribute( "index", xsct_non_negative_integer, "The index, in Rosetta numbering, of a residue on which this perturber should act." );
		XMLSchemaComplexTypeGenerator add_pert_subelement_add_residue_ct_gen;
		add_pert_subelement_add_residue_ct_gen.element_name( "AddResidue" )
			.complex_type_naming_func( gen_kic_add_pert_subelement_mangler )
			.add_attributes( add_pert_subelement_add_residue_attrs )
			.description( "The \"AddResidue\" sub-sub-tag in an \"AddPerturber\" sub-tag permits a residue to be defined on which the perturber may act.  Use multiple \"AddResidue\" tags to define more than one residue." )
			.write_complex_type_to_schema( xsd );

		//   2. AddAtoms
		AttributeList add_pert_subelement_add_atoms_attrs;
		add_pert_subelement_add_atoms_attrs
			+ Attr::required_attribute( "res1", xsct_non_negative_integer, "The residue (specified using Rosetta numbering) containing the first atom." )
			+ Attr::required_attribute( "atom1", xs_string, "The first atom, specified as an atom name string." )
			+ Attr( "res2", xsct_non_negative_integer, "The residue (specified using Rosetta numbering) containing the second atom." )
			+ Attr( "atom2", xs_string, "The second atom, specified as an atom name string." )
			+ Attr( "res3", xsct_non_negative_integer, "The residue (specified using Rosetta numbering) containing the third atom." )
			+ Attr( "atom3", xs_string, "The third atom, specified as an atom name string." )
			+ Attr( "res4", xsct_non_negative_integer, "The residue (specified using Rosetta numbering) containing the fourth atom." )
			+ Attr( "atom4", xs_string, "The fourth atom, specified as an atom name string." );
		XMLSchemaComplexTypeGenerator add_pert_subelement_add_atoms_ct_gen;
		add_pert_subelement_add_atoms_ct_gen.element_name( "AddAtoms" )
			.complex_type_naming_func( gen_kic_add_pert_subelement_mangler )
			.add_attributes( add_pert_subelement_add_atoms_attrs )
			.description( "The \"AddAtoms\" sub-sub-tag in an \"AddPerturber\" sub-tag permits one or more atoms to be defined, on which the perturber may act.  A particular type of perturber may take a list of one, two, three, or four atoms, depending on what it acts on." )
			.write_complex_type_to_schema( xsd );

		//    3. AddValue
		AttributeList add_pert_subelement_add_value_attrs;
		add_pert_subelement_add_value_attrs
			+ Attr::required_attribute( "value", xsct_real, "The value to pass to the perturber." );

		XMLSchemaComplexTypeGenerator add_pert_subelement_add_value_ct_gen;
		add_pert_subelement_add_value_ct_gen.element_name( "AddValue" )
			.complex_type_naming_func( gen_kic_add_pert_subelement_mangler )
			.add_attributes( add_pert_subelement_add_value_attrs )
			.description( "The \"AddValue\" sub-sub-tag in an \"AddPerturber\" sub-tag permits a value to be set for those perturbers that take a value.  (An example of this is the \"set_dihedral\" perturber, which takes a value for the dihedral angle to set.)  Use multiple \"AddValue\" tags to specify more than one value." )
			.write_complex_type_to_schema( xsd );

		// Now the simple subelement list for the AddPerturber ct
		XMLSchemaSimpleSubelementList add_pert_ssl;
		add_pert_ssl
			.add_already_defined_subelement( "AddResidue", gen_kic_add_pert_subelement_mangler )
			.add_already_defined_subelement( "AddAtoms", gen_kic_add_pert_subelement_mangler )
			.add_already_defined_subelement( "AddValue", gen_kic_add_pert_subelement_mangler );
		XMLSchemaComplexTypeGenerator add_pert_ct_gen;
		add_pert_ct_gen.element_name( "AddPerturber" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( add_perturber_attributes )
			.set_subelements_repeatable( add_pert_ssl )
			.description( "The \"AddPerturber\" sub-tag is used to add and configure GeneralizedKICperturbers.  GeneralizedKICperturbers alter the chain to be closed in some way prior to kinematic closure, allowing definition of a desired conformation or conformational sampling. They can only act on the chain to be closed, and have no effect on tail residues or on any other part of the input structure. Perturbers are applied in the order that they are defined. Different perturbers may alter the same degrees of freedom, sequentially." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// AddFilter subelement
		AttributeList add_filter_attributes;
		add_filter_attributes
			+ Attr::required_attribute( "type", "genkic_filter_name", "The name of the GeneralizedKICfilter to apply (e.g. \"backbone_bin\", \"rama_prepro_check\", \"atom_pair_distance\", etc.).  See the Rosetta help wiki for a full list of avaliable filters." )
			+ Attr( "residue", xsct_non_negative_integer, "A parameter required by several filters, including \"backbone_bin\", \"alpha_aa_rama_check\", and \"rama_prepro_check\": the index, in Rosetta numbering, of the residue that we're examining with this filter." )
			+ Attr::attribute_w_default( "bin_params_file", xs_string, "A parameter specific for the \"backbone_bin\" filter: a filename for a file defining mainchain torsion bins.", "ABBA" )
			+ Attr( "bin", xs_string, "A parameter specific for the \"backbone_bin\" filter: the mainchain torsion bin in which an amino acid residue must lie.  If the residue is not in this bin, the solution is discarded." )
			+ Attr::attribute_w_default( "rama_cutoff_energy", xsct_real, "A parameter specific for the \"alpha_aa_rama_check\" and \"rama_prepro_check\" filters: the value for the rama or rama_prepro score terms above which solutions are discarded.", "0.3" );

		// This one has its own subelements!
		//   1. AddFilterParameterReal
		AttributeList add_filter_subelement_add_filter_param_real_attrs;
		add_filter_subelement_add_filter_param_real_attrs

			+ Attr::required_attribute( "name", xs_string, "The name of the real-valued paramter.  (For example, \"distance\" for the \"atom_pair_distance\" filter.)" )
			+ Attr::required_attribute( "value", xsct_real, "The value of the real-valued parameter." );

		XMLSchemaComplexTypeGenerator add_filter_subelement_add_filter_param_real_ct_gen;
		add_filter_subelement_add_filter_param_real_ct_gen.element_name( "AddFilterParameterReal" )
			.complex_type_naming_func( gen_kic_add_filter_subelement_mangler )
			.add_attributes( add_filter_subelement_add_filter_param_real_attrs )
			.description( "An \"AddFilterParameterReal\" sub-sub-tag in an \"AddFilter\" tag sets a real-valued setting for a filter.  For example, the \"atom_pair_distance\" filter takes a real-valued parameter, called \"distance\", for the cutoff separation between two atoms used to discard solutions." )
			.write_complex_type_to_schema( xsd );

		//   2. AddFilterParameterInteger
		AttributeList add_filter_subelement_add_filter_param_int_attrs;
		add_filter_subelement_add_filter_param_int_attrs
			+ Attr::required_attribute( "name", xs_string, "The name of the integer parameter.  (For example, \"res1\" or \"res2\" for the \"atom_pair_distance\" filter.)" )
			+ Attr::required_attribute( "value", xsct_non_negative_integer, "The value of the integer parameter." );
		XMLSchemaComplexTypeGenerator add_filter_subelement_add_filter_param_int_ct_gen;
		add_filter_subelement_add_filter_param_int_ct_gen.element_name( "AddFilterParameterInteger" )
			.complex_type_naming_func( gen_kic_add_filter_subelement_mangler )
			.add_attributes( add_filter_subelement_add_filter_param_int_attrs )
			.description( "An \"AddFilterParameterInteger\" sub-sub-tag in an \"AddFilter\" tag sets an integer setting for a filter.  For example, the \"atom_pair_distance\" filter takes two integer parameters, called \"res1\" and \"res2\", for the indices of the two residues containing the two atoms whose separation is to be measured." )
			.write_complex_type_to_schema( xsd );

		//   3. AddFilterParameterBoolean
		AttributeList add_filter_subelement_add_filter_param_bool_attrs;
		add_filter_subelement_add_filter_param_bool_attrs
			+ Attr::required_attribute( "name", xs_string, "The name of the Boolean paramter.  (For example, \"greater_than\" for the \"atom_pair_distance\" filter.)" )
			+ Attr::required_attribute( "value", xsct_rosetta_bool, "The value (\"true\" or \"false\") of the Boolean parameter." );
		XMLSchemaComplexTypeGenerator add_filter_subelement_add_filter_param_bool_ct_gen;
		add_filter_subelement_add_filter_param_bool_ct_gen.element_name( "AddFilterParameterBoolean" )
			.complex_type_naming_func( gen_kic_add_filter_subelement_mangler )
			.add_attributes( add_filter_subelement_add_filter_param_bool_attrs )
			.description( "An \"AddFilterParameterBoolean\" sub-sub-tag in an \"AddFilter\" tag sets a Boolean setting for a filter.  For example, the \"atom_pair_distance\" filter has a Boolean parameter called \"greater_than\" which, when set to \"true\" selects for interatomic separation greater than a cutoff rather than less than the cutoff." )
			.write_complex_type_to_schema( xsd );

		//   4. AddFilterParameterString
		AttributeList add_filter_subelement_add_filter_param_string_attrs;
		add_filter_subelement_add_filter_param_string_attrs
			+ Attr::required_attribute( "name", xs_string, "The name of the string paramter.  (For example, \"atom1\" or \"atom2\" for the \"atom_pair_distance\" filter.)" )
			+ Attr::required_attribute( "value", xs_string, "The value of the string parameter." );
		XMLSchemaComplexTypeGenerator add_filter_subelement_add_filter_param_string_ct_gen;
		add_filter_subelement_add_filter_param_string_ct_gen.element_name( "AddFilterParameterString" )
			.complex_type_naming_func( gen_kic_add_filter_subelement_mangler )
			.add_attributes( add_filter_subelement_add_filter_param_string_attrs )
			.description( "An \"AddFilterParameterString\" sub-sub-tag in an \"AddFilter\" tag sets a string setting for a filter.  For example, the \"atom_pair_distance\" filter takes two string parameters, called \"atom1\" and \"atom2\", for the names of the two atoms whose separation is to be measured." )
			.write_complex_type_to_schema( xsd );

		// Now the simple subelement list for the AddFilter ct
		XMLSchemaSimpleSubelementList add_filter_ssl;
		add_filter_ssl
			.add_already_defined_subelement( "AddFilterParameterReal", gen_kic_add_filter_subelement_mangler )
			.add_already_defined_subelement( "AddFilterParameterInteger", gen_kic_add_filter_subelement_mangler )
			.add_already_defined_subelement( "AddFilterParameterBoolean", gen_kic_add_filter_subelement_mangler )
			.add_already_defined_subelement( "AddFilterParameterString", gen_kic_add_filter_subelement_mangler );
		XMLSchemaComplexTypeGenerator add_filter_ct_gen;
		add_filter_ct_gen.element_name( "AddFilter" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( add_filter_attributes )
			.set_subelements_repeatable( add_filter_ssl )
			.description( "An \"AddFilter\" sub-tag adds a GeneralizedKICfilter, which discards solutions that do not pass certain criteria.  Generally, GeneralizedKICfilters are written to be considerably faster than conventional Rosetta filters, since they usually do not act on full poses.  They are applied before pre-selection movers." )
			.write_complex_type_to_schema( xsd );
	}

	{
		// AddAtomPairDistanceFilter
		AttributeList atpair_filter_attrs;
		atpair_filter_attrs
			+ Attr::required_attribute( "res1", xsct_non_negative_integer, "The residue (Rosetta index) containing the first atom in the pair whose separation we will be measuring." )
			+ Attr::required_attribute( "atom1", xs_string, "The first atom in the pair whose separation we will be measuring, specified as an atom name string." )
			+ Attr::required_attribute( "res2", xsct_non_negative_integer, "The residue (Rosetta index) containing the second atom in the pair whose separation we will be measuring." )
			+ Attr::required_attribute( "atom2", xs_string, "The second atom in the pair whose separation we will be measuring, specified as an atom name string." )
			+ Attr::required_attribute( "distance", xsct_real, "The cutoff inter-atomic distance for discarding solutions." )
			+ Attr( "greater_than", xsct_rosetta_bool, "If \"greater_than\" is set to true, the filter will discard any solution for which the atoms are not separated by AT LEAST the cutoff distance.  If it is false (the default), then any solution for which the atoms are separated by more than the cutoff distance is discarded." );

		XMLSchemaComplexTypeGenerator atpair_filter_ct_gen;
		atpair_filter_ct_gen.element_name( "AddAtomPairDistanceFilter" )
			.complex_type_naming_func( & gen_kic_subelement_mangler )
			.add_attributes( atpair_filter_attrs )
			.description( "The \"AddAtomPairDistanceFilter\" sub-tag is a shorthand for adding an \"atom_pair_distance\" GeneralizedKICFilter.  It discards solutions if pairs of atoms are separated by more than a cutoff distance." )
			.write_complex_type_to_schema( xsd );
	}

	// OK! Now for the GeneralizedKIC complex type itself!
	XMLSchemaSimpleSubelementList genkic_ssl;
	genkic_ssl.add_already_defined_subelement( "AddResidue", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "AddTailResidue", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "SetPivots", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "SampleCisPeptideBond", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "CloseBond", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "AddPerturber", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "AddFilter", & gen_kic_subelement_mangler )
		.add_already_defined_subelement( "AddAtomPairDistanceFilter", & gen_kic_subelement_mangler );
	XMLSchemaComplexTypeGenerator genkic_ct_gen;
	genkic_ct_gen.element_name( mover_name() )
		.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.set_subelements_repeatable( genkic_ssl )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.description( "The GeneralizedKIC mover allows a user to define a chain of atoms whose conformation is to be sampled or altered subject to the constraint that it remain closed, or that it be closed if it starts out open.  The mover has myriad options for controlling sampling and selection of a solution, and can invoke other movers to apply to each closed solution that it finds.  It has been written to be fully general, and can operate on chains of atoms that pass through backbones, side-chains, ligands, etc.  No assumptions about alpha-amino acid backbones are hard-coded, so GeneralizedKIC should be fully compatible with any non-canonical entity." )
		.write_complex_type_to_schema( xsd );
}


std::string GeneralizedKICCreator::keyname() const {
	return GeneralizedKIC::mover_name();
}

protocols::moves::MoverOP
GeneralizedKICCreator::create_mover() const {
	return protocols::moves::MoverOP( new GeneralizedKIC );
}

void GeneralizedKICCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GeneralizedKIC::provide_xml_schema( xsd );
}



} //namespace generalized_kinematic_closure
} //namespace protocols
