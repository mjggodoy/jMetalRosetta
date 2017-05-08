// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::stepwise::modeler;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.options.StepWiseMonteCarloOptions" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace options {

//Constructor
StepWiseMonteCarloOptions::StepWiseMonteCarloOptions():
	verbose_scores_( false ),
	integration_test_mode_( false ),
	force_centroid_interaction_( true ),
	sampler_max_centroid_distance_( 0.0 ),
	use_phenix_geo_( true ),
	skip_deletions_( false ),
	erraser_( true ),
	cycles_( 500 ),
	add_proposal_density_factor_( 1.0 ),
	minimize_single_res_frequency_( 0.0 ),
	just_min_after_mutation_frequency_( 0.5 ),
	temperature_( 1.0 ),
	allow_split_off_( true ),
	virtual_sugar_keep_base_fixed_( false ),
	virtual_sugar_do_minimize_( false /*true*/ ),
	make_movie_( false ),
	sampler_perform_phosphate_pack_( true ),
	force_phosphate_instantiation_( true ),
	virtualize_packable_moieties_in_screening_pose_( true ),
	rebuild_bulge_mode_( false ),
	tether_jump_( true ),
	local_redock_only_( true ),
	skip_coord_constraints_( false ),
	filter_native_big_bins_( false ),
	allow_virtual_o2prime_hydrogens_( false ),
	allow_virtual_side_chains_( false ),
	n_sample_( 18 ),
	protein_prepack_( false ),
	o2prime_legacy_mode_( false ),
	recover_low_( false ),
	enumerate_( false ),
	preminimize_( false ),
	skip_preminimize_( false ),
	new_move_selector_( false ),
	test_all_moves_( false ),
	save_times_( false ),
	use_precomputed_library_( true ),
	minimize_after_delete_( true ),
	use_first_jump_for_submotif_( false ),
	checkpoint_( false ),
	checkpointing_frequency_( 0 )
{
	StepWiseBasicOptions::initialize_variables();
	set_silent_file( "default.out" );
}

//Destructor
StepWiseMonteCarloOptions::~StepWiseMonteCarloOptions()
{}

/// @brief copy constructor
StepWiseMonteCarloOptions::StepWiseMonteCarloOptions( StepWiseMonteCarloOptions const & src ) :
	ResourceOptions( src ),
	StepWiseBasicOptions( src ),
	StepWiseMoveSelectorOptions( src )
{
	*this = src;
}

/// @brief clone the options
StepWiseMonteCarloOptionsOP
StepWiseMonteCarloOptions::clone() const
{
	return StepWiseMonteCarloOptionsOP( new StepWiseMonteCarloOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
StepWiseMonteCarloOptions::initialize_from_command_line() {

	StepWiseBasicOptions::initialize_from_command_line();
	StepWiseMoveSelectorOptions::initialize_from_command_line();

	set_silent_file( option[ OptionKeys::out::file::silent ]() );
	set_verbose_scores( option[ OptionKeys::stepwise::monte_carlo::verbose_scores ]() );
	set_integration_test_mode( option[ OptionKeys::stepwise::rna::integration_test ]() );
	set_skip_deletions( option[ OptionKeys::stepwise::monte_carlo::skip_deletions ]() );
	set_num_pose_minimize( option[ OptionKeys::stepwise::num_pose_minimize ]() );
	set_erraser( option[ OptionKeys::stepwise::rna::erraser ]() );
	set_cycles( option[ OptionKeys::stepwise::monte_carlo::cycles ]() );
	set_add_proposal_density_factor( option[ OptionKeys::stepwise::monte_carlo::add_proposal_density_factor ]() );
	set_minimize_single_res_frequency( option[ OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency ]() );
	set_just_min_after_mutation_frequency( option[ OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency ]() );

	set_allow_split_off( option[ OptionKeys::stepwise::monte_carlo::allow_split_off ]() );
	set_temperature( option[ OptionKeys::stepwise::monte_carlo::temperature ]() );
	set_virtual_sugar_keep_base_fixed( option[ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed ]() );
	if ( option[ OptionKeys::stepwise::rna::virtual_sugar_do_minimize ].user() ) {
		set_virtual_sugar_do_minimize( option[ OptionKeys::stepwise::rna::virtual_sugar_do_minimize ]() );
	}
	force_centroid_interaction_ = true; // note default is different from stepwise enumeration
	if ( option[ OptionKeys::stepwise::rna::force_centroid_interaction ].user() ) set_force_centroid_interaction( option[ OptionKeys::stepwise::rna::force_centroid_interaction ]() );
	sampler_max_centroid_distance_ = option[ OptionKeys::stepwise::rna::sampler_max_centroid_distance ]();
	set_use_phenix_geo(  option[ OptionKeys::rna::corrected_geo ] );
	set_make_movie( option[ OptionKeys::stepwise::monte_carlo::make_movie ]() );
	sampler_perform_phosphate_pack_ = option[ OptionKeys::stepwise::rna::sampler_perform_phosphate_pack ]();
	force_phosphate_instantiation_ = option[ OptionKeys::stepwise::rna::force_phosphate_instantiation ]();
	rebuild_bulge_mode_ = option[ OptionKeys::stepwise::rna::rebuild_bulge_mode ]();
	tether_jump_ = option[ OptionKeys::stepwise::rna::tether_jump ]();
	o2prime_legacy_mode_ = option[ OptionKeys::stepwise::rna::o2prime_legacy_mode ]();
	recover_low_ = option[ OptionKeys::stepwise::monte_carlo::recover_low ]();
	enumerate_ = option[ OptionKeys::stepwise::enumerate ]();
	preminimize_ = option[ OptionKeys::stepwise::preminimize ]();
	skip_preminimize_ = option[ OptionKeys::stepwise::skip_preminimize ]();
	new_move_selector_ = option[ OptionKeys::stepwise::new_move_selector ]();
	test_all_moves_ = option[ OptionKeys::stepwise::test_all_moves ]();
	save_times_ = option[ OptionKeys::out::save_times ]();
	use_precomputed_library_ = option[ OptionKeys::stepwise::monte_carlo::use_precomputed_library ]();
	local_redock_only_ = option[ OptionKeys::stepwise::monte_carlo::local_redock_only ]();
	skip_coord_constraints_ = option[ OptionKeys::stepwise::protein::skip_coord_constraints ]();
	filter_native_big_bins_ = option[ OptionKeys::stepwise::protein::filter_native_big_bins ]();
	allow_virtual_o2prime_hydrogens_ = option[ OptionKeys::stepwise::rna::allow_virtual_o2prime_hydrogens ]();
	allow_virtual_side_chains_ = option[ OptionKeys::stepwise::protein::allow_virtual_side_chains ]();
	n_sample_ = option[ OptionKeys::stepwise::protein::n_sample ]();
	protein_prepack_ = option[ OptionKeys::stepwise::protein::protein_prepack ]();
	virtualize_packable_moieties_in_screening_pose_ = option[ OptionKeys::stepwise::virtualize_packable_moieties_in_screening_pose ]();
	use_first_jump_for_submotif_ = option[ OptionKeys::stepwise::monte_carlo::use_first_jump_for_submotif ]();
	checkpoint_ = option[ OptionKeys::stepwise::monte_carlo::checkpointing_frequency ]() != 0;
	checkpointing_frequency_ = option[ OptionKeys::stepwise::monte_carlo::checkpointing_frequency ]();

	if ( test_all_moves_ ) {
		set_num_random_samples( 0 );
		set_minimize_after_delete( false );
	}
}


//////////////////////////////////////////////////////////////////////////////////////
// Wait, would be smarter to just take advantage of clone(), right?
protocols::stepwise::modeler::options::StepWiseModelerOptionsOP
StepWiseMonteCarloOptions::setup_modeler_options() const{

	using namespace protocols::stepwise::options;
	using namespace protocols::stepwise::modeler::options;
	StepWiseModelerOptionsOP options( new StepWiseModelerOptions );

	// copy over StepWiseBasicOptions -- both StepWiseMonteCarloOptions and StepWiseModelerOptions inherit from it.
	StepWiseBasicOptions & options_basic_new( *options );
	StepWiseBasicOptions const & options_basic_this( *this );
	options_basic_new = options_basic_this;

	// general
	options->set_num_pose_minimize( 0 );
	options->set_num_random_samples( num_random_samples() );
	if ( num_pose_minimize() > 0 ) options->set_num_pose_minimize( num_pose_minimize() );
	options->set_rmsd_screen( rmsd_screen() );
	options->set_VDW_rep_screen_info( VDW_rep_screen_info() ); // caleb when you see a conflict, just get rid of this line.
	options->set_atr_rep_screen( atr_rep_screen() );

	// Might be smarter to inherit the following options from the same options parent class,
	// instead of setting manually here, where we could forget something.
	options->set_choose_random( !enumerate_ );
	options->set_virtualize_packable_moieties_in_screening_pose( virtualize_packable_moieties_in_screening_pose() );
	if ( preminimize_ ) options->set_use_packer_instead_of_rotamer_trials( true );  // to match SWA for proteins.
	if ( enumerate_ && !erraser() ) {
		options->set_output_minimized_pose_list( true ); // to match legacy SWA (turn off for erraser)
	}
	// protein-specific
	options->set_skip_coord_constraints( skip_coord_constraints() );
	options->set_filter_native_big_bins( filter_native_big_bins() );
	options->set_allow_virtual_o2prime_hydrogens( allow_virtual_o2prime_hydrogens() ); // RNA
	options->set_allow_virtual_side_chains( allow_virtual_side_chains() ); // protein
	options->set_n_sample( n_sample() );
	options->set_prepack( protein_prepack() );
	options->set_expand_loop_takeoff( true );

	// rna-specific
	options->set_integration_test_mode( integration_test_mode() );
	options->set_force_centroid_interaction( force_centroid_interaction() );
	options->set_sampler_max_centroid_distance( sampler_max_centroid_distance() );
	options->set_use_phenix_geo( use_phenix_geo() );
	options->set_kic_modeler_if_relevant( erraser() );
	options->set_virtual_sugar_keep_base_fixed( virtual_sugar_keep_base_fixed() );
	options->set_virtual_sugar_do_minimize( virtual_sugar_do_minimize() );
	if ( rebuild_bulge_mode_ ) options->set_force_centroid_interaction( false );
	options->set_tether_jump( tether_jump() );
	options->set_o2prime_legacy_mode( o2prime_legacy_mode() );
	options->set_sampler_perform_phosphate_pack( sampler_perform_phosphate_pack() );
	options->set_force_phosphate_instantiation( force_phosphate_instantiation() );
	options->set_allow_rebuild_bulge_mode( false );
	options->set_virtual_sugar_do_screens( false );

	return options;
}


} //options
} //monte_carlo
} //stepwise
} //protocols
