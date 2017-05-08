// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.cc
/// @brief Application-level code for the simple_cycpep_predict app.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifdef BOINC
#include <utility/io/izstream.hh>
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Package Headers
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/select/residue_selector/PhiSelector.hh>
#include <core/select/residue_selector/BinSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/aa_composition/AddCompositionConstraintMover.hh>
#include <protocols/aa_composition/ClearCompositionConstraintsMover.hh>
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.hh>
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/denovo_design/movers/FastDesign.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

//Constraints
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//Disulfides
#include <protocols/cyclic_peptide/TryDisulfPermutations.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>

//N-methylation
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

//TBMB
#include <protocols/cyclic_peptide/ThreefoldLinkerMover.hh>
#include <protocols/cyclic_peptide/threefold_linker/TBMB_Helper.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//numeric headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/angle.functions.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <cstdio>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication" );

/// @brief Register the set of options that this application uses (for the help menu).
///
void
protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sequence_file                        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cyclic_permutations                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_rama_filter                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rama_cutoff                          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_final_hbonds                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::default_rama_sampling_table          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rama_sampling_table_by_res           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::checkpoint_file                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rand_checkpoint_file                 );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_disulfides                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_prerelax               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_postrelax              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedral_perturbation );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::filter_oversaturated_hbond_acceptors );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::hbond_acceptor_energy_cutoff         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::design_peptide                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::prohibit_D_at_negative_phi           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::prohibit_L_at_positive_phi           );
	option.add_relevant( basic::options::OptionKeys::score::aa_composition_setup_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::do_not_count_adjacent_res_hbonds     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::n_methyl_positions                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_positions                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_TBMB_filters                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_sidechain_distance_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_constraints_energy_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_perturbation        );
#ifdef USEMPI //Options that are only needed in the MPI version:
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_processes_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_batchsize_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_sort_by                          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_choose_highest                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_output_fraction                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_stop_after_time                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_lambda );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_kbt );
#endif
	return;
}


namespace protocols {
namespace cyclic_peptide_predict {


/// @brief Constructor
/// @details If allow_file_read is true, initialization triggers reads
/// from the filesystem.
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication(
	bool const allow_file_read
) :
	my_rank_(0),
	already_completed_job_count_(0),
	scorefxn_(),
	suppress_checkpoints_(false),
	silent_out_(false),
	silentlist_out_(false),
	silentlist_(nullptr),
	summarylist_(nullptr),
	native_pose_(),
	out_filename_("S_"),
	out_scorefilename_("default.sc"),
	sequence_file_(""),
	sequence_string_(""),
	sequence_length_(0),
	genkic_closure_attempts_(1),
	genkic_min_solution_count_(1),
	cyclic_permutations_(true),
	use_rama_filter_(true),
	rama_cutoff_(0.3),
	high_hbond_weight_multiplier_(10),
	min_genkic_hbonds_(3.0),
	min_final_hbonds_(0.0),
	total_energy_cutoff_(0.0),
	use_total_energy_cutoff_(false),
	hbond_energy_cutoff_(-0.25),
	fast_relax_rounds_(3),
	count_sc_hbonds_(false),
	native_exists_(false),
	native_filename_(""),
	nstruct_(1),
	checkpoint_job_identifier_(""),
	checkpoint_filename_("checkpoint.txt"),
	default_rama_table_type_( core::scoring::unknown_ramatable_type ),
	rama_table_type_by_res_(),
	rand_checkpoint_file_("rng.state.gz"),
	try_all_disulfides_(false),
	disulf_energy_cutoff_prerelax_(15.0),
	disulf_energy_cutoff_postrelax_(0.5),
	user_set_alpha_dihedrals_(),
	user_set_dihedral_perturbation_(0.0),
	filter_oversaturated_hbond_acceptors_(true),
	oversaturated_hbond_cutoff_energy_(-0.1),
	sample_cis_pro_(true),
	sample_cis_pro_frequency_(0.3),
	design_peptide_(false),
	design_filename_(""),
	prevent_design_file_read_( !allow_file_read ),
	allowed_canonicals_by_position_(),
	allowed_noncanonicals_by_position_(),
	prohibit_D_at_negative_phi_(true),
	prohibit_L_at_positive_phi_(true),
	use_aa_comp_(false),
	L_alpha_comp_file_exists_(false),
	D_alpha_comp_file_exists_(false),
	L_beta_comp_file_exists_(false),
	D_beta_comp_file_exists_(false),
	comp_file_contents_L_alpha_(""),
	comp_file_contents_D_alpha_(""),
	comp_file_contents_L_beta_(""),
	comp_file_contents_D_beta_(""),
	abba_bins_(""),
	do_not_count_adjacent_res_hbonds_(true),
	angle_relax_rounds_(0),
	angle_length_relax_rounds_(0),
	cartesian_relax_rounds_(0),
	use_rama_prepro_for_sampling_(true),
	n_methyl_positions_(),
	tbmb_positions_(),
	link_all_cys_with_tbmb_(false),
	use_tbmb_filters_(true),
	tbmb_sidechain_distance_filter_multiplier_(1.0),
	tbmb_constraints_energy_filter_multiplier_(1.0),
	required_symmetry_repeats_(1),
	required_symmetry_mirroring_(false),
	required_symmetry_angle_threshold_(10.0),
	required_symmetry_perturbation_(0.0)
	//TODO -- initialize variables here.
{
	initialize_from_options();
}


/// @brief Explicit virtual destructor.
///
SimpleCycpepPredictApplication::~SimpleCycpepPredictApplication() = default;


/// @brief Explicit copy constructor.
///
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication( SimpleCycpepPredictApplication const &src ) :
	my_rank_(src.my_rank_),
	already_completed_job_count_( src.already_completed_job_count_ ),
	scorefxn_(), //Cloned below
	suppress_checkpoints_(src.suppress_checkpoints_),
	silent_out_(src.silent_out_),
	silentlist_out_(src.silentlist_out_),
	silentlist_(src.silentlist_),
	summarylist_(src.summarylist_),
	native_pose_(src.native_pose_),
	out_filename_(src.out_filename_),
	out_scorefilename_(src.out_scorefilename_),
	sequence_file_(src.sequence_file_),
	sequence_string_(src.sequence_string_),
	sequence_length_(src.sequence_length_),
	genkic_closure_attempts_(src.genkic_closure_attempts_),
	genkic_min_solution_count_(src.genkic_min_solution_count_),
	cyclic_permutations_(src.cyclic_permutations_),
	use_rama_filter_(src.use_rama_filter_),
	rama_cutoff_(src.rama_cutoff_),
	high_hbond_weight_multiplier_(src.high_hbond_weight_multiplier_),
	min_genkic_hbonds_(src.min_genkic_hbonds_),
	min_final_hbonds_(src.min_final_hbonds_),
	total_energy_cutoff_(src.total_energy_cutoff_),
	use_total_energy_cutoff_(src.use_total_energy_cutoff_),
	hbond_energy_cutoff_(src.hbond_energy_cutoff_),
	fast_relax_rounds_(src.fast_relax_rounds_),
	count_sc_hbonds_(src.count_sc_hbonds_),
	native_exists_(src.native_exists_),
	native_filename_(src.native_filename_),
	nstruct_(src.nstruct_),
	checkpoint_job_identifier_(src.checkpoint_job_identifier_),
	checkpoint_filename_(src.checkpoint_filename_),
	default_rama_table_type_( src.default_rama_table_type_ ),
	rama_table_type_by_res_( src.rama_table_type_by_res_ ),
	rand_checkpoint_file_(src.rand_checkpoint_file_),
	try_all_disulfides_(src.try_all_disulfides_),
	disulf_energy_cutoff_prerelax_(src.disulf_energy_cutoff_prerelax_),
	disulf_energy_cutoff_postrelax_(src.disulf_energy_cutoff_postrelax_),
	user_set_alpha_dihedrals_(src.user_set_alpha_dihedrals_),
	user_set_dihedral_perturbation_(src.user_set_dihedral_perturbation_),
	filter_oversaturated_hbond_acceptors_(src.filter_oversaturated_hbond_acceptors_),
	oversaturated_hbond_cutoff_energy_(src.oversaturated_hbond_cutoff_energy_),
	sample_cis_pro_(src.sample_cis_pro_),
	sample_cis_pro_frequency_(src.sample_cis_pro_frequency_),
	design_peptide_(src.design_peptide_),
	design_filename_(src.design_filename_),
	prevent_design_file_read_(src.prevent_design_file_read_),
	allowed_canonicals_by_position_(src.allowed_canonicals_by_position_),
	allowed_noncanonicals_by_position_(src.allowed_noncanonicals_by_position_),
	prohibit_D_at_negative_phi_(src.prohibit_D_at_negative_phi_),
	prohibit_L_at_positive_phi_(src.prohibit_L_at_positive_phi_),
	use_aa_comp_(src.use_aa_comp_),
	L_alpha_comp_file_exists_(src.L_alpha_comp_file_exists_),
	D_alpha_comp_file_exists_(src.D_alpha_comp_file_exists_),
	L_beta_comp_file_exists_(src.L_beta_comp_file_exists_),
	D_beta_comp_file_exists_(src.D_beta_comp_file_exists_),
	comp_file_contents_L_alpha_(src.comp_file_contents_L_alpha_),
	comp_file_contents_D_alpha_(src.comp_file_contents_D_alpha_),
	comp_file_contents_L_beta_(src.comp_file_contents_L_beta_),
	comp_file_contents_D_beta_(src.comp_file_contents_D_beta_),
	abba_bins_(src.abba_bins_),
	do_not_count_adjacent_res_hbonds_(src.do_not_count_adjacent_res_hbonds_),
	angle_relax_rounds_(src.angle_relax_rounds_),
	angle_length_relax_rounds_(src.angle_length_relax_rounds_),
	cartesian_relax_rounds_(src.cartesian_relax_rounds_),
	use_rama_prepro_for_sampling_(src.use_rama_prepro_for_sampling_),
	n_methyl_positions_(src.n_methyl_positions_),
	tbmb_positions_(src.tbmb_positions_),
	link_all_cys_with_tbmb_(src.link_all_cys_with_tbmb_),
	use_tbmb_filters_(src.use_tbmb_filters_),
	tbmb_sidechain_distance_filter_multiplier_(src.tbmb_sidechain_distance_filter_multiplier_),
	tbmb_constraints_energy_filter_multiplier_(src.tbmb_constraints_energy_filter_multiplier_),
	required_symmetry_repeats_(src.required_symmetry_repeats_),
	required_symmetry_mirroring_(src.required_symmetry_mirroring_),
	required_symmetry_angle_threshold_(src.required_symmetry_angle_threshold_),
	required_symmetry_perturbation_(src.required_symmetry_perturbation_)
	//TODO -- copy variables here.
{
	if ( src.scorefxn_ ) scorefxn_ = (src.scorefxn_)->clone();
}

/// @brief Initialize the application.
/// @details Initializes using the option system.
void
SimpleCycpepPredictApplication::initialize_from_options(
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::io;

	//Initial checks:
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user(), "Error in simple_cycpep_predict app: the user MUST provide a sequence file using the \"-cyclic_peptide:sequence_file\" flag." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() >= 0, "Error in simple_cycpep_predict app: the number of GeneralizedKIC closure attempts (\"-cyclic_peptide:genkic_closure_attempts\" flag) cannot be negative.  (Note also that setting this to zero is risky, since GenKIC will continue to seek solutions until the minimum number of solutions is reached.)" );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() >= 0, "Error in simple_cycpep_predict app: the minimum number of GenKIC solutions (\"-cyclic_peptide:genkic_min_solution_count\" flag) cannot be negative.  (Note also that setting this to zero means no minimum.)" );
	runtime_assert_string_msg( !(option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() == 0 && option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() == 0), "Error in simple_cycpep_predict app: both the \"-cyclic_peptide:genkic_closure_attempts\" and \"-cyclic_peptide:genkic_min_solution_count\" flags were set to zero.  This would result in GenKIC looping infinitely." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds during GenKIC steps (\"-cyclic_peptide:min_genkic_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds after relaxation steps (\"-cyclic_peptide:min_final_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( !( option[out::file::silent].user() && option[out::file::o].user() ), "Error in simple_cycpep_predict app: either silent file output (\"-out:file:silent\" flag) or PDB output (\"-out:file:o\") output may be used, but not both." );

	//Copy options to private member variables:
	sequence_file_ = option[basic::options::OptionKeys::cyclic_peptide::sequence_file]();
	genkic_closure_attempts_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() );
	genkic_min_solution_count_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() );
	cyclic_permutations_ = option[basic::options::OptionKeys::cyclic_peptide::cyclic_permutations]();
	use_rama_filter_ = option[basic::options::OptionKeys::cyclic_peptide::use_rama_filter]();
	rama_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::rama_cutoff]() );
	high_hbond_weight_multiplier_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier]() );
	min_genkic_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() );
	min_final_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() );
	if ( option[basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff].user() ) {
		set_total_energy_cutoff( option[basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff]() );
	}
	hbond_energy_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff]() );
	fast_relax_rounds_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds]() );
	count_sc_hbonds_ = option[basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds]();
	try_all_disulfides_ = option[basic::options::OptionKeys::cyclic_peptide::require_disulfides].user();
	disulf_energy_cutoff_prerelax_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_prerelax]() );
	disulf_energy_cutoff_postrelax_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_postrelax]() );

	//Get the scorefunction:
	if ( !prevent_design_file_read_ ) scorefxn_ = core::scoring::get_score_function(); //Reads from file.  Don't use in MPI mode.

	//Read in the comp files, if any (in non-MPI mode):
	if ( !prevent_design_file_read_ ) {
		if ( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_alpha_, option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file](), false );
			L_alpha_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_alpha_, option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file](), false );
			D_alpha_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_beta_, option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file](), false );
			L_beta_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_beta_, option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file](), false );
			D_beta_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file]() << "." << std::endl;
		}
		read_file_into_string( abba_bins_, "protocol_data/generalizedKIC/bin_params/ABBA.bin_params", true );
		TR << "Loaded ABBA.bin_params." << std::endl;
	}

	//Get the native, if it exists:
	if ( option[in::file::native].user() ) {
		native_exists_ = true;
		native_filename_ = option[in::file::native]();
	} else {
		native_exists_ = false;
		native_filename_ = "";
	}

	//Set up output file names:
	out_filename_ = "S_";
	out_scorefilename_ = "default.sc";
	if ( option[out::file::silent].user() ) {
		out_filename_=option[out::file::silent]();
		silent_out_=true;
		option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary"); //Force binary file output.
	} else if ( option[out::file::o].user() ) {
		out_filename_=option[out::file::o]();
		silent_out_=false;
	}
	if ( option[out::file::scorefile].user() ) { out_scorefilename_=option[out::file::scorefile](); }

	//Figure out number of structures to try to generate:
	if ( option[out::nstruct].user() ) {
		runtime_assert_string_msg( option[out::nstruct]() > 0, "Error in simple_cycpep_predict app: the \"-out:nstruct\" flag's value cannot be less than 1." );
		nstruct_ = static_cast<core::Size>(option[out::nstruct]());
	} else { nstruct_ = 1; }

	checkpoint_job_identifier_ = option[basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier]();
	checkpoint_filename_ = option[basic::options::OptionKeys::cyclic_peptide::checkpoint_file]();
	rand_checkpoint_file_ = option[basic::options::OptionKeys::cyclic_peptide::rand_checkpoint_file]();

	//Figure out what custom Ramachandran tables we're using, if any:
	set_default_rama_table_type( option[basic::options::OptionKeys::cyclic_peptide::default_rama_sampling_table]() );
	set_rama_table_type_by_res( option[basic::options::OptionKeys::cyclic_peptide::rama_sampling_table_by_res]() );

	//If the user has set particular dihedrals, read them in
	if ( option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals].user() ) {
		utility::vector1 < core::Real > const user_set_dihedrals( option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals]() );
		runtime_assert_string_msg( user_set_dihedrals.size() % 4 == 0, "Error in simple_cycpep_predict app: the \"-user_set_alpha_dihedrals\" option must be followed by one or more groups of four numbers.  Each group must consist of a sequence position, then phi/psi/omega values." );
		for ( core::Size i=1, imax=user_set_dihedrals.size(); i<=imax; i+=4 ) {
			core::Size const seqpos( static_cast< core::Size >( user_set_dihedrals[i] ) );
			runtime_assert_string_msg( user_set_alpha_dihedrals_.count(seqpos) == 0, "Error in simple_cycpep_predict app: a residue index was specified more than once with the \"-user_set_alpha_dihedrals\" option." );
			utility::vector1< core::Real > phipsiomegavect;
			phipsiomegavect.reserve(3);
			phipsiomegavect.push_back( user_set_dihedrals[i+1] );
			phipsiomegavect.push_back( user_set_dihedrals[i+2] );
			phipsiomegavect.push_back( user_set_dihedrals[i+3] );
			user_set_alpha_dihedrals_[seqpos] = phipsiomegavect;
		}
	}
	user_set_dihedral_perturbation_ = option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedral_perturbation]();

	filter_oversaturated_hbond_acceptors_ = option[basic::options::OptionKeys::cyclic_peptide::filter_oversaturated_hbond_acceptors]();
	oversaturated_hbond_cutoff_energy_ = option[basic::options::OptionKeys::cyclic_peptide::hbond_acceptor_energy_cutoff]();

	//Options related to sampling cis prolines:
	if ( option[basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency].user() ) { //Turn on cis proline sampling iff the user specifies it.
		set_sample_cis_pro_frequency( option[basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency]() );
	}

	//Options related to design:
	design_peptide_ = option[basic::options::OptionKeys::cyclic_peptide::design_peptide]();
	if ( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position].user() ) {
		design_filename_ = option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position]();
	} else {
		prevent_design_file_read_ = true;
	}
	prohibit_D_at_negative_phi_ = option[basic::options::OptionKeys::cyclic_peptide::prohibit_D_at_negative_phi]();
	prohibit_L_at_positive_phi_ = option[basic::options::OptionKeys::cyclic_peptide::prohibit_L_at_positive_phi]();

	//Read design options from file:
	if ( !prevent_design_file_read_ ) {
		read_peptide_design_file( design_filename_, allowed_canonicals_by_position_, allowed_noncanonicals_by_position_ );
		prevent_design_file_read_ = true; //Prevents re-read if reinitialized.
	}

	// If we're designing and the user has specified a comp file, turn on aa_composition for design steps.
	if ( option[basic::options::OptionKeys::cyclic_peptide::design_peptide]() &&
			( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ||
			option[basic::options::OptionKeys::score::aa_composition_setup_file].user()
			)
			) {
		TR << "One or more .comp files have been provided.  The app will turn on the aa_composition score term during design steps." << std::endl;
		use_aa_comp_ = true;
	}

	do_not_count_adjacent_res_hbonds_ = option[basic::options::OptionKeys::cyclic_peptide::do_not_count_adjacent_res_hbonds]();

	if ( option[basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds].user() ) {
		set_angle_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds]() ) );
	}
	if ( option[basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds].user() ) {
		set_angle_length_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds]() ) );
	}
	if ( option[basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds].user() ) {
		set_cartesian_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds]() ) );
	}

	if ( option[basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling].user() ) {
		set_use_rama_prepro_for_sampling( !option[basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling]() );
	}

	//Store the N-methylated positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions].user() ) {
		core::Size const npos( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]().size() );
		n_methyl_positions_.resize( npos );
		for ( core::Size i=1; i<=npos; ++i ) {
			runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]()[i] > 0, "Error in simple_cycpep_predict app: The N-methylated positions must all have indices greater than zero." );
			n_methyl_positions_[i] = static_cast< core::Size >( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]()[i] );
		}
	}

	//Store the TBMB positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions].user() ) {
		runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB].user(), "Error in simple_cycpep_predict application: The \"-TBMB_positions\" flag and the \"-link_all_cys_with_TBMB\" flag cannot be used together." );
		core::Size const ntbmbres(option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions]().size());
		runtime_assert_string_msg( ntbmbres > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TBMB_positions\" commandline option must be followed by a list of residues to link with 1,3,5-tris(bromomethyl)benzene." );
		runtime_assert_string_msg( ntbmbres % 3 == 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TBMB_positions\" commandline option must be followed by a list of residues, where the number of residues in the list is a multiple of three.  Groups of three residues will be linked with 1,3,5-tris(bromomethyl)benzene." );
		core::Size count(0);
		tbmb_positions_.resize(ntbmbres / 3);
		for ( core::Size i=1; i<=ntbmbres / 3; ++i ) {
			utility::vector1 <core::Size> innervect(3);
			for ( core::Size j=1; j<=3; ++j ) {
				++count;
				innervect[j] = option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions]()[count];
			}
			tbmb_positions_[i] = innervect;
		}
	}
	use_tbmb_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_TBMB_filters]();
	tbmb_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TBMB_sidechain_distance_filter_multiplier]();
	tbmb_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TBMB_constraints_energy_filter_multiplier]();
	link_all_cys_with_tbmb_ = option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB]();

	//Options for symmetric sampling:
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:require_symmetry_repeats\" flag must be provided with a positive value." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold]() > 0.0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:require_symmetry_angle_threshold\" flag must be provided with a positive value." );
	if ( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring]() ) {
		runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats].user() && option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() > 1 && option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() % 2 == 0,
			"Error in simple_cycpep_predict application: If the \"-cyclic_peptide:require_symmetry_mirroring\" option is used, then the \"-cyclic_peptide:require_symmetry_repeats\" option must be provided, must be set greater than 1, and must be set to a value divisible by 2." );
	}
	required_symmetry_repeats_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]();
	required_symmetry_mirroring_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring]();
	required_symmetry_angle_threshold_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold]();
	required_symmetry_perturbation_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_perturbation]();

	return;
} //initialize_from_options()

/// @brief Sets the default scorefunction to use.
/// @details The scorefunction is cloned.  The high-hbond version is constructed
/// from this one.  If necessary, the aa_composition score term will be turned on
/// in that one; it needn't be turned on in this one.
void
SimpleCycpepPredictApplication::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	debug_assert( sfxn_in );
	scorefxn_ = sfxn_in->clone();
}

/// @brief Allows external code to provide a native, so that the SimpleCycpepPredictApplication doesn't have to read
/// directly from disk.
void
SimpleCycpepPredictApplication::set_native(
	core::pose::PoseCOP native
) {
	runtime_assert(native); //Can't be NULL.
	native_exists_=true;
	core::pose::PoseOP native_pose_copy( native->clone() );
	set_up_native( native_pose_copy, 0);
	native_pose_ = native_pose_copy;
}

/// @brief Allows external code to provide a sequence, so that the SimpleCycpepPredictApplication doesn't have to read
/// directly from disk.
void
SimpleCycpepPredictApplication::set_sequence(
	std::string const &seq
) {
	runtime_assert( seq != "" ); //Can't be empty string.
	sequence_string_ = seq;
}

/// @brief Allows external code to set the allowed residues by position, so that this needn't be read directly
/// from disk.
void
SimpleCycpepPredictApplication::set_allowed_residues_by_position (
	std::map< core::Size, utility::vector1< std::string > > const &allowed_canonicals,
	std::map< core::Size, utility::vector1< std::string > > const &allowed_noncanonicals
) {
	allowed_canonicals_by_position_ = allowed_canonicals;
	allowed_noncanonicals_by_position_ = allowed_noncanonicals;
}

/// @brief Allows external code to specify that output should be appended to a list of SilentStructureOPs, so that the
/// SimpleCycpepPredictApplication doesn't have to write directly to disk.
void
SimpleCycpepPredictApplication::set_silentstructure_outputlist(
	utility::vector1 < core::io::silent::SilentStructOP > * silentlist,
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > * summarylist
) {
	runtime_assert( silentlist && summarylist );
	silentlist_ = silentlist;
	summarylist_ = summarylist;
	silentlist_out_ = true;
	silent_out_ = false;
}

/// @brief Allows external code to suppress checkpointing, to prevent direct file I/O from disk.
/// @details Useful on Blue Gene.
void
SimpleCycpepPredictApplication::set_suppress_checkpoints(
	bool const suppress_checkpoints
) {
	suppress_checkpoints_ = suppress_checkpoints;
}

/// @brief If called by MPI code, the rank of the current process can be stored here.
/// @details Used for output of job summaries.
void
SimpleCycpepPredictApplication::set_my_rank(
	int const rank_in
) {
	my_rank_ = rank_in;
}

/// @brief Set the number of jobs that this process has already completed.
///
void
SimpleCycpepPredictApplication::set_already_completed_job_count(
	core::Size const count_in
) {
	already_completed_job_count_ = count_in;
}

/// @brief Allows external code to override the number of structures that this should generate (otherwise
/// set by options system.
void
SimpleCycpepPredictApplication::set_nstruct(
	core::Size const nstruct_in
) {
	runtime_assert( nstruct_in > 0 );
	nstruct_ = nstruct_in;
}

/// @brief Allows external code to set the file contents for the L-alpha aa_composition file.
///
void
SimpleCycpepPredictApplication::set_L_alpha_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_L_alpha_ = contents_in;
	L_alpha_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the D-alpha aa_composition file.
///
void
SimpleCycpepPredictApplication::set_D_alpha_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_D_alpha_ = contents_in;
	D_alpha_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the L-beta aa_composition file.
///
void
SimpleCycpepPredictApplication::set_L_beta_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_L_beta_ = contents_in;
	L_beta_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the D-beta aa_composition file.
///
void
SimpleCycpepPredictApplication::set_D_beta_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_D_beta_ = contents_in;
	D_beta_comp_file_exists_=true;
}

/// @brief Allows external code to set the ABBA bin parameters without having to read the
/// bin params file directly from disk.
void
SimpleCycpepPredictApplication::set_abba_bins_binfile_contents(
	std::string const &contents_in
) {
	abba_bins_ = contents_in;
}

/// @brief Set the frequency with which we sample cis proline.
/// @details Implicitly sets sample_cis_pro_ to "true" if freq_in is not 0.0, "false" if it is.
void
SimpleCycpepPredictApplication::set_sample_cis_pro_frequency(
	core::Real const &freq_in
) {
	runtime_assert_string_msg( 0.0 <= freq_in && freq_in <= 1.0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_sample_cis_pro_frequency(): The frequency must be between 0 and 1." );
	sample_cis_pro_frequency_ = freq_in;
	if ( freq_in > 0.0 ) {
		sample_cis_pro_ = true;
	} else {
		sample_cis_pro_ = false;
		TR << "Disabling cis-proline sampling." << std::endl;
	}
}

/// @brief Set cis proline sampling OFF.
///
void
SimpleCycpepPredictApplication::disable_cis_pro_sampling() {
	sample_cis_pro_ = false;
	sample_cis_pro_frequency_ = 0.0;
}

/// @brief Set the total energy cutoff.
/// @details Also sets use_total_energy_cutoff_ to 'true'.
void
SimpleCycpepPredictApplication::set_total_energy_cutoff(
	core::Real const &value_in
) {
	total_energy_cutoff_ = value_in;
	use_total_energy_cutoff_ = true;
}

/// @brief Sets use_total_energy_cutoff_ to 'false'.
///
void
SimpleCycpepPredictApplication::disable_total_energy_cutoff() {
	use_total_energy_cutoff_ = false;
}

/// @brief Set the number of rounds of relaxation with flexible
/// bond angles.
void
SimpleCycpepPredictApplication::set_angle_relax_rounds(
	core::Size const rounds_in
) {
	angle_relax_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of relaxation with flexible
/// bond angles and bond lengths.
void
SimpleCycpepPredictApplication::set_angle_length_relax_rounds(
	core::Size const rounds_in
) {
	angle_length_relax_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of Cartesian relaxation.
///
void
SimpleCycpepPredictApplication::set_cartesian_relax_rounds(
	core::Size const rounds_in
) {
	cartesian_relax_rounds_ = rounds_in;
}

/// @brief Set whether we're using RamaPrePro tables for sampling.
/// @details Setting this to "false" lets us use classic rama tables.  True by default.
void
SimpleCycpepPredictApplication::set_use_rama_prepro_for_sampling(
	bool const setting
) {
	use_rama_prepro_for_sampling_ = setting;
}



/// @brief Actually run the application.
/// @details The initialize_from_options() function must be called before calling this.  (Called by default constructor.)
void
SimpleCycpepPredictApplication::run() const {

	//Get the scorefunction:
	debug_assert( scorefxn_ );
	core::scoring::ScoreFunctionOP sfxn_default( scorefxn_ );
	//Create a scorefunction variant with upweighted backbone hbond terms:
	core::scoring::ScoreFunctionOP sfxn_highhbond( sfxn_default->clone() );
	sfxn_highhbond->set_weight( core::scoring::hbond_lr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_lr_bb) ); //Upweight the long-range backbone hbonds
	sfxn_highhbond->set_weight( core::scoring::hbond_sr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_sr_bb) ); //Upweight the short-range backbone hbonds
	//Turn on aa_composition in this scorefunction if we're doing design and .comp files have been provided.
	if ( sfxn_highhbond->get_weight( core::scoring::aa_composition ) == 0.0 ) { sfxn_highhbond->set_weight( core::scoring::aa_composition, 1.0 ); }
	//Create variants of the above two scorefunctions with constraint weights turned on:
	core::scoring::ScoreFunctionOP sfxn_default_cst( sfxn_default->clone() ); //Will NOT have aa_compostion turned on.
	core::scoring::ScoreFunctionOP sfxn_highhbond_cst( sfxn_highhbond->clone() ); //Will have aa_composition turned on if we're doing design and .comp files are provided.
	if ( sfxn_default->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
	core::scoring::ScoreFunctionOP sfxn_highhbond_cst_cart, sfxn_default_cst_cart;
	if ( angle_relax_rounds() > 0 || angle_length_relax_rounds() > 0 || cartesian_relax_rounds() > 0 ) {
		sfxn_highhbond_cst_cart = sfxn_highhbond_cst->clone();
		debug_assert( sfxn_highhbond_cst_cart );
		if ( sfxn_highhbond_cst_cart->get_weight( core::scoring::cart_bonded ) == 0.0 ) {
			TR << "Activating cart_bonded (with weight=0.5) for high-hbonds scorefunction Cartesian variant, and de-activating pro_close term." << std::endl;
			sfxn_highhbond_cst_cart->set_weight( core::scoring::cart_bonded, 0.5 );
			sfxn_highhbond_cst_cart->set_weight( core::scoring::pro_close, 0.0 );
		}
		sfxn_default_cst_cart = sfxn_default_cst->clone();
		debug_assert( sfxn_default_cst_cart );
		if ( sfxn_default_cst_cart->get_weight( core::scoring::cart_bonded ) == 0.0 ) {
			TR << "Activating cart_bonded (with weight=0.5) for default scorefunction Cartesian variant, and de-activating pro_close term." << std::endl;
			sfxn_default_cst_cart->set_weight( core::scoring::cart_bonded, 0.5 );
			sfxn_default_cst_cart->set_weight( core::scoring::pro_close, 0.0 );
		}
	}

	//Get the sequence that we're considering:
	utility::vector1 < std::string > resnames;
	read_sequence( sequence_file_, resnames );
	sequence_length_ = resnames.size(); //Store the number of residues in the sequence, excluding crosslinkers.

	//Check that, if we're enforcing symmetry, the sequence length is an integer multiple of the number of symmetry repeats:
	if ( required_symmetry_repeats_ > 1 ) {
		runtime_assert_string_msg(
			sequence_length_ % required_symmetry_repeats_ == 0,
			"Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::run(): Symmetry has been specified by the user, but the number of residues in the peptide is not an integral multiple of the number of symmetry repeats."
		);
	}

	//Check that, if the link_all_cys_with_TBMB flag is used, the sequence has exactly three cysteine residues:
	if ( link_all_cys_with_tbmb_ ) {
		debug_assert( tbmb_positions_.size() == 0 ); //Should be true
		utility::vector1 < core::Size > cys_positions;
		for ( core::Size i=1; i<=sequence_length(); ++i ) {
			if ( !resnames[i].compare("CYS") || !resnames[i].compare("DCYS") ) cys_positions.push_back(i);
		}
		runtime_assert_string_msg( cys_positions.size() == 3, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::run(): The \"-cyclic_peptide:link_all_cys_with_TBMB\" flag was used, but the sequence does not contain exactly three CYS/DCYS residues." );
		tbmb_positions_.push_back( cys_positions );
	}

	//Get the native sequence that we will compare to.
	core::pose::PoseOP native_pose;
	if ( native_exists_ ) {
		if ( native_pose_ ) {
			native_pose=native_pose_->clone();
		} else {
			native_pose=core::pose::PoseOP(new core::pose::Pose);
			TR << "Importing native structure from " << native_filename_ << "." << std::endl;
			import_and_set_up_native ( native_filename_, native_pose, resnames.size() );
		}
#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( *native_pose );
#endif
	} else {
		TR << "No native structure specified by the user.  No RMSD values will be calculated." << std::endl;
	}

	//Set up a filter for total number of hbonds:
	protocols::filters::CombinedFilterOP total_hbond( new protocols::filters::CombinedFilter );
	set_up_hbond_filter( total_hbond, resnames.size(), sfxn_default, static_cast<core::Real>( min_genkic_hbonds_ ) );

	//Get the checkpoint information:
	core::Size success_count(0);
	core::Size curstruct(0);
	initialize_checkpointing( curstruct, success_count );

	//EVERYTHING ABOVE THIS POINT IS DONE ONCE PER PROGRAM EXECUTION.
	++curstruct;
	for ( core::Size irepeat=curstruct, irepeat_max=nstruct_; irepeat<=irepeat_max; ++irepeat ) { //Loop nstruct times
#ifdef BOINC
		{ //Increment the model count for BOINC.
			protocols::boinc::BoincSharedMemory* shmem = protocols::boinc::Boinc::get_shmem();
			shmem->model_count = shmem->model_count + 1;
		}
#endif

		//Cyclic permutation of sequence.
		core::Size cyclic_offset(0);
		utility::vector1 < std::string > resnames_copy;
		if ( cyclic_permutations_ ) {
			cyclic_offset = do_cyclic_permutation( resnames, resnames_copy );
		} else {
			resnames_copy = resnames;
		}
		runtime_assert(cyclic_offset < resnames_copy.size() ); //Should be true.

		//Create the pose:
		core::pose::PoseOP pose( new core::pose::Pose );
		build_polymer(pose, resnames_copy);

		//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
		protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
		set_up_termini_mover( termini, pose );
		termini->apply(*pose);

		//Add cyclic constraints:
		add_cyclic_constraints(pose);

		//Set all omega values to 180 and randomize mainchain torsions:
		set_mainchain_torsions(pose, cyclic_offset);

		//Add N-methylation:
		add_n_methylation( pose, cyclic_offset );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( *pose );
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( *pose );
#endif

		//Do the kinematic closure:
		bool const success( genkic_close(pose, sfxn_highhbond_cst, sfxn_highhbond_cst_cart, sfxn_default, total_hbond, cyclic_offset) );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
#endif

		if ( !success ) {
			TR << "Closure failed.";
			if ( irepeat < irepeat_max ) {
				TR << "  Continuing to next job." << std::endl;
				checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
				//Increment total jobs and check whether it's time to quit.
				if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
			} else {
				TR << std::endl;
			}
			TR.flush();
			continue;
		}

		//If we reach here, then closure was successful.  Time to relax the pose.

		TR << "Closure successful." << std::endl;

		if ( fast_relax_rounds_ > 0 ) do_final_fastrelax( pose, sfxn_default_cst, fast_relax_rounds_, false, false, false );
		if ( angle_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, angle_relax_rounds(), true, false, false );
		if ( angle_length_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, angle_length_relax_rounds(), true, true, false );
		if ( cartesian_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, cartesian_relax_rounds(), false, false, true );
		if ( angle_relax_rounds() > 0 || angle_length_relax_rounds() > 0 || cartesian_relax_rounds() > 0 ) {
			do_final_fastrelax( pose, sfxn_default_cst, 1, false, false, false ); //Do one more round of regular FastRelax if we've done any Cartesian, just to make sure we're in a pro_close minimum.
		}

		//If we're filtering by symmetry, do so here a final time:
		if ( required_symmetry_repeats_ > 1 ) {
			protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter3( new protocols::cyclic_peptide::CycpepSymmetryFilter );
			symmfilter3->set_symm_repeats( required_symmetry_repeats_ );
			symmfilter3->set_mirror_symm( required_symmetry_mirroring_ );
			symmfilter3->set_angle_threshold( required_symmetry_angle_threshold_ );
			core::select::residue_selector::ResidueIndexSelectorOP iselector( new core::select::residue_selector::ResidueIndexSelector );
			std::stringstream pep_indices("");
			pep_indices << "1-" << sequence_length();
			iselector->set_index( pep_indices.str() );
			symmfilter3->set_selector( iselector );
			if ( !symmfilter3->apply( *pose ) ) {
				TR << "Final symmetry filter passes.  This peptide has c" << required_symmetry_repeats_ << (required_symmetry_mirroring_ ? "/m " : " " ) << "symmetry." << std::endl;
			} else {
				TR << "Final symmetry filter failed.  This peptide lost c" << required_symmetry_repeats_ << (required_symmetry_mirroring_ ? "/m " : " " ) << "symmetry during the final relaxation." << std::endl;
				if ( irepeat < irepeat_max ) {
					TR << "Continuing to next job." << std::endl;
					checkpoint( irepeat, success_count ) ;
#ifdef BOINC
					//Increment total jobs and check whether it's time to quit.
					if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
				}
				continue;
			}
		}

		//pose->dump_pdb( "TEMP.pdb" ); //DELETE ME!!!

		//Undo the cyclic permutation in anticipation of re-aligning to the native:
		depermute( pose, cyclic_offset );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
#endif

		core::Real native_rmsd(0.0);

		//Score the pose before output:
		(*sfxn_default)(*pose);

		//Filter based on total energy:
		if ( use_total_energy_cutoff_ && pose->energies().total_energy() > total_energy_cutoff_ ) {
			TR << "Total final pose energy is " << pose->energies().total_energy() << ", which is greater than the cutoff of " << total_energy_cutoff_ << ".  Failing job." << std::endl;
			TR.flush();
			checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
			//Increment total jobs and check whether it's time to quit.
			if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
			continue;
		}

		//Re-filter based on number of Hbonds (using option[min_final_hbonds]()):
		core::Real const final_hbonds( total_hbond->compute( *pose ) );
		if ( final_hbonds > -1.0*min_final_hbonds_ ) {
			TR << "Final hbond count is " << -1.0*final_hbonds << ", which is less than the minimum.  Failing job." << std::endl;
			TR.flush();
			checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
			//Increment total jobs and check whether it's time to quit.
			if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
			continue;
		}

		++success_count; //Increment the count of number of successes.

		if ( native_pose ) {
			native_rmsd = align_and_calculate_rmsd(pose, native_pose);
		}

		core::Size const cis_peptide_bonds( count_cis_peptide_bonds( pose ) );

		TR << "Result\tRMSD\tEnergy\tHbonds\tCisPepBonds" << std::endl;
		TR << irepeat << "\t";
		if ( native_pose ) { TR << native_rmsd; }
		else { TR << "--"; }
		TR << "\t" << pose->energies().total_energy() << "\t" << -1.0*final_hbonds << "\t" << cis_peptide_bonds << std::endl;

		if ( silent_out_ || silentlist_out_ ) { //Writing directly to silent file or to a list of silent file data OPs
			core::io::silent::SilentFileOptions opts;
			core::io::silent::SilentStructOP ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts ) );
			char tag[512];
			if ( my_rank_ > 0 ) {
				sprintf(tag, "result_proc%04lu_%04lu", static_cast<unsigned long>(my_rank_), static_cast<unsigned long>(irepeat+already_completed_job_count_) );
			} else {
				sprintf(tag, "result_%04lu", static_cast<unsigned long>(irepeat) );
			}
			ss->fill_struct( *pose, std::string(tag) );
			if ( native_pose ) ss->add_energy( "RMSD", native_rmsd ); //Add the RMSD to the energy to be written out in the silent file.
			ss->add_energy( "HBOND_COUNT", -1.0*final_hbonds ); //Add the hbond count to be written out in the silent file.
			ss->add_energy( "CIS_PEPTIDE_BOND_COUNT", cis_peptide_bonds ); //Add the cis-peptide bond count to be written out in the silent file.
#ifdef BOINC_GRAPHICS
			protocols::boinc::Boinc::update_graphics_current( *pose );
			protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
			protocols::boinc::Boinc::update_graphics_last_accepted( *pose, pose->energies().total_energy() );
			protocols::boinc::Boinc::update_graphics_low_energy( *pose, pose->energies().total_energy() );
#endif
			if ( silent_out_ ) {
				core::io::silent::SilentFileOptions opts;
				core::io::silent::SilentFileDataOP silent_file (new core::io::silent::SilentFileData( opts ) );
				silent_file->set_filename( out_filename_ );
				silent_file->write_silent_struct( *ss, out_filename_ );
			}
			if ( silentlist_out_ ) {
				silentlist_->push_back(ss);
				core::Size curjob( summarylist_->size() + 1 );
				summarylist_->push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary( my_rank_, curjob, pose->energies().total_energy(), (native_pose ? native_rmsd : 0), static_cast< core::Size >( std::round(-1.0*final_hbonds) ), cis_peptide_bonds ) ) );
			}
		} else { //if pdb output
			char outstring[512];
			sprintf(outstring, "%s%04lu.pdb", out_filename_.c_str(), static_cast<unsigned long>(irepeat) );
			pose->dump_scored_pdb( std::string(outstring), *sfxn_default );
		}

		TR.flush();
		checkpoint( irepeat, success_count ); //This job has been attempted and has succeeded; don't repeat it.

#ifdef BOINC
		//Increment total jobs and check whether it's time to quit.
		if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
	} //Looping through nstruct

	TR << nstruct_ << " jobs attempted.  " << success_count << " jobs returned solutions." << std::endl;
	TR.flush();

	end_checkpointing(); //Delete the checkpoint file at this point, since all jobs have completed.
	return;
}

/// @brief Count the number of cis-peptide bonds in the pose.
/// @details Counts as cis if in the range (-90,90].
core::Size
SimpleCycpepPredictApplication::count_cis_peptide_bonds(
	core::pose::PoseCOP pose
) const {
	core::Size count(0);
	for ( core::Size i=1, imax=sequence_length(); i<=imax; ++i ) {
		core::Real const omegaval( numeric::principal_angle_degrees( pose->omega(i) /*Should handle terminal peptide bonds.*/ ) );
		TR.Debug << "omega" << i << "=" << omegaval << std::endl;
		if ( omegaval <= 90 && omegaval > -90 ) ++count; //Count this as cis if in the interval (-90,90]
	}
	return count;
}

/// @brief Carry out the final FastRelax.
///
void
SimpleCycpepPredictApplication::do_final_fastrelax(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Size const relax_rounds,
	bool const angle_min,
	bool const length_min,
	bool const cartesian_min
) const {
	protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn, 1) );
	if ( cartesian_min ) {
		frlx->cartesian(true);
	} else {
		if ( angle_min ) {
			frlx->minimize_bond_angles(true);
		}
		if ( length_min ) {
			frlx->minimize_bond_lengths(true);
		}
	}

	//Mover to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP final_termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( final_termini, pose );

	(*sfxn)(*pose);
	core::Real cur_energy( pose->energies().total_energy() );
	for ( core::Size i=1; i<=relax_rounds; ++i ) {
		core::pose::PoseOP pose_copy( pose->clone() );
		if ( TR.visible() ) {
			TR << "Applying final FastRelax, round " << i;
			if ( cartesian_min ) TR << ", with Cartesian minimization.";
			else if ( angle_min && length_min ) TR << ", with bond angle and bond length minimization.";
			else if ( angle_min && !length_min ) TR << ", with bond angle minimization.";
			TR << "." << std::endl;
		}
		frlx->apply( *pose_copy );
		final_termini->apply( *pose_copy );
		(*sfxn)(*pose_copy);
		if ( pose_copy->energies().total_energy() < cur_energy ) {
			cur_energy = pose_copy->energies().total_energy();
			(*pose) = (*pose_copy);
		}
#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
#endif
	}
}

/// @brief Actually build the geometry that we'll be working with.
///
void
SimpleCycpepPredictApplication::build_polymer(
	core::pose::PoseOP pose,
	utility::vector1<std::string> const &restypes
) const {
	using namespace protocols::cyclic_peptide;
	core::Size const nres( restypes.size() );
	runtime_assert(restypes.size() >=4 );

	TR << "Building sequence ";
	for ( core::Size i=1; i<=nres; ++i ) {
		TR << restypes[i];
		if ( i<nres ) TR << " ";
	}
	TR << "." << std::endl;

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for ( core::Size i=1; i<=nres; ++i ) {
		stubmover.add_residue( "Append", restypes[i], 0, false, "", 1, 0, "" );
	}

	stubmover.apply(*pose);

	TR << "Build successful." << std::endl;

	return;
} //build_polymer()

/// @brief Add N-methylation.
/// @details Must be called after pose is cyclized.
void
SimpleCycpepPredictApplication::add_n_methylation(
	core::pose::PoseOP pose,
	core::Size const cyclic_offset
) const {
	if ( n_methyl_positions_.size() == 0 ) { return; } //Do nothing if there's no N-methylation.

	core::Size const nres(sequence_length());

	TR << "Adding N-methylation" << std::endl;
	protocols::simple_moves::ModifyVariantTypeMover add_nmethyl;
	add_nmethyl.set_additional_type_to_add("N_METHYLATION");
	core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
	for ( core::Size i=1, imax=n_methyl_positions_.size(); i<=imax; ++i ) {
		runtime_assert_string_msg( n_methyl_positions_[i] > 0  && n_methyl_positions_[i] <= nres, "Error in simple_cycpep_predict app: The N-methylation position indices must be within the pose!" );
		int permuted_position( static_cast<int>(n_methyl_positions_[i]) - static_cast<int>(cyclic_offset) );
		if ( permuted_position < 1 ) permuted_position += static_cast<int>(nres);
		selector->append_index( static_cast<core::Size>(permuted_position) );
	}
	add_nmethyl.set_residue_selector(selector);
	add_nmethyl.apply(*pose);

}

/// @brief Given the name of a Rama_Table_Type, set the default Rama_Table_Type.
/// @details Error if unknown type.
void
SimpleCycpepPredictApplication::set_default_rama_table_type(
	std::string const &type_name
) {
	if ( type_name=="" ) return; //Default case -- no default rama table provided.

	default_rama_table_type_ = get_rama_table_type_from_name( type_name );
	if ( TR.visible() ) {
		TR << "Set default Rama table type to " << type_name << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Given a string vector that we need to parse, populate the rama_table_type_by_res_ map.
/// @details The string vector must be of the format: [integer] [rama type name] [integer] [rama type name] etc.
/// Throws error if could not parse.
void
SimpleCycpepPredictApplication::set_rama_table_type_by_res(
	utility::vector1 <std::string> const &type_name_vector
) {
	core::Size const vectsize( type_name_vector.size() );
	if ( vectsize == 0 ) return;

	core::Size i(0);
	while ( i<vectsize ) {
		++i;
		std::stringstream ss(type_name_vector[i]);
		core::Size tempval(0);
		ss >> tempval;
		if ( ss.fail() ) {
			std::stringstream msg("");
			msg << "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_rama_table_type_by_res(): Could not interpret " << type_name_vector[i] << " as an integer.  The \"-rama_sampling_table_by_res\" flag must be given a series of values of the pattern [res_index] [rama_table_type] [res_index] [rama_table_type] etc." << std::endl;
			utility_exit_with_message(msg.str());
		}
		++i;
		if ( i > vectsize ) {
			std::stringstream msg("");
			msg << "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_rama_table_type_by_res(): A residue index was found with no corresponding Rama table type." << std::endl;
			utility_exit_with_message(msg.str());
		}
		std::string tempstring("");
		std::stringstream ss2(type_name_vector[i]);
		ss2 >> tempstring;
		rama_table_type_by_res_[tempval] = get_rama_table_type_from_name(tempstring);
		if ( TR.visible() ) TR << "Set custom Ramachandran table for residue " << tempval << " to " << tempstring << "." << std::endl;
	}

	if ( TR.visible() ) TR.flush();

	return;
}

/// @brief Given a Rama_Table_Type name, return the Rama_Table_Type, or an informative error message on failure.
///
core::scoring::Rama_Table_Type
SimpleCycpepPredictApplication::get_rama_table_type_from_name(
	std::string const &type_name
) const {
	//Get an instance of the Ramachandran object:
	core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

	//Get the type from the name:
	core::scoring::Rama_Table_Type const type ( rama.get_ramatable_type_by_name(type_name) );

	//Check that the type was parsed sensibly:
	if ( type == core::scoring::unknown_ramatable_type ) {
		std::stringstream err_msg("Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_rama_table_type_from_name():");
		err_msg << "  The provided custom Ramachandran table type (\"" << type_name << "\") is unknown.  Allowed options are: ";
		for ( core::Size i=1; i<static_cast<core::Size>(core::scoring::unknown_ramatable_type); ++i ) {
			err_msg << rama.get_ramatable_name_by_type( static_cast<core::scoring::Rama_Table_Type>(i) );
			if ( i < (static_cast<core::Size>(core::scoring::unknown_ramatable_type) - 1) ) {
				err_msg << ", ";
			}
		}
		err_msg << "." << std::endl;
		utility_exit_with_message( err_msg.str() );
	}

	return(type);
}


/// @brief Read a sequence (as a series of full names, separated by whitespace) and store
/// it in a string vector.
void
SimpleCycpepPredictApplication::read_sequence (
	std::string const &seqfile,
	utility::vector1 < std::string > &resnames
) const {
	using namespace utility::io;
	resnames.clear();

	utility::vector1< std::string > lines; //Storing all lines
	if ( sequence_string_ == "" ) {
		izstream infile;
		infile.open( seqfile );
		runtime_assert_string_msg( infile.good(), "Error in read_sequence() in app simple_cycpep_predict:  Unable to open sequence file for read!" );

		TR << "Opened " << seqfile << " for read." << std::endl;

		std::string curline(""); //Buffer for current line.

		//Read the file:
		while ( getline(infile, curline) ) {
			if ( curline.size() < 1 ) continue; //Ignore blank lines.
			lines.push_back( curline );
		}
		infile.close();
	} else {
		lines.push_back( sequence_string_ );
	}

	//Parse the lines:
	for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) { //Loop through all lines
		if ( TR.Debug.visible() ) TR.Debug << "Parsing \"" << lines[i] << "\"." << std::endl;
		std::istringstream curline(lines[i]);
		std::string oneword("");
		while ( !curline.eof() ) {
			curline >> oneword;
			resnames.push_back( oneword );
		}
	}

	if ( TR.visible() ) {
		TR << "Parsed the following sequence:" << std::endl;
		for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
			TR << resnames[i];
			if ( i<imax ) TR << ", ";
		}
		TR << "." << std::endl;
	}

	runtime_assert_string_msg( resnames.size() >= 4, "Error in simple_cycpcp_predict app read_sequence() function!  The minimum number of residues for a cyclic peptide is 4.  (GenKIC requires three residues, plus a fourth to serve as an anchor)." );

	return;
}


/// @brief Set up the DeclareBond mover used to connect the termini.
///
void
SimpleCycpepPredictApplication::set_up_termini_mover (
	protocols::cyclic_peptide::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native,
	core::Size const last_res /*=0*/
) const {
	core::Size const nres(sequence_length());

	core::Size const cterm( last_res == 0 ? nres : last_res );

	runtime_assert_string_msg(pose->residue(1).has_lower_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: residue 1 does not have a LOWER_CONNECT.");
	runtime_assert_string_msg(pose->residue(cterm).has_upper_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: the final residue does not have an UPPER_CONNECT.");
	std::string firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string lastatom( pose->residue(cterm).atom_name( pose->residue(cterm).upper_connect_atom() ) );

	if ( native ) {
		TR << "Setting up terminal bond for the native pose between residue 1, atom " << firstatom << " and residue " << cterm << ", atom " << lastatom << "." << std::endl;
	} else {
		TR << "Setting up terminal bond between residue 1, atom " << firstatom << " and residue " << cterm << ", atom " << lastatom << "." << std::endl;
	}

	termini->set( cterm, lastatom, 1, firstatom, false, false, 0, 0, false  );

	return;
}


/// @brief Takes a vector of residue names, chooses a random number for cyclic offset, and
/// does a cyclic permutation.
/// @details Returns the offset and stores the new string vector in resnames_copy.
core::Size
SimpleCycpepPredictApplication::do_cyclic_permutation (
	utility::vector1 <std::string> const &resnames,
	utility::vector1 <std::string> &resnames_copy
) const {
	core::Size const nname( resnames.size() );//Number of residue names
	core::Size const offset( static_cast<core::Size>(numeric::random::rg().random_range(0,resnames.size()-1)) );

	resnames_copy.clear();
	resnames_copy.resize(nname, "");
	core::Size counter(0);
	for ( core::Size i=offset+1; i<=nname; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}
	for ( core::Size i=1; i<=offset; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}

	TR << "Circularly shifted residue list by " << offset << ".  New list is: ";
	for ( core::Size i=1; i<=nname; ++i ) {
		TR << resnames_copy[i];
		if ( i<nname ) TR << ", ";
	}
	TR << std::endl;
	TR.flush();

	return offset;
}


/// @brief Imports the native pose and sets up a terminial peptide bond.
///
void
SimpleCycpepPredictApplication::import_and_set_up_native (
	std::string const &native_file,
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) const {
	core::import_pose::pose_from_file(*native_pose, native_file, core::import_pose::PDB_file);
	TR << "Improrting native structure from " << native_file << "." << std::endl;

	set_up_native( native_pose, expected_residue_count );

	return;
}

/// @brief Sets up a terminial peptide bond and does some checks.
///
void
SimpleCycpepPredictApplication::set_up_native (
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) const {
	// Count residues and find the last residue.
	core::Size last_res(0), res_count(0);
	for ( core::Size ir=1, irmax=native_pose->total_residue(); ir<=irmax; ++ir ) {
		if ( native_pose->residue_type(ir).is_alpha_aa() || native_pose->residue_type(ir).is_beta_aa() || native_pose->residue_type(ir).is_gamma_aa() ) {
			++res_count;
			last_res = ir;
		}
	}

	if ( expected_residue_count != 0 ) {
		runtime_assert_string_msg(
			res_count == expected_residue_count,
			"Error in simple_cycpep_predict app!  The imported native pose has a different number of residues than the sequence provided."
		);
	}

	TR << "Stripping termini from native structure." << std::endl;
	core::pose::remove_lower_terminus_type_from_pose_residue(*native_pose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(*native_pose, last_res);

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, native_pose, true, last_res );
	termini->apply(*native_pose);
}


/// @brief Function to add cyclic constraints to a pose.
///
void
SimpleCycpepPredictApplication::add_cyclic_constraints (
	core::pose::PoseOP pose
) const {
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up cyclic constraints." << std::endl;

	core::Size const nres(sequence_length());

	//The four atoms defining the peptide bond:
	AtomID const atom_a( pose->residue(nres).type().icoor(pose->residue(nres).upper_connect_atom()).stub_atom1().atomno(), nres );
	AtomID const atom_b( pose->residue(nres).upper_connect_atom(), nres );
	AtomID const atom_c( pose->residue(1).lower_connect_atom(), 1 );
	core::Size atom_d_index(0);
	for ( core::Size i=1, imax=pose->residue(1).n_mainchain_atoms(); i<=imax; ++i ) { //Find the atom index of the first mainchain atom with the lower_connect atom as a parent.
		if ( i == atom_c.atomno() ) continue;
		if ( pose->residue(1).type().icoor(i).stub_atom1().atomno() == atom_c.atomno() ) {
			atom_d_index=i;
			break;
		}
	}
	AtomID const atom_d( atom_d_index, 1 );

	TR << "The following four atoms define the terminal bond:" << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose->residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose->residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose->residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose->residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;

	{//Peptide bond length constraint:
		FuncOP harmfunc1( new HarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_LENGTH, 0.01) );
		ConstraintCOP distconst1( new AtomPairConstraint ( atom_b, atom_c, harmfunc1 ) );
		pose->add_constraint (distconst1);
	}

	{ //Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1( new CircularHarmonicFunc( numeric::constants::d::pi, 0.02) );
		ConstraintCOP dihedconst1( new DihedralConstraint ( atom_a, atom_b, atom_c, atom_d, circharmfunc1) );
		pose->add_constraint (dihedconst1);
	}

	{ //Peptide bond angle constraints:
		FuncOP circharmfunc2a( new CircularHarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_C_ANGLE, 0.02) );
		FuncOP circharmfunc2b( new CircularHarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_N_ANGLE, 0.02) );
		ConstraintCOP angleconst1( new AngleConstraint ( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( new AngleConstraint ( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose->add_constraint (angleconst1);
		pose->add_constraint (angleconst2);
	}

	TR << "Finished setting up constraints." << std::endl;

	return;
}


/// @brief Sets all omega values to 180, and randomizes mainchain torsions.
/// @details For alpha-amino acids, mainchain torsions are randomized by the Ramachandran plot.
/// For other residue types, just randomizes mainchain torsions other than peptide bonds.
void
SimpleCycpepPredictApplication::set_mainchain_torsions (
	core::pose::PoseOP pose,
	core::Size const cyclic_offset
) const {
	TR << "Randomizing mainchain torsions." << std::endl;
	core::Size const nres(sequence_length());
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues
		if ( pose->residue(i).type().is_alpha_aa() ) {
			core::scoring::Ramachandran & rama(core::scoring::ScoringManager::get_instance()->get_Ramachandran_nonconst() ); //Get the Rama scoring function; must be nonconst to allow lazy loading
			core::scoring::RamaPrePro const & ramaprepro( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
			core::Real phi(0.0), psi(0.0);
			//TR << "aa" << i << "=" << pose->residue_type(i).aa() << std::endl; //DELETE ME
			core::Size const cur_abs_pos( original_position( i, cyclic_offset, pose->size() ) );
			if ( custom_rama_table_defined( cur_abs_pos ) ) {
				rama.draw_random_phi_psi_from_extra_cdf( rama_table_type_by_res(cur_abs_pos), phi, psi);
			} else if ( default_rama_table_type() != core::scoring::unknown_ramatable_type ) {
				rama.draw_random_phi_psi_from_extra_cdf( default_rama_table_type(), phi, psi);
			} else {

				if ( use_rama_prepro_for_sampling() ) { //Using rama_prepro tables for sampling:
					utility::vector1< core::Real > rand_torsions;
					core::chemical::ResidueTypeCOP following_rsd( pose->residue_type(
						pose->residue(i).residue_connection_partner( pose->residue(i).upper_connect().index() )
						).get_self_ptr()
					);
					ramaprepro.random_mainchain_torsions( pose->conformation(), pose->residue_type(i).get_self_ptr(), following_rsd, rand_torsions );
					phi=rand_torsions[1]; psi=rand_torsions[2];
				} else { //Using classic rama tables for sampling:
					if ( pose->residue(i).backbone_aa() != core::chemical::aa_unk ) {
						rama.random_phipsi_from_rama(pose->residue(i).backbone_aa(), phi, psi);
					} else {
						rama.random_phipsi_from_rama( pose->residue_type(i).aa(), phi, psi);
					}
				}

			}
			pose->set_phi(i,phi);
			pose->set_psi(i,psi);
			if ( i!=nres ) pose->set_omega(i, 180.0);
		} else { //If this is not an alpha-amino acid:
			for ( core::Size j=1, jmax=pose->residue(i).mainchain_torsions().size(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
				if ( i==nres && j==jmax ) continue; //Skip the last mainchain torsion (not a DOF).
				core::Real setting(180.0);
				if ( j!=jmax ) {
					setting = numeric::random::rg().uniform()*360.0 - 180.0;
				}
				pose->set_torsion( core::id::TorsionID(i, core::id::BB, j), setting );
			}
		}
	}
	return;
}

/// @brief Set up the filters for the mainchain hydrogen bonds that will
/// be used to discard solutions with too little mainchain hydrogen bonding.
void
SimpleCycpepPredictApplication::set_up_hbond_filter(
	protocols::filters::CombinedFilterOP total_hbond,
	core::Size const nres,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const &min_hbonds
) const {
	total_hbond->set_threshold( -1.0 * min_hbonds );
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues and add hbond counters
		protocols::protein_interface_design::filters::HbondsToResidueFilterOP hbondfilt( new protocols::protein_interface_design::filters::HbondsToResidueFilter );
		hbondfilt->set_resnum(i);
		hbondfilt->set_sidechain( count_sc_hbonds_ );
		hbondfilt->set_energy_cutoff( hbond_energy_cutoff_ );
		hbondfilt->set_partners(0);
		hbondfilt->set_scorefxn( sfxn );

		//Add ResidueSelectors to ensure that hydrogen bonds with adjacent residues are not counted.
		if ( do_not_count_adjacent_res_hbonds_ ) {
			std::stringstream indices_string("");
			core::Size const avoid1( i - 1 > 0 ? i - 1 : nres );
			core::Size const avoid2( i + 1 <= nres ? i + 1 : 1 );
			bool first(true);
			for ( core::Size j=1; j<=nres; ++j ) {
				if ( j == i || j == avoid1 || j == avoid2 ) continue; //Skip the current residue and its adjacent residues.
				if ( first ) { first=false; }
				else { indices_string << ","; }
				indices_string << j;
			}
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( new core::select::residue_selector::ResidueIndexSelector );
			index_selector->set_index( indices_string.str() );
			hbondfilt->set_selector( index_selector );
		}

		total_hbond->add_filter( hbondfilt, -0.5, false );
	}
	return;
}


/// @brief Use GeneralizedKIC to close the pose.
///
bool
SimpleCycpepPredictApplication::genkic_close(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn_highhbond,
	core::scoring::ScoreFunctionOP sfxn_highhbond_cart,
	core::scoring::ScoreFunctionCOP sfxn_default,
	protocols::filters::CombinedFilterOP total_hbond,
	core::Size const cyclic_offset
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::selector;

	TR << "Performing GeneralizedKIC closure of loop." << std::endl;

	//Number of residues in the pose:
	core::Size const nres( sequence_length() );
	runtime_assert( nres >= 4 ); //Already checked at sequence load time, so should be true, but let's make sure.
	core::Size const res_per_symm_repeat(
		required_symmetry_repeats_ > 1 ? nres / required_symmetry_repeats_ : nres
	);

	//Randomly pick one of the middle residues to be the anchor residue (or one of the middle residues of the asymmetric unit if this peptide has symmetry):
	core::Size const anchor_res(
		required_symmetry_repeats_ > 1 ? numeric::random::rg().random_range(2, res_per_symm_repeat) : numeric::random::rg().random_range(2, nres-1)
	);
	core::Size const first_loop_res( anchor_res + 1 );
	core::Size const last_loop_res( anchor_res - 1 );

	//Randomly pick a residue to be the middle pivot residue.  Can't be first in loop, last in loop, or anchor res.
	core::Size middle_loop_res( numeric::random::rg().random_range(1, nres-3 ) );
	if ( middle_loop_res == last_loop_res ) { middle_loop_res += 3; }
	else if ( middle_loop_res == anchor_res ) { middle_loop_res +=2; }
	else if ( middle_loop_res == first_loop_res ) { middle_loop_res +=1; }
	if ( middle_loop_res > nres ) { middle_loop_res -= nres; }

	//Create the pre-selection mover and set options.
	protocols::rosetta_scripts::ParsedProtocolOP pp( new protocols::rosetta_scripts::ParsedProtocol );

	//Update O and H atoms at the cyclization point:
	protocols::cyclic_peptide::DeclareBondOP update_OH( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( update_OH, pose );
	pp->add_mover_filter_pair( update_OH, "Update_cyclization_point_polymer_dependent_atoms_1", nullptr );

	//Filter for total hydrogen bonds:
	if ( min_genkic_hbonds_ > 0.0 ) pp->add_mover_filter_pair( nullptr, "Total_Hbonds", total_hbond );

	//Filter out poses with oversaturated hydrogen bond acceptors.
	if ( filter_oversaturated_hbond_acceptors_ ) {
		protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterOP oversat1( new protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter );
		oversat1->set_scorefxn( sfxn_default );
		oversat1->set_hbond_energy_cutoff( oversaturated_hbond_cutoff_energy_ );
		pp->add_mover_filter_pair( nullptr, "Oversaturated_Hbond_Acceptors", oversat1 );
	}

	//If we're filtering by symmetry, do so here:
	if ( required_symmetry_repeats_ > 1 ) {
		protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter1( new protocols::cyclic_peptide::CycpepSymmetryFilter );
		symmfilter1->set_symm_repeats( required_symmetry_repeats_ );
		symmfilter1->set_mirror_symm( required_symmetry_mirroring_ );
		symmfilter1->set_angle_threshold( required_symmetry_angle_threshold_ );
		core::select::residue_selector::ResidueIndexSelectorOP iselector( new core::select::residue_selector::ResidueIndexSelector );
		std::stringstream pep_indices("");
		pep_indices << "1-" << sequence_length();
		iselector->set_index( pep_indices.str() );
		symmfilter1->set_selector( iselector );
		pp->add_mover_filter_pair(nullptr, "Cycpep_Symmetry_Filter_1", symmfilter1);
	}

	//If we're considering TBMB, add it here.
	if ( tbmb_positions_.size() > 0 ) {
		for ( core::Size i=1, imax=tbmb_positions_.size(); i<=imax; ++i ) { //Loop through all sets of triples of residues.
			debug_assert(tbmb_positions_[i].size() == 3); //Should always be true.
			std::stringstream cys_indices;
			for ( core::Size j=1; j<=3; ++j ) {
				cys_indices << current_position( tbmb_positions_[i][j], cyclic_offset, nres );
				if ( j<3 ) cys_indices << ",";
			}
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( new core::select::residue_selector::ResidueIndexSelector );
			index_selector->set_index( cys_indices.str() );
			protocols::cyclic_peptide::ThreefoldLinkerMoverOP threelinker( new protocols::cyclic_peptide::ThreefoldLinkerMover );
			threelinker->set_residue_selector(index_selector);
			threelinker->set_linker_name("TBMB");
			threelinker->set_behaviour( true, true, true, false, false );
			threelinker->set_filter_behaviour( use_tbmb_filters_, use_tbmb_filters_, false, 0.0, tbmb_sidechain_distance_filter_multiplier_, tbmb_constraints_energy_filter_multiplier_ );
			threelinker->set_scorefxn( sfxn_highhbond );
			threelinker->set_sidechain_frlx_rounds(3);
			//TODO -- set options for filtering.
			std::stringstream movername;
			movername << "TBMB_link_" << i;
			pp->add_mover_filter_pair( threelinker, movername.str(), nullptr );
		}
	}

	core::Size disulf_count(0);
	//If we're considering disulfides, add the TryDisulfPermutations mover and a filter to the ParsedProtocol:
	if ( try_all_disulfides_ ) {
		protocols::cyclic_peptide::TryDisulfPermutationsOP trydisulf( new protocols::cyclic_peptide::TryDisulfPermutations ); //Default settings should be fine.
		core::Size disulf_res_count(0);
		for ( core::Size ir=1, irmax=pose->size(); ir<=irmax; ++ir ) { if ( pose->residue(ir).type().get_disulfide_atom_name() != "NONE" ) ++disulf_res_count; } //Count disulfide-forming residues in the pose.
		disulf_count = disulf_res_count / 2; //Div operator -- gives correct number of disulfides even in odd disulfide-forming residue case.
		protocols::simple_filters::ScoreTypeFilterOP disulf_filter1( new protocols::simple_filters::ScoreTypeFilter( sfxn_highhbond, core::scoring::dslf_fa13, disulf_energy_cutoff_prerelax_ * static_cast<core::Real>(disulf_count) ) );
		pp->add_mover_filter_pair( trydisulf, "Try_Disulfide_Permutations", disulf_filter1 );
	}

	//Add the FastRelax with high hbond weight to the pre-selection parsed protocol.
	if ( design_peptide_ ) {
		if ( L_alpha_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP L_alpha_cst( new protocols::aa_composition::AddCompositionConstraintMover );
			core::select::residue_selector::BinSelectorOP select_L_alpha( new core::select::residue_selector::BinSelector );
			select_L_alpha->set_bin_name("A");
			select_L_alpha->initialize_from_file_contents_and_check( abba_bins_ );
			L_alpha_cst->create_constraint_from_file_contents( comp_file_contents_L_alpha_ );
			L_alpha_cst->add_residue_selector( select_L_alpha );
			pp->add_mover_filter_pair( L_alpha_cst, "Add_L_Alpha_AACompositionConstraints", nullptr );
		}
		if ( D_alpha_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP D_alpha_cst( new protocols::aa_composition::AddCompositionConstraintMover );
			core::select::residue_selector::BinSelectorOP select_D_alpha( new core::select::residue_selector::BinSelector );
			select_D_alpha->set_bin_name("Aprime");
			select_D_alpha->initialize_from_file_contents_and_check( abba_bins_ );
			D_alpha_cst->create_constraint_from_file_contents( comp_file_contents_D_alpha_ );
			D_alpha_cst->add_residue_selector( select_D_alpha );
			pp->add_mover_filter_pair( D_alpha_cst, "Add_D_Alpha_AACompositionConstraints", nullptr );
		}
		if ( L_beta_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP L_beta_cst( new protocols::aa_composition::AddCompositionConstraintMover );
			core::select::residue_selector::BinSelectorOP select_L_beta( new core::select::residue_selector::BinSelector );
			select_L_beta->set_bin_name("B");
			select_L_beta->initialize_from_file_contents_and_check( abba_bins_ );
			L_beta_cst->create_constraint_from_file_contents( comp_file_contents_L_beta_ );
			L_beta_cst->add_residue_selector( select_L_beta );
			pp->add_mover_filter_pair( L_beta_cst, "Add_L_Beta_AACompositionConstraints", nullptr );
		}
		if ( D_beta_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP D_beta_cst( new protocols::aa_composition::AddCompositionConstraintMover );
			core::select::residue_selector::BinSelectorOP select_D_beta( new core::select::residue_selector::BinSelector );
			select_D_beta->set_bin_name("Bprime");
			select_D_beta->initialize_from_file_contents_and_check( abba_bins_ );
			D_beta_cst->create_constraint_from_file_contents( comp_file_contents_D_beta_ );
			D_beta_cst->add_residue_selector( select_D_beta );
			pp->add_mover_filter_pair( D_beta_cst, "Add_D_Beta_AACompositionConstraints", nullptr );
		}

		if ( fast_relax_rounds_ > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes( new protocols::denovo_design::movers::FastDesign(sfxn_highhbond, fast_relax_rounds_) );
			set_up_design_taskoperations( fdes, cyclic_offset, pose->size() );
			pp->add_mover_filter_pair( fdes, "High_Hbond_FastDesign", nullptr );
		}
		if ( angle_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes2( new protocols::denovo_design::movers::FastDesign(sfxn_highhbond_cart, angle_relax_rounds()) );
			fdes2->minimize_bond_angles(true);
			set_up_design_taskoperations( fdes2, cyclic_offset, pose->size() );
			pp->add_mover_filter_pair( fdes2, "High_Hbond_FastDesign_angle_relax", nullptr );
		}
		if ( angle_length_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes3( new protocols::denovo_design::movers::FastDesign(sfxn_highhbond_cart, angle_length_relax_rounds()) );
			fdes3->minimize_bond_angles(true);
			fdes3->minimize_bond_lengths(true);
			set_up_design_taskoperations( fdes3, cyclic_offset, pose->size() );
			pp->add_mover_filter_pair( fdes3, "High_Hbond_FastDesign_angle_length_relax", nullptr );
		}
		if ( cartesian_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes4( new protocols::denovo_design::movers::FastDesign(sfxn_highhbond_cart, cartesian_relax_rounds()) );
			fdes4->cartesian(true);
			set_up_design_taskoperations( fdes4, cyclic_offset, pose->size() );
			pp->add_mover_filter_pair( fdes4, "High_Hbond_FastDesign_Cartesian_relax", nullptr );
		}

		protocols::aa_composition::ClearCompositionConstraintsMoverOP clear_aacomp_cst( new protocols::aa_composition::ClearCompositionConstraintsMover );
		pp->add_mover_filter_pair( clear_aacomp_cst, "Clear_AACompositionConstraints", nullptr );
	} else {
		if ( fast_relax_rounds_ > 0 ) {
			protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn_highhbond, fast_relax_rounds_) );
			pp->add_mover_filter_pair( frlx, "High_Hbond_FastRelax", nullptr );
		}
		if ( angle_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx2( new protocols::relax::FastRelax(sfxn_highhbond_cart, angle_relax_rounds() ) );
			frlx2->minimize_bond_angles(true);
			pp->add_mover_filter_pair( frlx2, "High_Hbond_FastRelax_angles", nullptr );
		}
		if ( angle_length_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx3( new protocols::relax::FastRelax(sfxn_highhbond_cart, angle_length_relax_rounds() ) );
			frlx3->minimize_bond_angles(true);
			pp->add_mover_filter_pair( frlx3, "High_Hbond_FastRelax_angles_bondlengths", nullptr );
		}
		if ( cartesian_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx4( new protocols::relax::FastRelax(sfxn_highhbond_cart, cartesian_relax_rounds() ) );
			frlx4->minimize_bond_angles(true);
			pp->add_mover_filter_pair( frlx4, "High_Hbond_FastRelax_Cartesian", nullptr );
		}
	}

	//Update O and H atoms at the cyclization point:
	pp->add_mover_filter_pair( update_OH, "Update_cyclization_point_polymer_dependent_atoms_2", nullptr );

	//Add more stringent disulfide filtering post-relax:
	if ( try_all_disulfides_ ) {
		protocols::simple_filters::ScoreTypeFilterOP disulf_filter2( new protocols::simple_filters::ScoreTypeFilter( sfxn_highhbond, core::scoring::dslf_fa13, disulf_energy_cutoff_postrelax_ * static_cast<core::Real>(disulf_count) ) );
		pp->add_mover_filter_pair( nullptr, "Postrelax_disulfide_filter", disulf_filter2 );
	}

	if ( filter_oversaturated_hbond_acceptors_ ) {
		protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterOP oversat2( new protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter );
		oversat2->set_scorefxn( sfxn_default );
		oversat2->set_hbond_energy_cutoff( oversaturated_hbond_cutoff_energy_ );
		pp->add_mover_filter_pair( nullptr, "Postrelax_Oversaturated_Hbond_Acceptors", oversat2 );
	}

	//If we're filtering by symmetry, do so again here:
	if ( required_symmetry_repeats_ > 1 ) {
		protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter1( new protocols::cyclic_peptide::CycpepSymmetryFilter );
		symmfilter1->set_symm_repeats( required_symmetry_repeats_ );
		symmfilter1->set_mirror_symm( required_symmetry_mirroring_ );
		symmfilter1->set_angle_threshold( required_symmetry_angle_threshold_ );
		core::select::residue_selector::ResidueIndexSelectorOP iselector( new core::select::residue_selector::ResidueIndexSelector );
		std::stringstream pep_indices("");
		pep_indices << "1-" << sequence_length();
		iselector->set_index( pep_indices.str() );
		symmfilter1->set_selector( iselector );
		pp->add_mover_filter_pair(nullptr, "Cycpep_Symmetry_Filter_2", symmfilter1);
	}

	//Create the mover and set options:
	GeneralizedKICOP genkic( new GeneralizedKIC );
	genkic->set_selector_type( lowest_energy_selector );
	genkic->set_closure_attempts( genkic_closure_attempts_ );
	genkic->set_min_solution_count( genkic_min_solution_count_ );
	genkic->set_selector_scorefunction( sfxn_highhbond );
	genkic->set_preselection_mover(pp);
	genkic->set_correct_polymer_dependent_atoms(true);

	//If we're using BOINC graphics, let the GenKIC mover update the graphics with a "ghost" of the current
	//conformation being sampled:
#ifdef BOINC_GRAPHICS
	genkic->set_attach_boinc_ghost_observer(true);
#endif


	//Define the loop residues:
	for ( core::Size i=first_loop_res; i<=nres; ++i ) { genkic->add_loop_residue(i); }
	for ( core::Size i=1; i<=last_loop_res; ++i ) { genkic->add_loop_residue(i); }

	//Set pivots:
	std::string at1(""), at2(""), at3("");
	if ( pose->residue(first_loop_res).type().is_alpha_aa() ) { at1="CA"; }
	else if ( pose->residue(first_loop_res).type().is_beta_aa() ) { at1="CM"; }
	else if ( pose->residue(first_loop_res).type().is_gamma_aa() ) { at1="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop start.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(middle_loop_res).type().is_alpha_aa() ) { at2="CA"; }
	else if ( pose->residue(middle_loop_res).type().is_beta_aa() ) { at2="CM"; }
	else if ( pose->residue(middle_loop_res).type().is_gamma_aa() ) { at2="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(last_loop_res).type().is_alpha_aa() ) { at3="CA"; }
	else if ( pose->residue(last_loop_res).type().is_beta_aa() ) { at3="CM"; }
	else if ( pose->residue(last_loop_res).type().is_gamma_aa() ) { at3="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	genkic->set_pivot_atoms( first_loop_res, at1, middle_loop_res, at2, last_loop_res, at3 );

	//Close the bond:
	std::string const firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string const lastatom( pose->residue(nres).atom_name( pose->residue(nres).upper_connect_atom() ) );
	genkic->close_bond( nres, lastatom, 1, firstatom, 0, "", 0, "", SimpleCycpepPredictApplication_PEPBOND_LENGTH, SimpleCycpepPredictApplication_PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, SimpleCycpepPredictApplication_PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );

	//Add perturbers:
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( i==anchor_res ) continue; //Can't perturb the anchor residue.
		if ( pose->residue(i).type().is_alpha_aa() ) {
			core::Size const res_in_original( original_position(i, cyclic_offset, pose->size() ) ); //Get the index of this position in the original pose (prior to any circular permutation).
			if ( user_set_alpha_dihedrals_.count(res_in_original) ) { //If this position is being set to a particular value...
				genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::set_dihedral );
				core::id::NamedAtomID Natom( "N", i );
				core::id::NamedAtomID CAatom( "CA", i );
				core::id::NamedAtomID Catom( "C", i );
				core::Size const nextres( i==nres ? 1 : i+1 );
				core::id::NamedAtomID Nnextatom( "N", nextres );
				utility::vector1< core::id::NamedAtomID > phivect; phivect.push_back( Natom ); phivect.push_back( CAatom );
				utility::vector1< core::id::NamedAtomID > psivect; psivect.push_back( CAatom ); psivect.push_back( Catom );
				utility::vector1< core::id::NamedAtomID > omegavect; omegavect.push_back( Catom ); omegavect.push_back( Nnextatom );
				genkic->add_atomset_to_perturber_atomset_list(phivect);
				genkic->add_atomset_to_perturber_atomset_list(psivect);
				if ( nextres != anchor_res ) genkic->add_atomset_to_perturber_atomset_list(omegavect);
				genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[1] );
				genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[2] );
				if ( nextres != anchor_res ) genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[3] );
				if ( user_set_dihedral_perturbation_ ) {
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::perturb_dihedral );
					genkic->add_atomset_to_perturber_atomset_list(phivect);
					genkic->add_atomset_to_perturber_atomset_list(psivect);
					//if ( nextres != anchor_res ) genkic->add_atomset_to_perturber_atomset_list(omegavect);
					genkic->add_value_to_perturber_value_list( user_set_dihedral_perturbation_ );
				}
			} else { //If this position is not set, randomize it.
				if ( custom_rama_table_defined( res_in_original ) ) { //If there is a custom rama table defined for sampling at this position, use it.
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
					genkic->add_residue_to_perturber_residue_list(i);
					genkic->set_perturber_custom_rama_table( rama_table_type_by_res( res_in_original ) );
				} else {
					if ( use_rama_prepro_for_sampling() && default_rama_table_type() == core::scoring::unknown_ramatable_type ) {
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_backbone_by_rama_prepro );
						genkic->add_residue_to_perturber_residue_list(i);
					} else {
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
						genkic->add_residue_to_perturber_residue_list(i);
						if ( default_rama_table_type() != core::scoring::unknown_ramatable_type ) {
							genkic->set_perturber_custom_rama_table( default_rama_table_type() );
						}
					}
				}
				if ( required_symmetry_repeats_ > 1 && i > res_per_symm_repeat ) { //This is a symmetry repeat
					core::Size res_to_copy( i % res_per_symm_repeat );
					if ( res_to_copy == 0 ) { res_to_copy = res_per_symm_repeat; }
					if ( required_symmetry_mirroring_ && ( (i-1) / res_per_symm_repeat ) % 2 == 1 ) {
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::mirror_backbone_dihedrals );
					} else {
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::copy_backbone_dihedrals );
					}
					genkic->add_residue_to_perturber_residue_list( res_to_copy );
					genkic->add_residue_to_perturber_residue_list( i );
					if ( required_symmetry_perturbation_ != 0.0 ) {
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::perturb_dihedral );
						core::id::NamedAtomID Natom( "N", i );
						core::id::NamedAtomID CAatom( "CA", i );
						core::id::NamedAtomID Catom( "C", i );
						utility::vector1< core::id::NamedAtomID > phivect; phivect.push_back( Natom ); phivect.push_back( CAatom );
						utility::vector1< core::id::NamedAtomID > psivect; psivect.push_back( CAatom ); psivect.push_back( Catom );
						genkic->add_atomset_to_perturber_atomset_list(phivect);
						genkic->add_atomset_to_perturber_atomset_list(psivect);
						//if ( nextres != anchor_res ) genkic->add_atomset_to_perturber_atomset_list(omegavect);
						genkic->add_value_to_perturber_value_list( required_symmetry_perturbation_ );
					}
				}
			}
		} else {
			//TODO Randomize mainchain torsions here for beta- and gamma-amino acids.
			utility_exit_with_message( "Handling of beta- and gamma-amino acids in setup of the genKIC perturber in the simple_cycpep_predict app has not yet been written.  TODO." );
		}
	}
	//Additional perturber: sampling cis proline.  Must be after the other perturbers.
	if ( sample_cis_pro() ) {
		for ( core::Size i=1; i<=nres; ++i ) {
			if ( i==1 && nres==anchor_res ) continue; //Can't perturb the anchor residue.
			if ( i-1==anchor_res ) continue; //Can't perturb the anchor residue.
			if ( pose->residue_type(i).aa() == core::chemical::aa_pro || pose->residue_type(i).aa() == core::chemical::aa_dpr || pose->residue_type(i).is_n_methylated() ) {
				genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::sample_cis_peptide_bond );
				genkic->add_value_to_perturber_value_list( sample_cis_pro_frequency() );
				genkic->add_residue_to_perturber_residue_list( i==1 ? nres : i-1 ); //The residue PRECEDING the proline is perturbed
			}
		}
	}

	//Add bump check filter:
	genkic->add_filter( protocols::generalized_kinematic_closure::filter::loop_bump_check );

	//Add rama check filters:
	if ( use_rama_filter() ) {
		for ( core::Size i=1; i<=nres; ++i ) {
			if ( i!=first_loop_res && i!=middle_loop_res && i!=last_loop_res ) continue; //Just filter the pivots.
			if ( use_rama_prepro_for_sampling() ) {
				genkic->add_filter( protocols::generalized_kinematic_closure::filter::rama_prepro_check );
				genkic->set_filter_resnum(i);
				genkic->set_filter_rama_cutoff_energy( rama_cutoff_ );
				if ( i==first_loop_res ) genkic->set_filter_attach_boinc_ghost_observer(true);
			} else {
				if ( pose->residue(i).type().is_alpha_aa() ) {
					genkic->add_filter( protocols::generalized_kinematic_closure::filter::alpha_aa_rama_check );
					genkic->set_filter_resnum(i);
					genkic->set_filter_rama_cutoff_energy( rama_cutoff_ );
					if ( i==first_loop_res ) genkic->set_filter_attach_boinc_ghost_observer(true);
				}
			}
		}
	}

	//Apply the mover:
	genkic->apply( *pose );

	return genkic->last_run_successful();
}

/// @brief Set up the TaskOperations that conrol the design process, given user inputs.
/// @details Default behaviour is designing all positions with L-canonicals and their
/// D-equivalents EXCEPT cys and met (and gly), unless the user overrides this.
void
SimpleCycpepPredictApplication::set_up_design_taskoperations(
	protocols::denovo_design::movers::FastDesignOP fdes,
	core::Size const cyclic_offset,
	core::Size const nres
) const {
	fdes->set_up_default_task_factory();

	std::stringstream L_resfile("");
	std::stringstream L_empty_resfile("");
	std::stringstream D_resfile("");
	L_resfile << "start" << std::endl;
	L_empty_resfile << "start" << std::endl;
	D_resfile << "start" << std::endl;

	for ( core::Size i=1; i<=nres; ++i ) {
		signed long int orig_res( static_cast<signed long int>(i) - static_cast<signed long int>( cyclic_offset ) - 1 );
		if ( orig_res < 1 ) orig_res += nres;
		debug_assert( orig_res >= 1 && orig_res <= static_cast<signed long int>(nres) ); //Should be true.
		L_empty_resfile << orig_res << " A EMPTY" << std::endl;
		if ( allowed_canonicals_by_position_.count( static_cast<core::Size>(orig_res) ) == 1 ) {
			if ( allowed_canonicals_by_position_.at( static_cast<core::Size>(orig_res) ).size() > 0 ) {
				L_resfile << orig_res << " A PIKAA " << get_oneletter_codes( allowed_canonicals_by_position_.at( static_cast<core::Size>(orig_res) ) ) << std::endl;
			} else {
				L_resfile << orig_res << " A EMPTY" << std::endl;
			}
		} else if ( allowed_canonicals_by_position_.count( 0 ) == 1 ) {
			if ( allowed_canonicals_by_position_ .at( 0 ).size() > 0 ) {
				L_resfile << orig_res << " A PIKAA " << get_oneletter_codes( allowed_canonicals_by_position_.at( 0 ) ) << std::endl;
			} else {
				L_resfile << orig_res << " A EMPTY" << std::endl;
			}
		} else {
			L_resfile << orig_res << " A PIKAA ADEFHIKLNPQRSTVWY" << std::endl;
		}
		if ( allowed_noncanonicals_by_position_.count( static_cast<core::Size>(orig_res) ) == 1 ) {
			if ( allowed_noncanonicals_by_position_.at( static_cast<core::Size>(orig_res) ).size() > 0 ) {
				D_resfile << orig_res << " A " << get_nc_threeletter_codes( allowed_noncanonicals_by_position_.at( static_cast<core::Size>(orig_res) ) ) << std::endl;
			}
		} else if ( allowed_noncanonicals_by_position_.count( 0 ) == 1 ) {
			if ( allowed_noncanonicals_by_position_.at( 0 ).size() > 0 ) {
				D_resfile << orig_res << " A " << get_nc_threeletter_codes( allowed_noncanonicals_by_position_.at( 0 ) ) << std::endl;
			}
		} else {
			D_resfile << orig_res << " A NC DAL NC DAS NC DGU NC DPH NC DHI NC DIL NC DLY NC DLE NC DAN NC DPR NC DGN NC DAR NC DSE NC DTH NC DVA NC DTR NC DTY" << std::endl;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "L-resfile:\n" << L_resfile.str() << std::endl;
		TR.Debug << "L-empty resfile:\n" << L_empty_resfile.str() << std::endl;
		TR.Debug << "D-resfile:\n" << D_resfile.str() << std::endl;
	}

	core::pack::task::operation::ReadResfileOP L_resfile_taskop( new core::pack::task::operation::ReadResfile );
	core::pack::task::operation::ReadResfileOP L_empty_resfile_taskop( new core::pack::task::operation::ReadResfile );
	core::pack::task::operation::ReadResfileOP D_resfile_taskop( new core::pack::task::operation::ReadResfile );

	L_resfile_taskop->set_cached_resfile( L_resfile.str() );
	D_resfile_taskop->set_cached_resfile( D_resfile.str() );
	if ( prohibit_L_at_positive_phi_ ) {
		core::select::residue_selector::PhiSelectorOP neg_phi_selector( new core::select::residue_selector::PhiSelector );
		neg_phi_selector->set_select_positive_phi(false);
		L_resfile_taskop->set_residue_selector( neg_phi_selector );

		L_empty_resfile_taskop->set_cached_resfile( L_empty_resfile.str() );
		core::select::residue_selector::NotResidueSelectorOP notselector( new core::select::residue_selector::NotResidueSelector );
		notselector->set_residue_selector( neg_phi_selector );
		L_empty_resfile_taskop->set_residue_selector( notselector );
	}
	if ( prohibit_D_at_negative_phi_ ) {
		core::select::residue_selector::PhiSelectorOP pos_phi_selector( new core::select::residue_selector::PhiSelector );
		pos_phi_selector->set_select_positive_phi(true);
		D_resfile_taskop->set_residue_selector( pos_phi_selector );
	}

	fdes->get_task_factory()->push_back( L_resfile_taskop );
	if ( prohibit_L_at_positive_phi_ ) fdes->get_task_factory()->push_back( L_empty_resfile_taskop );
	fdes->get_task_factory()->push_back( D_resfile_taskop );
}

/// @brief Given a vector of full residue names of canonical residues, give me a concatenated list of one-letter codes.
/// @details Does no checking for duplicates.
std::string
SimpleCycpepPredictApplication::get_oneletter_codes(
	utility::vector1< std::string > const &fullnames
) const {
	std::stringstream outstr;
	for ( core::Size i=1, imax=fullnames.size(); i<=imax; ++i ) {
		if ( fullnames[i] == "HIS_D" ) {
			outstr << "H";
		} else { //We can use the fact that the full name and the three-letter code are the same for canonicals:
			outstr << core::chemical::oneletter_code_from_aa( core::chemical::aa_from_name( fullnames[i] ) );
		}
	}

	return outstr.str();
}

/// @brief Given a vector of full residue names, give me a string of the form "NC <3-letter code> NC <3-letter code> NC <3-letter code> ..."
/// @details Does no checking for duplicates.  Will fail gracelessly with invalid names.
std::string
SimpleCycpepPredictApplication::get_nc_threeletter_codes(
	utility::vector1< std::string> const &fullnames
) const {
	using namespace core::chemical;
	std::stringstream outstr;
	ResidueTypeSetCOP restype_set( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	for ( core::Size i=1, imax=fullnames.size(); i<=imax; ++i ) {
		ResidueTypeCOP curtype( ResidueTypeFinder( *restype_set ).residue_base_name( fullnames[i] ).get_representative_type() );
		runtime_assert_string_msg( curtype, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_nc_threeletter_codes(): The name " + fullnames[i] + " corresponds to no known residue type." );
		outstr << "NC " << curtype->name3();
		if ( i < imax ) outstr << " ";
	}
	return outstr.str();
}


/// @brief Given a pose, store a list of the disulfides in the pose.
/// @details Clears the old_disulfides list and repopulates it.
void
SimpleCycpepPredictApplication::store_disulfides (
	core::pose::PoseCOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > &old_disulfides
) const {
	old_disulfides.clear();
	core::conformation::disulfide_bonds( pose->conformation(), old_disulfides );
	return;
}

/// @brief Given a pose and a list of the disulfides in the pose, break the disulfides.
///
void
SimpleCycpepPredictApplication::break_disulfides (
	core::pose::PoseOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
) const {
	for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
		core::conformation::break_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second );
	}
	return;
}

/// @brief Given a pose and a list of the disulfides that should be in the pose, form the disulfides.
///
void
SimpleCycpepPredictApplication::rebuild_disulfides (
	core::pose::PoseOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
) const {
	for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
		core::conformation::form_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second, true, false );
	}
	return;
}

/// @brief Given a list of old disulfide positions, generate a list of new disulfide positions based on the offset.
/// @details Replaces the new_disulfides list.
void
SimpleCycpepPredictApplication::depermute_disulfide_list(
	utility::vector1 < std::pair < core::Size, core::Size > > const &old_disulfides,
	utility::vector1 < std::pair < core::Size, core::Size > > &new_disulfides,
	core::Size const offset,
	core::Size const nres
) const {
	//TR << "Depermuting disulfide list." << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.

	new_disulfides.clear();
	for ( core::Size i=1, imax=old_disulfides.size(); i<=imax; ++i ) {
		//TR << "Old:\t" << old_disulfides[i].first << "\t" << old_disulfides[i].second << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.
		core::Size res1( old_disulfides[i].first + offset );
		if ( res1 > nres ) res1 -= nres;
		core::Size res2( old_disulfides[i].second + offset );
		if ( res2 > nres ) res2 -= nres;
		new_disulfides.push_back( std::pair< core::Size, core::Size >( res1, res2 ) );
		//TR << "New:\t" << new_disulfides[i].first << "\t" << new_disulfides[i].second << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.
	}

	return;
}

/// @brief Given a pose that has undergone an N-residue cyclic permutation, restore
/// the original pose, without the permutation.
void
SimpleCycpepPredictApplication::depermute (
	core::pose::PoseOP pose,
	core::Size const offset
) const {

	/// 1 2 3 4 5 6 7 8
	/// 2 3 4 5 6 7 8 1
	/// 3 4 5 6 7 8 1 2
	/// 4 5 6 7 8 1 2 3

	if ( offset==0 ) return; //Do nothing if the pose was not offset.

	core::Size const nres( sequence_length() );
	debug_assert(nres > offset);
	core::Size const old_first_res_index( nres-offset+1 );

	//TR << "nres=" << nres << " offset=" << offset << " old_first_res_index=" << old_first_res_index << std::endl; //DELETE ME

	//Store the old disulfides:
	utility::vector1 < std::pair < core::Size, core::Size > > old_disulfides;
	store_disulfides( pose, old_disulfides );
	//Break the old disulfides:
	break_disulfides( pose, old_disulfides );

	core::pose::PoseOP newpose( new core::pose::Pose );

	for ( core::Size ir=old_first_res_index; ir<=nres; ++ir ) {
		if ( ir == old_first_res_index ) {
			newpose->append_residue_by_jump( *(pose->residue(ir).clone()), 0, "", "", true );
		} else {
			newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
		}
	}

	for ( core::Size ir=1; ir<old_first_res_index; ++ir ) {
		newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
	}

	//Depermute the old disulfide list:
	utility::vector1 < std::pair < core::Size, core::Size > > new_disulfides;
	depermute_disulfide_list( old_disulfides, new_disulfides, offset, nres );

	//Re-form the disulfides:
	rebuild_disulfides( newpose, new_disulfides );

	//Re-append linker residues:
	if ( tbmb_positions_.size() > 0 ) {
		re_append_tbmb_residues( pose, newpose, offset );
	}

	//I don't bother to set up cyclic constraints, since we won't be doing any more minimization after calling this function.

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, newpose );
	termini->apply(*newpose);

	(*pose) = (*newpose);

	return;
}


/// @brief Align pose to native_pose, and return the RMSD between the two poses.
/// @details Assumes that the pose has already been de-permuted (i.e. the native and the pose line up).
/// Only uses alpha-amino acids for the alignment, currently.
core::Real
SimpleCycpepPredictApplication::align_and_calculate_rmsd(
	core::pose::PoseOP pose,
	core::pose::PoseCOP native_pose
) const {
	core::Size const nres( sequence_length() );
	core::Size res_counter(0); //Residue indices might not match between native pose and pose, due to linkers.

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, *pose, core::id::BOGUS_ATOM_ID);
	for ( core::Size ir=1, irmax=native_pose->total_residue(); ir<=irmax; ++ir ) {
		if ( !native_pose->residue_type(ir).is_alpha_aa() ) continue;
		++res_counter;
		runtime_assert_string_msg( res_counter <= nres, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): The native pose has more residues than the input sequence." );
		for ( core::Size ia=1, iamax=native_pose->residue(ir).type().first_sidechain_atom(); ia<iamax; ++ia ) { //Loop through all mainchain heavyatoms (including atoms coming off mainchain that are not sidechain atoms, like peptide "O").
			if ( native_pose->residue_type(ir).atom_is_hydrogen(ia) ) continue;
			//TR << "ir=" << ir << " ia=" << ia << " res_counter=" << res_counter << " native=" << native_pose->residue_type(ir).atom_name(ia) << " pred=" << pose->residue_type(res_counter).atom_name(ia) << std::endl; TR.flush(); //DELETE ME.
			runtime_assert_string_msg(
				!native_pose->residue_type(ir).atom_name(ia).compare( pose->residue_type(res_counter).atom_name(ia) ),
				"Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): Residue types or atom indices don't match between native and prediction."
			);
			amap[ core::id::AtomID(ia,res_counter) ] = core::id::AtomID(ia,ir);
			//TR << "Adding ia=" << ia << " ir=" << ir << " to map." << std::endl; //DELETE ME
		}
	}

	runtime_assert_string_msg( res_counter == nres, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): The native pose has fewer residues than the input sequence." );

	return core::scoring::superimpose_pose( *pose, *native_pose, amap ); //Superimpose the pose and return the RMSD.
}

/// @brief Create a new checkpoint file.
///
void
SimpleCycpepPredictApplication::new_checkpoint_file() const {
	using namespace utility::io;

	if ( !suppress_checkpoints_ ) {
		ozstream outfile;
		outfile.open( checkpoint_filename_ );
		runtime_assert(outfile.good());
		outfile << checkpoint_job_identifier_ << std::endl;
		outfile << "LAST\t0\tSUCCESS\t0" << std::endl;
		outfile.flush();
		outfile.close();
	}

	erase_random_seed_info();
	store_random_seed_info();

	return;
}


/// @brief Initialize checkpointing for this run.
/// @details  This function does several things.  First, it checks for an existing checkpoint
/// file.  If one exists, it checks whether the unique job name in the file matches the current
/// job.  If it does, then this job has already been attempted, and we're somewhere in the middle
/// of it.  The function reads the last attempt number and success count from the checkpoint
/// file, and returns these values.  Otherwise, it creates a new checkpoint file with the current
/// job name and returns (0,0).  If checkpointing is disabled, this function does nothing, and
/// returns (0,0).
/// @param[out] lastjob The index of the last job run.  Set to zero if checkpointing is disabled
/// or if we're creating a new checkpoint file (first job run).
/// @param[out] successes The number of successes so far.  Set to zero if checkpointing is
/// disabled or if we're creating a new checkpoint file (first job run).
void
SimpleCycpepPredictApplication::initialize_checkpointing(
	core::Size &lastjob,
	core::Size &successes
) const {
	using namespace utility::io;

	//If we're not using checkpointing, return 0,0:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		lastjob=0;
		successes=0;
		return;
	}

	//Check for a checkpoint file:
	izstream infile;
	infile.open( checkpoint_filename_ );
	if ( !infile.good() ) {
		//If the checkpoint file doesn't exist/isn't readable, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//If we've reached this point, then the checkpoint file IS good, and we need to read the job name and the last job ID/success count:
	std::string curline;
	infile.getline(curline);
	if ( infile.eof() || curline=="" || curline!=checkpoint_job_identifier_ ) {
		//If the checkpoint file isn't readable or is for a different job, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//Loop through the checkpoint file and read the job lines:
	while ( !infile.eof() ) {
		infile.getline(curline);
		std::istringstream ss(curline);
		ss >> curline;
		ss >> lastjob;
		ss >> curline;
		ss >> successes;
	}
	infile.close();

	get_random_seed_info();

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initialized job to " << lastjob << ", successes to " << successes << "." << std::endl;
		TR.Debug.flush();
	}

	return;
}

/// @brief Add a checkpoint to the checkpoint file.
/// @details  The checkpoint file must already exist.  Does nothing if checkpointing is disabled.
/// @param[in] curjob The index of the current job just run, for writing to the checkpoint file.
/// @param[in] successes The number of successes so far, for writing to the checkpoint file.
void
SimpleCycpepPredictApplication::checkpoint(
	core::Size const curjob,
	core::Size const successes
) const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		return;
	}

	ozstream outfile;
	outfile.open_append( checkpoint_filename_ );
	runtime_assert(outfile.good());
	outfile << "LAST\t" << curjob << "\tSUCCESS\t" << successes << std::endl;
	outfile.flush();
	outfile.close();

	store_random_seed_info();

#ifdef BOINC_GRAPHICS
	protocols::boinc::Boinc::update_pct_complete();
#endif

	return;

}

/// @brief End checkpointing and delete the checkpoint file.
/// @details Does nothing if checkpointing is disabled.
void
SimpleCycpepPredictApplication::end_checkpointing() const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		return;
	}
	runtime_assert( remove( checkpoint_filename_.c_str() ) == 0 );
	erase_random_seed_info();
	return;
}

/// @brief Restore the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::get_random_seed_info() const {
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	if ( utility::file::file_exists(rand_checkpoint_file_) ) {
		utility::io::izstream izs(rand_checkpoint_file_);
		numeric::random::rg().restoreState(izs);
		izs.close();
	}
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Store the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::store_random_seed_info() const {
	if ( suppress_checkpoints_ ) return; //Do nothing if we're not checkpointing
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	utility::io::ozstream ozs(rand_checkpoint_file_);
	numeric::random::rg().saveState(ozs);
	ozs.close();
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Erase the stored state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::erase_random_seed_info() const {
	if ( suppress_checkpoints_ ) return; //Do nothing if we're not checkpointing
	if ( utility::file::file_exists(rand_checkpoint_file_) ) {
		utility::file::file_delete(rand_checkpoint_file_);
	}
	return;
}

/// @brief Given a pose with TBMB in it and another pose without TBMB, copy the TBMB residues from the first to the second,
/// and add back covalent bonds.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back constraints.
void
SimpleCycpepPredictApplication::re_append_tbmb_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose,
	core::Size const offset
) const {
	debug_assert( pose->total_residue() > newpose->total_residue() );

	core::Size lastres(0);
	for ( core::Size i=1, imax=tbmb_positions_.size(); i<=imax; ++i ) { //Loop through all TBMBs
		//For each one, loop through the sequence and find the next TBMB.
		for ( core::Size j=lastres+1, jmax=pose->total_residue(); j<=jmax; ++j ) {
			if ( !pose->residue_type(j).name3().compare("TBM") ) {
				lastres = j;
				break;
			}
		}
		newpose->append_residue_by_jump( pose->residue(lastres), tbmb_positions_[i][1] ); //Jump from the first cys that links this TBMB to the TBMB
		protocols::cyclic_peptide::threefold_linker::TBMB_Helper helper;
		core::Size cys1, cys2, cys3;
		if ( offset >= tbmb_positions_[i][3] || offset < tbmb_positions_[i][1] ) {
			cys1 = tbmb_positions_[i][1]; cys2 = tbmb_positions_[i][2]; cys3 = tbmb_positions_[i][3];
		} else if ( offset >= tbmb_positions_[i][2] ) {
			cys1 = tbmb_positions_[i][3]; cys2 = tbmb_positions_[i][1]; cys3 = tbmb_positions_[i][2];
		} else {
			cys1 = tbmb_positions_[i][2]; cys2 = tbmb_positions_[i][3]; cys3 = tbmb_positions_[i][1];
		}
		helper.add_linker_bonds_asymmetric( *newpose, cys1, cys2, cys3, newpose->total_residue() );
	}
}


} //cyclic_peptide_predict
} //protocols
