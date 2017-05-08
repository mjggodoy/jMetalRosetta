// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/ThreefoldLinkerMover.cc
/// @brief This mover links three cysteine residues with a three-way cross-linker.  It adds the crosslinker,
/// sets up constraints, optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker
/// and the side-chains to which it connects), andthen optionally relaxes the whole structure.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/ThreefoldLinkerMover.hh>
#include <protocols/cyclic_peptide/ThreefoldLinkerMoverCreator.hh>
#include <protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.hh>
#include <protocols/cyclic_peptide/threefold_linker/TBMB_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>

// Protocols headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/mover_schemas.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.ThreefoldLinkerMover" );

namespace protocols {
namespace cyclic_peptide {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ThreefoldLinkerMover::ThreefoldLinkerMover():
	protocols::moves::Mover( ThreefoldLinkerMover::class_name() ),
	residue_selector_(), //Defaults to NULL pointer.
	linker_(no_crosslinker), //Defaults to no_crosslinker
	add_linker_(true),
	constrain_linker_(true),
	pack_and_minimize_linker_and_sidechains_(true),
	do_final_fastrelax_(false),
	threefold_symmetric_(false),
	sfxn_(),
	sidechain_frlx_rounds_(3),
	final_frlx_rounds_(3),
	filter_by_sidechain_distance_(true),
	filter_by_constraints_energy_(true),
	filter_by_total_score_(false),
	filter_by_total_score_cutoff_energy_(0.0),
	sidechain_distance_filter_multiplier_(1.0),
	constraints_energy_filter_multiplier_(1.0)
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
ThreefoldLinkerMover::ThreefoldLinkerMover( ThreefoldLinkerMover const & src ):
	protocols::moves::Mover( src ),
	residue_selector_(src.residue_selector_),
	linker_(src.linker_),
	add_linker_(src.add_linker_),
	constrain_linker_(src.constrain_linker_),
	pack_and_minimize_linker_and_sidechains_(src.pack_and_minimize_linker_and_sidechains_),
	do_final_fastrelax_(src.do_final_fastrelax_),
	threefold_symmetric_(src.threefold_symmetric_),
	sfxn_(src.sfxn_),
	sidechain_frlx_rounds_(src.sidechain_frlx_rounds_),
	final_frlx_rounds_(src.final_frlx_rounds_),
	filter_by_sidechain_distance_(src.filter_by_sidechain_distance_),
	filter_by_constraints_energy_(src.filter_by_constraints_energy_),
	filter_by_total_score_(src.filter_by_total_score_),
	filter_by_total_score_cutoff_energy_(src.filter_by_total_score_cutoff_energy_),
	sidechain_distance_filter_multiplier_(src.sidechain_distance_filter_multiplier_),
	constraints_energy_filter_multiplier_(src.constraints_energy_filter_multiplier_)
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
ThreefoldLinkerMover::~ThreefoldLinkerMover(){}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ThreefoldLinkerMover::apply( core::pose::Pose& pose){
	//Check for a residue selector, then apply it to the pose:
	runtime_assert_string_msg( residue_selector(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::apply():  A residue selector must be specified before calling this function." );
	runtime_assert_string_msg( scorefxn(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::apply():  A scorefunction must be specified before calling this function." );

	core::select::residue_selector::ResidueSubset const selection( residue_selector()->apply(pose) );

	//Check that we've selected exactly three residues:
	runtime_assert_string_msg( exactly_three_selected( selection ), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::apply():  The residue selector did not select exactly three residues." );

	//Create the helper, which has the functions that set up specific types of crosslinkers:
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperOP helper;
	switch( linker_ ) {
	case TBMB :
		helper = protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperOP( protocols::cyclic_peptide::threefold_linker::TBMB_HelperOP( new protocols::cyclic_peptide::threefold_linker::TBMB_Helper ) );
		break;
	default :
		utility_exit_with_message( "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::apply(): Invalid crosslinker specified." );
	}

	// Decide whether to call symmetric or asymmetric functions from here on:
	if ( core::conformation::symmetry::is_symmetric( pose.conformation() ) ) {
		symmetric_apply( pose, selection, helper );
	} else {
		asymmetric_apply( pose, selection, helper );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
////////////////////////////////////////////////////////////////////////////////
void
ThreefoldLinkerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

//////////////////////////////
/// RosettaScripts Support ///
//////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ThreefoldLinkerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	runtime_assert_string_msg( tag->hasOption("residue_selector"), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::parse_my_tag(): A residue selector MUST be supplied with the \"residue_selector\" option." );
	set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, datamap ) );
	runtime_assert_string_msg( residue_selector(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::parse_my_tag(): Setting residue selector failed." );

	runtime_assert_string_msg( tag->hasOption("linker_name"), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::parse_my_tag(): The name of the linker residue MUST be supplied with the \"linker_name\" option." );
	set_linker_name( tag->getOption<std::string>("linker_name") );

	set_behaviour(
		tag->getOption<bool>( "add_linker", add_linker() ),
		tag->getOption<bool>( "constrain_linker", constrain_linker() ),
		tag->getOption<bool>( "pack_and_minimize_linker_and_sidechains", pack_and_minimize_linker_and_sidechains() ),
		tag->getOption<bool>( "do_final_fastrelax", do_final_fastrelax() ),
		tag->getOption<bool>( "threefold_symmetric", threefold_symmetric() )
	);

	set_filter_behaviour(
		tag->getOption<bool>( "filter_by_sidechain_distance", filter_by_sidechain_distance() ),
		tag->getOption<bool>( "filter_by_constraints_energy", filter_by_constraints_energy() ),
		tag->getOption<bool>( "filter_by_final_energy", filter_by_total_score() ),
		tag->getOption<core::Real>( "final_energy_cutoff", filter_by_total_score_cutoff_energy() ),
		tag->getOption<core::Real>( "sidechain_distance_filter_multiplier", sidechain_distance_filter_multiplier() ),
		tag->getOption<core::Real>( "constraints_energy_filter_multiplier", constraints_energy_filter_multiplier() )
	);

	runtime_assert_string_msg( tag->hasOption("scorefxn"), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::parse_my_tag(): A scorefunction must be supplied with the \"scorefxn\" option." );
	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );

	set_sidechain_frlx_rounds( tag->getOption<core::Size>("sidechain_fastrelax_rounds", sidechain_frlx_rounds()) );
	set_final_frlx_rounds( tag->getOption<core::Size>("final_fastrelax_rounds", final_frlx_rounds()) );
}

/// @brief Provide information on what options are available in XML tag.
///
void
ThreefoldLinkerMover::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "A unique name for this instance of the ThreefoldLinkerMover." )
		+ XMLSchemaAttribute::required_attribute( "linker_name", xs_string, "The name of the type of linker to use (e.g. TBMB for 1,3,5-tris(bromomethyl)benzene).")
		+ XMLSchemaAttribute( "add_linker", xsct_rosetta_bool, "Should the linker geometry be added to the pose?  Default true." )
		+ XMLSchemaAttribute( "constrain_linker", xsct_rosetta_bool, "Should constraints for the linker be added to the pose?  Default true." )
		+ XMLSchemaAttribute( "pack_and_minimize_linker_and_sidechains", xsct_rosetta_bool, "Should the linker and the connecting sidechains be repacked, and should the jump to the linker, and the linker and connnecting side-chain degrees of torsional freedom, be energy-minimized?  Default true." )
		+ XMLSchemaAttribute( "sidechain_fastrelax_rounds", xs_integer, "The number of rounds of FastRelax to apply when packing and minimizing side-chains and the liker.  Default 3." )
		+ XMLSchemaAttribute( "do_final_fastrelax", xsct_rosetta_bool, "Should the whole pose be subjected to a FastRelax?  Default false." )
		+ XMLSchemaAttribute( "final_fastrelax_rounds", xs_integer, "The number of rounds of FastRelax to apply when relaxing the whole pose.  Default 3." )
		+ XMLSchemaAttribute( "threefold_symmetric", xsct_rosetta_bool, "Is this a threefold-symmetric pose being connected by a threefold-symmetric crosslinker?  Default false." )
		+ XMLSchemaAttribute( "filter_by_sidechain_distance", xsct_rosetta_bool, "Prior to adding the linker geometry, should this mover abort with failure status if the selected side-chains are too far apart to connect to the linker?  Default true." )
		+ XMLSchemaAttribute( "sidechain_distance_filter_multiplier", xsct_real, "This is a multiplier for the sidechain distance cutoff filter.  Higher values make the filter less stringent.  Default 1.0." )
		+ XMLSchemaAttribute( "filter_by_constraints_energy", xsct_rosetta_bool, "After adding the linker geometry, adding constraints, and repacking and minimizing the linker and the connecting side-chains, should ths mover abort with failure status if the constraints energy is too high (i.e. the energy-minimized linker geometry is bad)?  Default true." )
		+ XMLSchemaAttribute( "constraints_energy_filter_multiplier", xsct_real, "This is a multiplier for the constraints energy cutoff filter.  Higher values make the filter less stringent.  Default 1.0." )
		+ XMLSchemaAttribute( "filter_by_final_energy", xsct_rosetta_bool, "At the end of this protocol, should this mover exit with error status if the final energy is above a user-defined cutoff?  Default false." )
		+ XMLSchemaAttribute( "final_energy_cutoff", xsct_real, "If we are exiting with error status if the final energy is too high, this is the energy cutoff.  Default 0.0." )
		;
	//get attributes for parse residue selector
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "A previously-defined residue selector that has been set up to select exactly three residues." );
	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description_when_required(attlist, "scorefxn", "A scorefunction to use for packing, energy-minimization, and filtering.  If constraints are turned off in this score function, they will be turned on automatically at apply time." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds a three-way crosslinker linking three user-specified side-chains.", attlist );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
////////////////////////////////////////////////////////////////////////////////
moves::MoverOP
ThreefoldLinkerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ThreefoldLinkerMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ThreefoldLinkerMover::clone() const
{
	return protocols::moves::MoverOP( new ThreefoldLinkerMover( *this ) );
}

/// @brief Get the name of the Mover
std::string
ThreefoldLinkerMover::get_name() const
{
	return ThreefoldLinkerMover::class_name();
}

std::string
ThreefoldLinkerMover::class_name()
{
	return "ThreefoldLinkerMover";
}

/// @brief Returns the name of this Mover.
std::string
ThreefoldLinkerMover::mover_name() {
	return "ThreefoldLinkerMover";
}


/// @brief Set the residue selector to use.
///
void
ThreefoldLinkerMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in, "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_residue_selector(): A null pointer was passed to this function." );
	residue_selector_ = selector_in;
}

/// @brief Set the linker name.
///
void
ThreefoldLinkerMover::set_linker_name(
	std::string const &name_in
) {
	runtime_assert_string_msg( !name_in.empty(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_linker_name(): An empty string was passed to this function." );
	CrossLinker linker_in( get_crosslinker_enum( name_in ) );
	runtime_assert_string_msg( linker_in != unknown_crosslinker, "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_linker_name(): Could not interpret the linker name \"" + name_in + "\"." );
	linker_ = linker_in;
}

/// @brief Get the linker name.
///
std::string
ThreefoldLinkerMover::linker_name() const {
	return get_crosslinker_name( linker_ );
}

/// @brief Set the behaviour of this mover.
///
void
ThreefoldLinkerMover::set_behaviour(
	bool const add_linker,
	bool const constrain_linker,
	bool const pack_and_minimize_linker_and_sidechains,
	bool const do_final_fastrelax,
	bool const threefold_symmetric
) {
	add_linker_ = add_linker;
	constrain_linker_ = constrain_linker;
	pack_and_minimize_linker_and_sidechains_ = pack_and_minimize_linker_and_sidechains;
	do_final_fastrelax_ = do_final_fastrelax;
	threefold_symmetric_ = threefold_symmetric;
}

/// @brief Set the filtering behaviour of this mover.
///
void
ThreefoldLinkerMover::set_filter_behaviour(
	bool const filter_by_sidechain_distance,
	bool const filter_by_constraints_energy,
	bool const filter_by_total_score,
	core::Real const &filter_by_total_score_cutoff_energy,
	core::Real const &sidechain_distance_filter_multiplier,
	core::Real const &constraints_energy_filter_multiplier
) {
	filter_by_sidechain_distance_ = filter_by_sidechain_distance;
	filter_by_constraints_energy_ = filter_by_constraints_energy;
	filter_by_total_score_ = filter_by_total_score;
	filter_by_total_score_cutoff_energy_ = filter_by_total_score_cutoff_energy;
	sidechain_distance_filter_multiplier_ = sidechain_distance_filter_multiplier;
	constraints_energy_filter_multiplier_ = constraints_energy_filter_multiplier;
}

/// @brief Set the scorefunction to use for packing and minimization.
/// @details Cloned at apply time.  (That is, the scorefunction is shared until apply time).
void
ThreefoldLinkerMover::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	runtime_assert_string_msg( sfxn_in, "Error in Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_scorefxn(): A null pointer was passed to this function." );
	sfxn_ = sfxn_in;
}

/// @brief Get the scorefunction to use for packing and minimization.
///
core::scoring::ScoreFunctionCOP
ThreefoldLinkerMover::scorefxn() const {
	return sfxn_;
}

/// @brief Set the number of rounds of FastRelax to apply when minimizing the linker and the
/// side-chains that connect to it.
void
ThreefoldLinkerMover::set_sidechain_frlx_rounds(
	core::Size const rounds_in
) {
	runtime_assert_string_msg(rounds_in > 0, "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_sidechain_frlx_rounds(): The number of rounds must be greater than zero.");
	sidechain_frlx_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of FastRelax to apply at the end.
///
void
ThreefoldLinkerMover::set_final_frlx_rounds(
	core::Size const rounds_in
) {
	runtime_assert_string_msg(rounds_in > 0, "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::set_final_frlx_rounds(): The number of rounds must be greater than zero.");
	final_frlx_rounds_ = rounds_in;
}

std::ostream &
operator<<( std::ostream & os, ThreefoldLinkerMover const & mover )
{
	mover.show(os);
	return os;
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
ThreefoldLinkerMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ThreefoldLinkerMover );
}

std::string
ThreefoldLinkerMoverCreator::keyname() const
{
	return ThreefoldLinkerMover::class_name();
}

void
ThreefoldLinkerMoverCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	ThreefoldLinkerMover::provide_xml_schema( xsd );
}


///////////////////////
/// private methods ///
///////////////////////

/// @brief Apply the mover to a symmetric pose.
/// @details Requires threefold symmetry in the pose, and threefold_symmetric_ = true.
void
ThreefoldLinkerMover::symmetric_apply(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) {
	runtime_assert_string_msg(threefold_symmetric(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::symmetric_apply(): The function was called, but the mover's \"threefold_symmetric\" option is not set.");
	runtime_assert_string_msg( helper->selection_is_symmetric( selection, pose ), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::symmetric_apply(): The selector does not select three equivalent, symmetric residues." );

	core::pose::PoseOP pose_copy( pose.clone() );

	bool failed(false);
	if ( filter_by_sidechain_distance() ) {
		failed=filter_by_sidechain_distance_symmetric( *pose_copy, selection, helper );
		TR << "Symmetric sidechain distance filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
	}
	if ( !failed ) {
		if ( add_linker() ) {
			add_linker_symmetric( *pose_copy, selection, helper );
		}
		if ( constrain_linker() ) {
			add_linker_constraints_symmetric( *pose_copy, selection, helper, add_linker() );
		}
		if ( pack_and_minimize_linker_and_sidechains() ) {
			pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, false, true);
		}
		if ( filter_by_constraints_energy() ) {
			failed = filter_by_constraints_energy_symmetric( *pose_copy, selection, helper, add_linker() );
			TR << "Symmetric linker constraints filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
		}
		if ( !failed ) {
			if ( do_final_fastrelax() ) {
				pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, true, true);
				if ( filter_by_constraints_energy() ) {
					failed = filter_by_constraints_energy_symmetric( *pose_copy, selection, helper, add_linker() );
					TR << "Symmetric linker constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
			if ( !failed ) {
				if ( filter_by_total_score() ) {
					failed = filter_by_total_score( *pose_copy );
					TR << "Symmetric total score + constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
		}
	}

	if ( !failed ) {
		TR << "Symmetric ThreefoldLinkerMover reports SUCCESS.  Updating pose." << std::endl;
		pose = *pose_copy;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << "Symmetric ThreefoldLinkerMover reports FAILURE.  Returning input pose." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}

}

/// @brief Determine whether the residues to be crosslinked are too far apart.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too far apart), FALSE for success.
bool
ThreefoldLinkerMover::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	return helper->filter_by_sidechain_distance_symmetric( pose, selection, sidechain_distance_filter_multiplier() );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
ThreefoldLinkerMover::filter_by_constraints_energy_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper,
	bool const linker_was_added
) const {
	return helper->filter_by_constraints_energy_symmetric( pose, selection, linker_was_added, constraints_energy_filter_multiplier() );
}

/// @brief Given a selection of exactly three residues, add a crosslinker, align it crudely to the
/// selected residues, and set up covalent bonds.  This version is for symmetric poses.
void
ThreefoldLinkerMover::add_linker_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	helper->add_linker_symmetric(pose, selection);
}

/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.  This version is for symmetric poses.
void
ThreefoldLinkerMover::add_linker_constraints_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper,
	bool const linker_was_added
) const {
	helper->add_linker_constraints_symmetric( pose, selection, linker_was_added );
}

/// @brief Apply the mover to an asymmetric pose.
/// @details Requires and asymmetric pose, and threefold_symmetric_ = false.
void
ThreefoldLinkerMover::asymmetric_apply(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) {
	runtime_assert_string_msg(!threefold_symmetric(), "Error in protocols::cyclic_peptide::ThreefoldLinkerMover::asymmetric_apply(): The function was called, but the mover's \"threefold_symmetric\" option is set.");

	core::pose::PoseOP pose_copy( pose.clone() );

	bool failed(false);

	if ( filter_by_sidechain_distance() ) {
		failed = filter_by_sidechain_distance_asymmetric( *pose_copy, selection, helper );
		TR << "Sidechain distance filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
	}
	if ( !failed ) {
		if ( add_linker() ) {
			add_linker_asymmetric( *pose_copy, selection, helper );
		}
		if ( constrain_linker() ) {
			add_linker_constraints_asymmetric( *pose_copy, selection, helper );
		}
		if ( pack_and_minimize_linker_and_sidechains() ) {
			pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, false, false);
		}
		if ( filter_by_constraints_energy() ) {
			failed = filter_by_constraints_energy_asymmetric( *pose_copy, selection, helper );
			TR << "Linker constraints filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
		}
		if ( !failed ) {
			if ( do_final_fastrelax() ) {
				pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, true, false);
				if ( filter_by_constraints_energy() ) {
					failed = filter_by_constraints_energy_asymmetric( *pose_copy, selection, helper );
					TR << "Linker constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
			if ( !failed ) {
				if ( filter_by_total_score() ) {
					failed = filter_by_total_score( *pose_copy );
					TR << "Total score + constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
		}
	}

	if ( !failed ) {
		TR << "ThreefoldLinkerMover reports SUCCESS.  Updating pose." << std::endl;
		pose = *pose_copy;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << "ThreefoldLinkerMover reports FAILURE.  Returning input pose." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}
}

/// @brief Determine whether the residues to be crosslinked are too far apart.
/// @details Returns TRUE for failure (too far apart), FALSE for success.
bool
ThreefoldLinkerMover::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	return helper->filter_by_sidechain_distance_asymmetric( pose, selection, sidechain_distance_filter_multiplier() );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
ThreefoldLinkerMover::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	return helper->filter_by_constraints_energy_asymmetric( pose, selection, constraints_energy_filter_multiplier() );
}

/// @brief Determine whether the overall system has too high an overall score (including constraints) at the end of the protocol.
/// @details Returns TRUE for failure (too high an overall score) and FALSE for success.
bool
ThreefoldLinkerMover::filter_by_total_score(
	core::pose::Pose const &pose
) const {
	core::pose::Pose pose_copy(pose);

	core::scoring::ScoreFunctionOP sfxn( scorefxn()->clone() ); //Local copy of scorefunction.
	if ( sfxn->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::dihedral_constraint, 1.0); }

	(*sfxn)(pose_copy);
	core::Real const total_eng( pose_copy.energies().total_energy() );
	bool const failed( total_eng > filter_by_total_score_cutoff_energy() );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Total energy (including constraints) at end of ThreefoldLinkerMover protocol: " << total_eng << std::endl;
		TR.Debug << "Energy cutoff for ThreefoldLinkerMover protocol: " << filter_by_total_score_cutoff_energy() << std::endl;
		if ( failed ) { TR.Debug << "Energy filter FAILED." << std::endl; }
		else { TR.Debug << "Energy filter PASSED." << std::endl; }
	}

	return failed;
}

/// @brief Given a selection of exactly three residues, add a crosslinker, align it crudely to the
/// selected residues, and set up covalent bonds.
void
ThreefoldLinkerMover::add_linker_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	helper->add_linker_asymmetric(pose, selection);
}

/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
void
ThreefoldLinkerMover::add_linker_constraints_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper
) const {
	helper->add_linker_constraints_asymmetric( pose, selection );
}

/// @brief Repack and minimize the sidechains.
/// @details Also repacks and minimzes the linker, letting all jumps vary.
void
ThreefoldLinkerMover::pack_and_minimize_linker_and_sidechains(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper,
	bool const whole_structure,
	bool const symmetric
) const {

	core::scoring::ScoreFunctionOP sfxn( scorefxn()->clone() ); //Local copy of scorefunction.
	if ( sfxn->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::dihedral_constraint, 1.0); }

	protocols::relax::FastRelax frlx( sfxn, ( whole_structure ? final_frlx_rounds() : sidechain_frlx_rounds() ) );

	if ( !whole_structure ) {
		core::Size res1, res2, res3, linker_index1, linker_index2, linker_index3;

		helper->get_sidechain_indices( selection, res1, res2, res3 );
		if ( symmetric ) {
			if ( add_linker() ) {
				res2 += 1;
				res3 += 2;
			}
			helper->get_linker_indices_symmetric(pose, res1, res2, res3, linker_index1, linker_index2, linker_index3);
		} else {
			linker_index1 = helper->get_linker_index_asymmetric( pose, res1, res2, res3);
		}

		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_bb(false);
		movemap->set_chi(false);
		movemap->set_jump(true);
		movemap->set_chi( res1, true );
		movemap->set_chi( res2, true );
		movemap->set_chi( res3, true );
		movemap->set_chi( linker_index1, true );
		if ( symmetric ) {
			movemap->set_chi( linker_index2, true );
			movemap->set_chi( linker_index3, true );
		}

		frlx.set_movemap(movemap);
	}

	frlx.apply(pose);

}

/// @brief Given a selection from a ResidueSelector, check that exactly three residues have been selected.
/// @details Returns true if exactly three have been selected, false otherwise.
bool
ThreefoldLinkerMover::exactly_three_selected(
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size const nres( selection.size() );
	if ( nres< 3 ) return false;
	core::Size count(0);
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( selection[i] ) ++count;
	}
	return (count == 3);
}

/// @brief Given a CrossLinker enum, get its name.
///
std::string
ThreefoldLinkerMover::get_crosslinker_name(
	CrossLinker const crosslinker
) const {
	switch( crosslinker) {
	case no_crosslinker :
		return "no_crosslinker";
	case TBMB :
		return "TBMB";
	default :
		break;
	}
	return "unknown_crosslinker";
}

/// @brief Given a CrossLinker name, get its enum.
///
CrossLinker
ThreefoldLinkerMover::get_crosslinker_enum(
	std::string const &name
) const {
	for ( core::Size i(2); i < static_cast<core::Size>(end_of_crosslinker_list); ++i ) {
		if ( !name.compare( get_crosslinker_name( static_cast<CrossLinker>(i) ) ) ) {
			return static_cast<CrossLinker>(i);
		}
	}
	return unknown_crosslinker;
}

} //protocols
} //cyclic_peptide
