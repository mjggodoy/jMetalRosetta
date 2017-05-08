// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/ConstraintPreparer.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/ConstraintPreparer.hh>
#include <protocols/abinitio/abscript/ConstraintPreparerCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/EnvClaim.hh>

// Project headers
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <core/id/Exceptions.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/LoopsFileIO.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

//Req'd on WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.abscript.ConstraintPreparer", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP ConstraintPreparerCreator::keyname() const {
// XRW TEMP  return ConstraintPreparer::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ConstraintPreparerCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ConstraintPreparer );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ConstraintPreparer::mover_name() {
// XRW TEMP  return "ConstraintPreparer";
// XRW TEMP }

ConstraintPreparer::ConstraintPreparer():
	Parent(),
	combine_ratio_( 1 ),
	skip_redundant_( false ),
	skip_redundant_width_( 1 ),
	rand_drop_rate_( 0.0 ),
	reprepare_( false ),
	combine_exclude_res_(/* 0 */), //initialize vector of zero size, will get resized later.
	filename_(""),
	constraints_( /* NULL */ )
{}

claims::EnvClaims ConstraintPreparer::yield_claims( core::pose::Pose const&,
	basic::datacache::WriteableCacheableMapOP ) {
	claims::EnvClaims e;
	return e;
}

void ConstraintPreparer::prepare( core::pose::Pose& pose, core::Real ){
	using namespace core::scoring::constraints;

	if ( constraints_ && !reprepare_ ) {
		// If constraints isn't null, we already have constraints and don't need to reprepare
		// unless reprepare option is set.
		return;
	}

	try{
		// it's not great that this is re-loaded each time, but it's safe and doesn't get called that often.
		constraints_ = ConstraintIO::get_instance()->read_constraints( cst_file(), ConstraintSetOP( new ConstraintSet ), pose );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		throw utility::excn::EXCN_BadInput( get_name() + " encountered problem loading constraint file '"
			+ cst_file() + "' : " + e.msg() );
	}

// we don't know the correct size of this vector until we see the pose.
	combine_exclude_res_.resize( pose.size(), false );

	ConstraintCOPs added_constraints = constraints_->get_all_constraints();
	if ( skip_redundant() ) {
		skip_redundant_constraints( added_constraints, pose.size(), skip_redundant_width() );
	}
	if ( random_drop_rate() > 0.0 ) {
		drop_constraints( added_constraints, random_drop_rate() );
	}

	core::kinematics::ShortestPathInFoldTree sp( pose.fold_tree() );
	choose_effective_sequence_separation( sp, added_constraints );

	// if combine_ratio_ > 1 this will randomly combine constraints into multipletts with OR logic
	combine_constraints( added_constraints, combine_ratio(), combine_exclude_res_, sp );

	pose.add_constraints( added_constraints );

	if ( tr.Trace.visible() ) {
		pose.constraint_set()->show_definition( tr.Trace, pose );
	}

}

void ConstraintPreparer::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	skip_redundant_width( tag->getOption< core::Size >( "skip_redundant", 1 ) );
	if ( tag->hasOption( "skip_redundant" ) ) {
		skip_redundant( true );
		if ( option[ constraints::skip_redundant ].user() &&
				option[ constraints::skip_redundant ] == false ) {
			tr.Error << "Command line option constraints::skip_redundant conflicts "
				<< "with RosettaScripts initialization of ConstraintPreparer "
				<< tag->getOption< std::string >( "name", "null" ) << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( "ConstraintPreparer " + this->get_name() +
				"skip_redundant conflict" );
		} else if ( option[ constraints::skip_redundant_width ].user() &&
				(Size) option[ constraints::skip_redundant_width ] != skip_redundant_width() ) {
			tr.Error << "Command line option constraints::skip_redundant_width conflicts "
				<< "with RosettaScripts initialization of ConstraintPreparer "
				<< tag->getOption< std::string >( "name", "null" ) << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( "ConstraintPreparer " + this->get_name() +
				"skip_redundant_width conflict" );
		}
	} else {
		skip_redundant( false );
	}

	reprepare_ = tag->getOption< bool >( "reprepare", false );

	tr.Debug << "Set skip_redundant to " << skip_redundant() << " with width "
		<< skip_redundant_width() << std::endl;

	if ( tag->hasOption( "cst_file" ) ) {
		cst_file( tag->getOption< std::string >( "cst_file", "null" ) );
	} else {
		tr.Error << "ConstraintPreparers require a constraint file.  "
			<< "Use the 'cst_file' tag in the mover definition. " << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("ConstraintPreparer has no constraint file.");
	}

	combine_ratio( tag->getOption< core::Size >( "combine_ratio", 1 ) );
	if ( option[ constraints::combine ].user() ) {
		tr.Warning << "User-set option -constraints:combine being ignored" << std::endl;
	}

}

void ConstraintPreparer::cst_file( std::string const& filename ){
	constraints_ = NULL;
	filename_ = filename;
}

bool loop_stop_comp( loops::Loop const& lhs,
	loops::Loop const& rhs ){
	return ( lhs.stop() < rhs.stop() );
}

void ConstraintPreparer::combine_exclude_file( std::string const& filename ){
	std::ifstream is( filename.c_str() );

	if ( !is.good() ) {
		utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + filename + "'" );
	}

	loops::PoseNumberedLoopFileReader reader;
	reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
	loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(is, filename, false );
	loops::Loops rigid_core = loops::Loops( loops );

	core::Size maxres = std::max_element( rigid_core.v_begin(), rigid_core.v_end(), loop_stop_comp )->stop();

	//initially, the vector is only as large as the largest "true" set value
	combine_exclude_res_.resize( maxres, false );
	rigid_core.transfer_to_residue_vector( combine_exclude_res_, true );
}

// XRW TEMP std::string ConstraintPreparer::get_name() const {
// XRW TEMP  return "ConstraintPreparer";
// XRW TEMP }

std::string ConstraintPreparer::get_name() const {
	return mover_name();
}

std::string ConstraintPreparer::mover_name() {
	return "ConstraintPreparer";
}

void ConstraintPreparer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"skip_redundant", xsct_non_negative_integer,
		"Skip redundant constraints? Corresponds to skip_redundant_width option", "1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"name", xs_string,
		"Name identifying this ConstraintPreparer", "null");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"reprepare", xsct_rosetta_bool,
		"Should constraints be reprepared if some have already been set?", "false");
	attlist + XMLSchemaAttribute(
		"cst_file", xs_string,
		"Name of file containing constraints to apply to pose");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"combine_ratio", xsct_non_negative_integer,
		"If greater than 1, it will randomly combine constraints into multiplets of this size with OR logic", "1");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"StagePreparer that adds the specified constraints to the constraint set.",
		attlist );
}

std::string ConstraintPreparerCreator::keyname() const {
	return ConstraintPreparer::mover_name();
}

protocols::moves::MoverOP
ConstraintPreparerCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConstraintPreparer );
}

void ConstraintPreparerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConstraintPreparer::provide_xml_schema( xsd );
}



} // abscript
} // abinitio
} // protocols
