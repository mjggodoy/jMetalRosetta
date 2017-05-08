// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockRbTrialMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/MonteCarlo.hh>

// Package headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/symmetric_docking/SymFoldandDockCreators.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace symmetric_docking {

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.symmetry.SymFoldandDockRbTrialMover" );

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover() :
	Mover( "SymFoldandDockRbTrialMover" ) {
	init();
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn
) : Mover( "SymFoldandDockRbTrialMover" ) {
	init();
	scorefxn_ = scorefxn;
}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move
) : Mover( "SymFoldandDockRbTrialMover" )
{
	init();
	scorefxn_ = scorefxn;
	smooth_move_ = smooth_move;
}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move,
	core::Real rot_mag,
	core::Real trans_mag
) : Mover( "SymFoldandDockRbTrialMover" )
{
	init();
	scorefxn_ = scorefxn;
	smooth_move_ = smooth_move;
	rot_mag_ = rot_mag;
	trans_mag_= trans_mag;
}

void
SymFoldandDockRbTrialMover::init() {
	using namespace basic::options;

	// fpd: the decision on whether to use smooth or nonsmooth is deferred to runtime
	//      for simplicity only expose nonsmooth parameters to RS
	smooth_move_ = false; // don't expose to RS
	rot_mag_ = option[ OptionKeys::fold_and_dock::rb_rot_magnitude ]();
	trans_mag_ = option[ OptionKeys::fold_and_dock::rb_trans_magnitude ]();
	rot_mag_smooth_ = option[ OptionKeys::fold_and_dock::rot_mag_smooth ]();  // don't expose to RS
	trans_mag_smooth_ = option[ OptionKeys::fold_and_dock::trans_mag_smooth ]();  // don't expose to RS
	rigid_body_cycles_ = option[ OptionKeys::fold_and_dock::rigid_body_cycles ]();
	mc_filter_ = !option[ OptionKeys::fold_and_dock::rigid_body_disable_mc ]();
	rotate_anchor_to_x_ = option[ OptionKeys::fold_and_dock::rotate_anchor_to_x ]();
}

void
SymFoldandDockRbTrialMover::apply( core::pose::Pose & pose )
{
	protocols::simple_moves::symmetry::SetupForSymmetryMover setup;
	setup.apply( pose );

	using namespace core::conformation::symmetry;
	using namespace basic::options;
	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	TR.Debug << "Rb move applied..." << std::endl;

	// Docking options
	core::Real rot_mag_trial = smooth_move_ ? rot_mag_smooth_ : rot_mag_;
	core::Real trans_mag_trial = smooth_move_ ? trans_mag_smooth_ : trans_mag_;

	// Setup Monte Carlo object
	protocols::moves::MonteCarloOP monteCarlo_( new protocols::moves::MonteCarlo(pose, *scorefxn_, 2.0 ) );

	// set up mover for docking
	protocols::rigid::RigidBodyDofSeqPerturbMover rb_perturb =
		protocols::rigid::RigidBodyDofSeqPerturbMover( dofs , rot_mag_trial, trans_mag_trial );

	if ( rotate_anchor_to_x_ ) {
		core::pose::symmetry::rotate_anchor_to_x_axis( pose );
	}

	for ( Size i = 1; i <= rigid_body_cycles_; ++i ) {
		rb_perturb.apply( pose );
		if ( mc_filter_ ) monteCarlo_->boltzmann( pose );
	}
}

void
SymFoldandDockRbTrialMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	using namespace core::scoring;

	if ( tag->hasOption( "rot_mag" ) ) {
		rot_mag_ = tag->getOption<core::Real>( "rot_mag" );
	}
	if ( tag->hasOption( "trans_mag" ) ) {
		trans_mag_ = tag->getOption<core::Real>( "trans_mag" );
	}
	if ( tag->hasOption( "cycles" ) ) {
		rigid_body_cycles_ = tag->getOption<int>( "cycles" );
	}
	if ( tag->hasOption( "use_mc" ) ) {
		mc_filter_ = !tag->getOption<bool>( "use_mc" );
	}
	if ( tag->hasOption( "rotate_anchor_to_x" ) ) {
		rotate_anchor_to_x_ = tag->getOption<bool>( "rotate_anchor_to_x" );
	}
}

// XRW TEMP std::string
// XRW TEMP SymFoldandDockRbTrialMover::get_name() const {
// XRW TEMP  return "SymFoldandDockRbTrialMover";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP SymFoldandDockRbTrialMoverCreator::keyname() const {
// XRW TEMP  return SymFoldandDockRbTrialMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SymFoldandDockRbTrialMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SymFoldandDockRbTrialMover() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SymFoldandDockRbTrialMover::mover_name() {
// XRW TEMP  return "SymFoldandDockRbTrialMover";
// XRW TEMP }

std::string SymFoldandDockRbTrialMover::get_name() const {
	return mover_name();
}

std::string SymFoldandDockRbTrialMover::mover_name() {
	return "SymFoldandDockRbTrialMover";
}

void SymFoldandDockRbTrialMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "rot_mag", xs_decimal, "Rotational perturbation size." )
		+ XMLSchemaAttribute( "trans_mag", xs_decimal, "Translation perturbation size." )
		+ XMLSchemaAttribute( "cycles", xsct_non_negative_integer, "Number of cycles." )
		+ XMLSchemaAttribute( "use_mc", xsct_rosetta_bool, "Use Monte Carlo?" )
		+ XMLSchemaAttribute( "rotate_anchor_to_x", xsct_rosetta_bool, "Rotate the anchor residue to the x-axis before applying rigid body transformations." ) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Something else for setting up symmetric fold and dock.", attlist );
}

std::string SymFoldandDockRbTrialMoverCreator::keyname() const {
	return SymFoldandDockRbTrialMover::mover_name();
}

protocols::moves::MoverOP
SymFoldandDockRbTrialMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymFoldandDockRbTrialMover );
}

void SymFoldandDockRbTrialMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymFoldandDockRbTrialMover::provide_xml_schema( xsd );
}



} // symmetric_docking
} // protocols
