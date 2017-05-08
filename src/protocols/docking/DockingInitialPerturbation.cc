// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingInitialPerturbation.cc
/// @brief initial position functions
/// @details
///  This contains the functions that create initial positions for docking
///  You can either randomize partner 1 or partner 2, spin partner 2, or
///  perform a simple perturbation.
///  Also contains docking mcm protocol
/// @author Monica Berrondo

#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingInitialPerturbationCreator.hh> // zhe
#include <protocols/docking/metrics.hh>

// Rosetta Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/RigidBodyInfo.hh> // zhe
#include <protocols/docking/EllipsoidalRandomizationMover.hh> // NM
#include <protocols/scoring/Interface.hh>
#include <basic/datacache/DataMap.hh> // zhe

#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>

#include <core/kinematics/FoldTree.hh>

#include <basic/options/option.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

#include <basic/options/keys/mp.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers
#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>

#include <protocols/membrane/util.hh>

using basic::T;

// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#undef interface
#endif

using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockingInitialPerturbation" );
static core::Size trans ( 1 ), rot ( 2 );

using namespace core;

namespace protocols {
namespace docking {


// Creator part for DockingInitialPerturbation, used in scripts
// XRW TEMP std::string
// XRW TEMP DockingInitialPerturbationCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DockingInitialPerturbation::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DockingInitialPerturbationCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new DockingInitialPerturbation );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DockingInitialPerturbation::mover_name()
// XRW TEMP {
// XRW TEMP  return "DockingInitialPerturbation";
// XRW TEMP }


// initial perturbation on one of the partners
// the type of perturbation is defined in the options
// some of the options are randomize1 (randomizes the first partner)
// randomize2 (randomizes the second partner), dock_pert, and spin
//------------------------------------------------------------------------------
//
//     there are several ways to perturb the structure before beginning
//     the search; they are controlled through command-line flags
//
//     at the end, partners are slid into contact and scored


//constructors
DockingInitialPerturbation::DockingInitialPerturbation()
:
	Mover(),
	slide_( true ),
	rigid_body_info_( /* NULL */ )
{
	Mover::type( "DockingInitialPerturbation" );
	movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	multiple_jumps_ = false;
	init();
}

DockingInitialPerturbation::DockingInitialPerturbation(
	core::Size const rb_jump,
	bool const slide
) :
	Mover(),
	slide_(slide),
	rigid_body_info_( /* NULL */ )
{
	Mover::type( "DockingInitialPerturbation" );
	movable_jumps_ = utility::tools::make_vector1<core::Size>(rb_jump);
	multiple_jumps_ = false;
	init();
}

DockingInitialPerturbation::DockingInitialPerturbation(
	DockJumps const movable_jumps,
	bool const slide
) :
	Mover(),
	slide_( slide ),
	movable_jumps_( movable_jumps ),
	rigid_body_info_( /* NULL */ )
{
	Mover::type( "DockingInitialPerturbation" );
	if ( movable_jumps_.size() > 1 ) multiple_jumps_ = true;
	else multiple_jumps_ = false;
	init();
}

void
DockingInitialPerturbation::init()
{
	set_default();
	init_from_options();  // put this into apply in case scripts is used, then this will not be needed.
}

void
DockingInitialPerturbation::set_default()
{
	randomize1_ = false;
	randomize2_ = false;
	use_ellipsoidal_randomization_ = false;
	if_dock_pert_ = false;
	if_uniform_trans_ = false;
	spin_ = false;
	center_at_interface_ = false;
	// dock_pert_ = new utility::vector1< Real >(NULL);
	// uniform_trans_ = NULL;
	slide_axis_.zero();
	spin_center_.zero();
}

protocols::moves::MoverOP
DockingInitialPerturbation::clone() const {
	return protocols::moves::MoverOP( new DockingInitialPerturbation( *this ) );
}

//destructor
DockingInitialPerturbation::~DockingInitialPerturbation() = default;

// ALERT!
// register_options() and init_from_options() are not called anywhere yet!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!11

void
DockingInitialPerturbation::init_from_options()
{
	using namespace basic::options;
	TR << "Reading options..." << std::endl;

	if ( option[ OptionKeys::docking::randomize1 ].user() ) {
		set_randomize1(option[ OptionKeys::docking::randomize1 ]());
	}

	if ( option[ OptionKeys::docking::randomize2 ].user() ) {
		set_randomize2(option[ OptionKeys::docking::randomize2 ]());
	}

	if ( option[ OptionKeys::docking::use_ellipsoidal_randomization ].user() ) {
		set_use_ellipsoidal_randomization(option[ OptionKeys::docking::use_ellipsoidal_randomization ]());
	}

	if ( option[ OptionKeys::docking::dock_pert ].user() ) {
		set_dock_pert(option[ OptionKeys::docking::dock_pert ]());
	}

	if ( option[ OptionKeys::docking::uniform_trans ].user() ) {
		set_uniform_trans(option[ OptionKeys::docking::uniform_trans ]());
	}

	if ( option[ OptionKeys::docking::spin ].user() ) {
		set_spin(option[ OptionKeys::docking::spin ]());
	}

	if ( option[ OptionKeys::docking::tilt ].user() ) {
		set_tilt(option[ OptionKeys::docking::tilt ]());
		set_tilt1_center( option[ OptionKeys::docking::tilt1_center ]());
		set_tilt2_center( option[ OptionKeys::docking::tilt2_center ]());
	}



	if ( option[ OptionKeys::docking::center_at_interface ].user() ) {
		set_center(option[ OptionKeys::docking::center_at_interface ]());
	}
}

void
DockingInitialPerturbation::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::randomize1 );
	option.add_relevant( OptionKeys::docking::randomize2 );
	option.add_relevant( OptionKeys::docking::use_ellipsoidal_randomization );
	option.add_relevant( OptionKeys::docking::dock_pert );
	option.add_relevant( OptionKeys::docking::uniform_trans );
	option.add_relevant( OptionKeys::docking::spin );
	option.add_relevant( OptionKeys::docking::tilt );
	option.add_relevant( OptionKeys::docking::tilt1_center );
	option.add_relevant( OptionKeys::docking::tilt2_center );
	option.add_relevant( OptionKeys::docking::center_at_interface );
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief   Make starting perturbations for rigid body moves
///
/// @details    There are several ways to perturb the structure before beginning
///     the search; they are controlled through command-line flags
///     At the end, partners are slid into contact and scored (including
///     mc_reset).
///     Also, they are tested for passing the FAB filter.
///
///
/// @references see dock_structure or pose_docking from rosetta++
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
void DockingInitialPerturbation::apply( core::pose::Pose & pose )
{
	if ( rigid_body_info_ ) {
		movable_jumps_ = rigid_body_info_->movable_jumps();
		TR.Debug << "finished reading movable_jumps_ from RigidBodyInfo" << std::endl;
		if ( movable_jumps_.empty() ) {
			utility_exit_with_message( "DockSetupMover has to be applied before DockingInitialPerturbation !" );
		}
	}


	runtime_assert( !movable_jumps_.empty() );

	for ( DockJumps::const_iterator it=movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
		apply_body( pose, *it );
		TR.Debug <<"movable_jumps_ value in apply:" << *it << std::endl;
	}
}

void
DockingInitialPerturbation::apply_body(core::pose::Pose & pose, core::Size jump_number )
{
	using namespace moves;

	if ( randomize1_ ) {
		TR << "randomize1: true" << std::endl;
		if ( use_ellipsoidal_randomization_ ) {
			TR << "use_ellipsoidal_randomization: true" << std::endl;
			EllipsoidalRandomizationMover mover( jump_number, true );
			mover.apply( pose );
			slide_axis_ = mover.get_slide_axis();
			spin_center_ = mover.get_spin_center();
			// Partners need to start in contact for ellipsoidal randomization to work
			if ( randomize2_ && !pose.is_fullatom() ) {
				DockingSlideIntoContact slide( jump_number, slide_axis_ );
				slide.apply( pose );
				TR.Debug << "centroid mode, DockingSlideIntoContact applied" << std::endl;
			} else if ( randomize2_ ) {
				FaDockingSlideIntoContact slide( jump_number, slide_axis_ );
				slide.apply( pose );
				TR.Debug << "fa-standard mode, FaDockingSlideIntoContact applied" << std::endl;
			}
		} else {
			rigid::RigidBodyRandomizeMover mover( pose, jump_number, rigid::partner_upstream );
			mover.apply( pose );
		}
	}

	if ( randomize2_ ) {
		TR << "randomize2: true" << std::endl;
		if ( use_ellipsoidal_randomization_ ) {
			TR << "use_ellipsoidal_randomization: true" << std::endl;
			EllipsoidalRandomizationMover mover( jump_number, false );
			mover.apply( pose );
			slide_axis_ = mover.get_slide_axis();
			spin_center_ = mover.get_spin_center();
		} else {
			rigid::RigidBodyRandomizeMover mover( pose, jump_number, rigid::partner_downstream );
			mover.apply( pose );
		}
	}
	if ( if_dock_pert_ ) {
		// DO NOT supply default values for this option -- reasonable values differ for protein and ligand protocols.
		// Also, default values will cause a perturbation to *always* occur, even with no command line flags -- very surprising.
		// Adding defaults WILL BREAK existing protocols in unexpected ways.
		// Decided by Jeff, Monica, Ian, and Sid in March 2008.
		//
		// Jeff notes that eventually there should be 3 parameters, like Rosetta++:
		// rotation, normal translation, and perpendicular translation.
		TR << "dock_pert: true" << std::endl;
		/// read in dock_pert options from commandline.  the first value is the
		/// rotation magnitude and the second value is the translational value
		utility::vector1< Real > pert_mags = dock_pert_;
		//TR << "option[ docking::rotational ]()" << rot << "\n";
		//TR << "option[ docking::parallel ]()" << trans << "\n";
		TR << "option[ docking::dock_pert ]()" << pert_mags[rot] << ' ' << pert_mags[trans] << std::endl;
		rigid::RigidBodyPerturbMoverOP mover;
		if ( center_at_interface_ ) mover = rigid::RigidBodyPerturbMoverOP( new rigid::RigidBodyPerturbMover( jump_number, pert_mags[rot], pert_mags[trans], rigid::partner_downstream, true ) );
		else mover = rigid::RigidBodyPerturbMoverOP( new rigid::RigidBodyPerturbMover( jump_number, pert_mags[rot], pert_mags[trans] ) );
		mover->apply( pose );
	}
	if ( if_uniform_trans_ ) {
		core::Real mag( uniform_trans_ );
		TR << "uniform_trans: " << mag << std::endl;
		rigid::UniformSphereTransMover mover( jump_number, mag );
		mover.apply( pose );
	}
	if ( spin_ ) {
		TR << "axis_spin: true" << std::endl;
		rigid::RigidBodySpinMover mover( jump_number );
		if ( slide_axis_.length() != 0 ) {
			mover.spin_axis( slide_axis_ );
			mover.rot_center( spin_center_ );
		}
		mover.apply( pose );
	}
	// Tilt partners around sliding axis up to a user specified angle
	if ( tilt_.size() ) {
		if ( tilt_.size() != 2 ) {
			utility_exit_with_message( "Tilt option of DockingInitialPerturbation must be given exactly two parameters." );
		}
		TR << "axis_tilt: true" << std::endl;
		// process user specified rotation centers
		core::Size tilt1_center_resid_=0;
		core::Size tilt2_center_resid_=0;
		if ( tilt1_center_ != "" ) tilt1_center_resid_ = core::pose::parse_resnum( tilt1_center_, pose );
		if ( tilt2_center_ != "" ) tilt2_center_resid_ = core::pose::parse_resnum( tilt2_center_, pose );

		// instantiate mover and apply tilt
		rigid::RigidBodyTiltMover mover(jump_number,tilt_[1],tilt_[2],tilt1_center_resid_,tilt2_center_resid_);
		if ( slide_axis_.length() != 0 ) {
			mover.spin_axis( slide_axis_ );
		}
		mover.apply( pose );
	}
	// DO NOT do this for e.g. ligand docking
	if ( slide_ && !pose.is_fullatom() ) {
		TR << "sliding into contact for centroid mode" << std::endl;
		DockingSlideIntoContact slide( jump_number, slide_axis_ );
		slide.apply( pose );
		TR.Debug << "centroid mode, DockingSlideIntoContact applied" << std::endl;
	} else if ( slide_ ) {
		TR << "sliding into contact for full-atom mode" << std::endl;
		FaDockingSlideIntoContact slide( jump_number, slide_axis_ );
		slide.apply( pose );
		TR.Debug << "fa-standard mode, FaDockingSlideIntoContact applied" << std::endl;
	}


}
// XRW TEMP std::string
// XRW TEMP DockingInitialPerturbation::get_name() const {
// XRW TEMP  return "DockingInitialPerturbation";
// XRW TEMP }

void
DockingInitialPerturbation::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	if ( !data_map.has( "RigidBodyInfo", "docking_setup" ) ) {
		TR << "RigidBodyInfo not found in basic::datacache::DataMap" << std::endl;
		rigid_body_info_ = protocols::docking::RigidBodyInfoOP( new protocols::docking::RigidBodyInfo );
		data_map.add( "RigidBodyInfo", "docking_setup", rigid_body_info_ );
		//  throw utility::excn::EXCN_RosettaScriptsOption( "RigidBodyInfo not found in basic::datacache::DataMap, DockingInitialPerturbation can not be done, so exit here!" );
	} else {
		rigid_body_info_ = data_map.get_ptr< protocols::docking::RigidBodyInfo >( "RigidBodyInfo", "docking_setup" );
		TR.Debug << "get RigidBodyInfo pointer from basic::datacache::DataMap" << std::endl;
	}

	if ( tag->hasOption( "randomize1" ) ) {
		set_randomize1( tag->getOption<bool>( "randomize1" ) );
	}

	if ( tag->hasOption( "randomize2" ) ) {
		set_randomize2( tag->getOption<bool>( "randomize2" ) );
	}

	if ( tag->hasOption( "use_ellipsoidal_randomization" ) ) {
		set_use_ellipsoidal_randomization( tag->getOption<bool>( "use_ellipsoidal_randomization" ) );
	}

	if ( tag->hasOption( "dock_pert" ) && tag->getOption<bool>( "dock_pert" ) ) {
		dock_pert_.push_back( tag->getOption<core::Real>("trans" ) );
		dock_pert_.push_back( tag->getOption<core::Real>( "rot" ) );
		set_dock_pert( dock_pert_ );
	}

	if ( tag->hasOption( "uniform_trans" ) ) {
		set_uniform_trans( tag->getOption<core::Real>( "uniform_trans" ) );
	}

	if ( tag->hasOption( "spin" ) ) {
		set_spin( tag->getOption<bool>( "spin" ) );
	}

	if ( tag->hasOption( "center_at_interface" ) ) {
		set_center( tag->getOption<bool>( "center_at_interface" ) );
	}

	slide_ = tag->getOption<bool>( "slide", true );
}

std::string DockingInitialPerturbation::get_name() const {
	return mover_name();
}

std::string DockingInitialPerturbation::mover_name() {
	return "DockingInitialPerturbation";
}

void DockingInitialPerturbation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "randomize1", xsct_rosetta_bool, "Randomize the first docking partner" )
		+ XMLSchemaAttribute( "randomize2", xsct_rosetta_bool, "Randomize the second docking partner" )
		+ XMLSchemaAttribute( "use_ellipsoidal_randomization", xsct_rosetta_bool, "Use the EllipsoidalRandomizationMover instead of the RigidBodyRandomizeMover" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_pert", xsct_rosetta_bool, "Read in translational and rotational perturbations and apply to pose", "false" )
		+ XMLSchemaAttribute( "trans", xsct_real, "Translational perturbation to apply before docking" ) //required if dock_pert is true
		+ XMLSchemaAttribute( "rot", xsct_real, "Rotational perturbation to apply before docking" ) //required if dock_pert is true
		+ XMLSchemaAttribute( "uniform_trans", xsct_real, "Use the UniformSphereTransMover" )
		+ XMLSchemaAttribute( "spin", xsct_rosetta_bool, "Spin partner about its axis" )
		+ XMLSchemaAttribute( "center_at_interface", xsct_rosetta_bool, "Center the spin at the interface" )
		+ XMLSchemaAttribute( "slide", xsct_rosetta_bool, "Slide docking partners into contact" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform an initial perturbation to docking partners before docking", attlist );
}

std::string DockingInitialPerturbationCreator::keyname() const {
	return DockingInitialPerturbation::mover_name();
}

protocols::moves::MoverOP
DockingInitialPerturbationCreator::create_mover() const {
	return protocols::moves::MoverOP( new DockingInitialPerturbation );
}

void DockingInitialPerturbationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DockingInitialPerturbation::provide_xml_schema( xsd );
}



////////////////////////////////////////// DockingSlideIntoContact ////////////////////////////////

// default constructor
DockingSlideIntoContact::DockingSlideIntoContact() : Mover()
{
	using namespace core::scoring;
	Mover::type( "DockingSlideIntoContact" );
	rb_jump_ = 1;
	//slide_axis_(0);
	scorefxn_ = ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//constructor
DockingSlideIntoContact::DockingSlideIntoContact(
	core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump), slide_axis_(0.0)
{
	using namespace core::scoring;
	Mover::type( "DockingSlideIntoContact" );
	scorefxn_ = ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

DockingSlideIntoContact::DockingSlideIntoContact(
	core::Size const rb_jump,
	core::Vector const & slide_axis
): Mover(), rb_jump_(rb_jump), slide_axis_(slide_axis)
{
	using namespace core::scoring;
	Mover::type( "DockingSlideIntoContact" );
	scorefxn_ = ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//destructor
DockingSlideIntoContact::~DockingSlideIntoContact() = default;


void DockingSlideIntoContact::apply( core::pose::Pose & pose )
{
	using namespace moves;
	using namespace protocols::membrane;
	using namespace core::scoring;

	bool vary_stepsize( false );

	// for a membrane pose the translation axis should be in the membrane
	if ( pose.conformation().is_membrane() ) {

		TR << "Adjusting slide axis for membrane proteins" << std::endl;
		TR << "     slide axis before: " << slide_axis_.to_string() << std::endl;

		// get membrane axis from docking metrics function
		slide_axis_ = membrane_axis( pose, rb_jump_ );

		// get direction of axis by making a trial move and seeing whether the partners
		// move together or apart from each other
		// this is stupid, but I don't know how else to fix this

		// create EmbeddinDef objects
		protocols::membrane::geometry::EmbeddingDef emb_up, emb_down;
		update_partner_embeddings( pose, rb_jump_, emb_up, emb_down );
		TR << "center1: " << emb_up.center().to_string() << std::endl;
		TR << "center2: " << emb_down.center().to_string() << std::endl;

		// get distance between points
		core::Real dist1 = ( emb_down.center() - emb_up.center() ).length();

		// trial move
		rigid::RigidBodyTransMoverOP trial( new rigid::RigidBodyTransMover( slide_axis_, rb_jump_, vary_stepsize ) );
		trial->apply( pose );

		// get new distance between points
		update_partner_embeddings( pose, rb_jump_, emb_up, emb_down );
		core::Real dist2 = ( emb_down.center() - emb_up.center() ).length();
		TR << "dist1: " << dist1 << ", dist2: " << dist2 << std::endl;
		TR << "center1: " << emb_up.center().to_string() << std::endl;
		TR << "center2: " << emb_down.center().to_string() << std::endl;

		// negate or not, this is super counter intuitive!!!
		if ( dist2 < dist1 ) {
			TR << "slide axis negated..." << std::endl;
			slide_axis_.negate();
		}

		slide_axis_.normalize();
		TR << "     slide axis after: " << slide_axis_.to_string() << std::endl;

		// set variable stepsize for membrane proteins
		vary_stepsize = true;
	}

	rigid::RigidBodyTransMoverOP mover;

	if ( slide_axis_.length() != 0 ) {
		mover = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover( slide_axis_, rb_jump_, vary_stepsize ) );
	} else {
		mover = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover( pose, rb_jump_, vary_stepsize ) );
	}
	( *scorefxn_ )( pose );

	TR << "sliding into contact" << std::endl;
	TR << "Moving away" << std::endl;
	core::Size const counter_breakpoint( 500 );
	core::Size counter( 0 );

	// first try moving away from each other
	while ( pose.energies().total_energies()[ interchain_vdw ] > 0.1 && counter <= counter_breakpoint ) {
		mover->apply( pose );
		( *scorefxn_ )( pose );
		++counter;
	}
	if ( counter > counter_breakpoint ) {
		TR<<"failed moving away with original vector. Aborting DockingSlideIntoContact."<<std::endl;
		set_current_tag( "fail" );
		return;
	}
	counter = 0;

	// then try moving towards each other
	TR << "Moving together" << std::endl;
	mover->trans_axis().negate();
	while ( counter <= counter_breakpoint && pose.energies().total_energies()[ interchain_vdw ] < 0.1 ) {

		//  TR << "interchain_vdw: " << pose.energies().total_energies()[ interchain_vdw ] << std::endl;
		mover->apply( pose );
		( *scorefxn_ )( pose );
		++counter;
	}
	if ( counter > counter_breakpoint ) {
		TR<<"moving together failed. Aborting DockingSlideIntoContact."<<std::endl;
		set_current_tag( "fail" );
		return;
	}

	// if the stepsize was variable, do a few more tries with stepsize of 1, before
	// moving it back out
	counter = 0;
	if ( mover->step_size() > 1.0 ) {
		TR << "step size of 1.0..." << std::endl;
		TR << "interchain scores: " << pose.energies().total_energies()[ interchain_vdw ] << std::endl;
		while ( counter <= 10 && pose.energies().total_energies()[ interchain_vdw ] < 0.1 ) {

			TR << "moving partners together" << std::endl;
			TR << "interchain scores: " << pose.energies().total_energies()[ interchain_vdw ] << std::endl;
			mover->vary_stepsize( false );
			mover->step_size( 1.0 );
			mover->apply( pose );
			( *scorefxn_ )( pose );
			++counter;
		}
	}
	// move away again until just touching
	mover->trans_axis().negate();
	mover->apply( pose );
}

std::string
DockingSlideIntoContact::get_name() const {
	return "DockingSlideIntoContact";
}

void
DockingSlideIntoContact::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Jump number: " << get_jump_num() << std::endl;
}

std::ostream &operator<< ( std::ostream & os, DockingSlideIntoContact const & mover )
{
	mover.show(os);
	return os;
}

////////////////////////////////////////// FaDockingSlideIntoContact ////////////////////////////////

// default constructor
FaDockingSlideIntoContact::FaDockingSlideIntoContact()
{
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
	//slide_axis_(0.0);
}


//constructor
FaDockingSlideIntoContact::FaDockingSlideIntoContact(
	core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump), tolerance_(0.2), slide_axis_(0.0)
{
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

FaDockingSlideIntoContact::FaDockingSlideIntoContact( utility::vector1<core::Size> rb_jumps):
	Mover(), rb_jumps_(rb_jumps),tolerance_(0.2),slide_axis_(0.0){
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

FaDockingSlideIntoContact::FaDockingSlideIntoContact( core::Size const rb_jump, core::Vector const & slide_axis): Mover(), rb_jump_(rb_jump), tolerance_(0.2), slide_axis_(slide_axis)
{
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

//destructor
FaDockingSlideIntoContact::~FaDockingSlideIntoContact() = default;

void FaDockingSlideIntoContact::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	(*scorefxn_)( pose );
	core::Real const initial_fa_rep = pose.energies().total_energies()[ fa_rep ];
	bool are_touching = false;

	utility::vector1<rigid::RigidBodyTransMover> trans_movers;

	if ( slide_axis_.length() != 0 ) {
		trans_movers.push_back( rigid::RigidBodyTransMover( slide_axis_, rb_jump_ ));
	} else if ( rb_jumps_.size()<1 ) {
		trans_movers.push_back( rigid::RigidBodyTransMover( pose,rb_jump_ ));
	} else {
		for ( core::Size & rb_jump : rb_jumps_ ) {
			trans_movers.push_back( rigid::RigidBodyTransMover(pose, rb_jump ));
		}
	}

	utility::vector1< rigid::RigidBodyTransMover >::iterator const end(trans_movers.end());

	//int i=1;
	// Take 2A steps till clash, then back apart one step.  Now you're within 2A of touching.
	// Repeat with 1A steps, 0.5A steps, 0.25A steps, etc until you're as close are you want.
	for ( core::Real stepsize = 2.0; stepsize > tolerance_; stepsize /= 2.0 ) {
		for (
				auto trans_mover(trans_movers.begin());
				trans_mover != end;
				++trans_mover
				) {
			trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move together
			trans_mover->step_size(stepsize);
		}
		core::Size const counter_breakpoint( 500 );
		core::Size counter( 0 );
		do
		{
			for (
					auto trans_mover(trans_movers.begin());
					trans_mover != end;
					++trans_mover
					) {
				trans_mover->apply( pose );
			}
			(*scorefxn_)( pose );
			core::Real const push_together_fa_rep = pose.energies().total_energies()[ fa_rep ];
			//std::cout << "fa_rep = " << push_together_fa_rep << std::endl;
			are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
			//std::ostringstream s;
			//s << "snapshot" << i << ".pdb";
			//pose.dump_pdb(s.str());
			//i += 1;
			++counter;
		} while( counter <= counter_breakpoint && !are_touching );
		if ( counter > counter_breakpoint ) {
			TR<<"Failed Fadocking Slide Together. Aborting."<<std::endl;
			set_current_tag( "fail" );
		}
		for (
				auto trans_mover(trans_movers.begin());
				trans_mover != end;
				++trans_mover
				) {
			trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move apart
			trans_mover->apply( pose );
		}
	}
}

std::string
FaDockingSlideIntoContact::get_name() const {
	return "FaDockingSlideTogether";
}

void
FaDockingSlideIntoContact::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Jump number: " << get_jump_num() << "\nTolerance:   " << get_tolerance() << std::endl;
}

std::ostream &operator<< ( std::ostream &os, FaDockingSlideIntoContact const &fadock )
{
	fadock.show(os);
	return os;
}

////////////////////////////////////////////////////////////////////////////

/// @brief More general function for low-res or highres DockingSlideIntoContact
/// @details Both DockingSlideIntoContact and FaDockingSlideIntoContact have their
///  issues and are not as general as they could be; they should be a single
///  class, not two; since they are called quite often, I am rewriting the
///  class for now and will replace it's use later on
/// @details Contrary to the name, slides things apart first, then together.
/// OK for proteins, bad for ligands (because they may escape the pocket permanently).

/// @brief default constructor, default jump = 1, default scorefunction is lowres
SlideIntoContact::SlideIntoContact() : Mover(),
	jump_( 1 ),
	slide_axis_( 0 ),
	vary_stepsize_( false ),
	stepsize_( 1.0 ),
	move_apart_first_ ( true ),
	scorefxn_( nullptr ),
	scoretype_( ),
	threshold_( 1.0 ),
	starting_rep_( 0 )
{}

/// @brief constructor with arguments
SlideIntoContact::SlideIntoContact( core::Size const jump ) : Mover(),
	jump_( jump ),
	slide_axis_( 0 ),
	vary_stepsize_( false ),
	stepsize_( 1.0 ),
	move_apart_first_ ( true ),
	scorefxn_( nullptr ),
	scoretype_( ),
	threshold_( 1.0 ),
	starting_rep_( 0 )
{}

/// @brief destructor
SlideIntoContact::~SlideIntoContact() = default;

/// @brief apply
void SlideIntoContact::apply( core::pose::Pose & pose ) {

	using namespace moves;
	using namespace core::scoring;
	using namespace protocols::membrane;

	register_options();
	init_from_cmd();

	// use centroid scorefunction if not previously defined
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
		scoretype_ = score_type_from_name( "interchain_contact" );
	}

	rigid::RigidBodyTransMoverOP mover;

	// if slide axis known
	if ( slide_axis_.length() != 0 ) {
		mover = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover( slide_axis_, jump_, vary_stepsize_ ) );
	} else {
		// if slide axis unknown, get it from partner COMs
		mover = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover( pose, jump_, vary_stepsize_ ) );
	}

	// score the pose to update the coordinates and get a starting score for the repulsive score
	// ( *scorefxn_ )( pose );
	// core::Real starting_rep = pose.energies().total_energies()[ scoretype_ ];
	TR << "starting repulsion: " << starting_rep_ << std::endl;

	core::Size const counter_breakpoint( 200 );
	core::Size counter( 0 );

	// if first moving apart
	if ( move_apart_first_ == true ) {

		TR << "Moving away" << std::endl;

		// negating axis for moving away
		mover->trans_axis().negate();
		mover->step_size( 500 );
		mover->apply( pose );

		//  // first try moving away from each other
		//  while ( pose.energies().total_energies()[ scoretype_ ] - starting_rep_ > 10 && counter <= counter_breakpoint ) {
		//   mover->apply( pose );
		//   ( *scorefxn_ )( pose );
		//
		//   pose.energies().show_total_headers( TR );
		//   TR << std::endl;
		//   pose.energies().show_totals( TR );
		//   TR << std::endl;
		//
		//   ++counter;
		//  }
		//  if ( counter > counter_breakpoint ) {
		//   TR<<"failed moving away with original vector. Aborting SlideIntoContact."<<std::endl;
		//   set_current_tag( "fail" );
		//   return;
		//  }
		//  counter = 0;

		// then try moving towards each other, negating axis again to move together
		mover->trans_axis().negate();
	}

	TR << "Moving together" << std::endl;
	counter = 0;
	( *scorefxn_ )( pose );
	core::Real diff = pose.energies().total_energies()[ scoretype_ ] - starting_rep_;
	while ( counter <= counter_breakpoint && diff <= threshold_ ) {

		TR << name_from_score_type( scoretype_ ) << ": " << pose.energies().total_energies()[ scoretype_ ] << std::endl;

		mover->apply( pose );
		( *scorefxn_ )( pose );
		diff = pose.energies().total_energies()[ scoretype_ ] - starting_rep_;

		TR << "counter" << counter << " break " << counter_breakpoint << " diff: " << diff << ", threshold: " << threshold_ << std::endl;
		++counter;
	}
	if ( counter > counter_breakpoint ) {
		TR<<"moving together failed. Aborting SlideIntoContact."<<std::endl;
		set_current_tag( "fail" );
		return;
	}

	// if the stepsize was variable, do a few more tries with stepsize of 1, before
	// moving it back out
	counter = 0;
	if ( mover->step_size() > 1.0 ) {
		TR << "step size of 1.0..." << std::endl;
		mover->vary_stepsize( false );
		mover->step_size( 1.0 );

		( *scorefxn_ )( pose );
		diff = pose.energies().total_energies()[ scoretype_ ] - starting_rep_;

		TR << name_from_score_type( scoretype_ ) << ": " << pose.energies().total_energies()[ scoretype_ ] << std::endl;
		TR << "starting rep: " << starting_rep_ << std::endl;
		TR << "diff: " << diff << ", threshold: " << threshold_ << std::endl;

		while ( counter <= 50 && diff <= threshold_ ) {

			TR << "moving partners together" << std::endl;
			TR << name_from_score_type( scoretype_ ) << ": " << pose.energies().total_energies()[ scoretype_ ] << std::endl;

			mover->apply( pose );

			( *scorefxn_ )( pose );
			diff = pose.energies().total_energies()[ scoretype_ ] - starting_rep_;
			TR << "diff: " << diff << ", threshold: " << threshold_ << std::endl;
			++counter;
		}
	}
	// move away again until just touching
	// mover->trans_axis().negate();
	// mover->apply( pose );

} // apply

/// @brief show
void SlideIntoContact::show(std::ostream & output) const {
	output << "Jump number: " << get_jump_num() << std::endl;
}

/// @brief set slide axis in Mover
void SlideIntoContact::slide_axis( core::Vector const & slide_axis ) {
	slide_axis_ = slide_axis;
}

/// @brief vary the stepsize depending on the distance between the partners
void SlideIntoContact::vary_stepsize( bool yesno ) {
	vary_stepsize_ = yesno;
}

/// @brief set the stepsize to walk between the partners
void SlideIntoContact::stepsize( core::Real stepsize ) {
	stepsize_ = stepsize;
}

/// @brief move partners apart first
void SlideIntoContact::move_apart_first( bool yesno ) {
	move_apart_first_ = yesno;
}

/// @brief set starting repulsion
void SlideIntoContact::set_starting_rep( core::Real starting_rep ) {
	starting_rep_ = starting_rep;
}

/// @brief scorefunction and scoreterm used for evaluating closeness of partners
void SlideIntoContact::scorefunction( std::string sfxn_name, std::string scoreterm_for_sliding ) {
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( sfxn_name );
	scoretype_ = core::scoring::score_type_from_name( scoreterm_for_sliding );
}

/// @brief Get the name of this mover
std::string SlideIntoContact::get_name() const {
	return "SlideIntoContact";
}

/// @brief get the jump number between the partners
core::Size SlideIntoContact::get_jump_num() const {
	return jump_;
}

/// @brief get the stepsize to walk between the partners
core::Real SlideIntoContact::get_stepsize() const {
	return stepsize_;
}

/// @brief get the scorefunction name that is used to evaluate closeness btw partners
std::string SlideIntoContact::get_sfxn_name() const {
	return scorefxn_->get_name();
}

////////////////////////////////////////////////////////////////////////////////

/// @brief Register Options with JD2
void SlideIntoContact::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::dock::slide_threshold );

} // register options

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize Mover options from the comandline
void SlideIntoContact::init_from_cmd() {

	using namespace basic::options;

	// docking partners
	if ( option[ OptionKeys::mp::dock::slide_threshold ].user() ) {
		threshold_ = option[ OptionKeys::mp::dock::slide_threshold ]();
	}

} // init from cmd

void move_apart( core::pose::Pose & pose, int jump, core::Vector const & axis ) {

	using namespace protocols::rigid;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	TR << "move apart" << std::endl;

	RigidBodyTransMoverOP mover( new RigidBodyTransMover( axis, jump ) );
	mover->step_size( 500 );
	mover->apply( pose );

	// compute distance between partners
	EmbeddingDef emb_up, emb_down;
	update_partner_embeddings( pose, jump, emb_up, emb_down );
	core::Real dist1 = ( emb_down.center() - emb_up.center() ).length();
	TR << "distance between partners: " << dist1 << std::endl;

}

void move_together( core::pose::Pose & pose, int jump, core::scoring::ScoreFunctionOP sfxn ) {

	using namespace basic::options;
	using namespace core::pose;
	using namespace protocols::rigid;
	using namespace core::scoring::methods;
	using namespace protocols::scoring;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	// docking partners
	core::Real max_contacts( 10 );
	if ( option[ OptionKeys::mp::dock::slide_threshold ].user() ) {
		max_contacts = option[ OptionKeys::mp::dock::slide_threshold ]();
	}

	// split pose by jump
	Pose pose1, pose2;
	partition_pose_by_jump( pose, jump, pose1, pose2 );

	// compute radii of upstream and downstream partners
	RG_Energy_Fast rg = RG_Energy_Fast();
	core::Real rgyr1 = rg.calculate_rg_score( pose1 );
	core::Real rgyr2 = rg.calculate_rg_score( pose2 );

	// compute axis and distance between partners
	EmbeddingDef emb_up, emb_down;
	update_partner_embeddings( pose, jump, emb_up, emb_down );
	core::Vector axis = emb_up.center() - emb_down.center();
	core::Real dist1 = axis.length();
	TR << "distance between partners before big step: " << dist1 << std::endl;

	// initialize mover
	RigidBodyTransMoverOP mover( new RigidBodyTransMover( axis, jump ) );

	// compute first step size
	core::Real step1 = axis.length() - rgyr1 - rgyr2 - 10;
	TR << "axis length " << axis.length() << std::endl;
	TR << "rgyr1 " << rgyr1 << std::endl;
	TR << "rgyr2 " << rgyr2 << std::endl;
	TR << "step1: " << step1 << std::endl;
	TR << "jump " << jump << std::endl;
	mover->step_size( step1 );
	mover->apply( pose );

	// compute distance between partners
	update_partner_embeddings( pose, jump, emb_up, emb_down );
	dist1 = ( emb_down.center() - emb_up.center() ).length();
	TR << "distance between partners after big step: " << dist1 << std::endl;

	// small steps to move closer, compute nres in interface
	mover->step_size( 1.0 );
	( *sfxn )( pose );
	Interface interface = Interface( jump );
	interface.calculate( pose );

	core::Size counter(0);
	core::Size cnt_stop(200);
	while ( interface.interface_nres() <= max_contacts && counter <= cnt_stop ) {

		mover->apply( pose );
		( *sfxn )( pose );
		interface.calculate( pose );

		update_partner_embeddings( pose, jump, emb_up, emb_down );
		dist1 = ( emb_down.center() - emb_up.center() ).length();
		TR << "small step " << counter << " distance between partners: " << dist1 << std::endl;

		TR << "interface nres " << interface.interface_nres() << std::endl;
		counter++;
	}
}


} // namespace docking
} // namespace protocols
