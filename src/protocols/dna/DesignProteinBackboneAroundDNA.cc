// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DesignProteinBackboneAroundDNA.cc
/// @brief performs a round flexible backbone sampling/design.  (Interim(?) encapsulation of some loose code in a mover class.)
/// @author ashworth

#include <protocols/dna/DesignProteinBackboneAroundDNA.hh>
#include <protocols/dna/DesignProteinBackboneAroundDNACreator.hh>

#include <core/chemical/ResidueType.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/options/option.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/Tracer.hh>

#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/make_loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
using utility::vector1;

#include <sstream>
#include <set>
#include <iostream>

// option key includes

#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/keys/Key3Vector.hh>
#include <iterator>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace dna {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace pose;
using namespace scoring;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.DesignProteinBackboneAroundDNA", t_info );

// XRW TEMP std::string
// XRW TEMP DesignProteinBackboneAroundDNACreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DesignProteinBackboneAroundDNA::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DesignProteinBackboneAroundDNACreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DesignProteinBackboneAroundDNA );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DesignProteinBackboneAroundDNA::mover_name()
// XRW TEMP {
// XRW TEMP  return "DesignProteinBackboneAroundDNA";
// XRW TEMP }

DesignProteinBackboneAroundDNA::DesignProteinBackboneAroundDNA() :
	protocols::simple_moves::PackRotamersMover( DesignProteinBackboneAroundDNA::mover_name() ),
	type_("ccd"),
	gapspan_( option[ OptionKeys::loops::gapspan ]() ),
	spread_( option[ OptionKeys::loops::spread ]() ),
	cycles_outer_( option[ OptionKeys::run::cycles_outer ]() ),
	cycles_inner_( option[ OptionKeys::run::cycles_inner ]() ),
	repack_rate_( option[ OptionKeys::run::repack_rate ]() ),
	temp_initial_( option[ OptionKeys::MonteCarlo::temp_initial ]() ),
	temp_final_( option[ OptionKeys::MonteCarlo::temp_final ]() ),
	designable_second_shell_( option[ OptionKeys::dna::design::designable_second_shell ]() )
{}

DesignProteinBackboneAroundDNA::DesignProteinBackboneAroundDNA(
	std::string  type,
	ScoreFunctionCOP scorefxn
) :
	protocols::simple_moves::PackRotamersMover( DesignProteinBackboneAroundDNA::mover_name() ),
	type_(std::move(type)),
	gapspan_( option[ OptionKeys::loops::gapspan ]() ),
	spread_( option[ OptionKeys::loops::spread ]() ),
	cycles_outer_( option[ OptionKeys::run::cycles_outer ]() ),
	cycles_inner_( option[ OptionKeys::run::cycles_inner ]() ),
	repack_rate_( option[ OptionKeys::run::repack_rate ]() ),
	temp_initial_( option[ OptionKeys::MonteCarlo::temp_initial ]() ),
	temp_final_( option[ OptionKeys::MonteCarlo::temp_final ]() ),
	designable_second_shell_( option[ OptionKeys::dna::design::designable_second_shell ]() )
{
	score_function( scorefxn );
}

DesignProteinBackboneAroundDNA::~DesignProteinBackboneAroundDNA()= default;

void
DesignProteinBackboneAroundDNA::targeted_dna( DnaDesignDefOPs const & defs ) {
	targeted_dna_ = defs;
}

DnaDesignDefOPs const &
DesignProteinBackboneAroundDNA::targeted_dna() const { return targeted_dna_; }

void
DesignProteinBackboneAroundDNA::apply( Pose & pose )
{
	if ( ! task_factory() ) {
		TaskFactoryOP new_tf( new TaskFactory );
		new_tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
		if ( option[ OptionKeys::packing::resfile ].user() ) new_tf->push_back( TaskOperationCOP( new ReadResfile ) );
		RestrictDesignToProteinDNAInterfaceOP rest( new RestrictDesignToProteinDNAInterface );
		rest->set_base_only( option[ OptionKeys::dna::design::base_contacts_only ]() );
		rest->copy_targeted_dna( targeted_dna_ );
		new_tf->push_back( rest );
		task_factory( new_tf );
	}

	// make a list of moveable positions as indicated by the primary TaskFactory
	PackerTaskOP ptask = task_factory()->create_task_and_apply_taskoperations( pose );

	std::set< Size > packing_positions;
	vector1< Size > design_positions;
	for ( Size pos(1); pos < ptask->total_residue(); ++pos ) {
		if ( ! pose.residue_type(pos).is_protein() ) continue;
		if ( ptask->being_packed( pos ) ) packing_positions.insert( pos );
		if ( ptask->being_designed( pos ) ) design_positions.push_back( pos );
	}
	if ( design_positions.empty() ) {
		TR << "WARNING: no designable positions" << std::endl;
		if ( packing_positions.empty() ) {
			TR << "WARNING: no packable positions either" << std::endl;
		} else {
			TR << "using packable positions instead of designable positions "
				<< "to define movable backbone regions" << std::endl;
			design_positions.insert( design_positions.begin(),
				packing_positions.begin(), packing_positions.end() );
		}
	}

	// debug a.k.a. level 400 output
	std::copy( design_positions.begin(), design_positions.end(),
		std::ostream_iterator<Size>(TR(t_debug),",") );

	//this 'second-shell' TaskFactory will be used for packer operations during the backbone protocols
	TaskFactoryOP task_factory2( new TaskFactory );

	task_factory2->push_back( operation::TaskOperationCOP( new IncludeCurrent ) );
	// the following will disable packing outside of the neighborhood around design_positions

	toolbox::task_operations::RestrictToNeighborhoodOperationOP nbop( new toolbox::task_operations::RestrictToNeighborhoodOperation( packing_positions ) );
	task_factory2->push_back( nbop );

	// option-dependent designable second shell
	if ( ! designable_second_shell_ ) task_factory2->push_back( TaskOperationCOP( new RestrictToRepacking ) );

	// prevent DNA packing/designing
	task_factory2->push_back( TaskOperationCOP( new OperateOnCertainResidues( ResLvlTaskOperationOP(new RestrictToRepackingRLT), ResFilterOP(new ResidueHasProperty("DNA")) ) ) );

	// make loops
	loops::LoopsOP loops_to_move( new Loops() );

	loops::loops_around_residues( *loops_to_move, pose, design_positions, gapspan_, spread_ );
	if ( loops_to_move->size() == 0 ) {
		TR << "WARNING: no loop regions were defined, aborting backbone design" << std::endl;
		return;
	}
	set_loop_info( pose, *loops_to_move );

	// rewrite foldtree to fix loop termini (prevents propagation of movement beyond loop)
	kinematics::FoldTree ft_new, ft_orig( pose.fold_tree() );
	loops::fold_tree_from_loops( pose, *loops_to_move, ft_new );
	pose.fold_tree( ft_new );

	// finally call backbone movement protocol
	if ( type_ == "ccd" ) {
		ccd( pose, loops_to_move, task_factory2 ); // creates chainbreaks
	} else if ( type_ == "backrub" ) {
		backrub( pose, loops_to_move, task_factory2 ); // conservative
	} else {
		TR << "Unknown backbone movement type: " << type_ << '\n';
	}

	loops::remove_cutpoint_variants( pose ); // purpose, exactly?
	pose.fold_tree( ft_orig ); // restore original fold tree
}

// XRW TEMP std::string
// XRW TEMP DesignProteinBackboneAroundDNA::get_name() const {
// XRW TEMP  return DesignProteinBackboneAroundDNA::mover_name();
// XRW TEMP }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void DesignProteinBackboneAroundDNA::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("type") ) type_ = tag->getOption< std::string >("type");
	if ( tag->hasOption("gapspan") ) gapspan_ = tag->getOption< Size >("gapspan");
	if ( tag->hasOption("spread") ) spread_ = tag->getOption< Size >("spread");
	if ( tag->hasOption("cycles_outer") ) cycles_outer_ = tag->getOption<Size>("cycles_outer");
	if ( tag->hasOption("cycles_inner") ) cycles_inner_ = tag->getOption<Size>("cycles_inner");
	if ( tag->hasOption("repack_rate") ) repack_rate_ = tag->getOption<Size>("repack_rate");
	if ( tag->hasOption("temp_initial") ) temp_initial_ = tag->getOption<Real>("temp_initial");
	if ( tag->hasOption("temp_final") ) temp_final_ = tag->getOption<Real>("temp_final");
	if ( tag->hasOption("designable_second_shell") ) {
		designable_second_shell_ = tag->getOption<bool>("designable_second_shell");
	}
	// the following are calls to PackRotamersMover methods
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DesignProteinBackboneAroundDNA::fresh_instance() const
{
	return moves::MoverOP( new DesignProteinBackboneAroundDNA );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DesignProteinBackboneAroundDNA::clone() const
{
	return moves::MoverOP( new DesignProteinBackboneAroundDNA( *this ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private methods
void
DesignProteinBackboneAroundDNA::set_loop_info( Pose const & pose, Loops const & loops )
{
	// clear out (any) pre-existing string info (from previous apply calls)
	info().clear();
	// store information about the loops that were modeled
	for ( auto const & loop : loops ) {
		Size const start( loop.start() ), cut( loop.cut() ), stop( loop.stop() );
		std::ostringstream loopinfo;
		loopinfo << "REMARK loop: ";
		if ( pose.pdb_info() ) {
			loopinfo << pose.pdb_info()->chain( start ) << "-" << pose.pdb_info()->number( start ) << ">"
				<< pose.pdb_info()->chain( cut ) << "-" << pose.pdb_info()->number( cut ) << ">"
				<< pose.pdb_info()->chain( stop ) << "-" << pose.pdb_info()->number( stop );
		} else {
			loopinfo << pose.chain( start ) << "-" << start << ">"
				<< pose.chain( cut ) << "-" << cut << ">"
				<< pose.chain( stop ) << "-" << stop;
		}
		info().push_back( loopinfo.str() );
	}
}

void
DesignProteinBackboneAroundDNA::ccd(
	Pose & pose,
	loops::LoopsOP const loops,
	TaskFactoryCOP task_factory2
)
{
	loops::loop_mover::refine::LoopMover_Refine_CCDOP refine_ccd( new loops::loop_mover::refine::LoopMover_Refine_CCD( loops, score_function()->clone() ) );
	refine_ccd->outer_cycles( cycles_outer_ );
	refine_ccd->max_inner_cycles( cycles_inner_ );
	refine_ccd->repack_period( repack_rate_ );
	refine_ccd->temp_initial( temp_initial_ );
	refine_ccd->temp_final( temp_final_ );
	refine_ccd->set_native_pose( PoseCOP( PoseOP( new Pose( pose ) ) ) );
	refine_ccd->set_task_factory( task_factory2 );
	refine_ccd->apply( pose );
}

void
DesignProteinBackboneAroundDNA::backrub(
	Pose & pose,
	loops::LoopsOP const loops,
	TaskFactoryCOP task_factory2
)
{
	ScoreFunctionOP br_scorefxn( score_function()->clone() );
	br_scorefxn->set_weight( mm_bend, 0.5 );
	// pivot atoms default to "CA" so that non-protein atoms are not considered during backrub scoring
	methods::EnergyMethodOptions emo( br_scorefxn->energy_method_options() );
	emo.bond_angle_central_atoms_to_score( option[ OptionKeys::backrub::pivot_atoms ] );
	br_scorefxn->set_energy_method_options( emo );

	protocols::backrub::BackrubMover backrubmover;
	// read known and unknown optimization parameters from the database
	backrubmover.branchopt().read_database();

	// this mover appears to require a separate copy of the input pose
	PoseCOP input_pose( PoseOP( new Pose( pose ) ) );
	backrubmover.set_input_pose( input_pose ); // virtual funtion in Mover base class

	// set up backrub segments
	backrubmover.clear_segments();

	for ( auto const & loop : *loops ) {
		Size const start( loop.start() ), stop( loop.stop() );
		backrubmover.add_segment(
			id::AtomID( pose.residue_type( start ).atom_index(" CA "), start ),
			id::AtomID( pose.residue_type( stop ).atom_index(" CA "), stop )
		);
	}

	// read MonteCarlo / cycles options
	Real const dtemp( ( temp_final_ - temp_initial_ ) / cycles_outer_ );

	for ( Size i(1); i <= cycles_outer_; ++i ) {
		moves::MonteCarlo mc( pose, *br_scorefxn, temp_initial_ + dtemp*i );

		for ( Size j(1); j <= cycles_inner_; ++j ) {
			backrubmover.apply( pose );
			mc.boltzmann( pose );
			if ( ( j % repack_rate_ ) == 0 ) {
				PackerTaskOP ptask =
					task_factory2->create_task_and_apply_taskoperations( pose );
				pack_rotamers( pose, *br_scorefxn, ptask );
				mc.boltzmann( pose );
			}
		}
		mc.show_counters();
		pose = mc.last_accepted_pose();
	}
}

std::string DesignProteinBackboneAroundDNA::get_name() const {
	return mover_name();
}

std::string DesignProteinBackboneAroundDNA::mover_name() {
	return "DesignProteinBackboneAroundDNA";
}

void DesignProteinBackboneAroundDNA::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction type_enum;
	type_enum.name("bb_movement_types");
	type_enum.base_type(xs_string);
	type_enum.add_restriction(xsr_enumeration, "ccd");
	type_enum.add_restriction(xsr_enumeration, "backrub");
	xsd.add_top_level_element(type_enum);

	attlist + XMLSchemaAttribute(
		"type", "bb_movement_types",
		"Backbone movement type");

	attlist + XMLSchemaAttribute(
		"gapspan", xsct_non_negative_integer,
		"Number of residues allowed between dna contacting residues "
		"for loop building.");

	attlist + XMLSchemaAttribute(
		"spread", xsct_non_negative_integer,
		"Loop start/ends are extended by this number of residues.");

	attlist + XMLSchemaAttribute(
		"cycles_outer", xsct_non_negative_integer,
		"outer cycles of CCD/backrub (MC)");

	attlist + XMLSchemaAttribute(
		"cycles_inner", xsct_non_negative_integer,
		"inner cycles of CCD/backrub (MC)");

	attlist + XMLSchemaAttribute(
		"repack_rate", xsct_non_negative_integer,
		"repack rate of CCD/backrub (MC)");

	attlist + XMLSchemaAttribute(
		"temp_initial", xsct_real,
		"Initial temp of backrub (MC)");

	attlist + XMLSchemaAttribute(
		"temp_final", xsct_real,
		"Final temp of backrub (MC)");

	attlist + XMLSchemaAttribute(
		"designable_second_shell", xsct_rosetta_bool,
		"Allow design of residues in the neighborhood of the design shell");

	rosetta_scripts::attributes_for_parse_score_function(attlist);
	rosetta_scripts::attributes_for_parse_task_operations(attlist);

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Performs a round flexible backbone sampling/design and loop design around DNA.",
		attlist );
}

std::string DesignProteinBackboneAroundDNACreator::keyname() const {
	return DesignProteinBackboneAroundDNA::mover_name();
}

protocols::moves::MoverOP
DesignProteinBackboneAroundDNACreator::create_mover() const {
	return protocols::moves::MoverOP( new DesignProteinBackboneAroundDNA );
}

void DesignProteinBackboneAroundDNACreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DesignProteinBackboneAroundDNA::provide_xml_schema( xsd );
}


} // namespace dna
} // namespace protocols
