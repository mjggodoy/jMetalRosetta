// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/BestHotspotCstMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/BestHotspotCstMover.hh>
#include <protocols/protein_interface_design/movers/BestHotspotCstMoverCreator.hh>
#include <protocols/protein_interface_design/util.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/tag/Tag.hh>

#include <core/id/AtomID.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.BestHotspotCstMover" );

// XRW TEMP std::string
// XRW TEMP BestHotspotCstMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return BestHotspotCstMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BestHotspotCstMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BestHotspotCstMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BestHotspotCstMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "BestHotspotCst";
// XRW TEMP }


BestHotspotCstMover::BestHotspotCstMover() :
	protocols::moves::Mover( BestHotspotCstMover::mover_name() )
{}

BestHotspotCstMover::BestHotspotCstMover(
	protocols::hotspot_hashing::HotspotStubSetOP stub_set,
	core::Size const host_chain,
	core::Size const n_resi
) :
	protocols::moves::Mover( BestHotspotCstMover::mover_name() ),
	host_chain_( host_chain ),
	n_resi_( n_resi )
{
	if ( stub_set ) stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *stub_set ) );
}

BestHotspotCstMover::BestHotspotCstMover( BestHotspotCstMover const & init ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( init ),
	host_chain_(init.host_chain_), n_resi_( init.n_resi_ )
{
	if ( init.stub_set_ ) stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *init.stub_set_ ) );
}

BestHotspotCstMover::~BestHotspotCstMover() {}


void
BestHotspotCstMover::apply( pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;
	HotspotStubSetOP working_stub_set( new HotspotStubSet( *stub_set_ ) ); // stub_set_ needs to remain pristine for next rounds of execution.
	core::scoring::ScoreFunctionCOP scorefxn( get_score_function() ); // default scorefxn for several definitions within the function

	// find the residues that have the best constraint backbone_stub_constraint scores
	utility::vector1< core::Size > const best_cst_residues( best_bbcst_residues( pose, host_chain_, n_resi_ ) );

	// make a packer task containing only the best residues
	core::pack::task::PackerTaskOP packer_task = core::pack::task::TaskFactory::create_packer_task( pose );
	for ( core::Size i=1; i <= pose.size(); ++i ) {
		if ( ( find( best_cst_residues.begin(), best_cst_residues.end(), i ) != best_cst_residues.end() ) )  continue;
		else packer_task->nonconst_residue_task( i ).prevent_repacking();
	}

	// Assign a fixed residue (for the stub constraints)
	core::Size fixed_res(1);
	if ( host_chain_ == 1 ) fixed_res = pose.size();
	core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );

	// reapply cst's, but use the packer task to only count the best residues
	working_stub_set->remove_all_hotspot_constraints( pose );
	working_stub_set->add_hotspot_constraints_to_pose(
		pose,
		fixed_atom_id,
		packer_task,
		working_stub_set,
		cb_force_constant_,
		-0.5, //worst allowed stub bonus
		true, // apply_self_energies
		8.0, // bump_cutoff
		true ///apply_ambiguous constraints
	);

	TR<<"Reapplied constraints to residues ";
	for ( core::Size i=1; i<=best_cst_residues.size(); ++i ) {
		TR<< best_cst_residues[i] << " ";
	}
	TR << std::endl;
}

// XRW TEMP std::string
// XRW TEMP BestHotspotCstMover::get_name() const {
// XRW TEMP  return BestHotspotCstMover::mover_name();
// XRW TEMP }

void BestHotspotCstMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;

	core::pose::PDBPoseMap const pose_map( pose.pdb_info()->pdb2pose() );

	host_chain_ = tag->getOption<core::Size>( "chain_to_design", 2 );
	n_resi_ = tag->getOption<core::Size>( "best_n", 3 );
	cb_force_constant_ = tag->getOption< core::Real >( "cb_force", 1.0 );

	std::string const hs( "hotspot_stubset" );
	stub_set_ = data.get_ptr<HotspotStubSet>( "constraints", hs );

	TR<<"BestHotspotCst mover on chain "<<host_chain_<<" with cbeta force " << cb_force_constant_ << "\n";
}

std::string BestHotspotCstMover::get_name() const {
	return mover_name();
}

std::string BestHotspotCstMover::mover_name() {
	return "BestHotspotCst";
}

void BestHotspotCstMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_non_negative_integer, "Chain to design, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "best_n", xsct_non_negative_integer, "Number of best-constrained residues for detailed analysis", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Force on CB atoms applied in the hotspot constraints", "1.0" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string BestHotspotCstMoverCreator::keyname() const {
	return BestHotspotCstMover::mover_name();
}

protocols::moves::MoverOP
BestHotspotCstMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BestHotspotCstMover );
}

void BestHotspotCstMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BestHotspotCstMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
