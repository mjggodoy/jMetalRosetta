// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/membrane/MembranePositionFromTopologyMover.cc
///
/// @brief      Computes and sets the initial position of the membrane
/// @details Computes and sets the initial position of the membrane from
///    sequence or structure (can be specified by the user at construction
///    or as a setup cmd flag).
///    CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///
///    NOTE: Requires a membrane pose!
///    NOTE: sequence not yet implemented
///    Last Modified: 6/21/14
///
/// @author  Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc

// Unit Headers
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMoverCreator.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Residue.hh>
#include <protocols/membrane/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.MembranePositionMoverFromTopologyMover" );

namespace protocols {
namespace membrane {

/// @brief Defualt Constructor
/// @details Compute the embedding of the pose based on xyz coordinates
/// and spanning topology provided in MembraneInfo
MembranePositionFromTopologyMover::MembranePositionFromTopologyMover() :
	protocols::moves::Mover(),
	anchor_at_res1_( true )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover
MembranePositionFromTopologyMover::MembranePositionFromTopologyMover( MembranePositionFromTopologyMover const & ) = default;

/// @brief Destructor
MembranePositionFromTopologyMover::~MembranePositionFromTopologyMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MembranePositionFromTopologyMover::clone() const {
	return ( protocols::moves::MoverOP( new MembranePositionFromTopologyMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembranePositionFromTopologyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MembranePositionFromTopologyMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembranePositionFromTopologyMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	if ( tag->hasOption( "anchor_at_res1" ) ) {
		anchor_at_res1_ = tag->getOption< bool >("anchor_at_res1");
	}

}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MembranePositionFromTopologyMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MembranePositionFromTopologyMover );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP MembranePositionFromTopologyMoverCreator::keyname() const {
// XRW TEMP  return MembranePositionFromTopologyMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP MembranePositionFromTopologyMover::mover_name() {
// XRW TEMP  return "MembranePositionFromTopologyMover";
// XRW TEMP }

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Anchor membrane at residue 1, default is true
void MembranePositionFromTopologyMover::anchor_at_res1( bool truefalse ) {
	anchor_at_res1_ = truefalse;
}

/// @brief Update Membrane position in pose
/// @details Compute membrane posiiton based on structure
/// and then call pose.update_membrane_position() to update the membrane position
void
MembranePositionFromTopologyMover::apply( core::pose::Pose & pose ) {

	using namespace core;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	// Check pose is a membrane pose
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Warning! Pose is not a membrane pose. Cannot perform mpframework operation on a non membrane pose~" );
	}
	if ( pose.conformation().membrane_info()->spanning_topology()->nspans() == 0 ) {
		utility_exit_with_message("The SpanningTopology object in MembraneInfo is empty!" );
	}

	// remember foldtree
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// if membrane not at root, reorder foldtree to have pose TM COM at root
	if ( is_membrane_fixed( pose ) ) {
		TR << "Reordering foldtree:" << std::endl;
		core::kinematics::FoldTree ft = pose.fold_tree();
		core::Size anchor;

		// for simple reorder
		if ( anchor_at_res1_ == true ) {
			anchor = 1;
		} else {
			// resetting the foldtree to anchor at the pose TM COM
			anchor = create_membrane_foldtree_anchor_pose_tmcom( pose );
		}

		// reorder foldtree
		ft.reorder( anchor );
		pose.fold_tree( ft );
	}

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

	// Initialize starting vectors
	Vector normal( 0, 0, 0 );
	Vector center( 0, 0, 0 );

	// Compute position from pose
	TR << "Computing initial membrane position from structure..." << std::endl;
	compute_structure_based_embedding( pose, center, normal );

	// Update membrane position - shift normal along center
	pose.conformation().update_membrane_position( center, normal );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree after reset: Is membrane fixed? " << is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

} // apply

/// @brief Get the name of this mover
// XRW TEMP std::string
// XRW TEMP MembranePositionFromTopologyMover::get_name() const {
// XRW TEMP  return "MembranePositionFromTopologyMover";
// XRW TEMP }

std::string MembranePositionFromTopologyMover::get_name() const {
	return mover_name();
}

std::string MembranePositionFromTopologyMover::mover_name() {
	return "MembranePositionFromTopologyMover";
}

void MembranePositionFromTopologyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute("anchor_at_res1", xsct_rosetta_bool, "Should we do a simple reorder anchoring the membrane at residue 1?");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Computes and sets the initial position of the membrane", attlist );
}

std::string MembranePositionFromTopologyMoverCreator::keyname() const {
	return MembranePositionFromTopologyMover::mover_name();
}

protocols::moves::MoverOP
MembranePositionFromTopologyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MembranePositionFromTopologyMover );
}

void MembranePositionFromTopologyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MembranePositionFromTopologyMover::provide_xml_schema( xsd );
}


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc
