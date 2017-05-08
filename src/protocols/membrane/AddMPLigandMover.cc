// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/membrane/AddMPLigandMover.cc
///
/// @brief  Add "single" ligand to to membrane pose
/// @details  Accommodate membrane protein ligand in the membrane framework by
///    reorganizing the current foldtree. Resulting foldtree will
///    keep the membrane attached to the transmembrane COM and ligand to the
///    closest binding pocket residue, provided in the constructor.
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author JKLeman (julia.koehler1982@gmail.com)
/// #RosettaMPMover

// Unit Headers
#include <protocols/membrane/AddMPLigandMover.hh>
#include <protocols/membrane/AddMPLigandMoverCreator.hh>
#include <protocols/membrane/util.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.AddMPLigandMover" );

namespace protocols {
namespace membrane {

////////////////////
/// Constructors ///
////////////////////

/// @brief Add membrane protein ligand mover
/// @details Attach ligand downstream in the foldtree
/// for refinement at the last residue as a default.
/// DO NOT USE
AddMPLigandMover::AddMPLigandMover() :
	protocols::moves::Mover(),
	closest_rsd_( 0 ),
	ligand_seqpos_( 0 )
{}

/// @brief Add Membrane protein ligand mover (custom)
/// @details Attach ligand downstream in the foldtree of the
/// closest residue to the binding pocket
AddMPLigandMover::AddMPLigandMover( core::Size closest_rsd, core::Size ligand_seqpos ) :
	protocols::moves::Mover(),
	closest_rsd_( closest_rsd ),
	ligand_seqpos_( ligand_seqpos )
{}

/// @brief Copy Constructor
/// @details Mkae a deep copy of this mover
AddMPLigandMover::AddMPLigandMover( AddMPLigandMover const & ) = default;

/// @brief Destructor
AddMPLigandMover::~AddMPLigandMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
AddMPLigandMover::clone() const {
	return ( protocols::moves::MoverOP( new AddMPLigandMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
AddMPLigandMover::fresh_instance() const {
	return ( protocols::moves::MoverOP( new AddMPLigandMover() ) );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
AddMPLigandMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// Read in closest residue option
	if ( tag->hasOption( "closest_rsd" ) ) {
		closest_rsd_ = tag->getOption< Size >( "closest_rsd" );
	}

	// Read in sequence position of the ligand
	if ( tag->hasOption( "ligand_seqpos" ) ) {
		ligand_seqpos_ = tag->getOption< Size >( "ligand_seqpos" );
	}

}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddMPLigandMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddMPLigandMover );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP AddMPLigandMoverCreator::keyname() const {
// XRW TEMP  return AddMPLigandMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP AddMPLigandMover::mover_name() {
// XRW TEMP  return "AddMPLigandMover";
// XRW TEMP }

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Mover Apply Method
void
AddMPLigandMover::apply( core::pose::Pose & pose ) {

	using namespace core::kinematics;

	// Check the pose is a membrane framework pose
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot perform add ligand to membrane pose operation because this pose is not a membrane pose!" );
	}

	TR << "Closest Residue: " << closest_rsd_ << " Ligand Seqpos: " << ligand_seqpos_ << std::endl;
	TR << "WARNING: MAKE SURE YOUR LIGAND CHAIN IS DIFFERENT FROM ANY PROTEIN CHAINS!" << std::endl;

	// Check the closest rsd and ligand seqpos parameters are valid
	// Note - default constructor parameters will make at least one of these statements
	// fail so please don't use it!
	if ( closest_rsd_ <= 0 || closest_rsd_ > pose.size() ) {
		utility_exit_with_message( "User specified closest residue to ligand is out of bounds.");
	}

	if ( ligand_seqpos_ <= 0 || ligand_seqpos_ > pose.size() ) {
		utility_exit_with_message( "User specified ligand position in sequence is out of bounds.");
	}

	if ( closest_rsd_ == ligand_seqpos_ ) {
		utility_exit_with_message( "Cannot sepcify closest residue as the ligand. Self-attachment is not valid");
	}

	// Get the following parameters from the foldtree: Residue COM, MP rsd position
	core::Size mp_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );

	// get all pose chains
	utility::vector1< core::Size > chain_end_residues( chain_end_res( pose ) );

	// get chainid of user-defined ligand anchor point
	core::Size ligand_anchor_chain = static_cast< core::Size >( pose.chain( closest_rsd_ ) );
	TR << "ligand anchor chain " << ligand_anchor_chain << std::endl;

	// get the anchor residues as chain tm COMs
	utility::vector1< core::Size > anchors( get_anchor_points_for_tmcom( pose ) );

	// create an anchors vector that contains all jumps and cutpoints
	utility::vector1< core::Vector > anchors_vector;
	for ( core::Size i = 1; i <= chain_end_residues.size()-2; ++i ) {

		// membrane residue
		if ( i == 1 ) {
			core::Vector jump1( anchors[ 1 ], mp_rsd, mp_rsd-1 );
			anchors_vector.push_back( jump1 );
		}

		// ligand
		if ( i == ligand_anchor_chain ) {
			core::Vector jump2( closest_rsd_, ligand_seqpos_, ligand_seqpos_-1 );
			anchors_vector.push_back( jump2 );
		}

		// else
		if ( i != 1 && i != ligand_anchor_chain && i != static_cast< core::Size >( pose.chain( ligand_seqpos_ ) ) && i != mp_rsd ) {
			core::Vector jump( anchors[ 1 ], anchors[ i ], chain_end_residues[i-1] );
			anchors_vector.push_back( jump );
		}
	}

	// create membrane foldtree from anchors, membrane jump is automatically set to 1
	create_specific_membrane_foldtree( pose, anchors_vector );

	// reorder such that pose_tm_com is at root
	FoldTree ft = pose.fold_tree();
	ft.reorder( anchors[1] );
	pose.fold_tree( ft );
	pose.fold_tree().show( TR );

}

/// @brief Show the name of this mvoer
// XRW TEMP std::string
// XRW TEMP AddMPLigandMover::get_name() const {
// XRW TEMP  return "AddMPLigandMover";
// XRW TEMP }

std::string AddMPLigandMover::get_name() const {
	return mover_name();
}

std::string AddMPLigandMover::mover_name() {
	return "AddMPLigandMover";
}

void AddMPLigandMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "closest_rsd", xsct_positive_integer, "Index of closest residue to ligand" )
		+ XMLSchemaAttribute( "ligand_seqpos", xsct_positive_integer, "Rosetta sequence position of ligand. Note the ligand should be in a separate chain." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add a ligand to a membrane protein", attlist );
}

std::string AddMPLigandMoverCreator::keyname() const {
	return AddMPLigandMover::mover_name();
}

protocols::moves::MoverOP
AddMPLigandMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMPLigandMover );
}

void AddMPLigandMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddMPLigandMover::provide_xml_schema( xsd );
}


} // membrane
} // protocols


