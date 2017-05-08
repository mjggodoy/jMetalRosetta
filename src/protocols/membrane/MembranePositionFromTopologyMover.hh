// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/membrane/MembranePositionFromTopologyMover.hh
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

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh

// Unit Headers
#include <protocols/membrane/MembranePositionFromTopologyMover.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/types.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace membrane {

/// @brief Compute the initial position of the membrane based upon sequence
/// or structure
class MembranePositionFromTopologyMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Default Constructor
	/// @details Compute the embedding of the pose based on xyz coordinates
	/// and spanning topology provided in MembraneInfo
	MembranePositionFromTopologyMover();

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover
	MembranePositionFromTopologyMover( MembranePositionFromTopologyMover const & src );

	/// @brief Destructor
	~MembranePositionFromTopologyMover() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Update Membrane position in pose
	/// @details Compute membrane posiiton based on sequence or structure
	/// and then call pose.update_membrane_position() to update the membrane position
	void apply( Pose & pose ) override;

	/// @brief Get the name of this mover
	// XRW TEMP  std::string get_name() const override;

	/// @brief Anchor membrane at residue 1, default is true
	void anchor_at_res1( bool truefalse );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // data

	// Anchor membrane residue at residue 1
	bool anchor_at_res1_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh
