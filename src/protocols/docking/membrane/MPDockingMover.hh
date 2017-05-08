// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/MPDockingMover.fwd.hh
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (6/24/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMover_hh
#define INCLUDED_protocols_docking_membrane_MPDockingMover_hh

// Unit Headers
#include <protocols/docking/membrane/MPDockingMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.fwd.hh>
#include <protocols/docking/DockingProtocol.fwd.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/PyMOLMover.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
//#include <basic/options/option.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

namespace protocols {
namespace docking {
namespace membrane {

/// @brief Docks two proteins together in the membrane bilayer
/// @details Requires running mpdocking_setup first to create a single pose;
///    Should also run docking_prepack before

class MPDockingMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Docks two proteins with default normal=(0,0,1) and center=(0,0,0)
	MPDockingMover( bool lowres=true, bool highres=true );

	/// @brief Constructor with jump number to dock over
	/// @details Docks two proteins with default normal=(0,0,1) and center=(0,0,0)
	MPDockingMover( Size jump_num, bool lowres=true, bool highres=true );

	/// @brief Copy Constructor
	/// @details Create a deep copy of this mover
	MPDockingMover( MPDockingMover const & src );

	/// @brief Destructor
	virtual ~MPDockingMover();

public: // methods

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

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

	/// @brief Get the name of this Mover (MPDockingMover)
	// XRW TEMP  virtual std::string get_name() const;

	/// @brief Setters for the DockingProtocol
	/// @details These are set here because the mover doesn't inherit from
	///    the DockingProtocol
	void set_cycles_inner( Size cycles_inner );

	/// @brief Set outer cycles in DockingProtocol
	void set_cycles_outer( Size cycles_outer );

	/// @brief Add membrane components to the pose, then dock proteins along
	///   the flexible jump
	void apply( Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods

	// setup
	void set_defaults( const Pose & pose );

	// register options with JD2
	void register_options();

	// overwrite defaults with stuff from cmdline
	void init_from_cmd();

	// finalize setup
	void finalize_setup();

private: // data

	// add membrane mover
	protocols::membrane::AddMembraneMoverOP add_membrane_mover_;

	// docking protocol
	DockingProtocolOP docking_protocol_;

	// lowres and highres
	bool lowres_;
	bool highres_;

	// scorefunction
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	// Membrane Center/Normal pair used for setup
	core::Vector center_;
	core::Vector normal_;

	// jump number to dock over
	Size jump_num_;

	// native for RMSD calculation
	core::pose::PoseOP native_;

};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMover_hh
