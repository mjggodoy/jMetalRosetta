// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/SnugDock.hh
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_antibody_snugdock_SnugDock_HH
#define INCLUDED_protocols_antibody_snugdock_SnugDock_HH

// Unit headers
#include <protocols/antibody/snugdock/SnugDock.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>

// Package headers
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/RefineOneCDRLoop.fwd.hh>
#include <protocols/antibody/CDRsMinPackMin.fwd.hh>


// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace antibody {
namespace snugdock {

/// @brief MoveSet step of Snugdock as described in:
///
///Sircar A, Gray JJ (2010) SnugDock: Paratope Structural Optimization duringAntibody-Antigen Docking Compensates
/// for Errors in Antibody Homology Models.
/// PLoS Comput Biol 6(1): e1000644. doi:10.1371/journal.pcbi.1000644
///
/// @details This is the main high resolution step of the overall SnugDock protocol.  See the Show method for a complete description of the algorithm.
///
class SnugDock: public docking::DockingHighRes {
public: // boiler plate / virtuals
	// default constructor
	SnugDock();

	// copy constructor
	SnugDock( SnugDock const & rhs );

	// assignment operator
	SnugDock & operator=( SnugDock const & rhs );

	// destructor
	virtual ~SnugDock();

	void apply( Pose & ) override;
	// XRW TEMP  virtual std::string get_name() const;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	bool reinitialize_for_new_input() const override;

	/// @brief Associates relevant options with the SnugDock class
	static void register_options();

public:
	void debug(); //Enables debug of the mover components of Snugdock - switches the random mover to sequence mover.
	Size number_of_high_resolution_cycles() const;
	void number_of_high_resolution_cycles( Size const number_of_high_resolution_cycles );
	void set_antibody_info( AntibodyInfoOP antibody_info );

	// Accessors for auto kink-constraint options
	bool high_res_kink_constraint() const { return high_res_kink_constraint_; }
	void high_res_kink_constraint( bool const setting ) { high_res_kink_constraint_ = setting; }

	void show( std::ostream & out=std::cout ) const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	friend std::ostream & operator<<(std::ostream& out, SnugDock const & snugdock );

private: // methods
	void setup_objects( Pose const & pose );
	void set_default();
	void init();
	void init_for_equal_operator_and_copy_constructor( SnugDock & lhs, SnugDock const & rhs);
	void init_from_options();

private: // data
	AntibodyInfoOP antibody_info_;
	moves::MonteCarloOP mc_;

	// Movers
	moves::MoverContainerOP high_resolution_step_;
	std::string loop_refinement_method_;
	antibody::CDRsMinPackMinOP pre_minimization_;

	// H3 filter options
	bool h3_filter_;
	bool debug_;
	Size h3_filter_tolerance_;

	// kink constraint
	bool high_res_kink_constraint_;

private:
	Size number_of_high_resolution_cycles_;

}; // class SnugDock

} // namespace snugdock
} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_snugdock_SnugDock_HH
