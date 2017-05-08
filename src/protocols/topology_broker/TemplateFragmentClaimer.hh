// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_TemplateFragmentClaimer_hh
#define INCLUDED_protocols_topology_broker_TemplateFragmentClaimer_hh


// Unit Headers
#include <protocols/topology_broker/TemplateFragmentClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/abinitio/Templates.hh>
#include <utility/vector1.hh>


//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {

/// @brief hacky wrapper to keep the old Template code alive a bit longer
/// this claimer deals with the Jumpy part of the Templates.
class TemplateFragmentClaimer : public FragmentClaimer {
public:
	TemplateFragmentClaimer(); //for factory
	TemplateFragmentClaimer( std::string config_file, core::Size fragsize, weights::AbinitioMoverWeightOP weight = nullptr );

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new TemplateFragmentClaimer( *this ) );
	}

	void read_config_file( std::string const& file );

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "TemplateFragmentClaimer";
	}

protected:
	bool read_tag( std::string tag, std::istream& is ) override;
	void init_after_reading() override;
private:
	// info about homologues structures --- if available
	abinitio::TemplatesOP templates_;
	core::Size frag_size_;
	std::string config_file_;
}; //class TemplateFragmentClaimer

}
}

#endif
