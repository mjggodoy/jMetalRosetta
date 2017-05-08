// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief adapted from backrub mover in protein interface design
/// @author Sinisa Bjelic sinibjelic@gmail.com

#ifndef INCLUDED_protocols_enzdes_BackboneSampler_hh
#define INCLUDED_protocols_enzdes_BackboneSampler_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>

#include <utility/vector1.hh>

// C++ headers

namespace protocols {
namespace enzdes {

class BackboneSampler : public protocols::ligand_docking::LigandBaseProtocol
{
public:
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pose::Pose Pose;

public:
	BackboneSampler();

	BackboneSampler(
		ScoreFunctionCOP scorefxn,
		core::Size const bb_moves,
		core::Real const mc_kt
	);

	void
	apply( Pose & pose ) override;

	/*
	std::string get_name() const override;
	*/
	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override
	{
		return protocols::moves::MoverOP( new BackboneSampler );
	}

	~BackboneSampler() override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) override;

private:
	core::Size bb_moves_;
	core::Real mc_kt_;
	core::scoring::ScoreFunctionOP scorefxn_repack_;
};

} // enzdes
} // protocols

#endif

