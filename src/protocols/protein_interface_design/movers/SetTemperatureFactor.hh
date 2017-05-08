// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/SetTemperatureFactor.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_SetTemperatureFactor_hh
#define INCLUDED_protocols_protein_interface_design_movers_SetTemperatureFactor_hh
#include <protocols/protein_interface_design/movers/SetTemperatureFactor.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief Set the temperature (b-)factor column in the PDB according to som
/// filter's value. The filter needs to be ResId-compatible, i.e. to report
/// values on a per-residue basis.
class SetTemperatureFactor : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	SetTemperatureFactor();
	void apply( Pose & pose ) override;
	// XRW TEMP  virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SetTemperatureFactor ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~SetTemperatureFactor();
	void filter( protocols::filters::FilterOP filter );
	void scaling( core::Real const scaling );
	core::Real scaling() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	protocols::filters::FilterOP filter_;
	core::Real scaling_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_SetTemperatureFactor_HH*/
