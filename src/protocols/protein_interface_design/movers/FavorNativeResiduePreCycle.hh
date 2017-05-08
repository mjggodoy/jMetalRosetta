// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/FavorNativeResiduePreCycle.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_FavorNativeResiduePreCycle_hh
#define INCLUDED_protocols_protein_interface_design_movers_FavorNativeResiduePreCycle_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/design_utils.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class FavorNativeResiduePreCycle : public protocols::moves::Mover
{
public:
	FavorNativeResiduePreCycle( core::Real const bonus = 1.5 ) : protocols::moves::Mover( "favor_native_residue" ), bonus_( bonus ) {
		runtime_assert( bonus_ >= 0 );
	}
	void apply( core::pose::Pose & pose ) override {
		FavorNativeResidue fnr( pose, bonus_ );
	}
	// XRW TEMP  virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const override {
		return( protocols::moves::MoverOP( new FavorNativeResiduePreCycle( bonus_ ) ) );
	}
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new FavorNativeResiduePreCycle ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~FavorNativeResiduePreCycle();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real bonus_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_FavorNativeResiduePreCycle_HH*/
