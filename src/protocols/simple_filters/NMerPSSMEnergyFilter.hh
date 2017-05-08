// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NMerPSSMEnergyFilter.hh
/// @brief definition of filter class NMerPSSMEnergyFilter.
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_NMerPSSMEnergyFilter_hh
#define INCLUDED_protocols_simple_filters_NMerPSSMEnergyFilter_hh

//unit headers
#include <protocols/simple_filters/NMerPSSMEnergyFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/methods/NMerPSSMEnergy.hh>

namespace protocols {
namespace simple_filters {

class NMerPSSMEnergyFilter : public filters::Filter
{
public:
	//default ctor
	NMerPSSMEnergyFilter();
	//full ctor
	NMerPSSMEnergyFilter(
		core::Real const score_type_threshold,
		std::string string_resnums
	);
	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new NMerPSSMEnergyFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new NMerPSSMEnergyFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute_residue( core::pose::Pose const & pose, core::Size const seqpos ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	~NMerPSSMEnergyFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real score_type_threshold_;
	std::string string_resnums_;
	core::scoring::methods::NMerPSSMEnergy energy_method_;
};

}
}

#endif
