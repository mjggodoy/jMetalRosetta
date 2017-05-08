// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/IRmsdFilter.hh
/// @brief Calculate I_RMSD using the method from protocols::docking::metrics
/// @author Justin Porter (jrporter@jhu.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_IRmsdFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_IRmsdFilter_hh


//CHECK: inclusion ordering?
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/docking/types.hh>

#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class IRmsdFilter : public protocols::filters::Filter
{
public:
	IRmsdFilter();
	IRmsdFilter(protocols::docking::DockJumps const movable_jumps,
		core::Real const threshold,
		core::pose::PoseOP reference_pose
	);

	//@brief applies the filter to the input pose
	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override{
		return protocols::filters::FilterOP( new IRmsdFilter() );
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	//@brief computes I_rmsd by calling protocols::docking::metrics.cc
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~IRmsdFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real threshold_;
	core::pose::PoseOP reference_pose_;
	protocols::docking::DockJumps movable_jumps_;
	core::scoring::ScoreFunctionOP scorefxn_;

};

} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_IRmsdFilter_HH_
