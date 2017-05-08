// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/RelativeSegmentFilter.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_RelativeSegmentFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_RelativeSegmentFilter_hh


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

/// @brief returns the residues aligned to a segment on the input pdb to the source pdb
class RelativeSegmentFilter : public protocols::filters::Filter
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
public :
	RelativeSegmentFilter() : Filter( "RelativeSegment" ) {}
	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override {
		return protocols::filters::FilterOP( new RelativeSegmentFilter( *this ) );
	}
	protocols::filters::FilterOP fresh_instance() const override {
		return protocols::filters::FilterOP( new RelativeSegmentFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	virtual ~RelativeSegmentFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	std::string source_pose() const { return source_pose_; }
	void source_pose( std::string const & s ) { source_pose_ = s; }
	core::Size start_res() const { return start_res_; }
	void start_res( core::Size const s ) { start_res_ = s;}
	core::Size stop_res() const { return stop_res_; }
	void stop_res( core::Size const s ) { stop_res_ = s; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string source_pose_;
	core::Size start_res_, stop_res_;
};

}
} // protein_interface_design
} // devel


#endif /*INCLUDED_DOCK_DESIGN_FILTERS_H_*/
