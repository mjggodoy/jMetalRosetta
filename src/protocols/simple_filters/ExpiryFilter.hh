// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ExpiryFilter.hh
/// @brief Simple filter that tests whether a file exists. Useful to test whether we're recovering from a checkpoint
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_filters_ExpiryFilter_hh
#define INCLUDED_protocols_simple_filters_ExpiryFilter_hh

//unit headers
#include <protocols/simple_filters/ExpiryFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class ExpiryFilter : public filters::Filter
{
public:
	//default ctor
	ExpiryFilter();
	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new ExpiryFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new ExpiryFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const &pose ) const;
	~ExpiryFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	core::Size seconds() const;
	void seconds( core::Size const s );
	core::Size start_time() const; // the time at start

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size seconds_;
	core::Size start_time_;
};

}
}

#endif
