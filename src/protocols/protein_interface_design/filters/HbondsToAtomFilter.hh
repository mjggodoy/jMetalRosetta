// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/HbondsToAtomFilter.hh
/// @brief definition of filter classes for iterations of docking/design.
/// @author Lei Shi (shileiustc @gmail.com)

#ifndef INCLUDED_protocols_protein_interface_design_filters_HbondsToAtomFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_HbondsToAtomFilter_hh


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


// C++ headers

// Unit headers
//#include <basic/datacache/DataMap.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

/// @brief returns true if the number of hbonding partners to a particular residue exceeds a certain value
/// This filter is useful in conjunction with DesignMinimizeHbonds class
class HbondsToAtomFilter : public protocols::filters::Filter
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
public :
	HbondsToAtomFilter() : Filter( "HbondsToAtom" ) {}
	HbondsToAtomFilter( Size const resnum, Size const partners, Real const energy_cutoff=-0.5,
		bool const backbone=false, bool const sidechain=true ) : Filter( "HbondsToAtom" ) {
		resnum_ = resnum; partners_ = partners; energy_cutoff_ = energy_cutoff; backbone_ = backbone;
		sidechain_ = sidechain;
		runtime_assert( backbone_ || sidechain_ );
		runtime_assert( partners_ );
		runtime_assert( energy_cutoff_ <= 0 );
	}
	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override {
		return protocols::filters::FilterOP( new HbondsToAtomFilter( *this ) );
	}
	protocols::filters::FilterOP fresh_instance() const override {
		return protocols::filters::FilterOP( new HbondsToAtomFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~HbondsToAtomFilter();

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
	Size resnum_, partners_;
	Real energy_cutoff_;
	std::string atomdesg_;
	bool backbone_, sidechain_, bb_bb_;
};

}
} // protein_interface_design
} // devel


#endif /*INCLUDED_DOCK_DESIGN_FILTERS_H_*/
