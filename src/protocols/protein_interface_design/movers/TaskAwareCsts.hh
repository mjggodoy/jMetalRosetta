// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/TaskAwareCsts.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_TaskAwareCsts_hh
#define INCLUDED_protocols_protein_interface_design_movers_TaskAwareCsts_hh
#include <protocols/protein_interface_design/movers/TaskAwareCsts.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief applies csts (currently only coordinate csts) to every designable position in pose according to taskoperations
class TaskAwareCsts : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	TaskAwareCsts();
	void apply( Pose & pose ) override;
	// XRW TEMP  virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new TaskAwareCsts ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~TaskAwareCsts();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	std::string cst_type() const { return cst_type_; }
	void cst_type( std::string const & type ) { cst_type_ = type; }

	std::string anchor_resnum() const { return anchor_resnum_; }
	void anchor_resnum( std::string const & c ) { anchor_resnum_ = c; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::pack::task::TaskFactoryOP task_factory_; /// dflt NULL (design all); designable residues will have the constraints applied to them.
	std::string cst_type_; //dflt coordinate;
	std::string anchor_resnum_; // dflt ""; if "", use the first residue as the anchor point
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_TaskAwareCsts_HH*/
