// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_Rotates_hh
#define INCLUDED_protocols_ligand_docking_Rotates_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/Rotates.fwd.hh>
#include <protocols/ligand_docking/Rotate.fwd.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Project Headers
#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class Rotates: public protocols::moves::Mover
{
public:
	Rotates();
	~Rotates() override;
	Rotates(Rotates const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	RotateOPs rotates_;
};


} //namespace ligand_docking
} //namespace protocols

#endif
