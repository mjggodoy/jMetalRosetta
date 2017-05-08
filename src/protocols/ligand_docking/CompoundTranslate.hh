// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/CompoundTranslate.hh
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_CompoundTranslate_hh
#define INCLUDED_protocols_ligand_docking_CompoundTranslate_hh

// Unit Headers
#include <protocols/ligand_docking/Translate.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class CompoundTranslate : public protocols::moves::Mover
{
public:
	CompoundTranslate();
	~CompoundTranslate() override;
	CompoundTranslate(CompoundTranslate const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
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
	TranslateOPs translates_;
	bool randomize_order_;
	bool allow_overlap_;
};

} //namespace ligand_docking
} //namespace protocols

#endif
