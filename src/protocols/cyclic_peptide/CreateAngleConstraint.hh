// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @file   protocols/cyclic_peptide/CreateAngleConstraint.hh
/// @brief  Add angle constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author modified by Parisa Hosseinzadeh (parisah@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CreateAngleConstraint_hh
#define INCLUDED_protocols_cyclic_peptide_CreateAngleConstraint_hh
#include <protocols/cyclic_peptide/CreateAngleConstraint.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace cyclic_peptide {

class CreateAngleConstraint : public moves::Mover
{
public:
	CreateAngleConstraint();
	~CreateAngleConstraint() override;

	void apply( core::pose::Pose &pose ) override;
	// XRW TEMP  std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	virtual void
	set (utility::vector1<Size>,utility::vector1<std::string>,utility::vector1<Size>, utility::vector1<std::string>, utility::vector1<Size>,utility::vector1<std::string>, utility::vector1<std::string>);

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1<Size> res_center_;
	utility::vector1<std::string> atom_center_;
	utility::vector1<Size> res1_;
	utility::vector1<std::string> atom1_;
	utility::vector1<Size> res2_;
	utility::vector1<std::string> atom2_;
	utility::vector1<std::string> cst_func_;
};

} // moves
} // protocols

#endif //INCLUDED_protocols_cyclic_peptide_CreateAngleConstraint_hh
