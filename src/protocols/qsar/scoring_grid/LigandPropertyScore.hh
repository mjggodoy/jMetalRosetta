// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/LigandPropertyScore.hh
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_qsar_scoring_grid_LigandPropertyScore_HH
#define INCLUDED_protocols_qsar_scoring_grid_LigandPropertyScore_HH

#include <protocols/qsar/scoring_grid/ConstantScoreBase.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class LigandPropertyScore : public ConstantScoreBase
{
public:
	LigandPropertyScore();

	virtual ~LigandPropertyScore() {}

	virtual void parse_my_tag(utility::tag::TagCOP tag);
	/// @brief return the current score of an UltraLightResidue using the current grid
	virtual core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map);

	/// @brief return the current score of a residue using the current grid
	virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map);

	/// @brief get the type of the grid
	virtual std::string get_type();

	/// @brief Serialize the GridBase object into a json_spirit Value
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json spirit Value into a GridBase object
	virtual void deserialize(utility::json_spirit::mObject data);

	static std::string grid_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string parameter_tag_;

};

}
}
}

#endif
