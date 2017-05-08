// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/residue_selector/AntibodyRegionSelector.hh
/// @brief  Creator for AntibodyRegionSelector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_residue_selector_AntibodyRegionSelectorCreator_HH
#define INCLUDED_protocols_antibody_residue_selector_AntibodyRegionSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace antibody {
namespace residue_selector {


class AntibodyRegionSelectorCreator : public core::select::residue_selector::ResidueSelectorCreator {
public:
	virtual core::select::residue_selector::ResidueSelectorOP create_residue_selector() const;

	virtual std::string keyname() const;

	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};


} //protocols
} //antibody
} //residue_selector


#endif //INCLUDED_protocols/antibody/residue_selector_AntibodyRegionSelector.hh

