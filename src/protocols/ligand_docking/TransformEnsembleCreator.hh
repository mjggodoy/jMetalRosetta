// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file protocols/ligand_docking/TransformCreator.hh
///@brief This class will create instances of Mover Transform for the MoverFactory
///@author Sam DeLuca

#ifndef INCLUDED_protocols_ligand_docking_TransformCreatorEnsemble_hh
#define INCLUDED_protocols_ligand_docking_TransformCreatorEnsemble_hh

#include <protocols/moves/Mover.fwd.hh>

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace ligand_docking {

class TransformEnsembleCreator : public protocols::moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}

#endif

