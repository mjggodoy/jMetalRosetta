// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/AdditionalOutputWrapperCreator.hh
/// @brief This class will create instances of Mover AdditionalOutputWrapper for the MoverFactory
/// @author Luki Goldschmidt

#ifndef INCLUDED_protocols_rosetta_scripts_AdditionalOutputWrapperCreator_hh
#define INCLUDED_protocols_rosetta_scripts_AdditionalOutputWrapperCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace rosetta_scripts {

class AdditionalOutputWrapperCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static std::string mover_name();
};

}
}

#endif

