// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CSTalosIO.hh
/// @brief Class that reads and writes chemical shifts in TALOS format
/// @author Nikolas Sgourakis sgourn@u.w.edu

#ifndef INCLUDED_protocols_frag_picker_JCouplingIO_hh
#define INCLUDED_protocols_frag_picker_JCouplingIO_hh

// utility headers
#include <core/types.hh>

#include <string>

// boost headers

#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {

class JCouplingIO {
public:

	JCouplingIO() {}

	JCouplingIO(std::string file_name) {
		read(file_name);
	}

	void read(std::string const&);

	std::pair< core::Real, core::Real > get_data( core::Size const res_num, bool & has_data );

	utility::vector1< core::Real > get_parameters();

	core::Size get_length() {
		return sequence_length_;
	}

private:
	utility::vector1< utility::vector1< core::Real > > data_;
	core::Real A_, B_, C_, THETA_;
	core::Size sequence_length_;

};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_JCouplingIO_HH */

