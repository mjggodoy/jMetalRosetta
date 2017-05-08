// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Parser.fwd.hh
/// @brief  header file for Parser class, part of August 2008 job distributor
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_Parser_fwd_hh
#define INCLUDED_protocols_jd2_Parser_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class Parser;
typedef utility::pointer::shared_ptr< Parser > ParserOP;
typedef utility::pointer::shared_ptr< Parser const > ParserCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_Parser_FWD_HH
