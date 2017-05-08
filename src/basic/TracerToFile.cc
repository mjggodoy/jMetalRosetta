// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file TracerToFile.hh
/// @brief Class for a tracer that writes all output to a file.
/// @author Matt O'Meara (mattjomeara@gmail.com)


#include <basic/TracerToFile.hh>
#include <cassert>          // for assert
#include <iosfwd>           // for ofstream, string, fstream
#include <ostream>          // for basic_ostream, operator<<, endl
#include <string>           // for operator+
#include <utility/exit.hh>  // for utility_exit_with_message


namespace basic {

TracerToFile::TracerToFile( std::string const & file_name )
{
	file_.open( file_name.c_str(), std::fstream::out );
	if ( !file_.good() || !file_.is_open() ) {
		utility_exit_with_message( "cannot open tracer file "+file_name );
	}
}

TracerToFile::~TracerToFile()
{
	(*this) << std::endl;
	file_.close();
}


void TracerToFile::t_flush( std::string const & s )
{
	assert( file_.is_open() );
	file_ << s;
}

} // namepsace corenested classes c+++
