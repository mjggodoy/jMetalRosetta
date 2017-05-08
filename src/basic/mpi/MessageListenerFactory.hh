// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MessageListenerFactor.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_basic_mpi_MessageListenerFactory_HH
#define INCLUDED_basic_mpi_MessageListenerFactory_HH

#include <basic/mpi/MessageListener.fwd.hh>

// C++ headers
#include <map>

// Utility headers
#include <utility/SingletonBase.hh>

namespace basic {
namespace mpi {

class MessageListenerFactory : public utility::SingletonBase< MessageListenerFactory >
{
public:
	friend class utility::SingletonBase< MessageListenerFactory >;

public:
	MessageListenerOP get_listener( listener_tags tag );

private:

	MessageListenerFactory();
	MessageListenerFactory(MessageListenerFactory const &);
	MessageListenerFactory const & operator = (MessageListenerFactory const &);

private:

	std::map< listener_tags, MessageListenerOP > listeners_;

};

} //namespace
} //namespace
#endif
