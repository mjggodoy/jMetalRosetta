// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/deallocation/DeallocationMessage.hh
/// @brief  The definition for class protocols::jd3::DeallocationMessage
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_deallocation_DeallocationMessage_HH
#define INCLUDED_protocols_jd3_deallocation_DeallocationMessage_HH

// Unit headers
#include <protocols/jd3/deallocation/DeallocationMessage.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace deallocation {

enum deallocation_msg_type : short
{
	unassigned_deallocation_msg,
	input_pose_deallocation_msg,
	string_string_pair_msg
};

/// @brief %DeallocationMessage class holds the output that's generated by a Job over the course
/// of its execution.  The %DeallocationMessage is handed by the JobQueen to the JobOutputWriter
/// objects, each of which have the opportunity to pull data out of the DeallocationMessage
/// class.
class DeallocationMessage : public utility::pointer::ReferenceCount
{
public:

	DeallocationMessage();
	DeallocationMessage( deallocation_msg_type msg_type );
	virtual ~DeallocationMessage();

	deallocation_msg_type
	deallocation_type() const;

	void deallocation_type( deallocation_msg_type setting );
private:
	deallocation_msg_type type_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // DeallocationMessage

} // namespace deallocation
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_deallocation_DeallocationMessage )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_deallocation_DeallocationMessage_HH