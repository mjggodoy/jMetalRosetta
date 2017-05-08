// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobResult.hh
/// @brief  The definition for class protocols::jd3::JobResult
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobResult_HH
#define INCLUDED_protocols_jd3_JobResult_HH

// Unit headers
#include <protocols/jd3/JobResult.fwd.hh>

// Package headers
#include <protocols/jd3/Job.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//C++ headers
#include <string>
#include <list>
#include <map>
#include <set>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @brief %JobResult class holds the output that's generated by a Job over the course
/// of its execution.  The %JobResult is handed by the JobQueen to the JobOutputWriter
/// objects, each of which have the opportunity to pull data out of the JobResult
/// class.
class JobResult : public utility::pointer::ReferenceCount
{
public:

	JobResult();
	~JobResult() override;

	virtual JobStatus status() const = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // JobResult

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_JobResult )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_JobResult_HH
