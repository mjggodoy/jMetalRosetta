// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_trajectory_DbTrajectoryWriter_hh
#define INCLUDED_protocols_trajectory_DbTrajectoryWriter_hh

// Unit Headers
#include <protocols/trajectory/DbTrajectoryWriter.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/noncopyable.hpp>

// C++ headers
#include <string>

namespace protocols {
namespace trajectory {

class DbTrajectoryWriter : private boost::noncopyable {

public:

	/// @brief Constructor.  The given pose will be recorded as the first frame
	/// in the trajectory.  Some performance settings can optionally be given.
	DbTrajectoryWriter(
		core::Size job_id, core::pose::Pose const & pose, core::Size frequency=1, core::Size cache_limit=200);

	/// @brief Specify that a pose should be saved once every N iterations.
	void set_frequency(core::Size setting);

	/// @brief Specify the number of poses that can be cached.
	void set_cache_limit(core::Size setting);

	/// @brief Add the given pose to the trajectory being saved.
	void update(core::pose::Pose const & pose);

	/// @brief Make sure any cached poses have been saved.
	void finalize() const;

private:

	/// @brief Generate the table schemas and write them to the database.
	void write_schema_to_db() const;

	/// @brief Write any cached poses into the database, then clear the cache.
	void write_cache_to_db() const;

private:

	core::Size job_id_;
	core::Size iteration_;
	core::Size frequency_;
	core::Size cache_limit_;

	struct Frame { core::Size iteration; core::pose::Pose pose; };
	mutable utility::vector1<Frame> frame_cache_;

};

} // trajectory namespace
} // protocols namespace

#endif
