// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Job classes
/// @author Steven Lewis smlewi@gmail.com

// Unit headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

// Project headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputterObserver.hh>

#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/sys_util.hh>

// C++ headers
#include <string>
#include <ctime>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.Job" );

namespace protocols {
namespace jd2 {


////////////////////////////Job/////////////////////////////
Job::Job( InnerJobOP inner_job, core::Size nstruct_index )
: inner_job_(std::move(inner_job)),
	nstruct_index_(nstruct_index),
	status_prefix_( "" ),
	long_strings_(),
	string_string_pairs_(),
	string_real_pairs_(),
	completed_(false),
	can_be_deleted_(false),
	start_time_(0),
	end_time_(0),
	timestamp_("")
{
}

/// @brief Copy constructor
///
Job::Job( Job const &src )
: utility::pointer::ReferenceCount(),
	inner_job_(src.inner_job_->clone()),
	nstruct_index_(src.nstruct_index_),
	status_prefix_( src.status_prefix_ ),
	long_strings_(src.long_strings_),
	string_string_pairs_(src.string_string_pairs_),
	string_real_pairs_(src.string_real_pairs_),
	completed_(src.completed_),
	can_be_deleted_(src.can_be_deleted_),
	start_time_(src.start_time_),
	timestamp_(src.timestamp_)
{
	//TR.Trace << "Using Job (base class) for JobDistributor" << std::endl;
}

Job::~Job()= default;



/// @brief returns a copy of this object whose "output fields" are zeroed out.  Used by the JobDistributor in cases where the job fails and must be retried to prevent accumulation of Job state after a failure.  This implementation was chosen over a clear_all_output function to prevent mover A from deleting mover B's hard work!  You probably should not be trying to call this function.
JobOP Job::copy_without_output() const{
	return JobOP( new Job(inner_job_, nstruct_index_) );
}

/// @brief Return an owning pointer to a copy of this object.
///
JobOP Job::clone() const
{
	return JobOP( new Job( *this ) );
}

InnerJobCOP Job::inner_job() const { return inner_job_; }

std::string const & Job::input_tag() const {
	return inner_job_->input_tag();
}

InnerJobOP Job::inner_job_nonconst() { return inner_job_; }

core::Size Job::nstruct_index() const { return nstruct_index_; }

core::Size Job::nstruct_max() const {
	return inner_job_->nstruct_max();
}

//functions for loading output info into the job
/// @brief add an output string
void Job::add_string( std::string const & string_in ){
	long_strings_.push_back(string_in);
}

/// @brief adds output strings
void Job::add_strings( Strings const & strings )
{
	long_strings_.insert( long_strings_.end(), strings.begin(), strings.end() );
}

/// @brief add a string/string pair
void Job::add_string_string_pair( std::string const & string1, std::string const & string2 ){
	string_string_pairs_[ string1 ] = string2;
}

/// @brief add a string/real pair
void Job::add_string_real_pair( std::string const & string_in, core::Real const real_in ){
	string_real_pairs_[ string_in ] = real_in;
}

/// @brief Delete the output strings, string/string pairs, and string/real pairs.
///
void Job::clear_output() {
	long_strings_.clear();
	string_string_pairs_.clear();
	string_real_pairs_.clear();
	return;
}

/// @brief return a COP to the input pose
core::pose::PoseCOP Job::get_pose() const {
	// if pose is loaded into job-object this pose has precedence.
	if ( inner_job_->get_pose() ) {
		return inner_job_->get_pose();
	} else { //if not ask job-inputter
		protocols::jd2::JobDistributor* jd
			= protocols::jd2::JobDistributor::get_instance();
		core::pose::PoseOP aPose( new core::pose::Pose );
		// in the following call we use a copy of the job-output. the idea is that if the inputter has not stored the
		// pose in the Job-Object already he will not do it now... in principle one could also use a const_cast which has
		// basically the same effect -- the copy copies a pointer and an integer...
		if ( jd->job_inputter() ) jd->job_inputter()->pose_from_job( *aPose, copy_without_output() /*or const_cast(this) */ );
		return aPose;
	}
	utility_exit_with_message( "Programming error: you asked for Job::get_pose() but there is neither a job_inputter nor a pose loaded into the Job-Object ");
	return nullptr;
}


/// @brief in-place copy of input pose
void Job::get_pose( core::pose::Pose& pose ) const {
	// if pose is loaded into job-object this pose has precedence.
	if ( inner_job_->get_pose() ) {
		pose = *( inner_job_->get_pose() );
		return;
	} else { //if not ask Job-inputter
		protocols::jd2::JobDistributor* jd
			= protocols::jd2::JobDistributor::get_instance();
		if ( jd->job_inputter() ) {
			// in the following call we use a copy of the job-output. the idea is that if the inputter has not stored the
			// pose in the Job-Object already he will not do it now... in principle one could also use a const_cast which has
			// basically the same effect -- the copy copies a pointer and an integer...
			jd->job_inputter()->pose_from_job( pose, copy_without_output() );
		}
		return;
	}
	utility_exit_with_message( "Programming error: you asked for Job::get_pose() but there is neither a job_inputter nor a pose loaded into the job-object");
}


//functions for returning output info from the job.  You get iterators so that this interface can stay constant as the underlying implementation changes
Job::Strings::const_iterator Job::output_strings_begin() const
{ return long_strings_.begin(); }

Job::Strings::const_iterator Job::output_strings_end() const
{ return long_strings_.end(); }

Job::StringStringPairs::const_iterator Job::output_string_string_pairs_begin() const
{ return string_string_pairs_.begin(); }

Job::StringStringPairs::const_iterator Job::output_string_string_pairs_end() const
{ return string_string_pairs_.end(); }

Job::StringRealPairs::const_iterator Job::output_string_real_pairs_begin() const
{ return string_real_pairs_.begin(); }

Job::StringRealPairs::const_iterator Job::output_string_real_pairs_end() const
{ return string_real_pairs_.end(); }


bool Job::bad() const {
	return inner_job_->bad();
}

void Job::set_bad(bool value)  {
	inner_job_->set_bad( value );
	if ( inner_job_->bad() ) can_be_deleted_=true; //If this job is bad, it can be deleted to free up memory.  If not, it may or may not be deletable, so don't touch the value of can_be_deleted_.
	return;
}

core::Size Job::start_time() const {
	return start_time_;
}

core::Size Job::end_time() const {
	return end_time_;
}

core::Size Job::elapsed_time() const {
	if ( start_time_ == 0 ) {
		// Haven't started job yet.
		return 0;
	} else if ( end_time_ == 0 ) {
		// In the middle of job
		return time(nullptr) - start_time_;
	}
	// job completed
	runtime_assert( end_time_ >= start_time_ );
	return end_time_ - start_time_;
}

std::string Job::timestamp() const {
	return timestamp_;
}

void Job::start_timing() {
	start_time_ = time(nullptr);
	end_time_ = 0;
	timestamp_ = utility::timestamp_short();
}

void Job::end_timing() {
	debug_assert( start_time_ != 0 );
	end_time_ = time(nullptr);
}

void Job::add_output_observer( JobOutputterObserverAP an_observer ) {
	output_observers_.insert( an_observer );
}


void Job::remove_output_observer( JobOutputterObserverAP an_observer ) {
	output_observers_.erase( an_observer );
}

void Job::call_output_observers( core::pose::Pose const & pose )
{
	for ( auto const & output_observer : output_observers_ ) {
		JobOutputterObserverOP observer = output_observer.lock();
		if ( observer ) {
			observer->add_values_to_job( pose, *this );
		}
	}
}


JobCOP const JD2_BOGUS_JOB(  JobOP( new Job( InnerJobOP( new InnerJob("EMPTY_JOB_use_jd2", 0)), 0 ) ) );

bool
operator==(
	Job const & a,
	Job const & b
) {
	return
		*(a.inner_job_) == *(b.inner_job_) &&
		a.nstruct_index_ == b.nstruct_index_ &&
		a.status_prefix_ == b.status_prefix_ &&
		a.long_strings_ == b.long_strings_ &&
		a.string_string_pairs_ == b.string_string_pairs_ &&
		a.string_real_pairs_ == b.string_real_pairs_ &&
		a.completed_ == b.completed_ &&
		a.can_be_deleted_ == b.can_be_deleted_;
}


bool
operator!=(
	Job const & a,
	Job const & b
) {
	return !(a == b);
}


void
Job::show(
	std::ostream & out
) const {
	out
		<< "Inner Job:";
	if ( inner_job_ ) {
		out << std::endl << *inner_job_ << std::endl;
	} else {
		out << " NULL" << std::endl;
	}

	out
		<< "nstruct index: " << nstruct_index_ << std::endl
		<< "status_prefix: " << status_prefix_ << std::endl
		<< "long_strings:" << std::endl;
	for ( auto const & long_string : long_strings_ ) {
		out << long_string << std::endl;
	}

	out
		<< "completed: " << completed_ << std::endl
		<< "can be deleted: " << (can_be_deleted_ ? "true" : "false") << std::endl
		<< "String -> String Pairs:" << std::endl;
	for ( auto const & string_string_pair : string_string_pairs_ ) {
		out
			<< "\t" << string_string_pair.first << ": " << string_string_pair.second << std::endl;
	}

	out << "String -> Real Pairs:" << std::endl;
	for ( auto const & string_real_pair : string_real_pairs_ ) {
		out
			<< "\t" << string_real_pair.first << ": " << string_real_pair.second << std::endl;
	}
}

std::ostream &
operator<< (
	std::ostream & out,
	const Job & job
) {
	job.show(out);
	return out;
}


} // jd2
} // protocols
