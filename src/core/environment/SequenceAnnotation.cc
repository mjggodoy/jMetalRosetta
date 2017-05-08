// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/moves/ClaimingMover.cc
/// @author Justin Porter

// Unit Headers
#include <core/environment/SequenceAnnotation.hh>

// Package headers

// Project headers
#include <utility/excn/Exceptions.hh>

#include <utility/string_util.hh>

// tracer
#include <basic/Tracer.hh>


// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "core.environment.SequenceAnnotation", basic::t_info );

namespace core {
namespace environment {

typedef utility::vector1< core::Size > vector1_size;

class LengthChecker {
public:
	LengthChecker( core::Size l ): l_( l ) {}
	void operator() ( core::Size i ) {
		if ( i > l_ || i < 1 ) {
			throw utility::excn::EXCN_RangeError( "Residue "+utility::to_string( i )+
				" is out of the range of the SequenceAnnotation object." );
		}
	}
private:
	core::Size l_;
};

SequenceAnnotation::SequenceAnnotation( core::Size length ):
	ReferenceCount(),
	length_( length ),
	pose_to_local_numbers_( length )
{
	utility::vector1< core::Size > base ( length );
	for ( core::Size i = 1; i <= length; ++i ) {
		base[i] = i;
	}
	_add_seq_label( "BASE", base );
}

void SequenceAnnotation::add_seq_label( std::string const& label, std::vector< core::Size > const& members ){
	_add_seq_label( label, vector1_size( members.begin(), members.end() ) );
}

void SequenceAnnotation::add_seq_label( std::string const& label, vector1_size const& members ){
	_add_seq_label( label, members );
}

void SequenceAnnotation::_add_seq_label( std::string const& label, vector1_size members ){
	std::sort( members.begin(), members.end() );

	//verify that all members are within the sequence range.
	std::for_each( members.begin(), members.end(), LengthChecker( length_ ) );

	if ( label_to_pose_numbers_.find( label ) != label_to_pose_numbers_.end() ) {
		if ( members != resolve_seq( label ) ) {
			std::ostringstream ss;
			ss << "The sequence label " << label << " already exists in the SequenceAnnotation object." << std::endl
				<< " Existing: " << resolve_seq( label ) << std::endl << "Conflicting: " << members;
			throw utility::excn::EXCN_KeyError( ss.str() );
		} else {
			//bail because nothing needs to be done--this exact label/member combination is already in there.
			return;
		}
	} else {
		label_to_pose_numbers_[label] = vector1_size( members.begin(), members.end() );
	}

	for ( Size i = 1; i <= members.size(); ++i ) {
		Size pose_number = members.at(i);
		pose_to_local_numbers_.at(pose_number).insert( std::make_pair( label, i ) );
	}
}

void SequenceAnnotation::add_jump_label( std::string const& label, core::Size i ) {
	if ( jump_label_to_number_.find(label) != jump_label_to_number_.end() &&
			jump_label_to_number_[ label ] != i ) {
		throw utility::excn::EXCN_KeyError( "The jump key '"+label+"' already exists in the SequenceAnnotation object." );
	}

	//jump_number_to_label_[ i ] = label;
	jump_label_to_number_[ label ] = i;
}

void SequenceAnnotation::rm_seq_label( std::string const& label ){
	if ( label_to_pose_numbers_.find( label ) == label_to_pose_numbers_.end() ) {
		throw utility::excn::EXCN_KeyError("The sequence key '"+label+"' does not exists in the SequenceAnnotation object" );
	}
	vector1_size& pose_numbers = label_to_pose_numbers_.find( label )->second;

	for ( Size pose_num = 1; pose_num <= pose_numbers.size(); ++pose_num ) {
		debug_assert( pose_num <= length_ );
		pose_to_local_numbers_[pose_num].erase( label );
	}

	label_to_pose_numbers_.erase( label );
}

void SequenceAnnotation::append_seq( std::string const& label ) {
	if ( label_to_pose_numbers_.find( label ) != label_to_pose_numbers_.end() ) {
		throw utility::excn::EXCN_KeyError("The sequence key '"+label+"' already exists in the SequenceAnnotation object" );
	}
	// TODO: Allow more than one residue to be appended as part of a label
	length_ += 1;
	label_to_pose_numbers_[ label ] = vector1_size( 1, length() );
	label_to_pose_numbers_[ "BASE" ].push_back( length() );
}

core::Size SequenceAnnotation::resolve_seq( LocalPosition const& local ) const {
	auto it = label_to_pose_numbers_.find( local.label() );

	if ( it == label_to_pose_numbers_.end() ) {
		std::ostringstream msg;
		msg << "SequenceAnnotation could not find local position " << local << "." << std::endl;
		throw utility::excn::EXCN_RangeError( msg.str() );
	} if ( local.position() < 1 ||
			local.position() > it->second.size() ) {
		std::ostringstream msg;
		msg << "SequenceAnnotation cannot resolve the local position " << local
			<< " because it is out of range." << std::endl;
		throw utility::excn::EXCN_RangeError( msg.str() );
	}

	return it->second.at( local.position() );
}

core::Size SequenceAnnotation::resolve_jump( std::string const& label ) const {

	auto l = jump_label_to_number_.find( label );

	if ( l == jump_label_to_number_.end() ) {
		std::ostringstream msg;
		msg << "SequenceAnnotation could not find jump with label "
			<< label << "." << std::endl;
		throw utility::excn::EXCN_RangeError( msg.str() );
	}

	return l->second;
}

bool SequenceAnnotation::has_seq_label( std::string const& label ) const {
	return label_to_pose_numbers_.find( label ) != label_to_pose_numbers_.end();
}

utility::vector1< core::Size > const& SequenceAnnotation::resolve_seq( std::string const& label ) const {
	if ( label_to_pose_numbers_.find( label ) == label_to_pose_numbers_.end() ) {
		throw utility::excn::EXCN_KeyError( "The jump key '"+label+
			"' already exists in the SequenceAnnotation object." );
	}
	return label_to_pose_numbers_.find( label )->second;
}

core::Size SequenceAnnotation::length( std::string const& label ) const {
	return label_to_pose_numbers_.find( label )->second.size();
}

} // environment
} // core
