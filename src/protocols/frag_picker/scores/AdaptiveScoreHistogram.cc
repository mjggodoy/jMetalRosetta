// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AdaptiveScoreHistogram.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// package headers
#include <basic/Tracer.hh>
#include <protocols/frag_picker/scores/AdaptiveScoreHistogram.hh>

#include <cmath>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using core::Real;
using core::Size;

static THREAD_LOCAL basic::Tracer trAdaptiveScoreHistogram(
	"fragment.picking.scores.AdaptiveScoreHistogram");

AdaptiveScoreHistogram::AdaptiveScoreHistogram(core::Real bin_size,core::Real initial_max_score) {
	bin_size_ = bin_size;
	is_up_to_date_ = true;
	core::Size new_size = (core::Size)(initial_max_score / bin_size_);
	data_.resize(new_size);
}

AdaptiveScoreHistogram::~AdaptiveScoreHistogram() {}

void AdaptiveScoreHistogram::insert(core::Real score) {

	core::Size bin_id = (core::Size)(score / bin_size_);
	if ( data_.size() <= bin_id ) {
		data_.resize(bin_id+1);
	}
	data_[bin_id+1]++;
	is_up_to_date_ = false;
}

core::Size AdaptiveScoreHistogram::sum() {

	core::Size sum = 0;
	for ( core::Size i=1; i<=data_.size(); i++ ) {
		sum += data_[i];
	}

	return sum;
}

core::Real AdaptiveScoreHistogram::p_value(core::Real score) {

	if ( cumulative_sums_.size() < data_.size() ) {
		cumulative_sums_.resize( data_.size() );
		is_up_to_date_ = false;
	}

	if ( !is_up_to_date_ ) {
		cumulative_sums_[1] = data_[1];
		for ( core::Size i=2; i<=data_.size(); i++ ) {
			cumulative_sums_[i] = cumulative_sums_[i-1] + data_[i];
		}
	}

	core::Size bin_id = (core::Size)(score / bin_size_);
	if ( bin_id >= cumulative_sums_.size() ) {
		return 0.0;
	} else {
		return log(cumulative_sums_[bin_id + 1] / (core::Real) cumulative_sums_[cumulative_sums_.size()]);
	}
}

}
} // frag_picker
} // protocols
