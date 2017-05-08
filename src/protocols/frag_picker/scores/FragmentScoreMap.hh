// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentScoreMap.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentScoreMap_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentScoreMap_hh

// package headers
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief holds all small scores (score components) for a given fragment
/// @details each scoring method puts its own result into the map. The total score
/// is a dot product of the vector from this map and a vector of weights
class FragmentScoreMap: public utility::pointer::ReferenceCount {

public:

	/// @brief creates a new map for a given number of components
	/// @details usually the new map should be created by FragmentScoreManager::create_empty_map()
	inline FragmentScoreMap(core::Size n_componens) {
		small_scores_.resize(n_componens);
		was_modified_ = false;
		recent_total_ = 0.0;
		recent_quota_score_ = 999.999;
	}

	/// @brief creates a new deep copy of this map
	FragmentScoreMapOP clone() const {
		FragmentScoreMapOP f( new FragmentScoreMap(small_scores_.size()) );
		for ( core::Size i = 1; i <= small_scores_.size(); i++ ) {
			f->small_scores_[i] = small_scores_[i];
		}
		f->was_modified_ = was_modified_;
		f->recent_total_ = recent_total_;

		return f;
	}

	inline core::Size size() const { return small_scores_.size(); }

	inline core::Real at(core::Size score_index) const { return small_scores_[score_index]; };

	/// @brief returns the vector of score components
	inline utility::vector1<core::Real>& get_score_components() const {
		return const_cast<utility::vector1<core::Real> &> (small_scores_);
	}

	/// @brief sets a new score value for a given component
	inline void set_score_component(core::Real score_value, core::Size component_id) {
		small_scores_[component_id] = score_value;
		was_modified_ = true;
	}

	/// @brief returns the total score that has been evaluated recently
	/// @details the method does not compute anything. You must check if the map has been
	/// modified after the last total score calculation. If so, you must recompute the total again
	inline core::Real get_most_recent_total_score() {
		return recent_total_;
	}

	/// @brief If the map has been modified, the total score it holds is not up-to-date and must be recalculated
	inline bool was_modified() {
		return was_modified_;
	}

	inline core::Real get_quota_score() { return recent_quota_score_; }

	inline void set_quota_score(core::Real quota_score) { recent_quota_score_ = quota_score; }

private:
	utility::vector1<core::Real> small_scores_;
	core::Real recent_total_;
	core::Real recent_quota_score_;
	bool was_modified_;
	friend class FragmentScoreManager;
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_FragmentScoreMap_HH */
