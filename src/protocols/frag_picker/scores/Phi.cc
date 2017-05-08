// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/Phi.cc
/// @brief  a base class for fragment scoring
/// @author David E Kim

#include <protocols/frag_picker/scores/Phi.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

//static const core::Real PHI_MIN_CONF = 0.5;

void Phi::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool Phi::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	return score( fragment, scores );

}

bool Phi::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {
	core::Real totalScore = 0;
	core::Size conf_positions = 0;
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		// skip low confidence positions
		core::Size qindex = i + f->get_first_index_in_query() - 1;
		VallResidueOP r = f->get_residue(i);
		//if (query_phi_prediction_conf_[qindex] < PHI_MIN_CONF) continue;
		// skip first residue in query and vall chunk
		if ( (i == 1 && ( qindex == 1 || f->get_first_index_in_vall() <= 1)) ||
				r->dssp_phi() == 360.0 || query_phi_prediction_[qindex] == 360.0 ) continue;
		// difference / 180
		totalScore += fabs( (query_phi_prediction_[qindex] - r->dssp_phi()) / 180.0 );
		conf_positions++;
	}
	totalScore /= (core::Real) conf_positions;
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}


} // scores
} // frag_picker
} // protocols


