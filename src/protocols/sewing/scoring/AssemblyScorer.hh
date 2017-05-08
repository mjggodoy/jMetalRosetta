// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AssemblyScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_AssemblyScorer_hh
#define INCLUDED_protocols_sewing_scoring_AssemblyScorer_hh

//Unit headers
#include <protocols/sewing/scoring/AssemblyScorer.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class AssemblyScorer : public utility::pointer::ReferenceCount {

public:

	///@brief default construct
	AssemblyScorer(){}

	virtual ~AssemblyScorer(){}

	virtual
	core::Real
	score(
		AssemblyCOP assembly
	) = 0;

private:

};

class AssemblyScoreFunction : public utility::pointer::ReferenceCount {

public:

	typedef std::map< std::string,  std::pair<core::Real, AssemblyScorerOP> > ScorerMap;

	///@brief return the weighted sum of all scorer scores
	core::Real
	score(
		AssemblyCOP assembly
	) const {
		core::Real score = 0.0;
		if ( assembly->total_residue() == 0 ) { return score; }
		ScorerMap::const_iterator it = scorers_.begin();
		ScorerMap::const_iterator it_end = scorers_.end();
		for ( ; it != it_end; ++it ) {
			score += it->second.first * it->second.second->score(assembly);
		}
		return score;
	}

	void
	add_scorer(
		std::string score_name,
		core::Real weight,
		AssemblyScorerOP scorer
	) {
		scorers_[score_name] = std::make_pair(weight, scorer);
	}

	utility::vector1< std::pair<std::string, core::Real> >
	get_all_scores(
		AssemblyCOP assembly
	) const {
		utility::vector1< std::pair<std::string, core::Real> > scores;
		ScorerMap::const_iterator it = scorers_.begin();
		ScorerMap::const_iterator it_end = scorers_.end();
		for ( ; it != it_end; ++it ) {
			core::Real score = it->second.first * it->second.second->score(assembly);
			scores.push_back( std::make_pair(it->first, score) );
		}
		return scores;
	}

private:

	//weights and their respective scores
	std::map< std::string,  std::pair<core::Real, AssemblyScorerOP> > scorers_;

};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
