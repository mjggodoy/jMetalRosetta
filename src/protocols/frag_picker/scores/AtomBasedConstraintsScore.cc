// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AtomBasedConstraintsScore.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>

#include <basic/Tracer.hh>
#include <core/chemical/ResidueType.hh>
#include <numeric/xyzVector.hh>

#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <protocols/frag_picker/scores/AtomBasedConstraintsScore.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

static THREAD_LOCAL basic::Tracer trAtomBasedConstraintsScore(
	"fragment.picking.scores.AtomBasedConstraintsScore");

/// @param priority - the priority for this scoring method. The lower the priority, the later the score will be evaluated
/// Because a fragment may be discarded when a score is too low, the most accurate and meaningful scores should have the highest priority
/// @param lowest_acceptable_value - a fragment for which this score is below a certain threshold will be discarded
/// @param query_size - the number of residues in the query sequence
/// @param constrainable_atoms - a vector of strings providing names of constrained atoms.
/// On every do_cahing() event these and only these atoms will be cached from a chunk's pose
/// @param score_name - name assigned to this scoring term; this string must show up in scores config file if the score is to be evaluated during picking
AtomBasedConstraintsScore::AtomBasedConstraintsScore(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, core::Size query_size, utility::vector1<
	std::string> constrainable_atoms, std::string score_name) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, score_name) {

	query_size_ = query_size;
	for ( core::Size i = 1; i < constrainable_atoms.size(); ++i ) {
		constrainable_atoms_.insert(std::pair<std::string, core::Size>(
			constrainable_atoms[i], i));
	}
}

/// @param priority - the priority for this scoring method. The lower the priority, the later the score will be evaluated
/// Because a fragment may be discarded when a score is too low, the most accurate and meaningful scores should have the highest priority
/// @param lowest_acceptable_value - a fragment for which this score is below a certain threshold will be discarded
/// @param query_size - the number of residues in the query sequence
/// @param score_name - name assigned to this scoring term; this string must show up in scores config file if the score is to be evaluated during picking
AtomBasedConstraintsScore::AtomBasedConstraintsScore(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, core::Size query_size, std::string score_name) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, score_name) {

	query_size_ = query_size;
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("N", 1));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("CA", 2));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("C", 3));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("O", 4));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("CB", 5));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("H", 6));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("HA", 7));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("HB", 8));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("1HB", 9));
	constrainable_atoms_.insert(std::pair<std::string, core::Size>("2HB", 10));
}

std::string AtomBasedConstraintsScore::get_constrained_atom_name(core::Size atom_id) {

	std::map<std::string, core::Size>::iterator it, end ;
	for ( it=constrainable_atoms_.begin(), end = constrainable_atoms_.end(); it != end; ++it ) {
		if ( (it)->second == atom_id ) {
			return (it)->first;
		}
	}

	return 0;
}

void AtomBasedConstraintsScore::do_caching(VallChunkOP chunk) {

	trAtomBasedConstraintsScore.Debug << "caching backbone atoms for "
		<< chunk->get_pdb_id() << " of size " << chunk->size() << std::endl;
	core::pose::PoseOP pose = chunk->get_pose();

	//  pose->dump_pdb("dump-"+chunk->get_pdb_id()+".pdb");

	numeric::xyzVector<core::Real> empty_one;
	for ( core::Size i = 1; i <= chunk->size(); ++i ) {
		utility::vector1<bool> flag_row(constrainable_atoms_.size());
		utility::vector1<numeric::xyzVector<core::Real> > row(
			constrainable_atoms_.size());
		for ( core::Size j = 1; j < constrainable_atoms_.size(); ++j ) {
			flag_row[j] = false;
			row[j] = empty_one;
		}
		chunk_atoms_xyz_.push_back(row);
		atom_flags_.push_back(flag_row);

		std::map<std::string, core::Size>::iterator it;
		core::chemical::ResidueType const & ith_res_type = pose->residue_type(i);
		for ( it = constrainable_atoms_.begin(); it
				!= constrainable_atoms_.end(); ++it ) {
			if ( ! ith_res_type.has( it->first ) ) {
				continue;
			}
			core::id::NamedAtomID idAtom(it->first, i);
			numeric::xyzVector<core::Real> xyz = pose->xyz(idAtom);
			chunk_atoms_xyz_[i][it->second] = xyz;
			atom_flags_[i][it->second] = true;
		}
	}
}

void AtomBasedConstraintsScore::clean_up() {
	chunk_atoms_xyz_.erase(chunk_atoms_xyz_.begin(), chunk_atoms_xyz_.end());
	atom_flags_.erase(atom_flags_.begin(), atom_flags_.end());
}

}
} // frag_picker
} // protocols
