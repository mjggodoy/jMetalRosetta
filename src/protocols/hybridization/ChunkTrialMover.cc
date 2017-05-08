// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Wrapper for InsertChunkMover. It can take a random template and steal coordinates of all chunks or a random one
/// @details
/// @author Yifan Song

#include <protocols/hybridization/ChunkTrialMover.hh>
#include <protocols/hybridization/util.hh>

#include <core/pose/PDBInfo.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.hybridization.ChunkTrialMover" );

namespace protocols {
namespace hybridization {

using namespace core;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::loops;
using namespace numeric::model_quality;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;

ChunkTrialMover::ChunkTrialMover(
	utility::vector1 < core::pose::PoseCOP > const & template_poses,
	utility::vector1 < protocols::loops::Loops > const & template_chunks,
	bool random_template,
	AlignOption align_option,
	utility::vector1<bool> sampling_chunk_in ) :
	template_poses_(template_poses),
	template_chunks_(template_chunks),
	random_template_(random_template),
	align_option_(align_option),
	align_chunk_(),
	max_registry_shift_global_(0)
{
	moves::Mover::type( "ChunkTrialMover" );

	Size count = 0;
	for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
		if ( template_chunks_[i_template].size() != 0 ) ++count;

	}
	if ( count == 0 ) {
		has_valid_moves_ = false;
		return;
	}
	has_valid_moves_ = true;

	// sequence mapping
	sequence_alignments_.clear();
	for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
		std::map <core::Size, core::Size> sequence_alignment;
		sequence_alignment.clear();
		get_alignment_from_template(template_poses_[i_template], sequence_alignment);
		sequence_alignments_.push_back(sequence_alignment);
	}

	// template coverage
	highest_tmpl_resnum_ = 0;
	for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
		for ( core::Size ires=1; ires<=template_poses_[i_template]->size(); ++ires ) {
			if ( template_poses_[i_template]->pdb_info()->number(ires) > static_cast<int>(highest_tmpl_resnum_) ) {
				highest_tmpl_resnum_ = template_poses_[i_template]->pdb_info()->number(ires);
			}
		}
	}
	residue_covered_by_template_.resize(highest_tmpl_resnum_, 0);
	for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
		for ( core::Size ires=1; ires<=template_poses_[i_template]->size(); ++ires ) {
			core::Size ires_template = template_poses_[i_template]->pdb_info()->number(ires);
			residue_covered_by_template_[ires_template] += 1;
		}
	}

	// allowed positions
	sampling_chunk_ = sampling_chunk_in;
}


void
ChunkTrialMover::get_alignment_from_template(
	core::pose::PoseCOP template_pose,
	std::map <core::Size, core::Size> & seqpos_alignment
) {
	// specific to this case, alignment comes from residue number
	for ( core::Size ires=1; ires<=template_pose->size(); ++ires ) {
		TR.Debug << "Sequence aln: " << template_pose->pdb_info()->number(ires) << " " << ires << std::endl;
		seqpos_alignment[template_pose->pdb_info()->number(ires)] = ires;
	}
}

void ChunkTrialMover::pick_random_template()
{
	assert(template_poses_.size() != 0);
	set_template(0);
	int ntrials=100;
	while ( !template_number() && --ntrials>0 ) {
		set_template( numeric::random::rg().random_range(1, template_poses_.size()) );
		if ( template_chunks_[template_number()].size() == 0 || ignore_template_indices_.count(template_number()) ) set_template(0);
	}
	if ( ntrials == 0 ) {
		utility_exit_with_message( "Fatal error in ChunkTrialMover::pick_random_template()");
	}
}

void ChunkTrialMover::set_template(core::Size const template_number) {
	template_number_ = template_number;
}

core::Size ChunkTrialMover::template_number() {
	return template_number_;
}

void ChunkTrialMover::set_max_registry_shift( core::Size max_registry_shift_in ) {
	max_registry_shift_global_ = max_registry_shift_in;
}


void
ChunkTrialMover::pick_random_chunk(core::pose::Pose & pose) {
	int ntrials=500;

	bool chosen_good_jump=false;

	while ( !chosen_good_jump && ntrials>0 ) {
		--ntrials;
		jump_number_ = numeric::random::rg().random_range(1, pose.num_jump());
		core::Size jump_residue_pose = pose.fold_tree().downstream_jump_residue(jump_number_);

		if ( pose.residue(jump_residue_pose).aa() == core::chemical::aa_vrt ) {
			continue;
		}

		// make sure downstream residues are in templates (avoid symmetry problems)
		std::list < core::Size > downstream_residues = downstream_residues_from_jump(pose, jump_number_);
		for ( core::Size & downstream_residue : downstream_residues ) {
			if ( downstream_residue < highest_tmpl_resnum_ ) {
				chosen_good_jump = true;
				break;
			}
		}
		if ( ! chosen_good_jump ) continue;

		// check if we're allowed to insert here
		chosen_good_jump = false;
		if ( sampling_chunk_.size() == 0 ) {
			chosen_good_jump = true;
		} else {
			for ( core::Size & downstream_residue : downstream_residues ) {
				if ( downstream_residue > sampling_chunk_.size() ) continue;
				if ( sampling_chunk_[downstream_residue]==true ) {
					chosen_good_jump=true;
					break;
				}
			}
		}
	}

	if ( ntrials == 0 ) {
		jump_number_ = pose.num_jump() + 1;
	}
}

Size ChunkTrialMover::trial_counter(Size ires) {
	return align_chunk_.trial_counter(ires);
}

void
ChunkTrialMover::apply(core::pose::Pose & pose) {
	// pick a random template
	if ( random_template_ ) {
		pick_random_template();
	}
	if ( ignore_template_indices_.count(template_number()) ) return;

	//TR << "templ number: " << template_number() << std::endl;
	align_chunk_.set_template(
		template_poses_[template_number()],
		template_number(),
		sequence_alignments_[template_number()]
	);

	// random chunk or loop all chunks
	if ( align_option_ == random_chunk ) {
		// pick a random jump
		pick_random_chunk(pose);
		if ( jump_number_ > pose.num_jump() ) {
			TR << "Warning! Fail to find a template chunk for sampling." << std::endl;
			return;
		}
		align_chunk_.set_aligned_chunk(pose, jump_number_, false);

		// apply alignment
		int registry_shift = numeric::random::rg().random_range(-max_registry_shift_global_, max_registry_shift_global_);
		align_chunk_.set_registry_shift(registry_shift);
		align_chunk_.apply(pose);
	} else {
		// loop over all jumps (we're initializing)
		for ( core::Size jump_number=1; jump_number<=pose.num_jump(); ++jump_number ) {
			// make sure that the downstream residues of this jump is within the template residues
			bool is_jump_affect_moveable_residue=false;
			std::list < core::Size > downstream_residues = downstream_residues_from_jump(pose, jump_number);
			for ( core::Size & downstream_residue : downstream_residues ) {
				if ( downstream_residue <= residue_covered_by_template_.size() ) {
					if ( residue_covered_by_template_[downstream_residue] > 0 ) {
						is_jump_affect_moveable_residue=true;
						break;
					}
				}
			}
			if ( is_jump_affect_moveable_residue ) {
				align_chunk_.set_aligned_chunk(pose, jump_number, true);
				align_chunk_.set_registry_shift(0);
				align_chunk_.apply(pose);
				if ( !align_chunk_.success() ) {
					TR.Debug << "Warning! This chunk might not be aligned, jump number " << jump_number << std::endl;
				}
			}
		}
	}
}

std::string ChunkTrialMover::get_name() const
{
	return "ChunkTrialMover";
}

} // hybridization
} // protocols

