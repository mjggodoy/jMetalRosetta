// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_CanonicalSamplingMover_HH
#define INCLUDED_protocols_canonical_sampling_CanonicalSamplingMover_HH


#include <core/pose/Pose.hh>

#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>


#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {

class CanonicalSamplingMover: public moves::Mover{

public:
	static void register_options();

	CanonicalSamplingMover();

	CanonicalSamplingMover(core::scoring::ScoreFunctionOP sfxn,
		protocols::canonical_sampling::mc_convergence_checks::Pool_RMSD_OP ptr,
		int ntrial);

	void add_mover(protocols::moves::MoverOP m,core::Real weight);

	std::string get_ABGEO_string( core::pose::Pose & pose, protocols::loops::Loops & loop);

	void ntrials(int ntrial);

	void set_defaults_from_cmdline();

	void set_temp(core::Real temperature);

	core::Real get_temp() {return temperature_;};

	void set_interval_pose_dump(int p_interval);

	void set_interval_data_dump(int d_interval);

	void set_scorefunction(core::scoring::ScoreFunctionOP sfxn);

	void detailed_balance(bool truefalse);
	bool detailed_balance() const { return detailed_balance_; };

	core::Real transition_threshold() const {return transition_threshold_; };

	void use_MPI_sync_pools(bool truefalse);
	bool use_MPI_sync_pools() const {return MPI_synchronize_pools_; };

	void use_MPI_bcast( bool truefalse );
	bool use_MPI_bcast() const { return MPI_bcast_; };

	void use_hierarchical_clustering( bool truefalse );
	bool use_hierarchical_clustering() const { return use_hierarchical_clustering_; };

	void output_only_cluster_transitions(bool truefalse);

	void set_poolrmsd(protocols::canonical_sampling::mc_convergence_checks::Pool_RMSD_OP ptr);

	void apply(core::pose::Pose & pose) override;
	std::string get_name() const override;

private:

	core::Real periodic_range( core::Real a, core::Real x);

	void dump_xtc_format_decoy(
		std::ostream& os,
		core::pose::Pose const& pose,
		loops::Loops const& loop_to_dump
	);


	void dump_decoy_or_score(
		std::ostream& os,
		core::pose::Pose const& pose,
		core::Size i_trial,
		std::string const& jobname,
		loops::Loops const& loop_to_dump,
		bool score_only = false
	);

	void setup_constraints(core::pose::Pose & pose);

	protocols::moves::MonteCarloOP mc_;
	core::scoring::ScoreFunctionOP sfxn_;
	moves::RandomMoverOP randmove_;
	protocols::canonical_sampling::mc_convergence_checks::Pool_RMSD_OP pool_rms_;

	core::Size interval_posedump_;
	core::Size interval_transitiondump_;
	core::Size ntrials_;
	bool detailed_balance_;
	bool MPI_synchronize_pools_;
	bool MPI_bcast_;
	bool use_hierarchical_clustering_;
	bool save_loops_only_;
	bool dump_loops_only_;
	bool output_only_cluster_transition_;
	core::Real transition_threshold_;
	core::Real temperature_;

	static bool options_registered_;
	bool ramp_temperature_;
	bool boinc_mode_;

};

}
}

#endif //  INCLUDED_protocols_canonical_sampling_CanonicalSamplingMover_HH
