// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/GenericMonteCarloMover.hh
/// @brief perform a given mover and sample structures by MonteCarlo
/// @details The "score" evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP you gave.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_GenericMonteCarloMover_hh
#define INCLUDED_protocols_simple_moves_GenericMonteCarloMover_hh

// Unit header
#include <protocols/simple_moves/GenericMonteCarloMover.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/moves/MonteCarloStatus.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered_map.hpp>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/MoverApplyingMover.hh>

namespace protocols {
namespace simple_moves {

/// @brief Trigger API definition
typedef boost::function<bool(core::Size,
	core::Size,
	const core::pose::Pose&,
	core::scoring::ScoreFunctionOP)> GenericMonteCarloMoverTrigger;


class GenericMonteCarloMover : public protocols::moves::MoverApplyingMover {
	typedef protocols::moves::MoverApplyingMover Super;
public:
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef protocols::rosetta_scripts::ParsedProtocolOP ParsedProtocolOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::Movers_map Movers_map;

	/// @brief default constructor
	GenericMonteCarloMover();

	/// @brief value constructor without score function
	GenericMonteCarloMover(
		Size const maxtrials,
		Size const max_accepted_trials,
		Size const task_scaling,
		MoverOP const & mover,
		Real const temperature = 0.0,
		String sample_type = "low",
		bool const drift = true
	);

	/// @brief value constructor with score function
	// Undefined, commenting out to fix PyRosetta build
	/* GenericMonteCarloMover(
	Size const maxtrials,
	MoverOP const & mover,
	ScoreFunctionOP const & sfxn,
	Real const temperature = 0.0,
	String const sample_type = "low",
	bool const drift = true
	); */


	/// @brief value constructor with task operation via TaskFactory
	GenericMonteCarloMover(
		Size const maxtrials,
		Size const max_accepted_trials,
		Size const task_scaling,
		MoverOP const & mover,
		TaskFactoryOP factory_in,
		Real const temperature = 0.0,
		String sample_type = "low",
		bool const drift = true
	);


	/// @brief destructor
	~GenericMonteCarloMover() override;

	/// @brief create copy constructor
	MoverOP clone() const override;

	/// @brief create this type of objectt
	MoverOP fresh_instance() const override;

	/// @brief initialize object used in constructor
	void initialize();

	/// @brief apply GenericMonteCarloMover (Mover)
	void apply( Pose & pose ) override;

	// XRW TEMP  String get_name() const override;

	/// @brief reset MC iterations, with pose used for the last & best structures
	void reset( Pose & pose );

	/// @brief return the simulation state to the lowest energy structure we've seen
	void recover_low( Pose & pose );

	/// @brief core of MC
	bool boltzmann( Pose & pose );

	/// @brief core of MC -- evaulates a pose based on the scores/filters + temperatures. random_num is a vector of random numbers between 0 and 1 with size equal to the number of MC criteria
	virtual bool boltzmann( Pose & pose, utility::vector1< core::Real > const & random_nums );

	/// @brief Does what the mover needs to do when a pose is accepted, given a pose and scores
	void
	accept( Pose & pose, utility::vector1< core::Real > const & provisional_scores, protocols::moves::MCA const mca_status );

	Size num_designable( Pose & pose, PackerTaskOP & task);

public: // accessor

	/// @brief return the last accepted pose
	PoseOP last_accepted_pose() const;

	/// @brief return the last accepted score
	Real last_accepted_score() const;

	/// @brief return the lowest score pose
	PoseOP lowest_score_pose() const;

	/// @brief return the lowest score
	Real lowest_score() const;

	/// @brief return the lowest score
	Real current_score() const;

	/// @brief return mc_accepted
	protocols::moves::MCA mc_accpeted() const;

	/// @brief Return the score function in use
	ScoreFunctionOP score_function() const {
		return scorefxn_;
	}

	/// @brief Adds a new trigger, returning its id.
	///
	/// Example:
	/// #include <boost/bind.hpp>
	/// #include <boost/function.hpp>
	///
	/// bool no_op(core::Size stage,
	///            core::Size num_stages,
	///            core::Size cycle,
	///            core::Size num_cycles,
	///            const core::pose::Pose&,
	///            core::scoring::ScoreFunctionOP) {}
	///
	/// Trigger callback = boost::bind(&no_op, STAGE, NUM_STAGES, _1, _2, _3, _4);
	/// Size trigger_id = add_trigger(callback);
	///
	/// The current stage and number of stages must be bound at creation time.
	/// This information provides the triggers with context about the current
	/// progress of the simulation as a whole.
	///
	/// If the trigger returns true, rescoring occurs.
	Size add_trigger(const GenericMonteCarloMoverTrigger& trigger);

	/// @brief Returns the number of triggers
	Size num_triggers() const;

	/// @brief Returns maximum number of trials
	Size maxtrials() const { return maxtrials_; }

	/// @brief Returns maximum number of accepted trials
	Size max_accepted_trials() const { return max_accepted_trials_; }

	/// @brief Returns the task scaling value
	Size task_scaling() const { return task_scaling_; }

public: // mutators
	/// @brief Removes the trigger with the specified id
	void remove_trigger(Size trigger_id);

	/// @brief set max trials of MC trials
	void set_maxtrials( Size const ntrial );

	/// @brief set max accepted trials of MC trials
	void set_max_accepted_trials( Size const n_max_accepted_trial );

	/// @brief set task multiplier to calculate trials from task
	void set_task_scaling( Size const scaling );

	/// @brief set mover
	void set_mover( MoverOP mover ) override;

	/// @brief Pose is evaluated by ScoreFunctionOP during MC trials
	void set_scorefxn( ScoreFunctionOP const & sfxn );

	/// @brief set temperature
	void set_temperature( Real const temp );

	/// @brief set sample type, max or min
	/// when sample_type == max, sample pose which have higher score
	/// when sample_type == min, sample pose which have lower score
	void set_sampletype( String const & type );

	/// @brief if drift=false, the pose is set back to the initial pose at each MC trial
	/// Of course, this is not MC sampling.
	void set_drift( bool const drift );

	/// @brief if preapply=true, auto-accept the first application of the submover,
	/// ignoring boltzman criteria.
	void set_preapply( bool const preapply=false );

	/// @brief if recover_low=true, after apply() the structure
	/// is the lowest energy structure, rather than the last accepted structure.
	void set_recover_low( bool const recover_low );

	/// @brief if boltz_rank=true, rank structures by the temperature-weighted
	/// sum of scores, rather than a single filter
	/// @details The score used here is the effective combined energy function
	/// that the Monte Carlo sampler is sampling over.
	void set_boltz_rank( bool const boltz_rank );

	/// @brief show scores of last_accepted_score and "best_score" ( = flip_sign_ * lowest_score )
	void show_scores( std::ostream & out ) const;

	/// @brief show counters of ntrial and acceptance ratio
	void show_counters( std::ostream & out ) const;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const &
	) override;

	//static void modify_attlist_and_ssl_for_children( utility::tag::AttributeList & attlist, utility::tag::XMLSchemaSimpleSubelementList & ssl, utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		Movers_map const &
	);

	void add_filter( FilterOP filter, bool const adaptive, Real const temp, String const & sample_type, bool rank_by=false);

	void stopping_condition( protocols::filters::FilterOP filter );
	protocols::filters::FilterOP stopping_condition() const;
	void adaptive_movers( bool const a ){ adaptive_movers_ = a; }
	bool adaptive_movers() const{ return adaptive_movers_; }
	void adaptation_period( core::Size const a ) { adaptation_period_ = a; }
	core::Size adaptation_period() const{ return adaptation_period_; }
	std::string saved_accept_file_name() const;
	void saved_accept_file_name( std::string const );
	std::string saved_trial_number_file() const;
	void saved_trial_number_file( std::string const );
	core::Size load_trial_number_from_checkpoint( core::pose::Pose & );
	void save_trial_number_to_checkpoint( core::Size const i ) const;
	void reset_baselines( bool const r ){ reset_baselines_ = r; };
	bool reset_baselines() const{ return reset_baselines_; }

	void set_keep_filters( bool const k ){ keep_filters_ = k; };
	bool keep_filters() const{ return keep_filters_; }

	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP task_factory() const{ return factory_; }

	Size trial_counter() const { return trial_counter_; }
	void trial_counter( core::Size const i ){ trial_counter_ = i; };

	//void set_trial_counter
	Size accept_counter() const { return accept_counter_; }
	utility::vector1< Real > temperatures() const { return temperatures_; }
	void temperatures( utility::vector1< Real > const & temps ) { temperatures_ = temps; }
	bool recover_low() const { return recover_low_; }
	MoverOP mover() const override { return mover_; }
	utility::vector1< FilterOP > const & filters() const { return filters_; }
	utility::vector1< bool > const & adaptive() const { return adaptive_; }


	utility::vector1< core::Size > const & num_rejections() const { return num_rejections_; }

	ScoreFunctionOP scorefxn() const { return scorefxn_; }
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< bool > > mover_stopping_condition() const { return mover_stopping_condition_; }
	bool preapply() const { return preapply_; }
	bool drift() const { return drift_; }
	bool boltz_rank() const { return boltz_rank_; }
	utility::vector1< String > const & sample_types() const { return sample_types_; }
	utility::vector1< core::Real > const & last_accepted_scores() const { return last_accepted_scores_; }
	void last_accepted_scores( utility::vector1< core::Real > const & scores ) { last_accepted_scores_ = scores; }
	utility::vector1< core::Real > const & last_tested_scores() const { return last_tested_scores_; }
	utility::vector1< core::Real > const & lowest_scores() const { return lowest_scores_; }
	void lowest_scores( utility::vector1< core::Real > const & scores ) { lowest_scores_ = scores; }
	void lowest_score( core::Real const score ) { lowest_score_ = score; }
	void last_accepted_score( core::Real const score ) { last_accepted_score_ = score; }
	void lowest_score_pose( core::pose::PoseOP pose );
	void last_accepted_pose( core::pose::PoseOP pose );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	define_composition_schema( utility::tag::XMLSchemaDefinition & );

protected:
	/// @brief evalute pose by ScoreFunctionOP or FilterOP
	Real scoring( Pose & pose );

	/// @brief Executes all triggers. The order of trigger execution is undefined.
	/// Do not assume, depend, or in any way rely on a particular ordering.
	void fire_all_triggers(
		Size cycle,
		Size num_cycles,
		const Pose& pose,
		ScoreFunctionOP scoring);

	/// @brief generate a set of random numbers to accompany the filters or scorefxn
	utility::vector1< core::Real >
	generate_random() const;

	std::string progress_file() const{ return progress_file_; }
	void progress_file( std::string const & s ){ progress_file_ = s; }

	/// @brief Sets mc_accpeted
	void set_mc_accepted( moves::MCA const m ){ mc_accepted_ = m; }

	int trial_counter_ = 0;

	/// @brief Count the number of rejections each filter resulted in.
	utility::vector1< core::Size > num_rejections_;

private:
	/// @brief max number of MC trials
	Size maxtrials_ = 10;

	/// @brief max number of accepted MC trials
	Size max_accepted_trials_ = 0;

	/// @brief number of designable positions
	Size number_designable_ = 0;

	/// @brief multiply the number_designable_ by task_scaling to get the number of trials
	Size task_scaling_ = 5;

	/// @brief mover
	MoverOP mover_;

	/// @brief task
	PackerTaskOP task_;

	/// @brief task factory
	TaskFactoryOP factory_;

	/// @brief Pose is evaluated by FilterOP which can do report_sm() during MC trials
	utility::vector1< FilterOP > filters_;
	utility::vector1< bool > adaptive_; // deflt true; are the filters adaptive or constant

	/// @brief acceptance criterion temperature, per temperature
	utility::vector1< Real > temperatures_;
	/// @brief low/high, dflt low
	utility::vector1< String > sample_types_;
	utility::vector1< Real > last_accepted_scores_;
	/// @brief best filter scores
	utility::vector1< Real > lowest_scores_;
	utility::vector1< Real > last_tested_scores_;

	/// @brief Pose is evaluated by ScoreFunctionOP during MC trials
	ScoreFunctionOP scorefxn_;

	/// @brief acceptance criterion temperature
	Real temperature_ = 0.0;// temperature for non-filters

	/// @brief set sample type, max or min
	/// when sample_type == "high", sample pose which have higher score
	/// when sample_type == "low", sample pose which have lower score
	String sample_type_ = "low";

	/// @brief if drift=false, the pose is set back to the initial pose at each MC trial
	/// Of course, this is not MC sampling.
	bool drift_ = true;

	/// @brief Should we apply (and accept) the first application of the mover
	/// regardless of boltzman criteria? (Defaults true for historical reasons.)
	bool preapply_ = true;

	/// @brief At the end of application, is the structure
	/// the last accepted structure (recover_low_==false)
	/// or is it the lowest energy structure (recover_low_=true)
	bool recover_low_ = true;

	/// @brief By which filter (by index) are poses ranked by for non-Boltzman purposes
	/// (e.g. recover_low())
	Size rank_by_filter_ = 1;

	/// @brief If boltz_rank_=true, instead of ranking by a single filter, rank by the
	/// temperature-weighted sum of all scores
	bool boltz_rank_ = false;

	/// @brief current
	Real current_score_ = 0.0;

	/// @brief accepted structure
	Real last_accepted_score_ = 0.0;

	/// @brief lowest energy structure we've seen
	Real lowest_score_ = 0.0;

	/// @brief accepted structure
	PoseOP last_accepted_pose_;

	/// @brief lowest energy structure we've seen
	PoseOP lowest_score_pose_;

	/// @brief result of the last call to boltzmann
	moves::MCA mc_accepted_ = moves::MCA_accepted_score_beat_last;

	/// @brief to change the sing of calculated "score"
	Real flip_sign_ = 1;

	/// @brief diagnostics
	int accept_counter_ = 0;
	Real energy_gap_counter_ = 0.0;

	/// @brief Next trigger identifier to be assigned
	Size next_trigger_id_ = 1;

	/// @brief Collection of function callbacks
	boost::unordered_map<Size, GenericMonteCarloMoverTrigger> triggers_;
	/// @brief dflt false_filter; use this to stop an MC trajectory before maxtrials_ (if filter evaluates to true)
	protocols::filters::FilterOP stopping_condition_;
	/// @brief dflt NULL; if the mover defined a stopping condition on the datamap then this assumes the mover's value. In this way, the Mover can tell GenericMC to stop execution, e.g., if it has iterated over all internal possibilities
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< bool > > mover_stopping_condition_;
	/// @brief  change the mover probabilities according to the accept rates?; only works if the mover is a ParsedProtocol type with mode=single_random
	bool adaptive_movers_ = false;
	/// @brief dflt max( 10, trials/10 ); only works with adaptive; how often should the run probabilities be adapted?
	core::Size adaptation_period_ = 0;
	/// @brief dflt ""; if a file name is specified, after each accept a pdb file is dumped to disk. This is useful for checkpointing
	std::string saved_accept_file_name_;
	/// @brief dflt ""; if specified checkpoints the current trial number and recovers from it
	std::string saved_trial_number_file_;
	/// @brief dflt NULL; this is used by the called movers to set a certain tag. If saved_accept_file_name_ is set, then at exit the tag coming from the chosen mover is written to disk as, <saved_accept_file_name>.mover_tag. To work, mover_tag_ must be exposed to the movers being called.
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< std::string > > mover_tag_;
	/// @brief dflt ""; the mover being called by GenericMC. Used to add values to the poses DataCache.
	std::string user_defined_mover_name_;
	/// @brief reset the filters' baseline at trial=1?
	bool reset_baselines_ = true;
	/// @brief avoid clearing filters if a score function is given. Used in the EvolutionDynamicsMover which inherits from genericMCmover.
	bool keep_filters_ = false;
	/// @brief a file name where data on each step is saved.
	std::string progress_file_;
};

} // namespace simple_moves
} // namespace protocols

#endif
