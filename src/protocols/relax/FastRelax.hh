// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/FastRelax.hh
/// @brief The FastRelax Protocol
/// @details
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_FastRelax_hh
#define INCLUDED_protocols_relax_FastRelax_hh

// Unit headers
#include <protocols/relax/FastRelax.fwd.hh>

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/checkpoint/CheckPointer.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace relax {


struct RelaxScriptCommand {
	RelaxScriptCommand(){
		command="";
		param1=0;
		param2=0;
		param3=0;
		param4=0;
		nparams=0;
	}
	std::string command;
	core::Real  param1;
	core::Real  param2;
	core::Real  param3;
	core::Real  param4;
	core::Size  nparams;
	utility::vector1< core::Real > params_vec;
};

class FastRelax : public RelaxProtocolBase {
public:

	/// @brief Initialize FastRelax with a specific script file, which encodes the script of steps
	///  the is to be applied
	FastRelax(
		core::Size                     standard_repeats = 0
	);

	/// @brief Initialize FastRelax using the default script with a varying number of rounds,
	///  defined by the standard_repeats repeats paramter. By default, 5.
	FastRelax(
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::Size                     standard_repeats = 0
	);

	/// @brief Initialize FastRelax with a specific script file, which encodes the script of steps
	///  the is to be applied
	FastRelax(
		core::scoring::ScoreFunctionOP scorefxn_in,
		std::string const &            script_file
	);

	/// @brief Initialize FastRelax with a specific script file, and allows fastrelax to only use a subset of the number of stages in the scriptfile.
	FastRelax(
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::Size                     standard_repeats,
		const std::string &          script_file
	);
	/// @brief Override the stored script with the default script for batchrelax
	///        Ignores '-default_repeat' value
	void set_script_to_batchrelax_default(
		core::Size                     standard_repeats = 0
	);

	/// @brief replace schedule by providing lines of commands explicitly
	void set_script_from_lines( std::vector< std::string > const & filelines );

	/// @brief Force us to batchrelax with nonideal geometry (using additional memory)
	void set_force_nonideal( bool val ) { force_nonideal_ = val; }

	/// @brief Use the dualspace (Dihedral + Cartesian) relax protocol (Default false)
	/// @details
	/// Sets to use the lbfgs_armijo_nonmonotone if true or FR default if false
	/// Recommended to set max_iter to 200.
	/// Recommended to set minimize_bond_angles to true as well.
	/// Requires scorefunction setup for non-ideal minimization.
	void dualspace( bool val );
	bool dualspace();

	/// @brief Return the number of repeats
	core::Size default_repeats() const { return default_repeats_; }

	/// @brief virtual constructor to allow derivation
	~FastRelax() override;

	/// @brief Parses the FastRelaxTags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/// @brief Initializes class using option system. This is called by the constructors
	void set_to_default();

	/// @brief Return the name of this mover.
	// XRW TEMP  std::string get_name() const override;

	/// @brief return a fresh instance of this class in an owning pointer
	protocols::moves::MoverOP clone() const override;

	/// @brief Apply the FastRelax. Overloaded apply function from mover base class.
	void apply( core::pose::Pose & pose ) override;

	/// @brief Batch Relax, a new even faster way to relax entire batches of structures.
	void batch_apply(
		std::vector < core::io::silent::SilentStructOP > &  input_structs,
		core::scoring::constraints::ConstraintSetOP input_csts = nullptr,
		core::Real decay_rate = 0.5 );

	/// @brief sets the enable_design option.
	///  Note - (Only used by ParseMyTag! Class will respect the passed TF otherwise, as it should.)
	///
	inline void set_enable_design( bool const val ) { enable_design_ = val; }

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
	complex_type_generator_for_fast_relax( utility::tag::XMLSchemaDefinition &);
protected:

	void cmd_accept_to_best(
		const core::scoring::ScoreFunctionOP local_scorefxn,
		core::pose::Pose &pose,
		core::pose::Pose &best_pose,
		const core::pose::Pose &start_pose,
		core::Real       &best_score,
		core::Size       &accept_count
	);

	void do_minimize(
		core::pose::Pose &pose,
		core::Real tolerance,
		core::kinematics::MoveMapOP local_movemap,
		core::scoring::ScoreFunctionOP local_scorefxn
	);

	/// @brief Use for constraint scaling -- sets the coordinate constraint weight given a scorefxn, energy map, and weight.
	/// TL - Derived classes which may want to scale more than coordinate constraints can override this without forking apply()
	virtual void
	set_constraint_weight( core::scoring::ScoreFunctionOP local_scorefxn,
		core::scoring::EnergyMap const & full_weights,
		core::Real const weight,
		core::pose::Pose & pose ) const;

	void do_md(
		core::pose::Pose &pose,
		core::Real nstep,
		core::Real temp0,
		core::kinematics::MoveMapOP local_movemap,
		core::scoring::ScoreFunctionOP local_scorefxn
	);

protected:
	bool script_file_specified_; // whether user specified script file

private:

	void read_script_file( const std::string &script_file, core::Size standard_repeats = 5  );

	void check_nonideal_mintype();

private:   // options

	/// @brief Number of repeats if using default script
	core::Size default_repeats_;


	/// @brief If true, delete virtual residues at the end of FastRelax for smooth rmsd calculation later
	bool delete_virtual_residues_after_FastRelax_;

	/// @brief Apply Ramady(tm) fix to residues with bad rama scores ?
	bool ramady_;

	/// @brief Number of bad ramas to fix per call to 'fix_worst_bad_ramas'
	core::Size ramady_num_rebuild_;

	/// @brief Reject ramady changes perturbing structure more than this amount
	core::Real ramady_rms_limit_;

	/// @brief Cutoff for calling a rama 'bad'
	core::Real ramady_cutoff_;

	/// @brief Force ramady to be run (normally skip rate of 10%)
	bool ramady_force_;

	/// @brief Allow Chi angles to move ?
	bool repack_;

	/// @brief Allow DNA to move?  default is False.  If True will setup DNA-specific relax settings.
	bool dna_move_;

	/// @brief Do only a few test_cycles ?
	bool test_cycles_;

	/// @brief [batch mode only] force structures to be nonideal?  Default uses -out:silent_struct_type to decide
	bool force_nonideal_;

	/// @brief Dump pdb after repack, min, or ramp_repack_min?
	bool dumpall_;

	/// @brief Use dualspace (Dihedral + Cartesian) Relax?
	bool dualspace_;

	/// @brief Enable design -- if false, a RestrictToRepacking task will be added to the TaskFactory
	bool enable_design_;

	/// @brief Quit after this many accepts ?// limits themaximum number of accepts, default is 1000000
	core::Size script_max_accept_;

	/// @brief Do a symmetric RMSD calculation rather then a monomeric one.
	bool symmetric_rmsd_;

private:   // other data
	protocols::checkpoint::CheckPointer checkpoints_;

	std::vector <RelaxScriptCommand> script_;
	utility::tag::TagCOP movemap_tag_; // this cannot be parsed before apply b/c the fold tree is likely to change during a run
};


}
} // protocols

#endif
