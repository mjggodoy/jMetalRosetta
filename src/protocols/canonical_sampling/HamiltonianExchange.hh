// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief Implement replica exchange in the MetropolisHastingsMover Framework.
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_HamiltonianExchange_hh
#define INCLUDED_protocols_canonical_sampling_HamiltonianExchange_hh

// Unit Headers
#include <protocols/canonical_sampling/HamiltonianExchange.fwd.hh>
#include <protocols/canonical_sampling/AsyncMPITemperingBase.hh>

// Module Headers
#include <protocols/canonical_sampling/BiasEnergy.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace canonical_sampling {

/// @details
class HamiltonianExchange : public protocols::canonical_sampling::AsyncMPITemperingBase {
	typedef AsyncMPITemperingBase Parent;
	// typedef utility::vector1< core::Size > GridCoord; // moving to base class TemperatureController
	typedef utility::vector1< GridCoord > Grid;

public:

	HamiltonianExchange();

	//important d'tor to delete some C-style arrays
	~HamiltonianExchange() override;

	HamiltonianExchange( HamiltonianExchange const& );

	HamiltonianExchange& operator=( HamiltonianExchange const& );


	void apply( core::pose::Pose& ) override {};


	// XRW TEMP  std::string
	// XRW TEMP  get_name() const override;

	protocols::moves::MoverOP
	clone() const override;


	protocols::moves::MoverOP
	fresh_instance() const override;


	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	/// @brief not possible for HamExchange -- exit with ERROR if called
	core::Real
	temperature_move( core::Real score ) override;

	// Undefined, commenting out to fix PyRosetta build  void update_bias_energies( int exchange_partner, bool is_master );

	/// @brief execute the temperatur move ( called by observer_after_metropolis )
	/// returns the current temperatur in kT.
	core::Real
	temperature_move( core::pose::Pose& pose ) override;

	/// @brief callback executed before any Monte Carlo trials
	void
	initialize_simulation(
		core::pose::Pose& pose,
		MetropolisHastingsMover const& metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	) override;


	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size level,
		core::Real temperature,
		core::Size cycle
	) override;

	void show( std::ostream& ) const override;

	void next_exchange_schedule();

	void
	set_monte_carlo(
		protocols::moves::MonteCarloOP monte_carlo
	) override;

	GridCoord
	level_2_grid_coord( core::Size level ) const override { // for higher dimension grids
		return exchange_grid_[level];
	}

	core::Size
	exchange_grid_dim() const override {
		return exchange_grid_dimension_;
	}

	core::Size
	nlevels_per_dim( core::Size ) const override;

protected:
	void set_defaults();
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options() override;

	/// @brief return to uninitialized status
	void clear();

	/// @brief initialize temperatures and weights from file, return false if IO error occurrs
	bool initialize_from_file( std::string const& filename ) override;

	/// @brief initialize exchange schedule from file, return false if IO error occurrs
	virtual bool initialize_exchange_schedule_from_file( std::string const& filename );


	core::Size next_exchange_level() const override;

private:
	void setup_exchange_schedule();
	/// @brief small helper function; compute unique key out of (z1, z2, ... , zN) excluding zD
	core::Size coord2key(
		GridCoord const& coord,
		GridCoord const& max_coord,
		Size exclude_dim = 0
	);
	/// ------------------ register cmdline options ---------------------------

private:
	static bool options_registered_;

public:
	static void register_options();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/// ---------------- member variables --------------------------

private:

	/// the different score-funcitons
	utility::vector1< core::scoring::ScoreFunctionOP > hamiltonians_;

	// exchange schedule determines which nodes try to swap at a given trial phase
	typedef utility::vector1<std::pair<int, int> > ExchangeSchedule;
	utility::vector0< ExchangeSchedule > exchange_schedules_;
	core::Size current_exchange_schedule_;

	// exchange happens between neighbouring cells in a D dimensional grid
	// in temperature replica-exchange usually D=1.
	Grid exchange_grid_;
	core::Size exchange_grid_dimension_;
	bool successfully_initialized_;
	GridCoord max_coord_;
	std::string exchange_schedule_file_;

	BiasEnergyOP bias_energy_;
}; //end HamiltonianExchange

/// @brief Test IO operator for debug and Python bindings
std::ostream& operator << ( std::ostream & os, HamiltonianExchange const& );

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_HamiltonianExchange_HH
