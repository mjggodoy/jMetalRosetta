// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/InterchainPairEnergy.hh
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///      and between the two interfaces need to be evaluated
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_scoring_methods_InterchainPairEnergy_hh
#define INCLUDED_protocols_scoring_methods_InterchainPairEnergy_hh

#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <protocols/scoring/InterchainPotential.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace protocols {
namespace scoring {
namespace methods {


class InterchainPairEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;
public:


	InterchainPairEnergy();


	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	virtual
	void
	setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::pose::Pose const &,// pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const;

	virtual
	void
	finalize_total_energy(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// emap
	) const;

	/// This method *should* admit to defining intraresidue energies
	virtual
	bool
	defines_intrares_energy( core::scoring::EnergyMap const & ) const { return false; }

	void
	eval_intrares_energy(
		core::conformation::Residue const &,
		core::pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const {}

	virtual
	core::Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	protocols::scoring::InterchainPotential const & interchain_potential_;
	core::scoring::EnvPairPotential const & env_potential_;

	virtual
	core::Size version() const;
};


}
}
}

#endif
