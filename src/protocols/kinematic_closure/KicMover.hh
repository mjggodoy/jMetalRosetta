// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_KicMover_HH
#define INCLUDED_protocols_kinematic_closure_KicMover_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/KicMover.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/LoopMover.hh>

// Utility headers
#include <boost/noncopyable.hpp>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace kinematic_closure {

/// @brief Find a new backbone conformation for some region of a protein.
///
/// @details The first step in the kinematic closure algorithm is to pick a
/// loop and three pivot residues.  Any residues in the loop that are not
/// pivots are called non-pivots.  The non-pivot backbone torsions are used to
/// make new conformations, while the pivot torsions are used to ensure that
/// the backbone stays closed.  Use set_loops() to specify loops to sample.
/// Use set_pivot_picker() to specify how the pivots should be chosen.  Use
/// add_perturber() to specify how the non-pivots should be sampled.  By
/// default, the algorithm will pick pivots randomly within the region being
/// sampled and will sample the non-pivot torsions from a rama distribution.
///
/// Given a set of pivot residues and nonpivot torsion angles, the algorithm
/// will find up to 16 possible solutions.  The set_solution_picker() method
/// allows you to control which solution, if any, is picked.  By default, the
/// first solution found which passes both a rama and a bump check is used.
///
/// @note The default kinematic closure algorithm samples both pivot and
/// nonpivot torsions from a rama distribution.  Because this algorithm is
/// often used in situations where the score function also contains a rama
/// term, the rama bias is usually double-counted.  This is bad, because it
/// means that backbone torsions are normally sampled too narrowly.  There are
/// two proper ways to deal with this, but they both have drawbacks.  The first
/// would be to sample from a uniform distribution and to let the score
/// function take care of preferring angles that fit the rama distribution.
/// Unfortunately, this would be much less efficient than the current approach,
/// primarily because the rama check is fast and filters out a lot of bad
/// solutions.  The second approach would be to sample from a rama distribution
/// and to remove the rama term from the score function.  This would be even
/// more efficient than the current approach, but it would be prone to bugs
/// because other parts of the code wouldn't expect the score function to be
/// changing on every invocation of KIC.
///
/// Once the algorithm has been setup using the helper methods described above,
/// apply() can be called to actually sample a new backbone conformation.

class KicMover : public protocols::loop_modeling::LoopMover {

public:

	/// @brief Default constructor.
	KicMover();

	/// @brief Default destructor.
	~KicMover() override;

	/// @brief Return the name of this mover.
	// XRW TEMP  std::string get_name() const override { return "KicMover"; }

public:

	/// @copydoc LoopMover::get_children_names
	void get_children_names(
		utility::vector1<std::string> & names, std::string indent="") const override;

	/// @brief Return the PivotPicker being used by this mover.
	pivot_pickers::PivotPickerOP get_pivot_picker();

	/// @brief Return the SolutionPicker being used by this mover.
	solution_pickers::SolutionPickerOP get_solution_picker();

	/// @brief Specify how the non-pivot torsions should be sampled.
	void add_perturber(perturbers::PerturberOP perturber);

	/// @brief Forget any perturbers that had been previously added.
	void clear_perturbers();

	/// @brief Specify how the pivot residues should be chosen.
	void set_pivot_picker(pivot_pickers::PivotPickerOP picker);

	/// @brief Specify how a solution should be chosen.
	void set_solution_picker(solution_pickers::SolutionPickerOP picker);

	/// @copydoc LoopMover::request_fold_tree
	protocols::loop_modeling::FoldTreeRequest request_fold_tree() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:

	/// @brief Sample a new backbone conformation for the given loop.
	bool do_apply(Pose & pose, Loop const & loop) override;

private:
	perturbers::PerturberSetOP perturbers_;
	pivot_pickers::PivotPickerOP pivot_picker_;
	solution_pickers::SolutionPickerOP solution_picker_;

};

}
}

#endif
