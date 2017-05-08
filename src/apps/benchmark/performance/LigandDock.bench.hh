// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/LigandDock.bench.hh
///
/// @brief Dock the ligand in the 7cpa complex.
/// Use all options (flexible ligand, flexible backbone)
/// @author Gordon Lemmon


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Auto-header: duplicate removed #include <core/import_pose/import_pose.hh>

//Auto Headers
#include <utility/vector1.hh>

class LigandDockBenchmark : public PerformanceBenchmark
{
public:
	LigandDockBenchmark(std::string name) : PerformanceBenchmark(name) {};

	core::pose::Pose ligand_dock_pose;

	virtual void setUp() {
		basic::options::option.load_options_from_file("ligand_dock/ligand_dock_flags.txt");

		std::string pdb_file_name= basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
		core::import_pose::pose_from_file(ligand_dock_pose, pdb_file_name, core::import_pose::PDB_file);
	};

	virtual void run(core::Real scaleFactor) {
		protocols::ligand_docking::LigandDockProtocol dock_protocol;
		core::Size reps( (core::Size)(1 * scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scale factor.
		for ( core::Size i=0; i< reps; i++ ) {
			dock_protocol.apply(ligand_dock_pose);
		}
	};

	virtual void tearDown() {};
};
