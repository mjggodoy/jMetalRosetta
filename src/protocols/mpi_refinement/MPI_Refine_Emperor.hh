// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mpi_refinement/MPI_Refine_Emperor.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_mpi_refinement_MPI_Refine_Emperor_hh
#define INCLUDED_protocols_mpi_refinement_MPI_Refine_Emperor_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/mpi_refinement/MPI_Refinement.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

#include <string>
#include <vector>

namespace protocols {
namespace mpi_refinement {

class MPI_Refine_Emperor: public MPI_Refinement {

public:

	MPI_Refine_Emperor():
		MPI_Refinement( 'E' )
	{
		set_defaults();
	}

	~MPI_Refine_Emperor() override= default;

	void set_defaults();

public:

	void go() override;

protected:

	void init() override;

	void process_inbound_wus() override;

	void process_outbound_wus() override;


	// adding arriving structures to library
	bool add_structures_to_library( protocols::wum::SilentStructStore &new_structs, std::string add_algorithm = "" ) override;

private:
	bool process_termination();

private:

	bool process_termination_;
	bool termination_broadcasted_;

	std::map< core::Size, bool > masters_done_;
	core::Size max_emperor_lib_round_;
	bool dump_rounds_;

	// For Emperor library report frequency
	core::Size n_accept_cummul_;
	core::Size n_addcall_;
	core::Size n_change_emperor_report_;
	core::Size n_addcall_emperor_report_;
	core::Size n_dump_;

};


} // namespace mpi_refinement
} // namespace protocols

#endif
