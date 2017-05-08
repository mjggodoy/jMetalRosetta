// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/PDBJobOutputter.hh
/// @brief  header file for PDBJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_PDBJobOutputter_hh
#define INCLUDED_protocols_jd2_PDBJobOutputter_hh

//unit headers
#include <protocols/jd2/PDBJobOutputter.fwd.hh>
#include <protocols/jd2/wwPDBJobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

/// @details this simplest implementation of JobOutputter outputs raw PDBs and associated files, uncompressed.
class PDBJobOutputter : public protocols::jd2::wwPDBJobOutputter
{
public:

	typedef protocols::jd2::wwPDBJobOutputter parent;

	PDBJobOutputter();

	~PDBJobOutputter() override;

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function takes a string and writes it to disk (separately from Tracer output).  This implementation writes a single file whose filename is based on the job and a user-specified extension (default .data)
	// virtual --> moved to FileJobOutputter
	// void file( JobCOP job, std::string const & data );

	//Exists in parent wwPDBJO
	/// @brief this function outputs the final result of a job.  This implementation will write a PDB file (plus scores).
	//virtual
	//void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag );

	//Exists in parent wwPDBJO
	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.  This implementation will write a PDB file (plus scores).
	//virtual
	//void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	//Exists in parent wwPDBJO
	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb with the job's name already in existence.
	//virtual
	//bool job_has_completed( JobCOP job );

	//Exists in parent wwPDBJO
	/// @brief this is the master function for determining the unique output identifier for a job
	//virtual
	//std::string output_name( JobCOP job );

protected:
	//Exists in parent wwPDBJO
	/// @brief this private function provides the extended name, not just the output name.  e.g output_name returns 1UBQ_0001, this returns 1UBQ_0001.pdb.  In this case the extension is .pdb

	//virtual
	//std::string extended_name( JobCOP job, std::string const & suffix = "" );

	////////////////////////////////////////score-related functions///////////////////////////////////

	/// @brief this function extracts the pose's scores for printing
	//virtual
	//std::string extended_name( JobCOP job, std::string const suffix = "" );

	////////////////////////////////////////score-related functions///////////////////////////////////

	//Exists in parent wwPDBJO
	/// @brief this function extracts the pose's scores and outputs them as a string to be packaged in an output structure.
	/// @details Refactored in the 2016 Chemical XRW (eXtreme Rosetta Workshop) by Vikram K. Mulligan (vmullig@uw.edu).
	/// @param[in] job Const-access owning pointer to the job from which the data will be extracted.
	/// @param[out] data_out A string in which the data will be stored, that can later be passed to whatever container wants it.
	//virtual
	//std::string extract_data_from_Job( JobCOP job );

protected:
	//////////////////////////////////////protected PDB output/////////////////////////////////////
	/// @brief This is the function actually different between mmCIF and PDB output, and unshared by the wwPDB parent class.  Here it causes a pdb file to be written.  Pure virtual in the base class.  filename is an optional label for the score data table, not an actual control.

	void dump_pose( JobCOP job, core::pose::Pose const & pose, utility::io::ozstream & out, std::string const &filename="" ) override;

	////////////////////////////////////////data////////////////////////////////////////
private:

	//std::string extension_;  //exists in parent class
	//std::string path_; //exists in parent class

}; // PDBJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_PDBJobOutputter_HH
