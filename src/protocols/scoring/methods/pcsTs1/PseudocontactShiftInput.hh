// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/PseudocontactShiftInput.hh
///
/// @brief Read input .npc input file
///
/// @details The following classes are responsable to read / parse the PCS input file (.npc format)
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcsTs1_PseudocontactShiftInput_hh
#define INCLUDED_protocols_scoring_methods_pcsTs1_PseudocontactShiftInput_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <string>
#include <map>

//Auto Headers
#include <utility/vector1_bool.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs1 {

///////////////////////////////////////////////////////////////////////////
/// @brief PCS_line_data_Ts1 class: hold a line of the input file information (.npc format)
/// One PCS_line_data_Ts1 per line in the input file
class PCS_line_data_Ts1 {

private:
	/// @brief Cannot use default constructor!
	PCS_line_data_Ts1();

public:
	~PCS_line_data_Ts1();

	PCS_line_data_Ts1(PCS_line_data_Ts1 const & other);

	PCS_line_data_Ts1 &
	operator=( PCS_line_data_Ts1 const & other );

	PCS_line_data_Ts1(
		core::Size const residue_num,
		std::string const atom_name,
		core::Real const PCS_experimental,
		core::Real const PCS_tolerance
	);

	core::Size
	residue_num() const;

	std::string
	atom_name() const;

	core::Real
	PCS_experimental() const;

	core::Real
	PCS_tolerance() const;

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_line_data_Ts1 &PCS_l_d);

private:
	core::Size const residue_num_;
	std::string const atom_name_;
	core::Real const PCS_experimental_;
	core::Real const PCS_tolerance_;
};


//////////////////////////////////////////////////////////
/// @brief PCS_file_data_Ts1 contain all the information of a .npc file
/// one per lanthanide.
class PCS_file_data_Ts1 {
private:
	std::string const filename_;
	utility::vector1<PCS_line_data_Ts1> PCS_data_line_all_;
	core::Real const weight_;

	void
	read_PCS_file();

public:
private:
	PCS_file_data_Ts1();
public:
	~PCS_file_data_Ts1();

	PCS_file_data_Ts1(PCS_file_data_Ts1 const & other);

	PCS_file_data_Ts1 &
	operator=( PCS_file_data_Ts1 const & other );

	PCS_file_data_Ts1(std::string const & filename, core::Real const my_weight );

	std::string
	get_filename() const;

	core::Real
	get_weight() const;

	utility::vector1<PCS_line_data_Ts1> &
	get_PCS_data_line_all_reference();

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_file_data_Ts1 &PCS_f_d);
};


//////////////////////////////////////////////////////////////
/// @brief PCS_data_input_Ts1 contain all the input information for the PCS.
/// This includes all the information from the .npc files
class PCS_data_input_Ts1 {
private:
	std::map< std::string, PCS_file_data_Ts1 > PCS_filename_and_data_;

public:
	PCS_data_input_Ts1();

	~PCS_data_input_Ts1();

	PCS_data_input_Ts1(PCS_data_input_Ts1 const & other);

	PCS_data_input_Ts1 &
	operator=( PCS_data_input_Ts1 const & other );

	PCS_data_input_Ts1(utility::vector1<std::string> const & filenames,  utility::vector1<core::Real> const & weight);

	std::map< std::string, PCS_file_data_Ts1 > &
	get_PCS_data_input_reference();

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_data_input_Ts1 &PCS_d_i);

};

class PCS_data_input_manager_Ts1 {

private:
	PCS_data_input_manager_Ts1();

	static PCS_data_input_manager_Ts1 * instance_;
	std::map<std::string, PCS_data_input_Ts1> file_2_data_map_;

public:
	static
	PCS_data_input_manager_Ts1 *
	get_instance();

	PCS_data_input_Ts1
	get_input_data(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & vec_weight);
};


}//namespace pcsTs1
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
