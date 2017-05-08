// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/StartFrom.hh
/// @brief  a mover to place a ligand at a defined position
/// @author Gordon Lemmon
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_ligand_docking_StartFrom_hh
#define INCLUDED_protocols_ligand_docking_StartFrom_hh

#include <protocols/ligand_docking/StartFrom.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class StartFrom : public protocols::moves::Mover
{
public:
	StartFrom();
	~StartFrom() override;
	StartFrom(StartFrom const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/// @brief Add the given coordinates as a valid starting position for the given pdb_tag
	void add_coords(core::Vector const & coords, std::string const & pdb_tag = "default");

	/// @brief Add the given coordinates as a valid starting position for the given structure hash
	void add_coords_hash(core::Vector const & coords, std::string const & hash_value);

	/// @brief Set which chains the mover operates on.
	void chain(utility::vector1<std::string> const & chains) { chains_ = chains; }

	/// @brief Add a chain to the mover targets
	void chain(std::string const & chain) { chains_.push_back(chain); }

	/// @brief Get which chains the mover operates on
	std::string chain()
	{
		std::string chains;
		for ( utility::vector1<std::string>::const_iterator i = chains_.begin(); i != chains_.end(); ++i ) {
			chains += *i;
		}
		return chains;
	}

	/// @brief Set if we should use the neighbor atom or the all-atom centroid to center the ligand
	void use_nbr(bool setting) { use_nbr_ = setting; }

	/// @brief Get whether we should use the neighbor atom or the all-atom centroid to center the ligand
	bool use_nbr() { return use_nbr_; }

	void apply(core::pose::Pose & pose) override;

	/// @brief Parse the json-format startfrom file
	void parse_startfrom_file(std::string const & filename);

	/// @brief Parse a PDB file, grabbing the positions from the heavy atom coordinates.
	/// If atom_name is not empty, only grab coordinates from the specified atom name.
	void parse_pdb_file(std::string const & filename, std::string const & atom_name = "", std::string const & tag="default");

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief The chain which to move
	utility::vector1<std::string> chains_;
	/// @brief If true, try to center the chain based on the neighbor atom of the first residue
	/// Otherwise, use the all-atom centroid of the chain.
	bool use_nbr_;

	/// @brief The possible starting positions, indexed by tag (or "default")
	std::map< std::string, utility::vector1<core::Vector> > starting_positions_;

	/// @brief The possible starting positions, indexed by hash
	std::map< std::string, utility::vector1<core::Vector>  > hash_starting_positions_;
};

} //namespace ligand_docking
} //namespace protocols

#endif
