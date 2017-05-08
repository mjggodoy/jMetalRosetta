// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/Transform.hh
/// @author Thomas Willcock and Darwin Fu
/// Adapted from code by Sam Deluca

#ifndef INCLUDED_protocols_ligand_docking_TransformEnsemble_hh
#define INCLUDED_protocols_ligand_docking_TransformEnsemble_hh

#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/TransformEnsemble.fwd.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/vector1.hh>
#include <vector>


namespace protocols {
namespace ligand_docking {


struct TransformEnsemble_info{ // including default values

public:
	//change first three into vectors (create a vector of chains)
	utility::vector1<std::string> chains;
	utility::vector1<core::Size> chain_ids;
	utility::vector1<core::Size> jump_ids;
	core::Real move_distance;
	core::Real box_size;
	core::Real angle;
	core::Size cycles;
	core::Real temperature;
	core::Size repeats;
	TransformEnsemble_info(): chains(), move_distance(0),box_size(0), angle(0), cycles(0),repeats(1){};

	TransformEnsemble_info(TransformEnsemble_info const & other) :
		chains(other.chains),
		chain_ids(other.chain_ids),
		jump_ids(other.jump_ids),
		move_distance(other.move_distance),
		box_size(other.box_size),
		angle(other.angle),
		cycles(other.cycles),
		temperature(other.temperature),
		repeats(other.repeats){}
};

class TransformEnsemble: public protocols::moves::Mover
{
public:
	TransformEnsemble();
	TransformEnsemble(
		utility::vector1<std::string> const & chains,
		core::Real const & box_size,
		core::Real const & move_distance,
		core::Real const & angle,
		core::Size const & cycles,
		core::Real const & temp
	);
	TransformEnsemble(TransformEnsemble const & other);
	~TransformEnsemble() override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose
	) override;
	core::Vector weighted_center(utility::vector1<core::conformation::UltraLightResidue> & residues);
	void apply(core::pose::Pose & pose) override;

	void translate_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues, core::Real distance);
	void transform_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues);
	void change_conformers(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, const utility::vector1<core::conformation::UltraLightResidue> & reference_residues);
	void change_conformer(core::conformation::UltraLightResidue & ligand_residue, const core::conformation::UltraLightResidue & reference_residue, core::Size resid);
	void dump_conformer(core::conformation::UltraLightResidue & residue, utility::io::ozstream & output);
	void print_xyz(core::Vector vector);
	bool monte_carlo(core::Real & current, core::Real & last);
	bool check_grid(qsar::scoring_grid::GridManager* grid, utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, core::Real distance = 0);

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

private:
	//qsar::scoring_grid::GridManagerOP grid_manager_;

	TransformEnsemble_info transform_info_;

	utility::vector1<core::conformation::UltraLightResidue> best_ligands_;
	utility::vector1<core::conformation::UltraLightResidue> ligand_residues_;
	utility::vector1<core::conformation::UltraLightResidue> reference_residues_;
	utility::vector1<core::conformation::UltraLightResidue> last_accepted_ligand_residues_;
	utility::vector1<core::conformation::UltraLightResidue> last_accepted_reference_residues_;

	std::map<core::Size, utility::vector1< core::conformation::ResidueOP > > ligand_conformers_;
	bool optimize_until_score_is_negative_;
	bool output_sampled_space_;
	bool optimize_until_ideal_;
	bool use_conformers_;


	std::string sampled_space_file_;
	core::Real initial_perturb_;


};

}
}

#endif /* TRANSFORM_HH_ */
