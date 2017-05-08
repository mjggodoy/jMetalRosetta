// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/PeptideStubMoverCreator.hh>

#include <core/id/AtomID.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/chemical/VariantType.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.PeptideStubMover" );

namespace protocols {
namespace cyclic_peptide {

PeptideStubMover::PeptideStubMover(){
	init();
}

PeptideStubMover::PeptideStubMover( PeptideStubMover const &src ) :
	Mover("PeptideStubMover"),
	reset_(src.reset_),
	update_pdb_numbering_(src.update_pdb_numbering_),
	stub_mode_(src.stub_mode_),
	stub_rsd_names_(src.stub_rsd_names_),
	stub_rsd_jumping_(src.stub_rsd_jumping_),
	stub_rsd_connecting_atom_(src.stub_rsd_connecting_atom_),
	stub_rsd_repeat_(src.stub_rsd_repeat_),
	stub_insert_pos_(src.stub_insert_pos_),
	stub_anchor_rsd_(src.stub_anchor_rsd_),
	stub_anchor_rsd_connecting_atom_(src.stub_anchor_rsd_connecting_atom_)
{
}

PeptideStubMover::~PeptideStubMover()= default;

void PeptideStubMover::init() {
	reset_ = false;
	update_pdb_numbering_ = true;
}

void PeptideStubMover::apply( core::pose::Pose & pose )
{
	core::chemical::ResidueTypeSetCOP standard_residues = pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );
	if ( reset_ ) {
		pose.clear();
	}

	using namespace core::chemical;

	if ( pose.size() != 0 ) {
		for ( core::Size ir=1; ir<=pose.size() ; ++ir ) {
			if ( pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, ir );
			}
			if ( pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, ir );
			}
			if ( pose.residue_type(ir).has_variant_type(UPPER_TERMINUS_VARIANT) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS_VARIANT, ir );
			}
			if ( pose.residue_type(ir).has_variant_type(LOWER_TERMINUS_VARIANT) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS_VARIANT, ir );
			}
		}
		pose.update_residue_neighbors();
	}
	for ( core::Size istub=1; istub<=stub_rsd_names_.size(); ++istub ) {
		core::conformation::ResidueOP new_rsd( nullptr );
		new_rsd = core::conformation::ResidueFactory::create_residue( standard_residues->name_map(stub_rsd_names_[istub]) );

		// first stub always starts with a jump
		if ( istub == 1 && pose.size() == 0 ) {
			runtime_assert_string_msg(stub_mode_[istub] == append, "Can only use append for the first residue");
			pose.append_residue_by_jump(*new_rsd, 1);
			for ( Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
				pose.append_residue_by_bond(*new_rsd, true);
			}
		} else {
			Size anchor_rsd(stub_anchor_rsd_[istub]);
			if ( anchor_rsd == 0 ) anchor_rsd = pose.size();

			if ( stub_rsd_jumping_[istub] ) {
				runtime_assert_string_msg(stub_mode_[istub] == append, "Can only use append for jumps");
				pose.append_residue_by_jump(*new_rsd, anchor_rsd);
				for ( Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
					pose.append_residue_by_bond(*new_rsd, true);
				}
			} else {
				core::Size connecting_id(0);
				if ( stub_rsd_connecting_atom_[istub] == "" ) {
					connecting_id = new_rsd->type().lower_connect_id();
				} else {
					core::Size atomid(new_rsd->atom_index(stub_rsd_connecting_atom_[istub]));
					connecting_id = new_rsd->type().residue_connection_id_for_atom(atomid);
					if ( connecting_id == 0 ) {
						TR << "Error! Residue " << stub_rsd_names_[istub] << " Atom " << stub_rsd_connecting_atom_[istub] << " cannot be connected" << std::endl;
					}
				}

				core::Size anchor_connecting_id(0);
				if ( stub_anchor_rsd_connecting_atom_[istub] == "" ) {
					anchor_connecting_id = pose.residue_type(anchor_rsd).upper_connect_id();
				} else {
					core::Size atomid(pose.residue_type(anchor_rsd).atom_index(stub_anchor_rsd_connecting_atom_[istub]));
					anchor_connecting_id = pose.residue_type(anchor_rsd).residue_connection_id_for_atom(atomid);

					if ( connecting_id == 0 ) {
						TR << "Error! Residue " << stub_rsd_names_[anchor_rsd] << " Atom " << stub_anchor_rsd_connecting_atom_[istub] << " cannot be connected" << std::endl;
					}
				}

				if ( stub_mode_[istub] == append ) {
					pose.append_residue_by_bond(*new_rsd, true, connecting_id, anchor_rsd, anchor_connecting_id, false);
					//rebuild the polymer bond dependent atoms:
					rebuild_atoms(pose, anchor_rsd);
				} else if ( stub_mode_[istub] == prepend ) {
					pose.prepend_polymer_residue_before_seqpos(*new_rsd, anchor_rsd, true);
					rebuild_atoms(pose, anchor_rsd);
				} else if ( stub_mode_[istub] == insert ) {
					if ( stub_insert_pos_[istub] != 0 ) {
						pose.insert_residue_by_bond(*new_rsd, stub_insert_pos_[istub], anchor_rsd, true, stub_anchor_rsd_connecting_atom_[istub], stub_rsd_connecting_atom_[istub], false, true);
						//pose.dump_pdb("test.pdb");
					} else {
						pose.append_polymer_residue_after_seqpos(*new_rsd, anchor_rsd, true);
					}
					rebuild_atoms(pose, anchor_rsd);
				}

				for ( Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
					if ( stub_rsd_connecting_atom_[istub] == "" && stub_anchor_rsd_connecting_atom_[istub] == "" ) {
						if ( stub_mode_[istub] == append ) {
							pose.append_residue_by_bond(*new_rsd, true);
						} else if ( stub_mode_[istub] == insert ) {
							if ( stub_insert_pos_[istub] == 0 ) {
								pose.append_polymer_residue_after_seqpos(*new_rsd, anchor_rsd + i_repeat - 1, true);
							}
						} else if ( stub_mode_[istub] == prepend ) {
							pose.prepend_polymer_residue_before_seqpos(*new_rsd, anchor_rsd, true);
						}
					} else {
						TR << "Cannot use repeat for non-canonical insertion" << std::endl;
						break;
					}
				}
			}

		}
		if ( TR.Debug.visible() ) TR.Debug << "Done appending " << stub_rsd_names_[istub] << std::endl;
	}

	//protocols::loops::add_cutpoint_variants( pose );
	for ( core::Size ir=1; ir<=pose.size() ; ++ir ) {
		if ( pose.residue_type(ir).lower_connect_id() != 0 ) {
			if ( pose.residue(ir).connected_residue_at_resconn(pose.residue_type(ir).lower_connect_id()) == 0 ) {
				if ( pose.residue(ir).is_protein() && !pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
					core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, ir);
				}
			}
		}
		if ( pose.residue_type(ir).upper_connect_id() != 0 ) {
			if ( pose.residue(ir).connected_residue_at_resconn(pose.residue_type(ir).upper_connect_id()) == 0 ) {
				if ( pose.residue(ir).is_protein() && !pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
					core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_UPPER, ir);
				}
			}
		}
	}

	//pose.dump_pdb("test.pdb");
	if ( TR.Debug.visible() ) {
		TR << pose.annotated_sequence() << std::endl;
		pose.fold_tree().show(TR.Debug);
	}
	for ( core::Size ires=1; ires<=pose.size(); ++ires ) {
		if ( pose.fold_tree().is_jump_point(ires) ) {
			TR.Debug << "jump: Residue " << ires << std::endl;
		}
		if ( pose.fold_tree().is_cutpoint(ires) ) {
			TR.Debug << "cut: Residue " << ires << std::endl;
		}

		for ( core::Size icon=1; icon<=pose.residue_type(ires).n_possible_residue_connections(); ++icon ) {
			TR.Debug << "connection: Residue " << ires << " Atom " << pose.residue_type(ires).residue_connection(icon).atomno() << " to residue " << pose.residue(ires).connected_residue_at_resconn(icon) << ", connect id:" << pose.residue(ires).connect_map(icon).connid() << std::endl;
		}
	}

	if ( update_pdb_numbering_ ) update_pdb_numbering(pose); //If residues have been added, they need PDB chain IDs and numbers to be updated.
}


/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
PeptideStubMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	using namespace core;
	if ( tag->hasOption( "reset" ) ) set_reset_mode( tag->getOption< bool >( "reset" ) );
	if ( tag->hasOption( "update_pdb_numbering" ) ) set_update_pdb_numbering_mode(tag->getOption< bool >( "update_pdb_numbering" ) );

	reset_mover_data();

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		add_residue(
			(*tag_it)->getName(),
			(*tag_it)->getOption<std::string>( "resname", "" ),
			(*tag_it)->getOption<Size>( "position", 0 ),
			(*tag_it)->getOption<bool>( "jump", false ),
			(*tag_it)->getOption<std::string>( "connecting_atom", "" ),
			(*tag_it)->getOption<Size>( "repeat", 1 ),
			(*tag_it)->getOption<core::Size>( "anchor_rsd", 0 ),
			(*tag_it)->getOption<std::string>( "anchor_atom", "" )
		);
		/*if ( (*tag_it)->getName() == "Append" ) {
		stub_mode_.push_back(append);
		}
		else if ( (*tag_it)->getName() == "Prepend" ) {
		stub_mode_.push_back(prepend);
		}
		else if ( (*tag_it)->getName() == "Insert" ) {
		stub_mode_.push_back(insert);
		}
		stub_rsd_names_.push_back( (*tag_it)->getOption<std::string>( "resname", "" ) );
		stub_insert_pos_.push_back( (*tag_it)->getOption<Size>( "position", 0 ) );
		stub_rsd_jumping_.push_back( (*tag_it)->getOption<bool>( "jump", false ) );
		stub_rsd_connecting_atom_.push_back( (*tag_it)->getOption<std::string>( "connecting_atom", "" ) );
		stub_rsd_repeat_.push_back( (*tag_it)->getOption<Size>( "repeat", 1 ) );
		stub_anchor_rsd_.push_back( (*tag_it)->getOption<core::Size>( "anchor_rsd", 0 ) );
		stub_anchor_rsd_connecting_atom_.push_back( (*tag_it)->getOption<std::string>( "anchor_atom", "" ) ); */
	}
}


/// @brief Adds a residue to the list of residues to be appended, prepended, or inserted.
void PeptideStubMover::add_residue(
	std::string const &stubmode,
	std::string const &resname,
	core::Size const position,
	bool const jumpmode,
	std::string const &connecting_atom,
	core::Size const repeat,
	core::Size const anchor_rsd,
	std::string const &anchor_atom
) {

	runtime_assert_string_msg( stubmode=="Append" || stubmode=="Prepend" || stubmode=="Insert", "In PeptideStubMover::add_residue(): unrecognized stub creation mode.  Mode must be one of \"Append\", \"Prepend\", or \"Insert\"." );

	if ( stubmode == "Append" ) {
		stub_mode_.push_back(append);
	} else if ( stubmode == "Prepend" ) {
		stub_mode_.push_back(prepend);
	} else if ( stubmode == "Insert" ) {
		stub_mode_.push_back(insert);
	}

	stub_rsd_names_.push_back( resname );
	stub_insert_pos_.push_back( position );
	stub_rsd_jumping_.push_back( jumpmode );
	stub_rsd_connecting_atom_.push_back( connecting_atom );
	stub_rsd_repeat_.push_back( repeat );
	stub_anchor_rsd_.push_back( anchor_rsd );
	stub_anchor_rsd_connecting_atom_.push_back( anchor_atom );

	return;
}


moves::MoverOP PeptideStubMover::clone() const { return moves::MoverOP( new PeptideStubMover( *this ) ); }
moves::MoverOP PeptideStubMover::fresh_instance() const { return moves::MoverOP( new PeptideStubMover ); }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PeptideStubMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PeptideStubMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PeptideStubMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PeptideStubMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PeptideStubMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "PeptideStubMover";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PeptideStubMover::get_name() const {
// XRW TEMP  return "PeptideStubMover";
// XRW TEMP }

//Private functions:

/// @brief Rebuilds all atoms that are dependent on bonds between residue_index and any other residues (including atoms on the other residues).
void PeptideStubMover::rebuild_atoms(
	core::pose::Pose &pose,
	core::Size const residue_index
) const {
	assert(residue_index <= pose.size() && residue_index > 0);

	core::Size const nresconn = pose.residue(residue_index).n_possible_residue_connections();
	if ( nresconn>0 ) {
		for ( core::Size ic=1; ic<=nresconn; ++ic ) {
			if ( !pose.residue(residue_index).connection_incomplete(ic) ) {
				core::Size const conn_at_index = pose.residue(residue_index).residue_connect_atom_index(ic); //The index of the connection atom
				for ( core::Size ia=1, iamax=pose.residue(residue_index).natoms(); ia<=iamax; ++ia ) {
					if ( pose.residue(residue_index).icoor(ia).stub_atom(1).atomno()==conn_at_index ) { //If this is a child of the connection atom
						//TR << "Rebuilding rsd " << residue_index << " atom " << ia << " (" << pose.residue(residue_index).atom_name(ia) << ")" << std::endl; TR.flush();
						pose.conformation().set_xyz( core::id::AtomID( ia, residue_index ), pose.residue(residue_index).icoor(ia).build(pose.residue(residue_index), pose.conformation()) );
					}
				}
				core::Size const other_rsd_index = pose.residue(residue_index).connect_map(ic).resid(); //The residue index of the other residue that this one is connected to at connection ID "ic"
				core::Size const other_rsd_connect_index = pose.residue(residue_index).connect_map(ic).connid(); //The connection index on the other residue that links to THIS residue.
				core::Size const other_at_index = pose.residue(other_rsd_index).residue_connect_atom_index( other_rsd_connect_index ); //The index of the atom on the other residue that links to THIS residue.
				for ( core::Size ia=1, iamax=pose.residue(other_rsd_index).natoms(); ia<=iamax; ++ia ) { //Loop through all atoms on the OTHER residue
					if ( pose.residue(other_rsd_index).icoor(ia).stub_atom(1).atomno()==other_at_index ) { //If this is an immediate child of the other residue's atom that connects to this atom, rebuild it.
						pose.conformation().set_xyz( core::id::AtomID(ia, other_rsd_index), pose.residue(other_rsd_index).icoor(ia).build( pose.residue(other_rsd_index), pose.conformation() ) );
					}
				}
			}
		}
	}

	return;
}


/// @brief Updates the PDB numbering (PDB number/chain ID) as residues are added.
void PeptideStubMover::update_pdb_numbering (
	core::pose::Pose &pose
) const {

	char const chainids []="0ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

	core::pose::PDBInfoOP pdbinfo = pose.pdb_info();
	if ( pdbinfo.get() != nullptr ) {
		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
			pdbinfo->set_resinfo(ir, (pose.chain(ir) <=62 ? chainids[ pose.chain(ir) ] : '0'), static_cast<int>(ir));
		}
	}

	return;
}

std::string PeptideStubMover::get_name() const {
	return mover_name();
}

std::string PeptideStubMover::mover_name() {
	return "PeptideStubMover";
}

void PeptideStubMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "reset", xsct_rosetta_bool, "Delete the input pose" )
		+ XMLSchemaAttribute( "update_pdb_numbering", xsct_rosetta_bool, "Update the PDB chain IDs and numbers for this pose" );
	//go through subelements
	AttributeList subelement_attributes; //all have same attributes
	subelement_attributes
		+ XMLSchemaAttribute( "resname", xs_string, "Name of stub residue" )
		+ XMLSchemaAttribute( "position", xsct_non_negative_integer, "Position to insert stub residue" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_rosetta_bool, "Append residue by jump?", "false" )
		+ XMLSchemaAttribute( "connecting_atom", xs_string, "Name of atom where residues should be connected" )
		+ XMLSchemaAttribute::attribute_w_default( "repeat", xsct_non_negative_integer, "Number of times to add this residue", "1" )
		+ XMLSchemaAttribute( "anchor_rsd", xsct_non_negative_integer, "Residue to which the stub residue should bond" )
		+ XMLSchemaAttribute( "anchor_atom", xs_string, "Atom to which the stub residue should bond" );
	XMLSchemaSimpleSubelementList subelements;
	//Append, Prepend, or Insert
	subelements
		.add_simple_subelement( "Append", subelement_attributes, "Specify a residue to append at a position" )
		.add_simple_subelement( "Prepend", subelement_attributes, "Specify a residue to prepend at a position" )
		.add_simple_subelement( "Insert", subelement_attributes, "Specify a residue to insert at a position" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Used to append/insert/prepend residues to a pose", attlist, subelements );
}

std::string PeptideStubMoverCreator::keyname() const {
	return PeptideStubMover::mover_name();
}

protocols::moves::MoverOP
PeptideStubMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PeptideStubMover );
}

void PeptideStubMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideStubMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
