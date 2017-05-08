// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/AddChainMover.hh>
#include <protocols/simple_moves/AddChainMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>


#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/xyzVector.hh>
#include <protocols/toolbox/superimpose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.AddChainMover" );

namespace protocols {
namespace simple_moves {

using basic::T;
using basic::Error;
using basic::Warning;

// XRW TEMP std::string
// XRW TEMP AddChainMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddChainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddChainMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddChainMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddChainMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddChain";
// XRW TEMP }

AddChainMover::AddChainMover()
: moves::Mover("AddChain"),
	fname_( "" ),
	new_chain_( true ),
	scorefxn_( /* NULL */ )
{
	random_access_ = false;
	swap_chain_number_ = 0;
}

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size const pos : positions ) {
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

void AddChainMover::add_new_chain( core::pose::Pose & pose ) const {// pose is passed by reference. So function does not return anything but modifies pose, and leave it modified. The function is const as it does not modify private member data.
	if ( swap_chain_number()!=0 ) {
		return;
	}
	TR<<"AddChainMover is adding a new chain to pose: "<<std::endl;

	using namespace core::pose;

	Pose new_pose;

	utility::vector1< std::string > const split_names( utility::string_split< std::string >( fname(), ',', std::string()) );
	TR<<"Found "<<split_names.size()<<" file names"<<std::endl;
	core::Size const random_num = (core::Size) (numeric::random::rg().uniform() * split_names.size()) + 1;
	std::string const curr_fname( split_names[ random_num ] );
	TR<<"choosing number: "<<random_num<<" "<<curr_fname<<std::endl;

	TR<<"Before addchain, total residues: "<<pose.size()<<std::endl;
	core::import_pose::pose_from_file( new_pose, curr_fname , core::import_pose::PDB_file);
	new_pose.conformation().detect_disulfides();
	(*scorefxn()) ( new_pose );

	append_pose_to_pose( pose, new_pose, new_chain() );
	pose.conformation().detect_disulfides();
	pose.update_residue_neighbors();
	pose.update_pose_chains_from_pdb_chains();
	(*scorefxn())( pose );
	pose.pdb_info( PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) ); //reinitialize the PDBInfo
	TR<<"After addchain, total residues: "<<pose.size()<<std::endl;

	core::pose::add_comment(pose,"AddedChainName ",curr_fname);
}

void AddChainMover::swap_chain( core::pose::Pose & pose ) const {
	if ( swap_chain_number()==0 ) {
		return;
	}
	TR<<"AddChainMover will swap chain: "<<swap_chain_number()<<std::endl;

	using namespace core::pose;
	using namespace protocols::toolbox;


	Pose new_chain;

	utility::vector1< std::string > const split_names( utility::string_split< std::string >( fname(), ',', std::string()) );
	TR<<"Found "<<split_names.size()<<" file names"<<std::endl;
	core::Size const random_num = (core::Size) (numeric::random::rg().uniform() * split_names.size()) + 1;
	std::string const curr_fname( split_names[ random_num ] );
	TR<<"choosing number: "<<random_num<<" "<<curr_fname<<std::endl;

	TR<<"Before addchain, total residues: "<<pose.size()<<std::endl;
	core::import_pose::pose_from_file( new_chain, curr_fname , core::import_pose::PDB_file);
	new_chain.conformation().detect_disulfides();
	(*scorefxn()) ( new_chain );

	// Here we have the new chain in new_chain. This needs to be aligned to chain 2 in pose, and current chain 2 should be deleted.
	utility::vector1< PoseOP > pose_chains( pose.split_by_chain() ); // splits pose into a vector of poses, with one chain per pose

	Pose template_pose;
	append_pose_to_pose( template_pose, *pose_chains[ swap_chain_number() ], true ); // the template chain is the chain that will be swapped

	// Now I should align new_chain to template_pose
	utility::vector1< core::Size > template_positions;
	for ( core::Size i = 1; i <= template_pose.size(); ++i ) {
		template_positions.push_back(i);
	}

	utility::vector1< core::Size > new_chain_positions;
	for ( core::Size i = 1; i <= new_chain.size(); ++i ) {
		new_chain_positions.push_back(i);
	}

	runtime_assert( new_chain_positions == template_positions ); // The two residue number vectors should be identical

	utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coords( new_chain, new_chain_positions ) ), ref_coords( Ca_coords( template_pose, template_positions ) );
	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;
	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );
	apply_superposition_transform( new_chain, rotation, to_init_center, to_fit_center );

	// new_chain is aligned to template_pose. Now everything that is not swap_chain_number in pose_chains should be combined with new_chain in pose.
	Pose new_pose;
	for ( core::Size i = 1; i<swap_chain_number(); ++i ) {
		append_pose_to_pose( new_pose, *pose_chains[i], true ); // true here means that it will be appended as new chain. I'm using a * here to dereference the OP.
	}
	append_pose_to_pose( new_pose, new_chain, true ); // Why new_chain()???
	for ( core::Size i = swap_chain_number() + 1; i <= pose_chains.size(); ++i ) {
		append_pose_to_pose( new_pose, *pose_chains[i], true );
	}

	// Now add the comments from pose to new_pose, and while doing this update the AddedChainName to the new chain name.
	std::map< std::string, std::string > const comments = core::pose::get_all_comments( pose );
	std::map< std::string, std::string>::const_iterator it;

	for ( it = comments.begin(); it != comments.end(); ++it ) {
		if ( it->first == "AddedChainName " ) { // update AddedChain name
			core::pose::add_comment(new_pose,"AddedChainName ",curr_fname);
		} else {
			core::pose::add_comment(new_pose,it->first,it->second);
		}
	}

	// Assign new_pose to pose and update score/disulfides/neighbors.
	pose = new_pose;
	pose.conformation().detect_disulfides();
	pose.update_residue_neighbors();
	(*scorefxn())( pose );
	pose.pdb_info( PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) ); //reinitialize the PDBInfo
	TR<<"After addchain, total residues: "<<pose.size()<<std::endl;
}

void
AddChainMover::apply( Pose & pose )
{
	//runtime_assert( !new_chain() != (swap_chain_number()!=0) ); // You cannot do both new chain and swap chain!
	//runtime_assert( ( new_chain() && swap_chain_number()==0) || (!new_chain() && swap_chain_number()!=0)); // You cannot do both new chain and swap chain and you have to do one of the two...
	add_new_chain(pose);
	swap_chain(pose);


}

// XRW TEMP std::string
// XRW TEMP AddChainMover::get_name() const {
// XRW TEMP  return AddChainMover::mover_name();
// XRW TEMP }

moves::MoverOP
AddChainMover::clone() const
{
	return moves::MoverOP( new AddChainMover( *this ) );
}

moves::MoverOP
AddChainMover::fresh_instance() const
{
	return moves::MoverOP( new AddChainMover );
}

void
AddChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	random_access( tag->getOption< bool >( "random_access", false ) );
	fname( tag->getOption< std::string >( "file_name" ) );
	if ( random_access() ) {
		utility::vector1< std::string > const split_names( utility::string_split< std::string >( fname(), ',', std::string()) );
		TR<<"Found "<<split_names.size()<<" file names"<<std::endl;
		//  core::Size const random_num = (core::Size) (numeric::random::rg().uniform() * split_names.size()) + 1;
		//  TR<<"choosing number: "<<random_num<<" "<<split_names[ random_num ]<<std::endl;
		//  fname( split_names[ random_num ] );
	}
	new_chain( tag->getOption< bool >( "new_chain", true ) );
	swap_chain_number( tag->getOption< core::Size >( "swap_chain_number", 0 ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	TR<<"AddChain sets fname: "<<fname()<<" new_chain: "<<new_chain()<<std::endl;
}

void AddChainMover::scorefxn( core::scoring::ScoreFunctionOP s ){ scorefxn_ = s; }
core::scoring::ScoreFunctionOP AddChainMover::scorefxn() const{ return scorefxn_; }

std::string AddChainMover::get_name() const {
	return mover_name();
}

std::string AddChainMover::mover_name() {
	return "AddChain";
}

void AddChainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	std::string const chain_warning("There is some interaction between swap_chain_number and new_chain; probably you can use only one.  ");

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "random_access", xsct_rosetta_bool, "if true randomly choose one file name from a list and work with that throughout the run.", "false" )
		+ XMLSchemaAttribute::required_attribute( "file_name", xs_string, "Either a path to the file to read chains from, or a comma-separated list of such if random_access is true." )
		+ XMLSchemaAttribute::attribute_w_default( "new_chain", xsct_rosetta_bool, chain_warning + "add as a new chain?", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "swap_chain_number", xsct_non_negative_integer, chain_warning + "swap chain with specified chain number", "0" ); //0 is valid as a flag value, so it really is a default and not required

	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Replace or add chains to a pose from other PDBs", attlist );
}

std::string AddChainMoverCreator::keyname() const {
	return AddChainMover::mover_name();
}

protocols::moves::MoverOP
AddChainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddChainMover );
}

void AddChainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddChainMover::provide_xml_schema( xsd );
}

} // simple_moves
} // protocols
