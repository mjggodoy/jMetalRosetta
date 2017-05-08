// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/RosettaScripts/util.cc
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu),
///     Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/task/PackerTask.hh>

// Project Headers
#include <core/types.hh>
#include <core/init/score_function_corrections.hh>
#include <protocols/filters/Filter.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <core/id/types.hh>
#include <core/chemical/VariantType.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sql_database/types.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.RosettaScripts.util" );

namespace protocols {
namespace rosetta_scripts {

using namespace core::scoring;
using namespace protocols::moves;
using namespace core;
using namespace std;
using utility::vector1;


/// @brief Return the number of the residue on source that is nearest to res on target. If the distance
/// is greater than 2.0 returns 0 to indicate error
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ ){
	//TR<<"looking for neiboughrs of: "<<source.pdb_info()->name()<< " and residue "<<res<<std::endl;
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	for ( core::Size i = 1; i <= source.size(); ++i ) {
		if ( source.residue( i ).is_ligand() ) continue;
		if ( chain && source.residue( i ).chain() != chain ) continue;
		core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
		if ( dist <= min_dist ) {
			min_dist = dist;
			nearest_res = i;
		}
	}
	if ( min_dist <= 3.0 ) return nearest_res;
	else return 0;
}

void
find_nearest_res(  core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size & target_res, core::Real & target_dist, core::Size const chain/*=0*/ ){
	target_res = 0; target_dist = 0.0;
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	for ( core::Size i = 1; i <= source.size(); ++i ) {
		if ( source.residue( i ).is_ligand() ) continue;
		if ( chain && source.residue( i ).chain() != chain ) continue;
		core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
		if ( dist <= min_dist ) {
			min_dist = dist;
			nearest_res = i;
		}
	}
	if ( min_dist <= 3.0 ) {
		target_res = nearest_res;
		target_dist = min_dist;
	}
}


utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/) {
	utility::vector1< core::Size > designable_vec, packable_vec, both;
	designable_vec.clear(); packable_vec.clear(); both.clear();
	core::pack::task::PackerTaskOP packer_task( tf->create_task_and_apply_taskoperations( pose ) );
	for ( core::Size resi=1; resi<=pose.size(); ++resi ) {
		if ( packer_task->being_designed( resi ) ) {
			designable_vec.push_back( resi );
		} else if ( packer_task->being_packed( resi ) ) {
			packable_vec.push_back( resi );
		}
	}
	if ( designable && packable ) {
		both.insert( both.begin(), designable_vec.begin(), designable_vec.end() );
		both.insert( both.end(), packable_vec.begin(), packable_vec.end() );
		return both;
	}
	if ( designable ) {
		return designable_vec;
	}
	return packable_vec;
}
/// @brief finds the nearest disulife to given residue on pose by searching both up and down stream to the given postion
core::Size
find_nearest_disulfide( core::pose::Pose const & pose, core::Size const res)
{
	core::Size disulfideN=0,disulfideC=0;
	for ( core::Size i = res; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
			disulfideC=i;
			break;
		}
	}
	// TR<<"C-ter disulfide: "<<disulfideC<<std::endl;
	for ( core::Size i = res ; i > 0; --i ) {
		if ( pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
			disulfideN=i;
			break;
		}
	}
	//TR<<"N-ter disulfide: "<<disulfideN<<std::endl;
	if ( (disulfideN==0)&&(disulfideC==0) ) {
		utility_exit_with_message("Could not find disulfides on: "+pose.pdb_info()->name());
	}
	if ( ((disulfideC-res)>(res-disulfideN))&&disulfideN!=0 ) {
		return disulfideN;
	}

	return disulfideC;
}


///////////////////////////////////////////////////////////
//////////////////// Task Operations //////////////////////


/// @details This is essentially a shameless copy of Justin's PackRotamersMover::parse_task_operations. In truth
/// DesignRepackMover should disappear into Justin's better organized class, but this will wait... (SJF)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP new_task_factory( new TaskFactory );
	if ( ! tag->hasOption("task_operations") ) return new_task_factory;
	TR<<"Object "<<tag->getOption< std::string >( "name", "no_name" )<<" reading the following task_operations: ";
	return( parse_task_operations( tag->getOption< std::string >( "task_operations"), data ) );
}

core::pack::task::TaskFactoryOP
parse_task_operations( std::string const & task_list, basic::datacache::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP new_task_factory( new TaskFactory );
	std::string const t_o_val( task_list );
	typedef utility::vector1< std::string > StringVec;
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations\n";
	for ( auto const & t_o_key : t_o_keys ) {
		if ( data.has( "task_operations", t_o_key ) ) {
			new_task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR<<t_o_key<<' ';
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	TR<<std::endl;
	return new_task_factory;
}

///option to add or refer to a Taskfactory through the datamap, similar to how to add/refer to movemap OPs (EMS)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap & data, core::pack::task::TaskFactoryOP & task_factory /*, bool const reset_taskfactory */)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	if ( tag->hasOption("task_factory" ) ) {
		std::string const name( tag->getOption< std::string >("task_factory") );
		TR <<"taskfacotry name: " << name << std::endl;

		if ( data.has( "TaskFactory", name ) ) {
			task_factory = data.get_ptr<TaskFactory>( "TaskFactory", name );
			TR<<"found helper task_factory: "<< name <<" for mover: "<<tag->getName()<< std::endl;
		} else { // ( !data.has( "TaskFactory", name ) ){
			std::string tf_string = "TaskFactory";
			task_factory = core::pack::task::TaskFactoryOP( new TaskFactory );
			data.add( tf_string , name , task_factory );
			TR<<"adding new TaskFactory to the datamap: "<< name  <<std::endl;
		}
	}

	if ( ! tag->hasOption("task_operations") ) return task_factory;

	std::string const t_o_val( tag->getOption<std::string>("task_operations") );
	typedef utility::vector1< std::string > StringVec;
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations to mover "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";

	for ( auto const & t_o_key : t_o_keys ) {

		if ( data.has( "task_operations", t_o_key ) ) {
			task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR<<t_o_key<<' ';
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	TR<<std::endl;
	return task_factory;
}

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	using core::pack::task::operation::TaskOperationOP;
	using core::pack::task::operation::TaskOperation;
	typedef std::string String;

	utility::vector1< TaskOperationOP > task_operations;
	String const t_o_val( tag->getOption<String>("task_operations", "" ) );
	if ( t_o_val != "" ) {
		utility::vector1< String > const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( auto const & t_o_key : t_o_keys ) {
			if ( data.has( "task_operations", t_o_key ) ) {
				task_operations.push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
			}
		}
	}
	return task_operations;
}

////////////////////////////////////////////////////////////////////
///  Get attributes ( i.e. options ) for movers to build xml schemas

/// @brief Appends the attributes read by parse_task_operations
void
attributes_for_parse_task_operations( utility::tag::AttributeList & attributes, std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute( "task_operations", utility::tag::xs_string,
		( description == ""  ? "A comma separated list of TaskOperations to use." : description )  );
}


/// @brief Appends the attributes read by parse_task_operations when handed a TaskFactory
void
attributes_for_parse_task_operations_w_factory( utility::tag::AttributeList & attributes, std::string const & used_for_descr )
{
	std::string to_descrip = "A comma separated list of TaskOperations to use";
	std::string tf_descrip = "A TaskFactory specification to use";
	if ( ! used_for_descr.empty() ) {
		to_descrip += " for " + used_for_descr;
		tf_descrip += " for " + used_for_descr;
	}

	attributes
		+ utility::tag::XMLSchemaAttribute( "task_operations", utility::tag::xs_string, to_descrip )
		+ utility::tag::XMLSchemaAttribute( "task_factory", utility::tag::xs_string , tf_descrip );
}

/////////////////////////////////////////////////////////////
//////////////////// Residue Selectors //////////////////////

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" option in tag
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
core::select::residue_selector::ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name )
{
	return core::select::residue_selector::parse_residue_selector( tag, data, option_name );
}

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data )
{
	return core::select::residue_selector::get_residue_selector( selector_name, data );
}


/// @brief Appends the attributes read by parse_residue_selector
void
attributes_for_parse_residue_selector( utility::tag::AttributeList & attributes, std::string const & description )
{
	return core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name(attributes, description );
}

///////////////////////////////////////////////////////////
//////////////////// ScoreFunction ////////////////////////

/// @details Utility function to find a scorefunction from
/// parser-provided data. This is essentially a shameless copy of
/// Justin's PackRotamersMover::parse_score_function.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	std::string const scorefxn_key( tag->getOption<std::string>(option_name, dflt_key) );
	if ( ! data.has( "scorefxns", scorefxn_key ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap. To add a score function to the data map, define a score function in the <SCOREFXNS/>.'");
	}

	try{
		core::init::check_score_function_sanity(scorefxn_key, false);
	} catch(utility::excn::EXCN_BadInput & e){
		stringstream err_msg;
		err_msg
			<< std::endl << std::endl
			<< "Error getting score function from";
		if ( tag->hasOption("name") ) {
			err_msg
				<< " the tag '" << tag->getOption<std::string>("name") << "' of type ";
		} else {
			err_msg
				<< " a tag of type ";
		}
		err_msg
			<< "'" << tag->getName() << "':" << std::endl;

		if ( !tag->hasOption(option_name) ) {
			err_msg
				<< "  No value for option '" << option_name << "' is provided, so it is using the default value of '" << dflt_key << "'" << std::endl;
		} else {
			err_msg
				<< "  The value for option '" << option_name << "' = '" << scorefxn_key << "'" << std::endl;
		}
		err_msg
			<< "ERROR MESSAGE: " << e << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption(err_msg.str());
	}
	return data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_key );
}



/// @details Utility function to find a scorefunction from
/// parser-provided data for the option 'scorefxn'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	return parse_score_function(tag, "scorefxn", data, dflt_key);
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name
) {
	return tag->getOption<std::string>(option_name, "commandline");
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag
) {
	return get_score_function_name(tag, "scorefxn");
}

///  Get attributes ( i.e. options ) for movers to build xml schemas
void
attributes_for_get_score_function_name( utility::tag::AttributeList & attributes )
{
	attributes_for_get_score_function_name_w_description( attributes, "" );
}

void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes,
	std::string const & option_name )
{
	attributes_for_get_score_function_name_w_description( attributes, option_name, "" );
}

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_get_score_function_name_w_description( attributes, "scorefxn", description);
}

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute( option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes )
{
	attributes_for_parse_score_function_w_description( attributes, "" );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes, std::string const & sfxn_option_name )
{
	attributes_for_parse_score_function_w_description( attributes, sfxn_option_name, "" );
}

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_parse_score_function_w_description( attributes, "scorefxn", description );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute(
		sfxn_option_name, utility::tag::xs_string ,
		( description == "" ? "Name of score function to use" : description ) );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument and description.
/// @details This version appends the attributes as required attributes.
/// @author Vikram K. Mulligan.
void
attributes_for_parse_score_function_w_description_when_required( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description
) {
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		sfxn_option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}


///////////////////////////////////////////////////////////
//////////////////// Reference Pose ///////////////////////

void
attributes_for_saved_reference_pose(
	utility::tag::AttributeList & attributes,
	std::string const & attribute_name)
{
	attributes_for_saved_reference_pose_w_description( attributes, "", attribute_name );
}

void
attributes_for_saved_reference_pose_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description,
	std::string const & attribute_name)
{
	attributes + utility::tag::XMLSchemaAttribute( attribute_name, utility::tag::xs_string,
		( description == "" ? "Name of reference pose to use" : description ) );
}

core::pose::PoseOP
saved_reference_pose( utility::tag::TagCOP const in_tag, basic::datacache::DataMap & data_map, std::string const & tag_name ){

	if ( in_tag->hasOption(tag_name) ) {
		core::pose::PoseOP refpose(nullptr);
		std::string refpose_name(in_tag->getOption<std::string>( tag_name) );
		TR<<"Loading PDB: "<<refpose_name<<std::endl;

		if ( !data_map.has("spm_ref_poses",refpose_name) ) {
			refpose = core::pose::PoseOP( new core::pose::Pose() );
			data_map.add("spm_ref_poses",refpose_name,refpose );
		} else refpose = data_map.get_ptr<core::pose::Pose>("spm_ref_poses",refpose_name );

		return refpose;
	} else std::cerr << "WARNING: saved_reference_pose function called even though tag has no " + tag_name + " entry. something's unclean somewhere." << std::endl;
	return nullptr;
}


/////////////////////////////////////////////////////////
//////////////////// MoveMap ////////////////////////////

/// @brief utility function for parse_movemap which goes over each of the tags in a movemap section
void
foreach_movemap_tag(
	utility::tag::TagCOP const in_tag,
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP mm
){
	using namespace core::kinematics;
	using namespace utility::tag;

	for ( TagCOP const tag : in_tag->getTags() ) {
		std::string const name( tag->getName() );
		runtime_assert( name == "Jump" || name == "Chain" || name == "Span" );
		if ( name == "Jump" ) {
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const setting( tag->getOption< bool >( "setting" ) );
			if ( num == 0 ) mm->set_jump( setting ); // set all jumps if number==0
			else mm->set_jump( num, setting );
		}
		if ( name == "Chain" ) {
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			core::Size const chain_begin( pose.conformation().chain_begin( num ) );
			core::Size const chain_end( pose.conformation().chain_end( num ) );
			for ( core::Size i( chain_begin ); i <= chain_end; ++i ) {
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
			bool const bondangle( tag->getOption< bool >( "bondangle", false ) );
			bool const bondlength( tag->getOption< bool >( "bondlength", false ) );
			if ( bondangle || bondlength ) {
				for ( core::Size i( chain_begin ); i <= chain_end; ++i ) {
					for ( core::Size j=1; j<=pose.residue_type(i).natoms(); ++j ) {
						mm->set( core::id::DOF_ID(core::id::AtomID(j,i), core::id::THETA ), bondangle );
						mm->set( core::id::DOF_ID(core::id::AtomID(j,i), core::id::D ), bondlength );
					}
				}
			}
		}
		if ( name == "Span" ) {
			core::Size const begin( tag->getOption< core::Size >( "begin" ) );
			core::Size const end( tag->getOption< core::Size >( "end" ) );
			runtime_assert( end >= begin );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			for ( core::Size i( begin ); i <= end; ++i ) {
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
			bool const bondangle( tag->getOption< bool >( "bondangle", false ) );
			bool const bondlength( tag->getOption< bool >( "bondlength", false ) );
			if ( bondangle || bondlength ) {
				for ( core::Size i( begin ); i <= end; ++i ) {
					for ( core::Size j=1; j<=pose.residue_type(i).natoms(); ++j ) {
						mm->set( core::id::DOF_ID(core::id::AtomID(j,i), core::id::THETA ), bondangle );
						mm->set( core::id::DOF_ID(core::id::AtomID(j,i), core::id::D ), bondlength );
					}
				}
			}
		}
	}//for tag
}

void
parse_movemap_tag(
	utility::tag::TagCOP const in_tag,
	core::kinematics::MoveMapOP mm )
{
	if ( in_tag->hasOption( "bb" ) ) mm->set_bb( in_tag->getOption< bool >( "bb" ) );
	if ( in_tag->hasOption( "chi" ) ) mm->set_chi( in_tag->getOption< bool >( "chi" ) );
	if ( in_tag->hasOption( "jump" ) ) mm->set_jump( in_tag->getOption< bool >( "jump" ) );
}

// @brief variant of parse_movemap that takes in a datamap and searches it for already existing movemaps
void
parse_movemap(
	utility::tag::TagCOP const in_tag,
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP & mm,
	basic::datacache::DataMap & data,
	bool const reset_map
){
	using utility::tag::TagCOP;
	using namespace core::kinematics;

	if ( in_tag == nullptr ) return;
	utility::vector1< TagCOP > const branch_tags( in_tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() =="MoveMap" ) {
			break;
		}
	}
	if ( reset_map ) {
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}
	if ( tag_it == branch_tags.end() ) return;

	if ( (*tag_it)->hasOption("name") ) {
		std::string const name( (*tag_it)->getOption< std::string >( "name" ) );
		if ( data.has( "movemaps", name ) ) {
			mm = data.get_ptr<MoveMap>( "movemaps", name );
			TR<<"Found movemap "<<name<<" on datamap"<<std::endl;
		} else {
			data.add( "movemaps", name, mm );
			TR<<"Adding movemap "<<name<<" to datamap"<<std::endl;
		}
	}
	parse_movemap_tag( *tag_it, mm );
	foreach_movemap_tag( *tag_it, pose, mm );
}

/// @details modifies an existing movemap according to tag
/// the movemap defaults to move all bb, chi, and jumps.
void
parse_movemap(
	utility::tag::TagCOP const in_tag,
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP mm,
	bool const reset_movemap)
{
	using utility::tag::TagCOP;
	using namespace core::kinematics;

	if ( in_tag == nullptr ) return;
	utility::vector1< TagCOP > const branch_tags( in_tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() =="MoveMap" ) {
			break;
		}
	}
	if ( reset_movemap ) {
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}
	if ( tag_it == branch_tags.end() ) return;
	if ( (*tag_it)->hasOption( "name" ) ) {
		TR<<"ERROR in "<<*tag_it<<'\n';
		throw utility::excn::EXCN_RosettaScriptsOption("Tag called with option name but this option is not available to this mover. Note that FastRelax cannot work with a prespecified movemap b/c its movemap is parsed at apply time. Sorry." );
	}

	parse_movemap_tag( *tag_it, mm );
	foreach_movemap_tag( *tag_it, pose, mm );
}

std::string move_map_tag_namer( std::string const & subelement_name ) { return "move_map_" + subelement_name + "_type"; }
std::string unnamed_move_map_ct_namer( std::string const & ) { return "unnamed_move_map_type"; }
std::string optionally_named_move_map_ct_namer( std::string const & ) { return "optionally_named_move_map_type"; }

void
common_movemap_complex_type_def( utility::tag::XMLSchemaComplexTypeGenerator & ct_gen )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList movemap_subelements;
	AttributeList jump_attributes;
	jump_attributes
		+ XMLSchemaAttribute::required_attribute( "number",  xsct_non_negative_integer , "Which jump number (in the FoldTree)" )
		+ XMLSchemaAttribute::required_attribute( "setting", xsct_rosetta_bool , "true for move, false for don't move" );
	movemap_subelements.add_simple_subelement( "Jump", jump_attributes , "jumps are the not-chemistry internal coordinate connections between separate parts of your pose");

	std::string const chi_desc("move sidechain chi torsions?");
	std::string const bb_desc("move backbone torsions?");

	AttributeList chain_attributes;
	chain_attributes
		+ XMLSchemaAttribute::required_attribute( "number", xsct_non_negative_integer , "which chain?" )
		+ XMLSchemaAttribute::required_attribute( "chi",    xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute::required_attribute( "bb",     xsct_rosetta_bool , bb_desc );
	movemap_subelements.add_simple_subelement( "Chain", chain_attributes , "this controls a kinematically contiguous chain (think protein chains)");

	AttributeList span_attributes;
	span_attributes
		+ XMLSchemaAttribute::required_attribute( "begin",     xsct_non_negative_integer , "beginning of span" )
		+ XMLSchemaAttribute::required_attribute( "end",       xsct_non_negative_integer , "end of span" )
		+ XMLSchemaAttribute::required_attribute( "chi",       xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute::required_attribute( "bb",        xsct_rosetta_bool , bb_desc )
		+ XMLSchemaAttribute( "bondangle", xsct_rosetta_bool , "move 3-body angles?" )
		+ XMLSchemaAttribute( "bondlength", xsct_rosetta_bool , "move 2-body lengths?" );
	movemap_subelements.add_simple_subelement( "Span", span_attributes , "XRW TO DO, probably a user-defined region of the Pose");
	movemap_subelements.complex_type_naming_func( & move_map_tag_namer );


	AttributeList movemap_tag_attributes;
	movemap_tag_attributes
		+ XMLSchemaAttribute( "bb", xsct_rosetta_bool , bb_desc )
		+ XMLSchemaAttribute( "chi", xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute( "jump", xsct_rosetta_bool , "move all jumps?" );

	ct_gen.element_name( "MoveMap" )
		.description( "MoveMaps dicate what degrees of freedom are mobile to other bits of code, especially minimization - they do NOT work with packing" )
		.add_attributes( movemap_tag_attributes )
		.set_subelements_repeatable( movemap_subelements );

}

/// @brief Edits the complex type for an object that parses a MoveMap subtag
/// and which calls parse_movemap (without passing a DataMap parameter)
void
append_subelement_for_parse_movemap(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::XMLSchemaSimpleSubelementList & subelements,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGenerator unnamed_move_map;
	common_movemap_complex_type_def( unnamed_move_map );
	unnamed_move_map
		.description( description == "" ? "MoveMap specification" : description )
		.complex_type_naming_func( & unnamed_move_map_ct_namer )
		.write_complex_type_to_schema( xsd );

	subelements.add_already_defined_subelement( "MoveMap", & unnamed_move_map_ct_namer );
}


void
append_subelement_for_parse_movemap_w_datamap(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::XMLSchemaSimpleSubelementList & subelements,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGenerator named_move_map;
	common_movemap_complex_type_def( named_move_map );
	named_move_map
		.description( description == "" ? "MoveMap specification" : description )
		.add_optional_name_attribute()
		.complex_type_naming_func( & optionally_named_move_map_ct_namer )
		.write_complex_type_to_schema( xsd );

	subelements.add_already_defined_subelement( "MoveMap", & optionally_named_move_map_ct_namer );
}


void
add_movemaps_to_datamap(
	utility::tag::TagCOP in_tag,
	core::pose::Pose const & pose,
	basic::datacache::DataMap & data,
	bool initialize_mm_as_true)
{
	using utility::tag::TagCOP;
	using namespace core::kinematics;

	if ( in_tag == nullptr ) return;
	utility::vector1< TagCOP > const branch_tags( in_tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for ( auto const & branch_tag : branch_tags ) {
		if ( branch_tag->getName() =="MoveMap" ) {

			if ( ! branch_tag->hasOption("name") ) continue;
			std::string const name( branch_tag->getOption< std::string >( "name" ) );
			if ( data.has("movemaps", name) ) continue;


			MoveMapOP mm( new MoveMap() );
			if ( initialize_mm_as_true ) {
				mm->set_bb( true );
				mm->set_chi( true );
				mm->set_jump( true );
			}
			foreach_movemap_tag( branch_tag, pose, mm );
			data.add("movemaps", name, mm);
		}

	}
}


/////////////////////////////////////////////////////////
//////////////////// Filter ////////////////////////////

protocols::filters::FilterOP
parse_filter( std::string const & filter_name, protocols::filters::Filters_map const & filters ){
	auto filter_it( filters.find( filter_name ) );
	if ( filter_it == filters.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Filter "+filter_name+" not found" );
	}
	return filter_it->second;
}



/////////////////////////////////////////////////////////
//////////////////// Mover //////////////////////////////

protocols::moves::MoverOP
parse_mover( std::string const & mover_name, protocols::moves::Movers_map const & movers ){
	auto mover_it( movers.find( mover_name ) );
	if ( mover_it == movers.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Mover "+mover_name+" not found" );
	}
	return mover_it->second;
}



/////////////////////////////////////////////////////////
//////////////////// XYZ Vector /////////////////////////

void
attributes_for_parse_xyz_vector( utility::tag::AttributeList & attlist ){
	using namespace utility::tag;
	attlist
		+ XMLSchemaAttribute::required_attribute( "x", xsct_real, "X coordinate for this vector" )
		+ XMLSchemaAttribute::required_attribute( "y", xsct_real, "Y coordinate for this vector" )
		+ XMLSchemaAttribute::required_attribute( "z", xsct_real, "Z coordinate for this vector" );
}


/// @brief utility function for parsing xyzVector
numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagCOP const xyz_vector_tag ){
	if ( ! xyz_vector_tag->hasOption("x") ) throw utility::excn::EXCN_RosettaScriptsOption("xyz_vector requires 'x' coordinates option");
	if ( ! xyz_vector_tag->hasOption("y") ) throw utility::excn::EXCN_RosettaScriptsOption("xyz_vector requires 'y' coordinates option");
	if ( ! xyz_vector_tag->hasOption("z") ) throw utility::excn::EXCN_RosettaScriptsOption("xyz_vector requires 'z' coordinates option");

	numeric::xyzVector< core::Real > xyz_v (
		xyz_vector_tag->getOption<core::Real>("x"),
		xyz_vector_tag->getOption<core::Real>("y"),
		xyz_vector_tag->getOption<core::Real>("z")
	);

	return xyz_v;

}


// void
// attributes_for_report_to_db( utility::tag::AttributeList & attlist, utility::tag::XMLSchemaDefinition & xsd)
// {
//  using namespace utility::tag;
//
//  XMLSchemaRestriction database_mode;
//  database_mode.name( "database_mode_string" );
//  database_mode.base_type( xs_string );
//  database_mode.add_restriction( xsr_enumeration, "sqlite3" );
//  database_mode.add_restriction( xsr_enumeration, "postgres" );
//  database_mode.add_restriction( xsr_enumeration, "mysql" );
//  xsd.add_top_level_element( database_mode );
//  XMLSchemaRestriction database_transaction_mode;
//  database_transaction_mode.name( "database_transaction_mode_string" );
//  database_transaction_mode.base_type( xs_string );
//  database_transaction_mode.add_restriction( xsr_enumeration, "none" );
//  database_transaction_mode.add_restriction( xsr_enumeration, "standard" );
//  database_transaction_mode.add_restriction( xsr_enumeration, "chunk" );
//  xsd.add_top_level_element( database_transaction_mode );
//
//
//
//  XMLSchemaRestriction mode;
//  mode.name( "report_to_db_relevant_residues_mode" );
//  mode.base_type( xs_string );
//  //can be "explicit" or "implicit" with any capitalization
//  mode.add_restriction( xsr_pattern, "[eEiI][xXmM][pP][lL][iI][cC][iI][tT]" );
//  xsd.add_top_level_element( mode );
//
//
//  attlist
//   + XMLSchemaAttribute( "resource_description", xs_string, "Resource description referring to the resource to be output" )
//   + XMLSchemaAttribute( "database_resource", xs_string, "Resource description referring to the output database" )
//   + XMLSchemaAttribute( "database_resource_tag", xs_string, "Tag referring to the output database" )
//   + XMLSchemaAttribute::attribute_w_default( "transaction_mode", "database_transaction_mode_string", "Transaction mode for database output", "standard" )
//   + XMLSchemaAttribute( "database_mode", "database_mode_string", "Which type of output database to use?" )
//   + XMLSchemaAttribute( "database_name", xs_string, "Name of output database" )
//   + XMLSchemaAttribute( "database_pq_schema", xs_string, "Schema name within the database" )
//   + XMLSchemaAttribute( "database_separate_db_per_mpi_process", xsct_rosetta_bool, "Use a separate database for each MPI process? Specific to sqlite3" )
//   + XMLSchemaAttribute( "database_partition", xs_integer, "Database partition to use" )
//   + XMLSchemaAttribute( "database_read_only", xsct_rosetta_bool, "Is the database read-only?" )
//   + XMLSchemaAttribute( "database_host", xs_string, "URI to the database server (postgres only)" )
//   + XMLSchemaAttribute( "database_user", xs_string, "Username for database access( postgres only)" )
//   + XMLSchemaAttribute( "database_password", xs_string, "Password for database access (postgres only)" )
//   + XMLSchemaAttribute( "database_port", xsct_non_negative_integer, "Port to use for access to database server (postgres only)" )
//   + XMLSchemaAttribute( "batch_description", xs_string, "Description of features database" )
//   + XMLSchemaAttribute( "protocol_id", xsct_non_negative_integer, "Manually controls protocol_id associated with this ReportToDB tag. Autoincrements by default." )
//   + XMLSchemaAttribute( "batch_id", xsct_non_negative_integer, "Manually controls the batch_id associated with this ReportToDB tag. Autoincrements by default." )
//   + XMLSchemaAttribute( "use_transactions", xsct_rosetta_bool, "Use transactions to group database i/o to be more efficient. Turning them off can help debugging." )
//   + XMLSchemaAttribute::attribute_w_default( "cache_size", xsct_non_negative_integer, "Specify the maximum number 1k pages to keep in memory before writing to disk.", "2000" )
//   + XMLSchemaAttribute::attribute_w_default( "remove_xray_virt", xsct_rosetta_bool, "Remove virtual residue attached during xray refine process", "false" )
//   + XMLSchemaAttribute::attribute_w_default( "relevant_residues_mode", "report_to_db_relevant_residues_mode", "Determine what features are reported given the relevant residues", "explicit" );
//
//  rosetta_scripts::attributes_for_parse_task_operations( attlist );
// }

} //RosettaScripts
} //protocols
