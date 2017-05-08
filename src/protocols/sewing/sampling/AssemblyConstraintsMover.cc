// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AssemblyConstraintsMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/AssemblyConstraintsMover.hh>
#include <protocols/sewing/sampling/AssemblyConstraintsMoverCreator.hh>

#include <protocols/sewing/sampling/ReadNativeRotamersFileCreator.hh>
#include <protocols/sewing/sampling/ReadRepeatNativeRotamersFileCreator.hh>

// Package Headers
#include <protocols/sewing/util/io.hh>

// Protocol Headers
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <protocols/toolbox/task_operations/LinkResidues.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.AssemblyConstraintsMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////// ReadNativeRotamersFile TaskOperation  Functions   ////////////////////
////////////////////////////////////////////////////////////////////////////////////

using namespace core::pack::task::operation;
using namespace utility::tag;

ReadNativeRotamersFile::ReadNativeRotamersFile() : parent()
{}

ReadNativeRotamersFile::~ReadNativeRotamersFile() {}

TaskOperationOP
ReadNativeRotamersFileCreator::create_task_operation() const {
	return TaskOperationOP( new ReadNativeRotamersFile );
}

std::string
ReadNativeRotamersFileCreator::keyname() const
{
	return ReadNativeRotamersFile::keyname();
}

void ReadNativeRotamersFileCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReadNativeRotamersFile::provide_xml_schema( xsd );
}

core::pack::task::operation::TaskOperationOP
ReadNativeRotamersFile::clone() const {
	return core::pack::task::operation::TaskOperationOP( new ReadNativeRotamersFile( *this ) );
}

void
ReadNativeRotamersFile::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const {
	NativeRotamersMap::const_iterator map_it = nat_ro_map_.begin();
	NativeRotamersMap::const_iterator map_it_end = nat_ro_map_.end();
	for ( ; map_it != map_it_end; ++map_it ) {
		core::Size seqpos = map_it->first;
		utility::vector1< std::pair<bool, core::conformation::ResidueOP> > residues = map_it->second;
		utility::vector1<core::conformation::ResidueOP> saved_residues;
		for ( core::Size i=1; i<=residues.size(); ++i ) {
			if ( residues[i].first ) {
				saved_residues.push_back(residues[i].second);
			}
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "Adding " << saved_residues.size() << " extra rotamers to position " << seqpos << std::endl;
		}
		TR << "Adding " << saved_residues.size() << " extra rotamers to position " << seqpos << std::endl;
		core::pack::rotamer_set::AddResiduesRotamerSetOperation const nat_ro_set(saved_residues);
		core::pack::task::operation::AppendResidueRotamerSet append_rotamers_task_op(seqpos, nat_ro_set.clone());
		append_rotamers_task_op.apply(pose, task);
	}
}

void
ReadNativeRotamersFile::native_rotamers_map(
	NativeRotamersMap const & nat_ro_map
){
	nat_ro_map_ = nat_ro_map;
}

NativeRotamersMap const &
ReadNativeRotamersFile::native_rotamers_map() const{
	return nat_ro_map_;
}

void
ReadNativeRotamersFile::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("filename") ) {
		std::string filename = tag->getOption<std::string>("filename");
		nat_ro_map_ = read_native_residue_file(filename);
	} else {
		utility_exit_with_message("You must give a filename option to the ReadNativeRotamersFile TaskOperation tag");
	}
}

// AMW: surely this is a good time to make something required...
void ReadNativeRotamersFile::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "filename", xs_string , "XRW TO DO" ) );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO: not sure if this is necessary bc going to be deprecated with sewing movers" );
}

////////////////////////////////////////////////////////////////////////////////////
////////// ReadRepeatNativeRotamersFile TaskOperation  Functions   /////////////////
////////////////////////////////////////////////////////////////////////////////////

ReadRepeatNativeRotamersFile::ReadRepeatNativeRotamersFile() :
	parent()
{}

ReadRepeatNativeRotamersFile::~ReadRepeatNativeRotamersFile() {}


core::pack::task::operation::TaskOperationOP
ReadRepeatNativeRotamersFile::clone() const {
	return core::pack::task::operation::TaskOperationOP( new ReadRepeatNativeRotamersFile( *this ) );
}

///@brief apply the normal ReadNativeRotamerFile operations, then add RotamerLinks for the repeating unit
void
ReadRepeatNativeRotamersFile::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const {
	//ReadNativeRotamersFile::apply(pose, task);

	//core::Size repeat_size = native_rotamers_map().size();
	std::string first10 = pose.sequence().substr(0, 9);
	TR << "First10: " << first10 << std::endl;
	core::Size repeat_size = pose.sequence().find(first10, 10);
	TR << "Applying repeat native rotamers file for " << repeat_size << std::endl;

	//Setup rotamer links
	core::pack::rotamer_set::RotamerLinksOP links ( new core::pack::rotamer_set::RotamerLinks );
	links->resize(pose.size());
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		core::Size repeat_num = 1;
		int equiv_position = i-(int)(repeat_num*repeat_size);
		while ( equiv_position > 0 ) {
			TR << "Setting " << i << " equivalent to " << equiv_position << std::endl;
			links->set_equiv(i, equiv_position);
			++repeat_num;
			equiv_position = i-(int)(repeat_num*repeat_size);
		}

		repeat_num = 0;
		equiv_position = i+(repeat_num*repeat_size);
		while ( equiv_position <= (int)pose.size() ) {
			TR << "Setting " << i << " equivalent to " << equiv_position << std::endl;
			links->set_equiv(i, equiv_position);
			++repeat_num;
			equiv_position = i+(repeat_num*repeat_size);
		}
	}

	core::pack::task::operation::SetRotamerLinksOP rotamer_links (
		new core::pack::task::operation::SetRotamerLinks() );
	rotamer_links->set_links(links);
	rotamer_links->apply(pose, task);
}

std::string ReadRepeatNativeRotamersFile::keyname() { return "ReadRepeatNativeRotamersFile"; }

void ReadRepeatNativeRotamersFile::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() , "XRW TO DO" );
}

core::pack::task::operation::TaskOperationOP
ReadRepeatNativeRotamersFileCreator::create_task_operation() const {
	return core::pack::task::operation::TaskOperationOP( new ReadRepeatNativeRotamersFile );
}

std::string ReadRepeatNativeRotamersFileCreator::keyname() const
{ return ReadRepeatNativeRotamersFile::keyname(); }

void ReadRepeatNativeRotamersFileCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReadRepeatNativeRotamersFile::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Boiler Plate Mover Code   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AssemblyConstraintsMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new AssemblyConstraintsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AssemblyConstraintsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AssemblyConstraintsMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AssemblyConstraintsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AssemblyConstraintsMover";
// XRW TEMP }

protocols::moves::MoverOP
AssemblyConstraintsMover::clone() const {
	return( protocols::moves::MoverOP( new AssemblyConstraintsMover( *this ) ) );
}

protocols::moves::MoverOP
AssemblyConstraintsMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AssemblyConstraintsMover );
}

// XRW TEMP std::string
// XRW TEMP AssemblyConstraintsMover::get_name() const {
// XRW TEMP  return "AssemblyConstraintsMover";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Mover  Functions   ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

AssemblyConstraintsMover::AssemblyConstraintsMover():
	neighbor_cutoff_(16),
	base_native_bonus_(1.0),
	base_native_pro_bonus_(1.0) // Doonam plans to use with 2.0, but for other users, default value is 1.0
{}

AssemblyConstraintsMover::AssemblyConstraintsMover(
	NativeRotamersMap const & nat_ro_map,
	core::Real neighbor_cutoff,
	core::Real base_native_bonus,
	core::Real base_native_pro_bonus
):
	nat_ro_map_(nat_ro_map),
	neighbor_cutoff_(neighbor_cutoff),
	base_native_bonus_(base_native_bonus),
	base_native_pro_bonus_(base_native_pro_bonus)
{}

void
AssemblyConstraintsMover::nat_ro_map(
	NativeRotamersMap const & nat_ro_map
) {
	nat_ro_map_ = nat_ro_map;
}

void
AssemblyConstraintsMover::neighbor_cutoff(
	core::Real neighbor_cutoff
) {
	neighbor_cutoff_ = neighbor_cutoff;
}

void
AssemblyConstraintsMover::base_native_bonus(
	core::Real base_native_bonus
) {
	base_native_bonus_ = base_native_bonus;
}

void
AssemblyConstraintsMover::apply( core::pose::Pose & pose ) {

	using namespace core::scoring::constraints;

	pose.update_residue_neighbors();
	NativeRotamersMap::const_iterator map_it = nat_ro_map_.begin();
	NativeRotamersMap::const_iterator map_it_end = nat_ro_map_.end();
	TR << "Adding ResidueType constraints for " << nat_ro_map_.size() << " residue positions" << std::endl;
	for ( ; map_it != map_it_end; ++map_it ) {

		core::Size seqpos = map_it->first;
		if ( seqpos > pose.size() ) {
			break;
		}
		utility::vector1<std::pair<bool, core::conformation::ResidueOP> > residues = map_it->second;
		std::set<std::string> favored_types;
		for ( core::Size i=1; i<=residues.size(); ++i ) {
			if ( pose.energies().tenA_neighbor_graph().get_node(seqpos)->num_neighbors_counting_self() >= neighbor_cutoff_ ) {
				if ( favored_types.find(residues[i].second->type().name3()) != favored_types.end() ) {
					continue;
				}
				TR << "Favoring " << residues[i].second->type().name3()<< " at position " << seqpos << std::endl;
				if ( residues[i].second->type().name3() != "PRO" ) {
					ResidueTypeConstraintOP matched_nat_res_constraint ( new ResidueTypeConstraint(seqpos, residues[i].second->type().name3(), residues[i].second->type().name3(), base_native_bonus_) );
					pose.add_constraint(matched_nat_res_constraint);
					//Need to update residue neighbors every time we add a constraint. As adding a constraint clears the energies object
					pose.update_residue_neighbors();
					favored_types.insert(residues[i].second->type().name3());
				} else {
					TR << "Favoring " << residues[i].second->type().name3()<< " at position " << seqpos << " with base_native_pro_bonus_ " << base_native_pro_bonus_ << std::endl;
					ResidueTypeConstraintOP matched_nat_res_constraint ( new ResidueTypeConstraint(seqpos, residues[i].second->type().name3(), residues[i].second->type().name3(), base_native_pro_bonus_) );
					pose.add_constraint(matched_nat_res_constraint);
					//Need to update residue neighbors every time we add a constraint. As adding a constraint clears the energies object
					pose.update_residue_neighbors();
					favored_types.insert(residues[i].second->type().name3());
				}
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::sewing::repeat ].value() ) {
		apply_repeat(pose);
	}

}

void
AssemblyConstraintsMover::apply_repeat(
	core::pose::Pose & pose
) {

	/* hack to use the first 10 amino acids to define the repeat*/
	std::string first10 = pose.sequence().substr(0, 9);
	TR << "First10: " << first10 << std::endl;
	core::Size repeat_size = pose.sequence().find(first10, 10);
	if ( pose.size() % repeat_size != 0 ) {
		TR.Warning << "Pose does not contain an exact number of repeats! Repeat size: "
			<< repeat_size << ". Pose size: " << pose.size() << std::endl;
	}
	TR << "Applying NCS Constraints to a repeating unit of size " << repeat_size << std::endl;

	//Setup non-crystallographic constraints
	protocols::simple_moves::symmetry::SetupNCSMoverOP ncs_mover (
		new protocols::simple_moves::symmetry::SetupNCSMover() );

	//Create the target
	std::stringstream source;
	core::Size source_start = 2;
	core::Size source_end = repeat_size-1;
	source << source_start << "-" << source_end;

	//for(core::Size i=1; i<=num_repeats_; ++i) {
	core::Size repeat_num = 0;
	while ( true ) {
		++repeat_num;
		std::stringstream target;
		core::Size target_start = repeat_size*(repeat_num)+2;
		core::Size target_end = repeat_size*(repeat_num+1)-1;

		//If this is the last repeat, we don't have a connecting loop
		if ( target_end > pose.size() ) {
			core::Size difference = target_end - pose.size();
			target_end = target_end - difference;
			source_end = source_end - difference;
			std::stringstream source_temp;
			source_temp << source_start << "-" << source_end;
			target << target_start << "-" << target_end;
			ncs_mover->add_group(source_temp.str(), target.str());
			break;
		} else {
			target << target_start << "-" << target_end;
			ncs_mover->add_group(source.str(), target.str());
		}
	}
	TR << "Found " << repeat_num << " repeats" << std::endl;
	ncs_mover->apply(pose);

}


void
AssemblyConstraintsMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	/////////Segments////////////
	if ( tag->hasOption("native_rotamers_file") ) {
		std::string native_rotamers_file = tag->getOption<std::string>("native_rotamers_file");
		nat_ro_map_ = read_native_residue_file(native_rotamers_file);
	} else {
		utility_exit_with_message("You must specify a path to the native_rotamers_file for the AssemblyConstraintsMover");
	}
	if ( tag->hasOption("native_bonus") ) {
		base_native_bonus_ = tag->getOption<core::Real>("native_bonus");
	}
	if ( tag->hasOption("native_pro_bonus") ) { //native_proline_bonus
		base_native_pro_bonus_ = tag->getOption<core::Real>("native_pro_bonus");
	}
}

std::string AssemblyConstraintsMover::get_name() const {
	return mover_name();
}

std::string AssemblyConstraintsMover::mover_name() {
	return "AssemblyConstraintsMover";
}

void AssemblyConstraintsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	// TO DO: perhaps this is not the right function to call? -- also, delete this comment
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string AssemblyConstraintsMoverCreator::keyname() const {
	return AssemblyConstraintsMover::mover_name();
}

protocols::moves::MoverOP
AssemblyConstraintsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AssemblyConstraintsMover );
}

void AssemblyConstraintsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AssemblyConstraintsMover::provide_xml_schema( xsd );
}




} //sewing
} //protocols
