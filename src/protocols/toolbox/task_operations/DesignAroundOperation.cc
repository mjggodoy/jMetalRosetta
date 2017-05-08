// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/DesignAroundOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperationCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

// Project Headers
#include <core/pose/Pose.hh>


// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
// Auto-header: duplicate removed #include <core/pack/task/operation/TaskOperations.hh>


// C++ Headers
#include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.DesignAroundOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
DesignAroundOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DesignAroundOperation );
}

void DesignAroundOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DesignAroundOperation::provide_xml_schema( xsd );
}

std::string DesignAroundOperationCreator::keyname() const
{
	return DesignAroundOperation::keyname();
}

DesignAroundOperation::DesignAroundOperation() :
	design_shell_( 8.0 ),
	repack_shell_( 8.0 ),
	allow_design_( true ),
	resnums_allow_design_( true ),
	string_resnums_( "" )
{
	resid_.clear();
}

DesignAroundOperation::~DesignAroundOperation() {}

core::pack::task::operation::TaskOperationOP DesignAroundOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DesignAroundOperation( *this ) );
}

/// @brief restricts to repacking all residues outside of design_shell_ around each residue
void
DesignAroundOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task::operation;


	runtime_assert( repack_shell() >= design_shell() );
	set< core::Size > focus_residues;// all of the residues that were input (notice that the method is const, so I can't change resid_)
	focus_residues.clear();
	focus_residues.insert( resid_.begin(), resid_.end() );
	set< core::Size > const res_vec( core::pose::get_resnum_list( string_resnums_, pose ) );
	focus_residues.insert( res_vec.begin(), res_vec.end() );

	utility::vector1< core::Size > packing_residues, prevent_repacking_residues;
	packing_residues.clear(); prevent_repacking_residues.clear();
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		bool allow_design_res( false );
		bool allow_packing( false );
		for ( core::Size const res : focus_residues ) { // don't change anything for focus residues
			if ( i == res ) {
				if ( resnums_allow_design() ) allow_design_res = true;
				break;
			}
		}
		//check if design res nbr
		if ( allow_design() && !allow_design_res ) {
			for ( core::Size const res : focus_residues ) {
				if ( i == res ) continue; //dont need to check if is nbr of self
				core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
				if ( distance <= design_shell() || distance <= 0.0001 /*if design_shell is specified as 0 ensure that focus residues are allowed to design*/ ) {
					allow_design_res = true;
					break;
				}// fi distance
			} //foreach res
		}
		if ( allow_design_res ) continue;
		for ( core::Size const res : focus_residues ) {
			core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
			if ( distance <= repack_shell() ) {
				allow_packing = true;
				break;
			}// fi distance
		} //foreach res
		if ( allow_packing ) packing_residues.push_back( i );
		else prevent_repacking_residues.push_back( i );
	}//for i
	///for some unfathomable reason OperateOnCertainResidues defaults to applying to all residues if none are defined, so you have to be careful here...
	OperateOnCertainResidues oocr_repacking, oocr_prevent_repacking;
	if ( packing_residues.size() ) {
		oocr_repacking.op( ResLvlTaskOperationCOP( new RestrictToRepackingRLT ) );
		oocr_repacking.residue_indices( packing_residues );
		oocr_repacking.apply( pose, task );
	}
	if ( prevent_repacking_residues.size() ) {
		oocr_prevent_repacking.op( ResLvlTaskOperationCOP( new PreventRepackingRLT ) );
		oocr_prevent_repacking.residue_indices( prevent_repacking_residues );
		oocr_prevent_repacking.apply( pose, task );
	}
}

void
DesignAroundOperation::design_shell( core::Real const radius )
{
	design_shell_ = radius;
	if ( radius >= repack_shell() ) {
		repack_shell_ = radius;
	}
}

void
DesignAroundOperation::include_residue( core::Size const resid )
{
	resid_.insert( resid );
}

void
DesignAroundOperation::parse_tag( TagCOP tag , DataMap & )
{
	string_resnums_ = tag->getOption< std::string >( "resnums" );// these are kept in memory until the pose is available (at apply time)
	design_shell( tag->getOption< core::Real >( "design_shell", 8.0 ) );
	allow_design( tag->getOption< bool >( "allow_design", 1 ) );
	resnums_allow_design( tag->getOption< bool >( "resnums_allow_design", 1 ) );
	repack_shell( tag->getOption< core::Real >( "repack_shell", 8.0 ));
	runtime_assert( design_shell() <= repack_shell() );
	TR<<"repack_shell = "<<repack_shell()<<" design shell = "<<design_shell()<<std::endl;
}

void DesignAroundOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute(
		"resnums", xs_string,
		"can be a list of pdb numbers, such as 291B,101A." )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_shell", xsct_real,
		"what sphere to design around the target residues. Must be at least as large as design_shell.",
		"8.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"allow_design", xsct_rosetta_bool,
		"allow design in the sphere, else restrict to repacking.",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"resnums_allow_design", xsct_rosetta_bool,
		"allow design in the resnums list, else restrict to repacking.",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"repack_shell", xsct_real,
		"what sphere to pack around the target residues. Must be at least as large as design_shell.",
		"8.0" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
