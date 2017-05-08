// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictIdentitiesOperation.cc
/// @brief Restricts specified residue types to only repack, no design.
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictIdentitiesOperation.hh>
#include <protocols/toolbox/task_operations/RestrictIdentitiesOperationCreator.hh>

// Project Headers

#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pose/symmetry/util.hh>


// C++ Headers
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

//Auto Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <boost/functional/hash.hpp>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.RestrictIdentitiesOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
RestrictIdentitiesOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictIdentitiesOperation );
}

void RestrictIdentitiesOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictIdentitiesOperation::provide_xml_schema( xsd );
}

std::string RestrictIdentitiesOperationCreator::keyname() const
{
	return RestrictIdentitiesOperation::keyname();
}

// @brief default constructor
RestrictIdentitiesOperation::RestrictIdentitiesOperation() {}

// @brief constructor with arguments
RestrictIdentitiesOperation::RestrictIdentitiesOperation( utility::vector1 < std::string > identities, bool prevent_repacking ) :
	identities_( identities ),
	prevent_repacking_( prevent_repacking )
{}

// @brief destructor
RestrictIdentitiesOperation::~RestrictIdentitiesOperation() {}

// @brief copy constructor
core::pack::task::operation::TaskOperationOP RestrictIdentitiesOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictIdentitiesOperation( *this ) );
}

// @brief getters
utility::vector1< std::string > RestrictIdentitiesOperation::identities() const { return identities_; }
bool RestrictIdentitiesOperation::prevent_repacking() const { return prevent_repacking_; }

// @brief setters
void RestrictIdentitiesOperation::identities( utility::vector1 < std::string > identities_vec )
{
	runtime_assert( identities_vec.size() != 0 );
	identities_.clear();
	for ( std::string const & item : identities_vec ) {
		identities_.push_back( item );
	}
}
void RestrictIdentitiesOperation::prevent_repacking( bool const prevent_repacking) { prevent_repacking_ = prevent_repacking; }

// @brief apply function
void
RestrictIdentitiesOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	runtime_assert( identities_.size() != 0 );

	for ( Size j=1; j<=identities_.size(); j++ ) {
		for ( Size i=1; i<=pose.size(); i++ ) {
			std::string aa_name = pose.residue(i).name3();
			if ( aa_name == identities_[j] ) {
				if ( prevent_repacking_ == 1 ) {
					task.nonconst_residue_task(i).prevent_repacking();
					TR.Debug << "preventing " << aa_name << i << " from repacking" << std::endl;
				} else {
					task.nonconst_residue_task(i).restrict_to_repacking();
					TR.Debug << "restricting " << aa_name << i << " to repacking" << std::endl;
				}
			}
		}
	}

}

// @brief parse xml
void
RestrictIdentitiesOperation::parse_tag( TagCOP tag , DataMap & )
{
	unparsed_identities_ = tag->getOption< std::string >( "identities" ) ;
	prevent_repacking( tag->getOption< bool >( "prevent_repacking", false ));
	if ( unparsed_identities_ != "" ) {

		utility::vector1< std::string > const ids( utility::string_split( unparsed_identities_ , ',' ) );
		identities_.clear();

		std::string action = (prevent_repacking_) ? "Preventing from repacking " : "Restricting to repacking " ;
		TR << action << "residues of type(s): ";
		for ( std::string const & id : ids ) {
			TR << id << " ";
			identities_.push_back( id );
		}
		TR.Debug << std::endl;
	}
}

void RestrictIdentitiesOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	// See, I think this has to be there. But it still checks for a non-empty.
	// AMW TODO check
	attributes
		+ XMLSchemaAttribute::required_attribute( "identities", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default(  "prevent_repacking", xsct_rosetta_bool, "XRW TO DO",  "false"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
