// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/OptCysHG.cc
/// @brief  run optH on non-disulfided bonded CYS only; meant to relieve
///         any clashes caused by swapping of CYD->CYS after calling
///         Conformation::detect_disulfides()
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pack/task/operation/OptCysHG.hh>
#include <core/pack/task/operation/OptCysHGCreator.hh>

// package headers
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// package headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief default constructor
OptCysHG::OptCysHG() :
	Super()
{}


/// @brief copy constructor
OptCysHG::OptCysHG( OptCysHG const & rval ) :
	Super( rval )
{}


/// @brief default destructor
OptCysHG::~OptCysHG() {}

/// @brief clone this object
OptCysHG::TaskOperationOP OptCysHG::clone() const {
	return OptCysHG::TaskOperationOP( new OptCysHG( *this ) );
}


/// @brief apply operations to PackerTask
void OptCysHG::apply( Pose const & pose, PackerTask & task ) const {
	using core::chemical::aa_cys;
	using core::chemical::DISULFIDE;

	OptH optH;

	// restrict to only non-disulfide bonded CYS
	for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
		if ( ! pose.residue_type( i ).is_sidechain_thiol() || pose.residue_type( i ).has_variant_type( DISULFIDE ) ) { //check both names to be double sure; used in fake Cys catalytic residues
			optH.disallow_resid( i );
		}
	}

	optH.apply( pose, task );
}


std::string OptCysHG::keyname() { return "OptCysHG"; }

void OptCysHG::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "run optH on non-disulfided bonded CYS only; meant to relieve any clashes caused by swapping of CYD to CYS after calling conformation::detect_disulfides" );
}

TaskOperationOP OptCysHGCreator::create_task_operation() const
{
	return TaskOperationOP( new OptCysHG );
}

std::string OptCysHGCreator::keyname() const { return OptCysHG::keyname(); }

void OptCysHGCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	OptCysHG::provide_xml_schema( xsd );
}


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core
