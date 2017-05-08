// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @file   protocols/cyclic_peptide/CreateAngleConstraint.cc
/// @brief  Add angle constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author modified by Parisa Hosseinzadeh (parisah@uw.edu)

#include <protocols/cyclic_peptide/CreateAngleConstraint.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraintCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.CreateAngleConstraint" );

namespace protocols {
namespace cyclic_peptide {

CreateAngleConstraint::CreateAngleConstraint() //:
{}
CreateAngleConstraint::~CreateAngleConstraint()= default;


void CreateAngleConstraint::set(utility::vector1<Size> res_center,
	utility::vector1<std::string> atom_center,
	utility::vector1<Size> res1,
	utility::vector1<std::string> atom1,
	utility::vector1<Size> res2,
	utility::vector1<std::string> atom2,
	utility::vector1<std::string> cst_func
)
{
	res_center_=res_center;
	atom_center_=atom_center;
	res1_=res1;
	atom1_=atom1;
	res2_=res2;
	atom2_=atom2;
	cst_func_=cst_func;
}

////////////////////////////////////////////////////
//////////        APPLY FUNCTION      //////////////
///////////////////////////////////////////////////
//Actual apply function

void CreateAngleConstraint::apply( core::pose::Pose & pose )
{
	for ( Size i_cst=1; i_cst<=cst_func_.size(); ++i_cst ) {
		if ( cst_func_[i_cst] == "" ) {}
		else {
			std::istringstream data(cst_func_[i_cst]);
			std::string func_type;
			data >> func_type;
			core::scoring::func::FuncFactory func_factory;
			core::scoring::func::FuncOP func = func_factory.new_func( func_type );
			func->read_data( data );
			Size atomno0 = pose.residue_type(res_center_[i_cst]).atom_index(atom_center_[i_cst]);
			Size atomno1 = pose.residue_type(res1_[i_cst]).atom_index(atom1_[i_cst]);
			Size atomno2 = pose.residue_type(res2_[i_cst]).atom_index(atom2_[i_cst]);
			pose.add_constraint(
				core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint( core::id::AtomID(atomno1,res1_[i_cst]),
				core::id::AtomID(atomno0,res_center_[i_cst]),
				core::id::AtomID(atomno2,res2_[i_cst]), func ) ) ) );
		}
	}
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
CreateAngleConstraint::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Add" ) {
			res_center_.push_back( (*tag_it)->getOption< Size >( "res_center" ) );
			atom_center_.push_back( (*tag_it)->getOption< std::string >( "atom_center" ) );
			res1_.push_back( (*tag_it)->getOption< Size >( "res1" ) );
			atom1_.push_back( (*tag_it)->getOption< std::string >( "atom1" ) );
			res2_.push_back( (*tag_it)->getOption< Size >( "res2" ) );
			atom2_.push_back( (*tag_it)->getOption< std::string >( "atom2" ) );
			cst_func_.push_back( (*tag_it)->getOption< std::string >( "cst_func", "" ) );
		}
	}
}

moves::MoverOP CreateAngleConstraint::clone() const { return moves::MoverOP( new CreateAngleConstraint( *this ) ); }
moves::MoverOP CreateAngleConstraint::fresh_instance() const { return moves::MoverOP( new CreateAngleConstraint ); }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CreateAngleConstraintCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CreateAngleConstraint );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateAngleConstraintCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CreateAngleConstraint::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateAngleConstraint::mover_name()
// XRW TEMP {
// XRW TEMP  return "CreateAngleConstraint";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateAngleConstraint::get_name() const {
// XRW TEMP  return "CreateAngleConstraint";
// XRW TEMP }

std::string CreateAngleConstraint::get_name() const {
	return mover_name();
}

std::string CreateAngleConstraint::mover_name() {
	return "CreateAngleConstraint";
}

void CreateAngleConstraint::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist; //no attributes for main tag
	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "res_center", xsct_non_negative_integer, "Residue at center of angle" )
		+ XMLSchemaAttribute::required_attribute( "atom_center", xs_string, "Atom at center of angle" )
		+ XMLSchemaAttribute::required_attribute( "res1", xsct_non_negative_integer, "Residue on one side of angle" )
		+ XMLSchemaAttribute::required_attribute( "atom1", xs_string, "Atom on one side of angle" )
		+ XMLSchemaAttribute::required_attribute( "res2", xsct_non_negative_integer, "Residue on other side of the angle" )
		+ XMLSchemaAttribute::required_attribute( "atom2", xs_string, "Atom on other side of the angle" )
		+ XMLSchemaAttribute::required_attribute( "cst_func", xs_string, "Function to use for this constraint" );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_simple_subelement( "Add", subelement_attlist, "Specifies an angle constraint to add" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Adds angle constraints to a pose", attlist, subelements );
}

std::string CreateAngleConstraintCreator::keyname() const {
	return CreateAngleConstraint::mover_name();
}

protocols::moves::MoverOP
CreateAngleConstraintCreator::create_mover() const {
	return protocols::moves::MoverOP( new CreateAngleConstraint );
}

void CreateAngleConstraintCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CreateAngleConstraint::provide_xml_schema( xsd );
}


} // moves
} // protocols
