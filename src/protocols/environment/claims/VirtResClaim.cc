// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file VirtResClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/VirtResClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/DofPassport.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/ClientMover.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

#include <core/kinematics/AtomTree.hh>

// Project Headers
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>

#include <core/pose/util.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// C++ headers

// option key includes

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.VirtResClaim", basic::t_info );

namespace protocols {
namespace environment {
namespace claims {

std::string const VRT_LABEL_OPTION = "vrt_name";
std::string const PARENT_LABEL_OPTION = "parent";

VirtResClaim::VirtResClaim( ClientMoverOP owner,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const& datamap ):
	EnvClaim( owner ),
	vrt_label_( tag->getOption< std::string >( VRT_LABEL_OPTION ) ),
	j_claim_( owner,
	tag->getOption( "jump_label", vrt_label()+"_jump" ),
	tag->getOption< std::string >( PARENT_LABEL_OPTION ),
	LocalPosition( vrt_label(), 1 ) ),
	xyz_claim_( owner, vrt_label() )
{
	if ( datamap.has( "ResidueSelector", j_claim_.pos2().label() ) ) {
		using core::select::residue_selector::ResidueSelector;
		this->queue_for_annotation( j_claim_.pos2().label(), datamap.get_ptr< ResidueSelector const >( "ResidueSelector", j_claim_.pos2().label() ) );
	}

	j_claim_.physical( true );
	j_claim_.strength( Parent::parse_ctrl_str( tag->getOption< std::string >( "jump_control_strength", "DOES_NOT_CONTROL" ) ),
		DOES_NOT_CONTROL );
	xyz_claim_.strength( MUST_CONTROL,
		DOES_NOT_CONTROL );
}


VirtResClaim::VirtResClaim( ClientMoverOP owner,
	LocalPosition parent,
	std::string const& jump_label,
	std::string const& vrt_label ):
	EnvClaim( owner ),
	vrt_label_( vrt_label ),
	j_claim_( owner,
	jump_label,
	parent,
	LocalPosition( vrt_label, 1 ) ),
	xyz_claim_( owner, vrt_label )
{
	j_claim_.physical( true );
	xyz_claim_.strength( MUST_CONTROL, MUST_CONTROL );
}

void VirtResClaim::yield_elements( FoldTreeSketch const&, ResidueElements& elements ) const{
	ResidueElement e;

	e.label = vrt_label();

	elements.push_back( e );
}

void VirtResClaim::yield_elements( FoldTreeSketch const& fts, JumpElements& elements ) const {
	j_claim_.yield_elements( fts, elements );
}

void VirtResClaim::yield_elements( FoldTreeSketch const& fts, CutElements& elements ) const {
	j_claim_.yield_elements( fts, elements );
}

void VirtResClaim::yield_elements( core::pose::Pose const& pose, DOFElements& elements ) const {
	xyz_claim_.yield_elements( pose, elements );

	// Not using the jump claim's yield_elements code; we want to claim *all* jumps associated with the VRT
	// TODO: make this optional

	core::Size vrt_pos = static_cast< ProtectedConformation const* >( &pose.conformation() )->annotations()->resolve_seq( vrt_label() )[1];
	core::kinematics::FoldTree const& ft = pose.conformation().fold_tree();

	for ( int j_num = 1; j_num <= (int) ft.num_jump(); ++j_num ) {
		if ( ft.upstream_jump_residue( j_num ) == (int) vrt_pos ||
				ft.downstream_jump_residue( j_num ) == (int) vrt_pos ) {
			for ( core::Size rb_i = core::id::RB1; rb_i <= core::id::RB6; ++rb_i ) {
				DOFElement e = Parent::wrap_dof_id( core::id::DOF_ID( pose.conformation().jump_atom_id( j_num ),
					core::id::DOF_Type( rb_i ) ) );

				e.i_str = MUST_CONTROL;
				e.c_str = MUST_CONTROL;

				elements.push_back( e );
			}
		}
	}
}

void VirtResClaim::strength( ControlStrength const& c_str, ControlStrength const& i_str ){
	jump().strength( c_str, i_str );
	xyz_claim_.strength( c_str, i_str );
}


std::string
VirtResClaim::class_name(){
	//Currently this is what tags look for instead of the actual class name; is this right?
	return "VirtResClaim";
}



void
VirtResClaim::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define control strings
	XMLSchemaRestriction restr;
	restr.name( "envclaim_ctrl_str" );
	restr.base_type( xs_string );
	restr.add_restriction( xsr_pattern, "((([Dd][Oo][Ee][Ss]_[Nn][Oo][Tt])|([Cc][Aa][Nn])|([Mm][Uu][Ss][Tt]))_[Cc][Oo][Nn][Tt][Rr][Oo][Ll])|([Ee][Xx][Cc][Ll][Uu][Ss][Ii][Vv][Ee])" );
	xsd.add_top_level_element( restr );

	AttributeList attlist;
	attlist

		+ XMLSchemaAttribute::required_attribute( "vrt_name", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "jump_label", xs_string, "Defaults to vrt_name with _jump appended" )
		+ XMLSchemaAttribute::required_attribute( "parent", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_control_strength", "envclaim_ctrl_str", "XRW TO DO", "DOES_NOT_CONTROL" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( class_name() )
		.complex_type_naming_func ( & EnvClaim::envclaim_ct_namer )
		.description( "XRW TO DO" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );

}



EnvClaimOP VirtResClaim::clone() const {
	return EnvClaimOP( new VirtResClaim( *this ) );
}

std::string const& VirtResClaim::vrt_label() const {
	return vrt_label_;
}

std::string VirtResClaim::type() const{
	return "VirtRes";
}

void VirtResClaim::show( std::ostream& os ) const {
	os << type() << " '" << vrt_label() << "with jump '"
		<< j_claim_.label() << "' owned by a " << owner()->get_name();
}

} //claims
} //environment
} //protocols
