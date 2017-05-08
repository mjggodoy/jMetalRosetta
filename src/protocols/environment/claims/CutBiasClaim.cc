// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CutBiasClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/CutBiasClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/ClientMover.hh>

// Project Headers
#include <core/fragment/SecondaryStructure.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// C++ headers

// option key includes


static THREAD_LOCAL basic::Tracer tr( "protocols.environment.CutBiasClaim", basic::t_info );

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;

CutBiasClaim::CutBiasClaim( ClientMoverOP owner,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const& datamap ):
	Parent( owner ),
	label_( tag->getOption< std::string >( "label" ) )
{
	using core::select::residue_selector::ResidueSelector;

	core::Real bias = tag->getOption< core::Real >( "bias" );

	std::pair< core::Size, core::Size > range;
	range.first = tag->getOption< core::Size >( "region_start" );
	range.second = tag->getOption< core::Size >( "region_end" );

	for ( Size i = range.first; i <= range.second; ++i ) {
		biases_[ LocalPosition( label_, i ) ] = bias;
	}

	if ( datamap.has( "ResidueSelector", label_ ) ) {
		this->queue_for_annotation( label_, datamap.get_ptr< ResidueSelector const >( "ResidueSelector", label_ ) );
	}
}

CutBiasClaim::CutBiasClaim( ClientMoverOP owner,
	std::string const& label,
	core::fragment::SecondaryStructure const& ss_in ):
	Parent( owner ),
	label_( label )
{
	ObjexxFCL::FArray1D_float const& loop_frac = ss_in.loop_fraction();

	for ( Size i = 1; i <= ss_in.total_residue(); ++i ) {
		biases_[ LocalPosition( label_, i ) ] = loop_frac( (int) i );
	}

}

CutBiasClaim::CutBiasClaim( ClientMoverOP owner,
	std::string const& label,
	std::map< LocalPosition, core::Real > const& biases ):
	Parent( owner ),
	label_( label ),
	biases_( biases )
{}

CutBiasClaim::CutBiasClaim( ClientMoverOP owner,
	std::string const& label,
	std::pair< core::Size, core::Size > const& range,
	core::Real bias ):
	Parent( owner ),
	label_( label )
{
	for ( Size i = range.first; i <= range.second; ++i ) {
		biases_[ LocalPosition( label_, i ) ] = bias;
	}
}

void CutBiasClaim::yield_elements( FoldTreeSketch const&, CutBiasElements& elements ) const {

	for ( std::map< LocalPosition, core::Real >::const_iterator it = biases_.begin();
			it != biases_.end(); ++it ) {
		CutBiasElement e;

		e.p = it->first;
		e.bias = it->second;

		elements.push_back( e );
	}

}

std::string
CutBiasClaim::class_name(){
	return "CutBiasClaim";
}



void
CutBiasClaim::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "label", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "bias", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "region_start", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "region_end", xsct_non_negative_integer, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( class_name() )
		.complex_type_naming_func ( & EnvClaim::envclaim_ct_namer )
		.description( "XRW TO DO" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );


}

EnvClaimOP CutBiasClaim::clone() const {
	return EnvClaimOP( new CutBiasClaim( *this ) );
}

std::string CutBiasClaim::type() const{
	return "CutBias";
}

void CutBiasClaim::show( std::ostream& os ) const {
	os << type() << " owned by a " << owner()->get_name();
}

} //claims
} //environment
} //protocols
