// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Class to store ingformation about symmetrical dofs
/// @file   core/conformation/symmetry/SymSlideInfo.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymSlideInfo.hh>

// Utility header
#include <algorithm>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace symmetry {

SymSlideInfo::SymSlideInfo()
: slide_type_( RANDOM ),
	score_criteria_( CEN_DOCK_SCORE ),
	SlideCriteriaVal_( "AUTOMATIC" )
{
}

SymSlideInfo::SymSlideInfo( SymSlideInfo const & src )
: slide_type_( src.slide_type_ ),
	score_criteria_(  src.score_criteria_ ),
	SlideCriteriaVal_(  src.SlideCriteriaVal_ ),
	slide_order_( src.slide_order_ )
{
}

SymSlideInfo &
SymSlideInfo::operator=( SymSlideInfo const & src ) {
	slide_type_ = src.slide_type_;
	score_criteria_ = src.score_criteria_;
	SlideCriteriaVal_ = src.SlideCriteriaVal_;
	slide_order_ = src.slide_order_;
	return *this;
}

SymSlideInfo::~SymSlideInfo() {}

// setter functions
void SymSlideInfo::set_slide_type( SlideType slide_type )
{
	slide_type_ = slide_type;
}

void SymSlideInfo::set_SlideCriteriaType( SlideCriteriaType score_criteria )
{
	score_criteria_ = score_criteria;
}

void SymSlideInfo::set_SlideCriteriaVal( std::string SlideCriteriaVal )
{
	SlideCriteriaVal_ = SlideCriteriaVal;
}

void SymSlideInfo::set_slide_order( std::vector<core::Size> slide_order )
{
	slide_order_ = slide_order;
}

// get functions
SlideType
SymSlideInfo::get_slide_type() const
{
	return slide_type_;
}

SlideCriteriaType
SymSlideInfo::get_SlideCriteriaType() const
{
	return score_criteria_;
}

std::string
SymSlideInfo::get_SlideCriteriaVal() const
{
	return SlideCriteriaVal_;
}

std::vector<core::Size>
SymSlideInfo::get_slide_order() const
{
	return slide_order_;
}


bool
operator==(
	SymSlideInfo const & a,
	SymSlideInfo const & b
) {
	return
		(a.slide_type_ == b.slide_type_) &&
		(a.score_criteria_ == b.score_criteria_) &&
		(a.SlideCriteriaVal_ == b.SlideCriteriaVal_) &&
		std::equal(
		a.slide_order_.begin(), a.slide_order_.end(), b.slide_order_.begin());
}

bool
operator!=(
	SymSlideInfo const & a,
	SymSlideInfo const & b
) {
	return !(a == b);
}


} // symmetry
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::SymSlideInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( slide_type_ ) ); // enum core::conformation::symmetry::SlideType
	arc( CEREAL_NVP( score_criteria_ ) ); // enum core::conformation::symmetry::SlideCriteriaType
	arc( CEREAL_NVP( SlideCriteriaVal_ ) ); // std::string
	arc( CEREAL_NVP( slide_order_ ) ); // std::vector<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::SymSlideInfo::load( Archive & arc ) {
	arc( slide_type_ ); // enum core::conformation::symmetry::SlideType
	arc( score_criteria_ ); // enum core::conformation::symmetry::SlideCriteriaType
	arc( SlideCriteriaVal_ ); // std::string
	arc( slide_order_ ); // std::vector<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::SymSlideInfo );
#endif // SERIALIZATION
