// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Phil Bradley
/// @details
///
///     The ResidueTypeSelector is an object the picks out a subset of ResidueTypes, via a
///     bool operator[](ResidueType const &) method.  It is implemented as a logical AND of individual constraints,
///     each of which typically has an OR structure.  The system allows NOT at the beginning.
///
///     AA aa1 aa2 aa3
///     VARIANT_TYPE type1 type2 type3
///     PROPERTY property1 property2
///
///     So, e.g., the lines
///
///     PROPERTY PROTEIN
///     AA PRO GLY
///     NAME3 HPR
///     NOT VARIANT_TYPE PHOSPHO TERMINUS
///
///     would define a selector that matched residues with property PROTEIN, with aa types
///     pro or gly, with a three-letter code of HPR and not of variant type PHOSPHO or TERMINUS
///
///     The individual constraints that make up the ResidueTypeSelector object are subclasses of
///     ResidueTypeSelectorSingle; ResidueTypeSelector has a vector1 of ResidueTypeSelectorSingleOP's


#ifndef INCLUDED_core_chemical_ResidueTypeSelector_hh
#define INCLUDED_core_chemical_ResidueTypeSelector_hh

// Unit headers
#include <core/chemical/ResidueTypeSelector.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/ResidueProperties.hh>

// Utility headers
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

////////////////////////////////////////////////////////////////////////////////////////
/// @brief  A base class for defining a ResidueTypeSelector by a single criterion
class ResidueTypeSelectorSingle : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ResidueTypeSelectorSingle() override;

	// constructor
	ResidueTypeSelectorSingle( bool const result ):
		desired_result_( result )
	{}

	// select positively or negatively
	bool
	desired_result() const
	{
		return desired_result_;
	}

	virtual
	bool
	operator[]( ResidueType const & rsd ) const = 0;

private:
	bool desired_result_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ResidueTypeSelectorSingle();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue belong to ANY of these AAs?
class Selector_AA : public ResidueTypeSelectorSingle {
public:
	Selector_AA(
		utility::vector1< AA > const & aas_in,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		aas_( aas_in )
	{}

	// select by AA type
	bool
	operator[]( ResidueType const & rsd ) const override {
		// left-hand side will be TRUE if rsd.aa() is present in our list of AA's
		//std::cout << "Selector_AA: " << rsd.aa() << ' ' << aas_.size() << ' ' << desired_result() << std::endl;
		return ( ( std::find( aas_.begin(), aas_.end(), rsd.aa() ) != aas_.end() ) == desired_result() );
	}

	// data
private:
	utility::vector1< AA > aas_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_AA();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Is a certain string in the command-line option -chemical:allow_patch present ?
/// this selector does actually not depend on the residuetype it is queried for
class Selector_CMDFLAG : public ResidueTypeSelectorSingle {
public:
	Selector_CMDFLAG(
		std::string const& flags_in,
		bool const result
	);

	// select by AA type
	bool
	operator[]( ResidueType const & ) const override {
		return b_flag_is_present_ == desired_result();
	}

	// data
private:
	bool b_flag_is_present_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_CMDFLAG();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue have to ANY of these three-letter codes?
class Selector_NAME3 : public ResidueTypeSelectorSingle {
public:

	Selector_NAME3(
		utility::vector1< std::string > const & name3s_in,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		name3s_( name3s_in )
	{}

	// select by three-letter code
	bool
	operator[]( ResidueType const & rsd ) const override {
		return (  ( std::find( name3s_.begin(), name3s_.end(), rsd.name3() ) != name3s_.end() ) == desired_result() );
	}

private:
	utility::vector1< std::string > name3s_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_NAME3();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue have ANY of these properties?
class Selector_PROPERTY : public ResidueTypeSelectorSingle {
public:
	Selector_PROPERTY(
		utility::vector1< std::string > const & properties_in,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		properties_( properties_in )
	{}

	// select by PROPERTY
	bool
	operator[]( ResidueType const & rsd ) const override {
		for ( auto const & propertie : properties_ ) {
			if ( rsd.has_property( propertie ) ) return desired_result();
		}
		return !desired_result();
	}

	// data
private:
	utility::vector1< std::string > properties_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_PROPERTY();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue have ANY of variant types?
class Selector_VARIANT_TYPE : public ResidueTypeSelectorSingle {
public:
	Selector_VARIANT_TYPE(
		utility::vector1< std::string > const & variants_in,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		variants_( variants_in )
	{}

	// select by VARIANT_TYPE
	bool
	operator[]( ResidueType const & rsd ) const override {
		for ( auto const & variant : variants_ ) {
			if ( rsd.has_variant_type( variant ) ) return desired_result();
		}
		return !desired_result();
	}

	// data
private:
	utility::vector1< std::string > variants_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_VARIANT_TYPE();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief   Does the main chain of this residue follow from the given position label?
/// @details By position, it is meant the single digit integer by which atoms are labeled to indicate their position.
/// For example, if an upper connection is from C8 of a generic residue, the position is 8.
/// @remarks This selector was added primarily for use with carbohydrate residues, which have a wide assortment of main-
/// chain designations.  To properly patch upper terminus variants, it is necessary to know which atoms need to be added
/// and, particularly, at which position to add them.  However, this selector could be adapted for any residue subclass
/// that contains variability in the main chain, provided the nomenclature is consistent and numerical.
/// See patches/carbohydrates/upper_terminus.txt for an example of use.
/// @author  Labonte
class Selector_UPPER_POSITION : public ResidueTypeSelectorSingle {
public:
	Selector_UPPER_POSITION(uint const position, bool const result) :
		ResidueTypeSelectorSingle(result),
		position_(position)
	{}

	/// @brief Select by the position at which the upper connection is attached.
	/// @note  This is currently written to work only with carbohydrates, since only carbohydrates store their main-
	/// chain position.  If needed for another residue subclass, one would need to modify this method.  OR, it would be
	/// better to create a datum for storing the main-chain position in ResidueType.
	bool
	operator[](ResidueType const & rsd) const override {
		if ( rsd.is_carbohydrate() ) {  // only option for now
			if ( rsd.carbohydrate_info()->mainchain_glycosidic_bond_acceptor() == position_ ) {
				return desired_result();
			}
		}
		return !desired_result();
	}

private:
	uint position_;  // This is the position label at which the upper connection is attached.
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_UPPER_POSITION();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue have ALL of the variant types and no more
class Selector_MATCH_VARIANTS : public ResidueTypeSelectorSingle {
public:
	Selector_MATCH_VARIANTS(
		utility::vector1< std::string > const & variants_in,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		variants_( variants_in )
	{}

	// select by VARIANT_TYPE
	bool
	operator[]( ResidueType const & rsd ) const override {
		for ( auto const & variant : variants_ ) {
			if ( !rsd.has_variant_type( variant ) ) return !desired_result(); // rsd is missing one of our variants
		}
		if ( rsd.properties().get_list_of_variants().size() == variants_.size() ) return desired_result();
		return !desired_result(); // residue has an extra variant
	}

	// data
private:
	utility::vector1< std::string > variants_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_MATCH_VARIANTS();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue have NO variant types?
class Selector_NO_VARIANTS : public ResidueTypeSelectorSingle {
public:
	Selector_NO_VARIANTS(
		bool const result
	):
		ResidueTypeSelectorSingle( result )
	{}

	// select by VARIANT_TYPE
	bool
	operator[]( ResidueType const & rsd ) const override {
		return ( rsd.properties().get_list_of_variants().empty() == desired_result() );
	}
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_NO_VARIANTS();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Does the residue belong to ANY of these one-letter codes?
class Selector_NAME1 : public ResidueTypeSelectorSingle {
public:
	Selector_NAME1(
		char const n,
		bool const result
	):
		ResidueTypeSelectorSingle( result ),
		name1_( n )
	{}

	// select by name1 type
	bool
	operator[]( ResidueType const & rsd ) const override {
		return ( ( rsd.name1() == name1_ ) == desired_result() );
	}

	// data
private:
	char name1_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Selector_NAME1();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief create a singe ResidueTypeSelector from an input line.
ResidueTypeSelectorSingleOP
residue_selector_single_from_line( std::string const & line );


////////////////////////////////////////////////////////////////////////////////////////
/// @brief A class picking out a subset of ResidueType by multiple criteria
class ResidueTypeSelector : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ResidueTypeSelector() override;

	// [] operator: selector[ResidueType] => yes or no
	bool
	operator[]( ResidueType const & rsd ) const
	{
		//std::cout << "ResidueTypeSelector::operator[] " << rsd.name() << ' ' << selectors_.size() << std::endl;
		for ( uint i=1, i_end = selectors_.size(); i<= i_end; ++i ) {
			if ( !( ( *selectors_[i] )[ rsd ] ) ) return false;
		}
		return true;
	}

	// add a new selector single
	ResidueTypeSelector & // allow chaining
	add_line( std::string const & line )
	{
		ResidueTypeSelectorSingleOP new_selector( residue_selector_single_from_line( line ) );
		if ( new_selector ) {
			//std::cout << "add_line: success: " << line << std::endl;
			selectors_.push_back( new_selector );
		} else {
			std::cout << "ResidueTypeSelector::add_line: bad line:" << line << std::endl;
		}
		return *this;
	}

	// reset
	ResidueTypeSelector & // allow chaining
	clear()
	{
		selectors_.clear();
		return *this;
	}


	ResidueTypeSelector & // allow chaining
	set_name1( char const n )
	{
		selectors_.push_back( utility::pointer::shared_ptr<class core::chemical::ResidueTypeSelectorSingle>( new Selector_NAME1( n, true ) ) );
		return *this;
	}


	ResidueTypeSelector & // allow chaining
	set_aa( AA const aa )
	{
		utility::vector1< AA > aas( 1, aa );
		selectors_.push_back( utility::pointer::shared_ptr<class core::chemical::ResidueTypeSelectorSingle>( new Selector_AA( aas, true ) ) );
		return *this;
	}


	ResidueTypeSelector & // allow chaining
	set_property( std::string const property )
	{
		utility::vector1< std::string > properties( 1, property );
		selectors_.push_back( utility::pointer::shared_ptr<class core::chemical::ResidueTypeSelectorSingle>( new Selector_PROPERTY( properties, true ) ) );
		return *this;
	}


	ResidueTypeSelector & // allow chaining
	exclude_variants()
	{
		selectors_.push_back( utility::pointer::shared_ptr<class core::chemical::ResidueTypeSelectorSingle>( new Selector_NO_VARIANTS( true ) ) );
		return *this;
	}


	ResidueTypeSelector & // allow chaining
	match_variants( ResidueType const & rsd_type_to_match )
	{
		selectors_.push_back( utility::pointer::shared_ptr<class core::chemical::ResidueTypeSelectorSingle>( new Selector_MATCH_VARIANTS( rsd_type_to_match.properties().get_list_of_variants(), true ) ) );
		return *this;
	}

	// data
private:
	// a vector of single ResidueTypeSelector
	utility::vector1< ResidueTypeSelectorSingleOP > selectors_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueTypeSelector )
#endif // SERIALIZATION


#endif
