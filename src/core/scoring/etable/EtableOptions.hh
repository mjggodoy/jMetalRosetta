// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_etable_EtableOptions_hh
#define INCLUDED_core_scoring_etable_EtableOptions_hh


#include <core/types.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

// C++ headers
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace etable {

class EtableOptions : public utility::pointer::ReferenceCount {

public:

	EtableOptions();
	EtableOptions( utility::options::OptionCollection const & options );

	~EtableOptions();

	EtableOptions( EtableOptions const & src );

	EtableOptions &
	operator=( EtableOptions const & src );

	friend
	bool
	operator < ( EtableOptions const & a, EtableOptions const & b );

	friend
	bool
	operator==( EtableOptions const & a, EtableOptions const & b );

	friend
	std::ostream &
	operator<< ( std::ostream & out, const EtableOptions & options );

	void
	show( std::ostream & out ) const;

	void
	parse_my_tag( utility::tag::TagCOP tag );

	static
	void
	append_schema_attributes( utility::tag::AttributeList & attributes );

	void initialize_from_options();
	void initialize_from_options( utility::options::OptionCollection const & options );

	static
	void
	list_options_read( utility::options::OptionKeyList & option_list );

public:

	std::string etable_type;
	bool analytic_etable_evaluation;
	Real max_dis;
	int  bins_per_A2;
	Real Wradius;
	Real lj_switch_dis2sigma;
	bool no_lk_polar_desolvation;
	bool proline_N_is_lk_nonpolar; // bazzoli: is the N atom of proline treated as non-polar for LK solvation?
	Real lj_hbond_OH_donor_dis;
	Real lj_hbond_hdis;
	bool enlarge_h_lj_wdepth;
	bool fa_hatr;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // etable
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_etable_EtableOptions )
#endif // SERIALIZATION


#endif
