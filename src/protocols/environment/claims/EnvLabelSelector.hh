// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/environment/claims/EnvLabelSelector.hh
/// @brief  The EnvLabelSelector holds a set of local positions that are converted to seqids at apply.
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_claims_EnvLabelSelector_HH
#define INCLUDED_protocols_environment_claims_EnvLabelSelector_HH

// Unit headers
#include <protocols/environment/claims/EnvLabelSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Package headers
#include <core/environment/LocalPosition.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace environment {
namespace claims {

class EnvLabelSelector : public core::select::residue_selector::ResidueSelector {
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;
	typedef core::environment::LocalPositions LocalPositions;
	typedef core::environment::LocalPosition LocalPosition;

public:
	// derived from base class
	EnvLabelSelector();

	/// @brief Copy constructor
	///
	EnvLabelSelector( EnvLabelSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual core::select::residue_selector::ResidueSelectorOP clone() const;

	EnvLabelSelector( LocalPositions const& );

	EnvLabelSelector( LocalPosition const& );

	EnvLabelSelector( std::string const& label,
		std::pair< core::Size, core::Size > const& range );

	virtual ~EnvLabelSelector();

	virtual
	core::select::residue_selector::ResidueSubset
	apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	void set_local_positions( LocalPositions const& );

	LocalPositions const& local_positions() const{ return positions_; }

	void add_position( LocalPosition const& p ){ positions_.push_back( utility::pointer::shared_ptr<class core::environment::LocalPosition>( new LocalPosition( p ) ) ); }

	virtual
	std::string
	get_name() const;

	static std::string class_name();

private: // data members
	LocalPositions positions_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //namespace claims
} //namespace environment
} //namespace protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_environment_claims_EnvLabelSelector )
#endif // SERIALIZATION


#endif
