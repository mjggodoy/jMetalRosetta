// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/ConstraintsEnergyContainer.cc
/// @brief  Constraints Energy Container class implementation
/// @author Spencer Bliven <blivens@u.washington.edu>

// Unit headers
#include <core/scoring/disulfides/DisulfideMatchingEnergyContainer.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// STL Headers
#include <utility/assert.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Core serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace disulfides {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.disulfides.DisulfideMatchingEnergyContainer" );

/// @brief constructor
DisulfideMatchingNeighborIterator::DisulfideMatchingNeighborIterator(
	DisulfideMatchingEnergyContainer * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

/// @brief constructor, default to no disulfide bond
DisulfideMatchingNeighborIterator::DisulfideMatchingNeighborIterator(
	DisulfideMatchingEnergyContainer * owner
)
:
	owner_( owner ),
	focused_residue_( DisulfideMatchingEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( DisulfideMatchingEnergyContainer::NO_DISULFIDE )
{}

DisulfideMatchingNeighborIterator::~DisulfideMatchingNeighborIterator()
{}

/// @brief Assignment
ResidueNeighborIterator &
DisulfideMatchingNeighborIterator::operator = ( ResidueNeighborIterator const & rhs)
{
	debug_assert( &(dynamic_cast< DisulfideMatchingNeighborIterator const & > ( rhs )) );
	DisulfideMatchingNeighborIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;
}

/// @note Incrementing an iterator in a list with exactly one element moves that
/// iterator off the end of the list.
ResidueNeighborIterator const &
DisulfideMatchingNeighborIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	focused_residue_ = DisulfideMatchingEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = DisulfideMatchingEnergyContainer::NO_DISULFIDE;
	return *this;
}

bool
DisulfideMatchingNeighborIterator::operator == ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfideMatchingNeighborIterator const & > ( rhs )) );
	DisulfideMatchingNeighborIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}

bool
DisulfideMatchingNeighborIterator::operator != ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfideMatchingNeighborIterator const & > ( rhs )) );
	DisulfideMatchingNeighborIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}

/// @brief Get the higher-numbered residue for this disulfide bond
Size
DisulfideMatchingNeighborIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

/// @brief Get the lower-numbered residue for this disulfide bond
Size
DisulfideMatchingNeighborIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}

/// @brief Which residue are we looking for disulfide bonds to?
Size
DisulfideMatchingNeighborIterator::residue_iterated_on() const
{
	return focused_residue_;
}

/// @brief Which residue is disulfide bonded to the current residue?
Size
DisulfideMatchingNeighborIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}

/// @brief Save the specified energies for this disulfide to the
/// DisulfideMatchingEnergyContainer associated with this iterator.
void
DisulfideMatchingNeighborIterator::save_energy( EnergyMap const & emap )
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->save_energy( disulfide_index_, emap );
}

/// @brief Get the energies for the current disulfide bond from the
/// DisulfideMatchingEnergyContainer associated with this iterator.
void
DisulfideMatchingNeighborIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}

/// @brief Add some energies to the totals already in DisulfideMatchingEnergyContainer
void
DisulfideMatchingNeighborIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

void DisulfideMatchingNeighborIterator::mark_energy_computed()
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_computed( disulfide_index_ );
}

void DisulfideMatchingNeighborIterator::mark_energy_uncomputed()
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_uncomputed( disulfide_index_ );
}


bool
DisulfideMatchingNeighborIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}

////////////////////////////////////////////////////////////////////////
///// Disulfide Residue Neighbor Constant Iterator class implementation
////////////////////////////////////////////////////////////////////////

DisulfideMatchingNeighborConstIterator::DisulfideMatchingNeighborConstIterator(
	DisulfideMatchingEnergyContainer const * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

DisulfideMatchingNeighborConstIterator::DisulfideMatchingNeighborConstIterator(
	DisulfideMatchingEnergyContainer const * owner
) :
	owner_( owner ),
	focused_residue_( DisulfideMatchingEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( DisulfideMatchingEnergyContainer::NO_DISULFIDE )
{}

DisulfideMatchingNeighborConstIterator::~DisulfideMatchingNeighborConstIterator()
{}

ResidueNeighborConstIterator &
DisulfideMatchingNeighborConstIterator::operator = ( ResidueNeighborConstIterator const & rhs )
{
	debug_assert( &(dynamic_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs )) );
	DisulfideMatchingNeighborConstIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;

}

ResidueNeighborConstIterator const &
DisulfideMatchingNeighborConstIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	focused_residue_ = DisulfideMatchingEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = DisulfideMatchingEnergyContainer::NO_DISULFIDE;
	return *this;
}

/// @brief returns true if the two edge-list iterators are equal
bool
DisulfideMatchingNeighborConstIterator::operator == ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs )) );
	DisulfideMatchingNeighborConstIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}


/// @brief returns true if the two edge-list iterators are not equal
bool
DisulfideMatchingNeighborConstIterator::operator != ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs )) );
	DisulfideMatchingNeighborConstIterator const & drni_rhs = static_cast< DisulfideMatchingNeighborConstIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}

Size
DisulfideMatchingNeighborConstIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

Size
DisulfideMatchingNeighborConstIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}


Size
DisulfideMatchingNeighborConstIterator::residue_iterated_on() const
{
	return focused_residue_;
}

Size
DisulfideMatchingNeighborConstIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}

/// @brief overwrites the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not zero out the other positions in the emap.
void
DisulfideMatchingNeighborConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}

/// @brief accumulates the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not touch the other positions in the emap.
void
DisulfideMatchingNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

bool
DisulfideMatchingNeighborConstIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != DisulfideMatchingEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}


/////////////////////////////////////////////////////
/// Disulfide Energy Container Class Implementation
/////////////////////////////////////////////////////

Size const DisulfideMatchingEnergyContainer::NO_DISULFIDE( 0 );


DisulfideMatchingEnergyContainer::DisulfideMatchingEnergyContainer()
{}

bool
DisulfideMatchingEnergyContainer::empty() const
{
	return num_disulfides() == 0;
}


DisulfideMatchingEnergyContainer::DisulfideMatchingEnergyContainer( pose::Pose const & pose )
{
	find_disulfides( pose );
}

void
DisulfideMatchingEnergyContainer::update( pose::Pose const & pose )
{
	if ( disulfides_changed( pose ) ) find_disulfides( pose );
}

DisulfideMatchingEnergyContainer::~DisulfideMatchingEnergyContainer()
{}

LREnergyContainerOP
DisulfideMatchingEnergyContainer::clone() const
{
	DisulfideMatchingEnergyContainerOP dec( new DisulfideMatchingEnergyContainer );
	if ( !empty() ) {
		dec->disulfide_atom_indices_ = disulfide_atom_indices_;
		dec->disulfide_residue_types_ = disulfide_residue_types_;
		dec->resid_2_disulfide_index_ = resid_2_disulfide_index_;
		dec->disulfide_partners_ =  disulfide_partners_;
		dec->disulfide_info_ = disulfide_info_;
	}
	return dec;
}

bool
DisulfideMatchingEnergyContainer::any_neighbors_for_residue( int resid ) const
{
	return (Size) resid <= resid_2_disulfide_index_.size() && resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

bool
DisulfideMatchingEnergyContainer::any_upper_neighbors_for_residue( int resid ) const
{
	if ( (Size) resid <= resid_2_disulfide_index_.size() ) {
		if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
			// disulfide_partners_ stores an ordered pair s.t. the first index is to the
			// lower residue and thes second index is to the upper residue; if we want to
			// know whether resid has an upper neighbor, what we need to know is if
			// it's the lower residue.
			return disulfide_partners_[ resid_2_disulfide_index_[ resid ] ].first == (Size) resid;
		}
	}
	return false;
}

ResidueNeighborConstIteratorOP
DisulfideMatchingEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );
	if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
DisulfideMatchingEnergyContainer::const_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this ) );
}

ResidueNeighborConstIteratorOP
DisulfideMatchingEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );

	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid >= 0 &&
			Size( resid ) <= resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
DisulfideMatchingEnergyContainer::const_upper_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new DisulfideMatchingNeighborConstIterator( this ) );
}

ResidueNeighborIteratorOP
DisulfideMatchingEnergyContainer::neighbor_iterator_begin( int resid )
{
	debug_assert( !empty() );
	if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this ) );
	}
}

ResidueNeighborIteratorOP
DisulfideMatchingEnergyContainer::neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this ) );
}

ResidueNeighborIteratorOP
DisulfideMatchingEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid >= 0 &&
			Size( resid ) <= resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this ) );
	}
}

ResidueNeighborIteratorOP
DisulfideMatchingEnergyContainer::upper_neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new DisulfideMatchingNeighborIterator( this ) );
}

bool
DisulfideMatchingEnergyContainer::disulfide_bonded( Size res1id, Size res2id ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ res1id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res2id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res1id ] == resid_2_disulfide_index_[ res2id ];
}

bool
DisulfideMatchingEnergyContainer::residue_forms_disulfide( Size resid ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

Size
DisulfideMatchingEnergyContainer::other_neighbor_id( Size resid ) const
{
	return other_neighbor_id( resid_2_disulfide_index_[ resid ], resid );
}


// Mutators
void
DisulfideMatchingEnergyContainer::save_energy( Size disulfide_index, EnergyMap const & emap )
{
	disulfide_info_[ disulfide_index ].first.dslfc_rot() = emap[ dslfc_rot ];
	disulfide_info_[ disulfide_index ].first.dslfc_trans() = emap[ dslfc_trans ];
	disulfide_info_[ disulfide_index ].first.dslfc_RT() = emap[ dslfc_RT ];
}

void
DisulfideMatchingEnergyContainer::mark_energy_computed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = true;
}

void
DisulfideMatchingEnergyContainer::mark_energy_uncomputed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = false;
}

// Accessors
Size DisulfideMatchingEnergyContainer::lower_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].first;
}

Size DisulfideMatchingEnergyContainer::upper_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].second;
}

Size DisulfideMatchingEnergyContainer::other_neighbor_id( Size disulfide_index, Size resid ) const
{
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_partners_[ disulfide_index ].second :
		disulfide_partners_[ disulfide_index ].first );
}

DisulfideAtomIndices const &
DisulfideMatchingEnergyContainer::disulfide_atom_indices( Size resid ) const
{
	Size const disulfide_index( resid_2_disulfide_index_[ resid ] );
	debug_assert( disulfide_index != NO_DISULFIDE );
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_atom_indices_[ disulfide_index ].first :
		disulfide_atom_indices_[ disulfide_index ].second );
}


DisulfideAtomIndices const &
DisulfideMatchingEnergyContainer::other_neighbor_atom_indices( Size resid ) const
{
	Size const disulfide_index( resid_2_disulfide_index_[ resid ] );
	debug_assert( disulfide_index != NO_DISULFIDE );
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_atom_indices_[ disulfide_index ].second :
		disulfide_atom_indices_[ disulfide_index ].first );
}


void DisulfideMatchingEnergyContainer::accumulate_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslfc_rot ] += disulfide_info_[ disulfide_index ].first.dslfc_rot();
	emap[ dslfc_trans  ] += disulfide_info_[ disulfide_index ].first.dslfc_trans();
	emap[ dslfc_RT     ] += disulfide_info_[ disulfide_index ].first.dslfc_RT();
}

void DisulfideMatchingEnergyContainer::retrieve_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslfc_rot ] = disulfide_info_[ disulfide_index ].first.dslfc_rot();
	emap[ dslfc_trans  ] = disulfide_info_[ disulfide_index ].first.dslfc_trans();
	emap[ dslfc_RT     ] = disulfide_info_[ disulfide_index ].first.dslfc_RT();
}

bool DisulfideMatchingEnergyContainer::energy_computed( Size disulfide_index ) const
{
	return disulfide_info_[ disulfide_index ].second;
}


void
DisulfideMatchingEnergyContainer::find_disulfides( pose::Pose const & pose )
{
	TR.Debug << "In find_disulfides():" << std::endl;

	disulfide_partners_.clear();
	disulfide_atom_indices_.clear();
	disulfide_info_.clear();
	resid_2_disulfide_index_.resize( pose.size() );
	disulfide_residue_types_.resize( pose.size() );
	std::fill( resid_2_disulfide_index_.begin(), resid_2_disulfide_index_.end(), NO_DISULFIDE );
	std::fill( disulfide_residue_types_.begin(), disulfide_residue_types_.end(), chemical::ResidueTypeCOP(0) );

	Size count_disulfides( 0 );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		conformation::Residue res = pose.residue( ii );
		if ( // res.aa() == chemical::aa_cys &&
				res.type().is_disulfide_bonded() && // redundant with atom name check that used to be here
				res.has_variant_type( chemical::DISULFIDE ) &&
				resid_2_disulfide_index_[ ii ] == NO_DISULFIDE ) {
			++count_disulfides;

			//Centroid models are bonded CEN to CEN, fullatom are bonded SG to SG.
			//This code originally forgot to put the whole SG bit in so maybe this is a hack. I don't care!
			// -rv
			Size ii_connect_atom = res.atom_index( res.type().get_disulfide_atom_name() );

			Size other_res_ii( 0 );
			for ( Size jj = 1; jj <= res.type().n_possible_residue_connections(); ++jj ) {
				if ( (Size) res.type().residue_connection( jj ).atomno() == ii_connect_atom ) {
					other_res_ii = res.connect_map( jj ).resid();
					break;
				}
			}
			if ( other_res_ii == 0 ) {
				TR.Error << "ERROR: Could not find disulfide partner for residue " << ii << std::endl;
				utility_exit();
			}
			debug_assert( other_res_ii > ii );
			//Can only bond residues of the same residue type set (eg centroid to centroid)
			debug_assert( pose.residue_type(other_res_ii).mode() ==
				pose.residue_type(ii).mode() );

			TR.Debug << "Found disulf between " << ii << " and " << other_res_ii << std::endl;

			resid_2_disulfide_index_[ ii ] = count_disulfides;
			resid_2_disulfide_index_[ other_res_ii ] = count_disulfides;
			disulfide_residue_types_[ ii ] = pose.residue_type( ii ).get_self_ptr();
			disulfide_residue_types_[ other_res_ii ] = pose.residue_type( other_res_ii ).get_self_ptr();
			disulfide_partners_.push_back( std::pair< Size, Size >( ii, other_res_ii ) );
			disulfide_atom_indices_.push_back( std::pair< DisulfideAtomIndices, DisulfideAtomIndices > (
				DisulfideAtomIndices( pose.residue(ii ) ), DisulfideAtomIndices( pose.residue( other_res_ii ) ) ));
			DisulfideMatchingEnergyComponents temp;
			disulfide_info_.push_back( std::pair< DisulfideMatchingEnergyComponents, bool > ( temp, false ) );

			debug_assert(! empty());
		}
	}
	TR.Debug << "Found " << num_disulfides() << " DS" << std::endl;
}

// we could do something like keep a flag for when minimization is occurring and assume that
// disulfide connectivity information does not change during the course of minimization...
bool
DisulfideMatchingEnergyContainer::disulfides_changed( pose::Pose const & pose )
{
	Size const total_residue( pose.size() );
	if ( resid_2_disulfide_index_.size() != total_residue ) return true;

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		if ( resid_2_disulfide_index_[ ii ] != NO_DISULFIDE ) {
			conformation::Residue res = pose.residue( ii );
			if ( // res.aa() != chemical::aa_cys ||
					!res.type().is_disulfide_bonded() ||
					disulfide_residue_types_[ ii ].get() != & (pose.residue_type( ii )) ||
					/// subsumed by residue type check ! pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) ||
					! pose.residue_type( ii ).has( "CEN" ) || // not centroid
					res.connect_map(
					res.type().residue_connection_id_for_atom(
					res.atom_index( "CEN" ) ) ).resid() !=
					other_neighbor_id( resid_2_disulfide_index_[ ii ], ii ) ) {
				return true;
			}
		} else if ( ( pose.residue( ii ).type().is_sidechain_thiol() || pose.residue( ii ).type().is_disulfide_bonded() ) &&
				pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) ) {
			return true;
		}
	}
	return false;
}

Size DisulfideMatchingEnergyContainer::num_disulfides() const
{
	return disulfide_partners_.size();
}


}
}
}



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::disulfides::DisulfideMatchingEnergyComponents::save( Archive & arc ) const {
	arc( CEREAL_NVP( dslfc_rot_ ) ); // Energy
	arc( CEREAL_NVP( dslfc_trans_ ) ); // Energy
	arc( CEREAL_NVP( dslfc_RT_ ) ); // Energy
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::disulfides::DisulfideMatchingEnergyComponents::load( Archive & arc ) {
	arc( dslfc_rot_ ); // Energy
	arc( dslfc_trans_ ); // Energy
	arc( dslfc_RT_ ); // Energy
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::disulfides::DisulfideMatchingEnergyComponents );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::disulfides::DisulfideMatchingEnergyContainer::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::LREnergyContainer >( this ) );
	arc( CEREAL_NVP( resid_2_disulfide_index_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( disulfide_residue_types_ ) );
	arc( CEREAL_NVP( disulfide_partners_ ) ); // utility::vector1<std::pair<Size, Size> >
	arc( CEREAL_NVP( disulfide_atom_indices_ ) ); // utility::vector1<std::pair<DisulfideAtomIndices, DisulfideAtomIndices> >
	arc( CEREAL_NVP( disulfide_info_ ) ); // utility::vector1<std::pair<DisulfideMatchingEnergyComponents, _Bool> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::disulfides::DisulfideMatchingEnergyContainer::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::LREnergyContainer >( this ) );
	arc( resid_2_disulfide_index_ ); // utility::vector1<Size>
	arc( disulfide_residue_types_ );
	arc( disulfide_partners_ ); // utility::vector1<std::pair<Size, Size> >
	arc( disulfide_atom_indices_ ); // utility::vector1<std::pair<DisulfideAtomIndices, DisulfideAtomIndices> >
	arc( disulfide_info_ ); // utility::vector1<std::pair<DisulfideMatchingEnergyComponents, _Bool> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::disulfides::DisulfideMatchingEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::disulfides::DisulfideMatchingEnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_disulfides_DisulfideMatchingEnergyContainer )
#endif // SERIALIZATION
