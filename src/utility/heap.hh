// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file heap.hh
/// @brief class definition for a heap object based on Charlie
/// Strauss's heap code ported over from rosetta++. Stores a sorted list of
/// integers based on floating-point values.
/// @author James Thompson

#ifndef INCLUDED_utility_heap_hh
#define INCLUDED_utility_heap_hh

// Unit headers
#include <utility/heap.fwd.hh>

// TEMP

// ObjexxFCL Headers
#include <utility/vector0.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace utility {

class heap : public utility::pointer::ReferenceCount {

public:

	/// @brief Create a heap with this number of items.
	heap( int max_items );

	virtual ~heap() ; // auto-removing definition from header{}

	/// @brief Inserts a value into the heap that is sorted by coval. The err status
	/// will be set to true if the heap is already at capacity; this function w
	/// not increase the size of the heap.
	void
	heap_insert( int val, float coval, bool & err );

	/// @brief Extracts the val,coval pair with the lowest coval from the heap.
	/// This modifies the heap, and the returned values are put into the arguments
	/// val and coval.
	void
	heap_extract( int & val, float & coval, bool & err );

	void
	heap_replace( int val, float coval );

	void
	reset_coval( int val, float coval );

	/// @brief returns the smallest covalue stored in the heap.
	float
	heap_head() const;

	/// @brief returns the item with the smallest covalue
	int
	head_item() const;

	/// @brief Return the colvaue for a particular position in the heap
	float
	coval( int index ) const;

	/// @brief Return the item index for a particular position in the heap
	int
	val( int index ) const;

	/// @brief Return the colvaue for a particular item; O(N)
	float
	coval_for_val( int val ) const;

	/// @brief Return the size of the heap.
	int
	size() const;

	/// @brief Return the capacity of the heap.
	int
	capacity() const;

private:
	void
	heap_init( int max_items );

	void
	heap_down( int index_in );

	void
	heap_up( int index_in );

	int heap_size() const;
	int heap_capacity() const;

	int & heap_size();
	int & heap_capacity();

	void
	decrease_coval( int index, float coval );

	void
	increase_coval( int index, float coval );

	int index_for_val( int val ) const;

private:
	utility::vector0< int > heap_;
	utility::vector0< float > coheap_;

}; // class heap

} // ns utility

#endif
