// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResidueLevelTask.hh
/// @brief  Task class to describe packer's behavior header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Steven Lewis

#ifndef INCLUDED_core_pack_task_ResidueLevelTask_hh
#define INCLUDED_core_pack_task_ResidueLevelTask_hh

// Unit Headers

// Package Headers
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerLinks.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.fwd.hh>

// Project Headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>

// STL Headers
#include <iosfwd>


#ifdef PYROSETTA
	#include <core/id/SequenceMapping.hh>
#endif


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

class ResidueLevelTask
{
public:
	typedef std::list< chemical::ResidueTypeCOP > ResidueTypeCOPList;
	typedef std::list< chemical::ResidueTypeCOP >::iterator ResidueTypeCOPListIter;
	typedef std::list< chemical::ResidueTypeCOP >::const_iterator ResidueTypeCOPListConstIter;
public:
	virtual ~ResidueLevelTask();

	virtual
	ExtraRotSample
	extrachi_sample_level(
		bool buried,
		int chi,
		chemical::ResidueType const & concrete_residue
	) const = 0;

	virtual void initialize_from_command_line() = 0;
	virtual void initialize_extra_rotamer_flags_from_command_line() = 0;

	virtual void initialize_from_options( utility::options::OptionCollection const & options ) = 0;
	virtual void initialize_extra_rotamer_flags_from_options( utility::options::OptionCollection const & options ) = 0;

	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

	virtual void or_include_current( bool include_current ) = 0;
	virtual bool include_current() const = 0;

	virtual void add_behavior( std::string const & behavior ) = 0;
	virtual bool has_behavior( std::string const & behavior ) const = 0;
	virtual bool has_behavior() const = 0;

	virtual void target_type( chemical::ResidueTypeCOP type ) = 0;
	virtual void target_type( chemical::AA aa ) = 0;
	virtual void target_type( std::string name ) = 0;

	virtual void or_adducts( bool setting ) = 0;
	virtual bool adducts() const = 0;

	virtual void or_ex1( bool ex1 ) = 0;
	virtual void or_ex2( bool ex2 ) = 0;
	virtual void or_ex3( bool ex3 ) = 0;
	virtual void or_ex4( bool ex4 ) = 0;

	virtual void or_ex1_sample_level( ExtraRotSample ex1_sample_level ) = 0;
	virtual void or_ex2_sample_level( ExtraRotSample ex2_sample_level ) = 0;
	virtual void or_ex3_sample_level( ExtraRotSample ex3_sample_level ) = 0;
	virtual void or_ex4_sample_level( ExtraRotSample ex4_sample_level ) = 0;

	virtual void or_ex1aro( bool ex1aro ) = 0;
	virtual void or_ex2aro( bool ex2aro ) = 0;
	virtual void or_ex1aro_exposed( bool ex1aro_exposed ) = 0;
	virtual void or_ex2aro_exposed( bool ex2aro_exposed ) = 0;

	virtual void or_ex1aro_sample_level( ExtraRotSample ex1aro_sample_level ) = 0;
	virtual void or_ex2aro_sample_level( ExtraRotSample ex2aro_only_sample_level ) = 0;
	virtual void or_ex1aro_exposed_sample_level( ExtraRotSample ex1aro_exposed_sample_level ) = 0;
	virtual void or_ex2aro_exposed_sample_level( ExtraRotSample ex2aro_exposed_sample_level ) = 0;

	virtual void or_exdna_sample_level( ExtraRotSample exdna_sample_level ) = 0;

	virtual void or_operate_on_ex1( bool operate ) = 0;
	virtual void or_operate_on_ex2( bool operate ) = 0;
	virtual void or_operate_on_ex3( bool operate ) = 0;
	virtual void or_operate_on_ex4( bool operate ) = 0;

	virtual bool ex1() const = 0;
	virtual bool ex2() const = 0;
	virtual bool ex3() const = 0;
	virtual bool ex4() const = 0;

	virtual ExtraRotSample ex1_sample_level() const = 0;
	virtual ExtraRotSample ex2_sample_level() const = 0;
	virtual ExtraRotSample ex3_sample_level() const = 0;
	virtual ExtraRotSample ex4_sample_level() const = 0;

	virtual bool ex1aro() const = 0;
	virtual bool ex2aro() const = 0;
	virtual bool ex1aro_exposed() const = 0;
	virtual bool ex2aro_exposed() const = 0;

	virtual ExtraRotSample ex1aro_sample_level() const = 0;
	virtual ExtraRotSample ex2aro_sample_level() const = 0;
	virtual ExtraRotSample ex1aro_exposed_sample_level() const = 0;
	virtual ExtraRotSample ex2aro_exposed_sample_level() const = 0;

	virtual ExtraRotSample exdna_sample_level() const = 0;

	virtual bool operate_on_ex1() const = 0;
	virtual bool operate_on_ex2() const = 0;
	virtual bool operate_on_ex3() const = 0;
	virtual bool operate_on_ex4() const = 0;

	virtual void sample_proton_chi( bool setting ) = 0;
	virtual bool sample_proton_chi() const = 0;

	virtual void or_optimize_h( bool setting ) = 0;
	virtual bool optimize_h() const = 0;
	virtual void or_preserve_c_beta( bool setting ) = 0;
	virtual bool preserve_c_beta() const = 0;
	virtual void or_flip_HNQ( bool setting ) = 0;
	virtual bool flip_HNQ() const = 0;
	virtual void or_fix_his_tautomer( bool setting ) = 0;
	virtual bool fix_his_tautomer() const = 0;

	virtual void or_include_virtual_side_chain( bool include_virtual_side_chain ) = 0;
	virtual bool include_virtual_side_chain() const = 0;

	virtual void and_extrachi_cutoff( Size num_neighbors_to_be_called_buried ) = 0;

	virtual Size extrachi_cutoff() const = 0;

	// remove all ResidueTypes from the list of allowed residue types
	virtual void prevent_repacking() = 0;

	// contract (and) the list of available aas for canonical aa's
	// if an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code.
	virtual void restrict_absent_canonical_aas( utility::vector1< bool > const & ) = 0;

	// if an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code.
	// If the mode tag is specified, use that it when outputing the task mode in the resfile format rather than PIKAA or NOTAA.
	virtual void restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas, std::string const & mode ) = 0;

	//Same behavior restrict_absent_canonical_aas except that it always allows the native aa at a position even if it is not included in the allowed residues
	virtual void restrict_nonnative_canonical_aas( utility::vector1< bool > const & allowed_aas) = 0;


	/// @brief disables designing to nucleic acid residues not in the passed list
	virtual void restrict_absent_nas( utility::vector1< chemical::AA > const & keep_nas ) = 0;

	// only let this residue repack -- prevent redesign
	virtual void restrict_to_repacking() = 0;

	/// @brief
	virtual bool is_original_type( chemical::ResidueTypeCOP type ) const = 0;

	/// @brief
	virtual chemical::ResidueTypeSetCOP get_original_residue_set() const = 0;

	/// @brief
	virtual chemical::AA const & get_original_residue() const = 0;

	/// @brief expand (or) the list of available residue types for non-cannonicals
	virtual void allow_noncanonical_aa(
		std::string const & interchangeability_group,
		chemical::ResidueTypeSet const & residue_set // who gives this reference to the rlt?; maybe rlt holds a rts (what is an rlt?)
	) = 0;

	/// @brief expand (or) the list of available residue types for non-cannonicals.  Assumes same restypeset as original residue
	virtual void allow_noncanonical_aa( std::string const & aaname ) = 0;

	/// @brief explicitly allow a NCAA; assumes same ResidueTypeSet as original_residue_type_
	virtual void allow_noncanonical_aa( chemical::AA aa ) = 0;

	/// @brief explicitly disallow all NCAAs
	virtual void disallow_noncanonical_aas() = 0;

	// expand (or) the list of rsdtypes by including types with this aa and which variant-match the original rsd
	virtual
	void
	allow_aa( chemical::AA const & aa ) = 0;

	virtual ResidueTypeCOPList const & allowed_residue_types() const = 0;
	virtual ResidueTypeCOPListConstIter allowed_residue_types_begin() const = 0;
	virtual ResidueTypeCOPListConstIter allowed_residue_types_end() const = 0;
	virtual chemical::ResidueTypeCOP target_type() const = 0;

	virtual void print_allowed_types( std::ostream & os ) const = 0;

	virtual bool being_designed() const = 0; // is this residue up for design?
	virtual bool being_packed() const = 0; // is this residue being modified at all by the packer

	/// @brief ONLY for the RESET command in resfiles: completely reset this position.
	/// @details This does several things.  It:
	/// - Removes all noncanonicals allowed at this position.
	/// - Resets the list of allowed canonicals to the 20 standard canonicals.
	/// - Resets the designability of this position (design allowed).
	/// - Resets the repacking of this position (repacking allowed).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	virtual void reset() = 0;

	virtual
	rotamer_set::RotamerOperations const &
	rotamer_operations() const = 0;

	virtual
	void
	append_rotamer_operation(
		rotamer_set::RotamerOperationOP rotop
	) = 0;

	virtual rotamer_set::RotSetOperationListIterator
	rotamer_set_operation_begin() const = 0;

	virtual rotamer_set::RotSetOperationListIterator
	rotamer_set_operation_end() const = 0;

	virtual
	void
	append_rotamerset_operation(
		rotamer_set::RotamerSetOperationOP rotsetop
	) = 0;

	virtual
	std::string
	command_string() const = 0;

	virtual rna::RNA_ResidueLevelTask const & rna_task() const = 0;

	virtual rna::RNA_ResidueLevelTask & nonconst_rna_task() = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace task
} //namespace pack
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_ResidueLevelTask )
#endif // SERIALIZATION


#endif
