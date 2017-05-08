// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EtableTrie.fwd.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_elec_electrie_ElecTrie_fwd_hh
#define INCLUDED_core_scoring_elec_electrie_ElecTrie_fwd_hh

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace elec {
namespace electrie {

typedef trie::RotamerTrieBaseOP ElecRotamerTrieOP;
typedef trie::RotamerTrieBaseCOP ElecRotamerTrieCOP;

} // hbtrie
} // hbonds
} // scoring
} // core

#endif
