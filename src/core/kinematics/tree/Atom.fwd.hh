// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/tree/Atom.fwd.hh
/// @brief  Kinematics Atom forward declarations header
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_tree_Atom_fwd_hh
#define INCLUDED_core_kinematics_tree_Atom_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace kinematics {
namespace tree {

// Forward
class Atom;
class JumpAtom;

// Types
typedef  utility::pointer::weak_ptr< Atom >  AtomAP;
typedef  utility::pointer::weak_ptr< Atom const >  AtomCAP;

typedef  utility::pointer::shared_ptr< Atom >  AtomOP;
typedef  utility::pointer::shared_ptr< Atom const >  AtomCOP;

typedef  utility::pointer::weak_ptr< JumpAtom >  JumpAtomAP;
typedef  utility::pointer::weak_ptr< JumpAtom const >  JumpAtomCAP;

typedef  utility::pointer::shared_ptr< JumpAtom >  JumpAtomOP;
typedef  utility::pointer::shared_ptr< JumpAtom const >  JumpAtomCOP;

}
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_Atom_FWD_HH
