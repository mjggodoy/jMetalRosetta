// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/sequence/BasicSequenceCreator.hh
/// @brief  Base class for BasicSequenceCreators for the BasicSequence load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_sequence_BasicSequenceCreators_hh
#define INCLUDED_core_sequence_BasicSequenceCreators_hh

#include <core/sequence/SequenceCreator.hh>


#include <core/types.hh>
#include <utility/vector1.hh>


namespace core {
namespace sequence {

class SimpleSequenceCreator : public SequenceCreator {
public:
	SimpleSequenceCreator();
	~SimpleSequenceCreator() override;

	SequenceOP create_sequence() const override;
	std::string keyname() const override;
};

class SequenceProfileCreator : public SequenceCreator {
public:
	SequenceProfileCreator();
	~SequenceProfileCreator() override;

	SequenceOP create_sequence() const override;
	std::string keyname() const override;
};

class SequenceCouplingCreator : public SequenceCreator {
public:
	SequenceCouplingCreator();
	~SequenceCouplingCreator() override;

	SequenceOP create_sequence() const override;
	std::string keyname() const override;
};

class CompositeSequenceCreator : public SequenceCreator {
public:
	CompositeSequenceCreator();
	~CompositeSequenceCreator() override;

	SequenceOP create_sequence() const override;
	std::string keyname() const override;
};

class ChemicalShiftSequenceCreator : public SequenceCreator {
public:
	ChemicalShiftSequenceCreator();
	~ChemicalShiftSequenceCreator() override;

	SequenceOP create_sequence() const override;
	std::string keyname() const override;
};

} //namespace sequence
} //namespace core

#endif
