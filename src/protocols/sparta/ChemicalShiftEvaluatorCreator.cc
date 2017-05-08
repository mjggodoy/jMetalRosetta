// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/evaluation/ChemicalShiftsEvaluatorCreator.hh
/// @brief  Header for ChemicalShiftsEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/sparta/ChemicalShiftEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/sparta/ChemicalShiftEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// due to template function
#include <core/io/silent/SilentStruct.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/vector0.hh>

//Auto Headers


static THREAD_LOCAL basic::Tracer tr( "protocols.sparta.ChemicalShiftEvaluatorCreator" );

namespace protocols {
namespace sparta {

ChemicalShiftEvaluatorCreator::~ChemicalShiftEvaluatorCreator() = default;

void ChemicalShiftEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::chemical_shifts );

}

void ChemicalShiftEvaluatorCreator::add_evaluators( protocols::evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::evaluation;


	if ( option[ OptionKeys::evaluation::chemical_shifts ].user() ) {
		typedef utility::vector1< std::string > CSVector;
		CSVector const& cs_shifts( option[ OptionKeys::evaluation::chemical_shifts ]() );

		for ( auto it=cs_shifts.begin(); it!=cs_shifts.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != cs_shifts.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
					"need to specify dupletss <cs_shifts> <column> with option -evaluation:chemical_shifts   last read: "+fname );
			}
			eval.add_evaluation( PoseEvaluatorOP( new ChemicalShiftEvaluator( column, fname ) ) );
		}
	}


}

std::string ChemicalShiftEvaluatorCreator::type_name() const {
	return "ChemicalShiftEvaluator";
}

} //namespace
} //namespace
