// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/ParsedProtocol.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Modified this to facilitate use of ParsedProtocols to combine movers and filters in code outside of a RosettaScripts context.

#ifndef INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH
#define INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputterObserver.hh>

#include <core/types.hh>
#include <utility>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/exit.hh>
#include <protocols/moves/ResId.hh>

// C++ headers
#include <string>

// Unit headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace rosetta_scripts {

class ParsedProtocol :
	public protocols::moves::Mover,
	public protocols::moves::ResId,
	public protocols::jd2::JobOutputterObserver
{
public:
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	class MoverFilterPair {
		// OFL want to have more state in the MoverFilterPair -- transform from std:pair into a class, but keep first and second members intact
		// JRP it would sure be nice if this operated like a real class. Calls to first.first suck to read.
	public:
		MoverFilterPair( protocols::moves::MoverOP mover,
			std::string const& mover_name,
			protocols::filters::FilterOP filter,
			bool report_filter_at_end = false ) :
			first( std::make_pair( mover, mover_name ) ),
			second(std::move( filter )),
			report_filter_at_end_( report_filter_at_end )
		{}
		protocols::filters::Filter const& filter() const { return *second; }

		std::pair< protocols::moves::MoverOP, std::string > first;
		protocols::filters::FilterOP second;
		bool report_filter_at_end_;
	};

	typedef utility::vector1< MoverFilterPair > MoverFilterVector;
	typedef MoverFilterVector::iterator iterator;
	typedef MoverFilterVector::const_iterator const_iterator;

public:
	ParsedProtocol();
	~ParsedProtocol() override;
	void apply( Pose & pose ) override;
	core::pose::PoseOP get_additional_output( ) override;
	// XRW TEMP  std::string get_name() const override;
	void final_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );
	core::scoring::ScoreFunctionCOP final_scorefxn() const;
	void final_score(core::pose::Pose & pose) const;
	void report_all( Pose const & pose ) const; // cycles over all filter->report methods to output their values to a common stream.
	void report_filters_to_pose( Pose & pose ); // as above but reports to pose DataCache
	void report_filters_to_job( Pose const & pose ) const;  // as above but reports to job object
	//as above but is called directly from JobOutputter via Observer pattern

	void add_values_to_job( Pose const & pose, protocols::jd2::Job & ) const override;


	// void report_all_sm( std::map< std::string, core::Real > & score_map, Pose const & pose ) const; // ditto, but outputs filter values into score_map object
	protocols::moves::MoverCOP get_mover( core::Size const mover_number ) const {
		runtime_assert( movers_.size() >= mover_number && mover_number > 0 );
		return( movers_[ mover_number ].first.first );
	}
	MoverFilterPair get_mover_filter_pair( core::Size const mover_number ) const {
		runtime_assert( movers_.size() >= mover_number && mover_number > 0 );
		return( movers_[ mover_number ] );
	}

	/// @brief Add a mover-filter pair.
	/// @details Indended for use OUTSIDE of a RosettaScripts context.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void add_mover_filter_pair(
		protocols::moves::MoverOP mover,
		std::string const &mover_name,
		protocols::filters::FilterOP filter,
		bool const report_filter_at_end=false
	);

	void set_resid( core::Size const resid ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new ParsedProtocol ); }
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override; // this is defined as public here, b/c I need to circumvent the name-check, since this is called both by the Movers section (as ParsedProtocol) and the PROTOCOLS section.
	void clear() { movers_.clear(); }
	std::string mode() const{ return mode_; }
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	void apply_probability( utility::vector1< core::Real > const & a );
	utility::vector1< core::Real > apply_probability();
	core::Size size() { return movers_.size(); }
	core::Size last_attempted_mover_idx() { return last_attempted_mover_idx_; }
	void last_attempted_mover_idx( core::Size const s ){ last_attempted_mover_idx_ = s;}
	bool report_call_order() const { return report_call_order_; }
	void report_call_order( bool const c ) { report_call_order_ = c; }
	std::string call_order() const{ return call_order_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void finish_protocol(Pose & pose);

	/// @brief apply the mover of the pair
	void apply_mover(Pose & pose, MoverFilterPair const & mover_pair);

	/// @brief apply the filter of the pair
	bool apply_filter(Pose & pose, MoverFilterPair const & mover_pair);
	void sequence_protocol(Pose & pose, MoverFilterVector::const_iterator mover_it_in);

	void random_order_protocol(Pose & pose);
	void random_single_protocol(Pose & pose);

private:

	MoverFilterVector movers_;
	core::scoring::ScoreFunctionCOP final_scorefxn_;
	std::string mode_;
	utility::vector1< core::Real > apply_probability_; // if mode_="single_random", assigns a probability of execution to each mover/filter pair. Defaults to equal probabilities to all.
	core::Size last_attempted_mover_idx_; //index to last attempted mover; useful for adaptive monte carlo
	bool report_call_order_; //dflt false; At the end of the run, write to out the sequence of mover/filter calls (good for stochastic application
	std::string call_order_; // saved call order, not writeable
	protocols::moves::MoverOP last_mover_;
	bool resume_support_;

};

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH
