// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/ConnectChainsMover.hh
/// @brief connects chains and returns the connection with lowest rmsd
/// @author TJ Brunette tjbrunette@gmail.com
///
#ifndef INCLUDED_protocols_pose_length_moves_ConnectChainsMover_hh
#define INCLUDED_protocols_pose_length_moves_ConnectChainsMover_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/ConnectChainsMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

//#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ctime>
#include <boost/range/algorithm/count.hpp>


namespace protocols {
namespace pose_length_moves {

struct Chain{
public:
	core::Real rmsd;
	core::pose::PoseOP poseOP;
	Chain(core::pose::PoseOP poseOP_i, core::Real rmsd_i){
		rmsd = rmsd_i;
		poseOP = poseOP_i;
	}
};

class ConnectChainsMover : public protocols::moves::Mover {


public:
	ConnectChainsMover();
	void parse_input(utility::vector1<std::string> & individual_chains,utility::vector1< utility::vector1 <std::string> > & chains_in_poses);
	std::map<std::string, Chain> generate_connected_chains(core::pose::Pose const pose,utility::vector1<std::string> individual_chains);
	void assemble_missing_chain(std::map<std::string, Chain> & connected_chains,std::string chain_assembled,std::string chain_remainder);
	void generate_best_final_pose(core::pose::Pose & pose,utility::vector1< utility::vector1 <std::string> > chains_in_poses,std::map<std::string, Chain> connected_chains);
	void apply( Pose & pose ) override;
	// XRW TEMP  virtual std::string get_name() const;
	moves::MoverOP clone() const override { return moves::MoverOP( new ConnectChainsMover( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	int resAdjustmentStartLow_;
	int resAdjustmentStartHigh_;
	int resAdjustmentStopLow_;
	int resAdjustmentStopHigh_;
	int resAdjustmentStartLow_sheet_;
	int resAdjustmentStartHigh_sheet_;
	int resAdjustmentStopLow_sheet_;
	int resAdjustmentStopHigh_sheet_;
	core::Size loopLengthRangeLow_;
	core::Size loopLengthRangeHigh_;
	core::Real rmsThreshold_;
	std::string output_chains_;

	core::indexed_structure_store::SSHashedFragmentStore * SSHashedFragmentStore_;
};

/////////////
// used to sort chains so the least elements are on top of this list
/////////////
class chain_lt : public std::binary_function<std::string, std::string, bool> {
public:
	bool operator()(std::string x, std::string y) {
		core::Size x_chain_length = boost::count(x,'+');
		core::Size y_chain_length = boost::count(y,'+');
		if ( x_chain_length < y_chain_length ) return true;
		else if ( x_chain_length > y_chain_length ) return false;
		else return(x<y);
	}
};

} // pose_length_moves
} // protocols

#endif
