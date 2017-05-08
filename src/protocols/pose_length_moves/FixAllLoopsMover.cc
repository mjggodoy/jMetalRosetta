// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/FixAllLoopsMover.fwd.hh
/// @brief connects chains using a very fast RMSD lookback. only works for chains <5 residues. Designed to make loops look within .4 RMSD to naturally occuring loops
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/pose_length_moves/FixAllLoopsMover.hh>
#include <protocols/pose_length_moves/FixAllLoopsMoverCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>

#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/pose_length_moves/NearNativeLoopCloser.hh>

// Core Headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <ctime>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.pose_length_moves.FixAllLoopsMover" );

namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;



FixAllLoopsMover::FixAllLoopsMover():moves::Mover("FixAllLoopsMover"){}

protocols::loops::Loops FixAllLoopsMover::get_loops(core::pose::Pose const & pose){
	protocols::loops::Loops pose_loops;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	string lastSecStruct = dssp_string.substr(0,1);
	Size startLoop = 0;
	Size endLoop = 0;
	if ( dssp_string.substr(0,1) == "L" ) {
		startLoop = 1;
	}
	for ( core::Size ii = 2; ii <= pose.size(); ++ii ) {
		if ( dssp_string.substr(ii-1,1) == "L" && lastSecStruct != "L" ) {
			startLoop = ii;
		}
		if ( dssp_string.substr(ii-1,1) != "L" && lastSecStruct == "L" ) {
			endLoop = ii-1;
			if ( (startLoop != 1) && (endLoop!=pose.size()) ) {
				pose_loops.add_loop(startLoop,endLoop);
			}
		}
		lastSecStruct = dssp_string.substr(ii-1,1);
	}
	return(pose_loops);
}

void FixAllLoopsMover::apply(core::pose::Pose & pose) {
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("***fixing loops on a symmetric structure would be a bad idea because the fold tree is not properly maintained");
	}
	protocols::loops::Loops pose_loops = get_loops(pose);
	bool failure = false;
	if ( firstResidue_==1 && lastResidue_>pose.size() ) {
		TR << "working on all residues" << std::endl;
	} else {
		TR << "working on residues" << firstResidue_ << " to " << lastResidue_ << std::endl;
	}
	TR << "loop RMSD before optimization" << std::endl;
	for ( Size ii=pose_loops.num_loop(); ii>=1; --ii ) {
		vector1<Size> resids;
		resids.push_back(pose_loops[ii].start()-3);
		resids.push_back(pose_loops[ii].start()-1);
		resids.push_back(pose_loops[ii].start()+1);
		resids.push_back(pose_loops[ii].start()+2);
		Real loop_rmsd = SSHashedFragmentStore_->max_rmsd_in_region(pose,resids);
		TR << "Loop" << pose_loops[ii].start() << "-" << pose_loops[ii].stop() << " rmsd:" << loop_rmsd << std::endl;
	}
	for ( Size ii=pose_loops.num_loop(); ii>=1; --ii ) {
		if ( pose_loops[ii].start()>=firstResidue_ && pose_loops[ii].stop()<=lastResidue_ ) {
			TR << "working on " << pose_loops[ii].start() << "-" <<  pose_loops[ii].stop() << std::endl;
			vector1<Size> resids;
			resids.push_back(pose_loops[ii].start()-3);
			resids.push_back(pose_loops[ii].start()-1);
			resids.push_back(pose_loops[ii].start()+1);
			resids.push_back(pose_loops[ii].start()+2);
			Real loop_rmsd = SSHashedFragmentStore_->max_rmsd_in_region(pose,resids);
			if ( loop_rmsd > rmsThreshold_ ) {
				NearNativeLoopCloserOP loopCloserOP(new NearNativeLoopCloser(resAdjustmentStartLow_,resAdjustmentStartHigh_,resAdjustmentStopLow_,resAdjustmentStopHigh_,resAdjustmentStartLow_sheet_,resAdjustmentStartHigh_sheet_,resAdjustmentStopLow_sheet_,resAdjustmentStopHigh_sheet_,loopLengthRangeLow_,loopLengthRangeHigh_,pose_loops[ii].start()-1,pose_loops[ii].stop()+1,'A','A',rmsThreshold_,max_vdw_change_,true,ideal_,true));
				loopCloserOP->apply(pose);
				if ( reject_failed_loops_ ) {
					if ( loopCloserOP->get_last_move_status()!=protocols::moves::MS_SUCCESS ) {
						failure = true;
						TR << "no closure below threshold found" << std::endl;
						set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
					}
				} else {
					if ( loopCloserOP->get_last_move_status()!=protocols::moves::MS_SUCCESS ) {
						TR << "no closure below threshold found but being ignored" << std::endl;
					}
				}
			}
		}
	}
	if ( !failure ) {
		set_last_move_status(protocols::moves::MS_SUCCESS);
		pose_loops = get_loops(pose);
		TR << "loop RMSD after optimization (note:Residue numbers may have changed)" << std::endl;
		for ( Size ii=pose_loops.num_loop(); ii>=1; --ii ) {
			vector1<Size> resids;
			resids.push_back(pose_loops[ii].start()-3);
			resids.push_back(pose_loops[ii].start()-1);
			resids.push_back(pose_loops[ii].start()+1);
			resids.push_back(pose_loops[ii].start()+2);
			Real loop_rmsd = SSHashedFragmentStore_->max_rmsd_in_region(pose,resids);
			TR << "Loop" << pose_loops[ii].start() << "-" << pose_loops[ii].stop() << " rmsd:" << loop_rmsd << std::endl;
		}
	}
	//  //time_t end_time = time(NULL);
	//  //std::cout << "total_time" << end_time-start_time_ << std::endl;
}

void
FixAllLoopsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	//start_time_ = time(NULL);
	std::string loopLengthRange( tag->getOption< std::string >( "loopLengthRange", "1,4") );
	rmsThreshold_ = tag->getOption< core::Real >( "RMSthreshold", 0.4 );
	std::string resAdjustmentRange1( tag->getOption< std::string >( "resAdjustmentRangeSide1", "-3,3") );
	std::string resAdjustmentRange2( tag->getOption< std::string >( "resAdjustmentRangeSide2","-3,3") );
	std::string resAdjustmentRange1_sheet( tag->getOption< std::string >( "resAdjustmentRangeSide1_sheet", "-1,1") );
	std::string resAdjustmentRange2_sheet( tag->getOption< std::string >( "resAdjustmentRangeSide2_sheet","-1,1") );
	std::string residueRange(tag->getOption<std::string> ( "residue_range","1,999999") );
	max_vdw_change_ = tag->getOption<core::Real>("max_vdw_change",10.0);
	ideal_ = tag->getOption<bool>("ideal",false);
	reject_failed_loops_ = tag->getOption<bool>("reject_failed_loops",true);//If one loop fails just skip that one.
	utility::vector1< std::string > resAdjustmentRange1_split( utility::string_split( resAdjustmentRange1 , ',' ) );
	utility::vector1< std::string > resAdjustmentRange2_split( utility::string_split( resAdjustmentRange2 , ',' ) );
	utility::vector1< std::string > resAdjustmentRange1_sheet_split( utility::string_split( resAdjustmentRange1_sheet , ',' ) );
	utility::vector1< std::string > resAdjustmentRange2_sheet_split( utility::string_split( resAdjustmentRange2_sheet , ',' ) );
	utility::vector1< std::string > loopLengthRange_split( utility::string_split( loopLengthRange , ',' ) );
	utility::vector1< std::string > residueRange_split( utility::string_split( residueRange , ',' ) );
	if ( resAdjustmentRange1_split.size()==2 ) {
		resAdjustmentStartLow_ = atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange1_split[2].c_str());
	}
	if ( resAdjustmentRange1_split.size()==1 ) {
		resAdjustmentStartLow_= atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_= atoi(resAdjustmentRange1_split[1].c_str());
	}
	if ( resAdjustmentRange2_split.size()==2 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStopHigh_ = atoi(resAdjustmentRange2_split[2].c_str());
	}
	if ( resAdjustmentRange2_split.size()==1 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange2_split[1].c_str());
	}
	if ( resAdjustmentRange1_sheet_split.size()==2 ) {
		resAdjustmentStartLow_sheet_ = atoi(resAdjustmentRange1_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_ = atoi(resAdjustmentRange1_sheet_split[2].c_str());
	}
	if ( resAdjustmentRange1_sheet_split.size()==1 ) {
		resAdjustmentStartLow_sheet_= atoi(resAdjustmentRange1_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_= atoi(resAdjustmentRange1_sheet_split[1].c_str());
	}
	if ( resAdjustmentRange2_sheet_split.size()==2 ) {
		resAdjustmentStopLow_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
		resAdjustmentStopHigh_sheet_ = atoi(resAdjustmentRange2_sheet_split[2].c_str());
	}
	if ( resAdjustmentRange2_sheet_split.size()==1 ) {
		resAdjustmentStopLow_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
	}
	if ( loopLengthRange_split.size()==2 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[2].c_str());
	}
	if ( loopLengthRange_split.size()==1 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[1].c_str());
	}
	if ( residueRange_split.size()==2 ) {
		firstResidue_ = atoi(residueRange_split[1].c_str());
		lastResidue_ = atoi(residueRange_split[2].c_str());
	}
	SSHashedFragmentStore_ = core::indexed_structure_store::SSHashedFragmentStore::get_instance();
	SSHashedFragmentStore_->set_threshold_distance(rmsThreshold_);
	SSHashedFragmentStore_->init_SS_stub_HashedFragmentStore();
	TR << "database loaded!!" << std::endl;
	std::cout << resAdjustmentStartLow_ <<"," << resAdjustmentStartHigh_ << ",:," << resAdjustmentStopLow_ << "," << resAdjustmentStopHigh_ << ",:," << loopLengthRangeLow_ <<"," << loopLengthRangeHigh_ << std::endl;
}

std::string FixAllLoopsMover::get_name() const {
	return mover_name();
}

std::string FixAllLoopsMover::mover_name() {
	return "FixAllLoopsMover";
}

void FixAllLoopsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"loopLengthRange", xs_string,
		"XSD_XRW: TO DO", "1,4");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"RMSthreshold", xsct_real,
		"XSD_XRW: TO DO", "0.4");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide1", xs_string,
		"XSD_XRW: TO DO", "-3,3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide2", xs_string,
		"XSD_XRW: TO DO", "-3,3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide1_sheet", xs_string,
		"XSD_XRW: TO DO", "-1,1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide2_sheet", xs_string,
		"XSD_XRW: TO DO", "-1,1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"residue_range", xs_string,
		"XSD_XRW: TO DO", "1,999999");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_vdw_change", xsct_real,
		"XSD_XRW: TO DO", "10.0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"ideal", xsct_rosetta_bool,
		"XSD_XRW: TO DO", "false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"reject_failed_loops", xsct_rosetta_bool,
		"XSD_XRW: TO DO", "true");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Connects chains using a very fast RMSD lookback. "
		"Only works for chains less than 5 residues. Designed to make loops look within .4 RMSD "
		"to naturally occuring loops",
		attlist );
}

std::string FixAllLoopsMoverCreator::keyname() const {
	return FixAllLoopsMover::mover_name();
}

protocols::moves::MoverOP
FixAllLoopsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FixAllLoopsMover );
}

void FixAllLoopsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FixAllLoopsMover::provide_xml_schema( xsd );
}


}//pose_length_moves
}//protocols
