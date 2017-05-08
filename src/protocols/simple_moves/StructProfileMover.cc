// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/StuctProfileMover.fwd.hh
/// @brief Quickly generates a structure profile.
///
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <protocols/simple_moves/StructProfileMover.hh>
#include <protocols/simple_moves/StructProfileMoverCreator.hh>
#include <protocols/moves/Mover.hh>

#include <basic/database/open.hh>

// Core Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
//
#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>


#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>

#include <core/types.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <map>
#include <set>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.StructProfileMover" );

namespace protocols {
namespace simple_moves {
using namespace ObjexxFCL::format;

using namespace core;
using namespace std;
using utility::vector1;
using namespace core::indexed_structure_store;
using core::pose::Pose;

StructProfileMover::StructProfileMover():moves::Mover("StructProfileMover"){
	aa_order_="ACDEFGHIKLMNPQRSTVWY";
	read_P_AA_SS_cen6();
}

StructProfileMover::StructProfileMover(Real rmsThreshold,Size consider_topN_frags, Real burialWt, bool only_loops , Real allowed_deviation, Real allowed_deviation_loops, bool eliminate_background, bool outputProfile, bool add_csts_to_pose, bool ignore_terminal_res):moves::Mover("StructProfileMover"){
	using namespace core::indexed_structure_store;
	aa_order_="ACDEFGHIKLMNPQRSTVWY";
	read_P_AA_SS_cen6();
	rmsThreshold_ = rmsThreshold;
	consider_topN_frags_ = consider_topN_frags;
	burialWt_ = burialWt;
	only_loops_ = only_loops;
	allowed_deviation_ = allowed_deviation;
	allowed_deviation_loops_ = allowed_deviation_loops;
	eliminate_background_ = eliminate_background;
	outputProfile_ = outputProfile;
	add_csts_to_pose_ = add_csts_to_pose;
	ignore_terminal_res_ = ignore_terminal_res;
	cenType_=6;
	SSHashedFragmentStore_ = SSHashedFragmentStore::get_instance();
	SSHashedFragmentStore_->set_threshold_distance(rmsThreshold_);
}

// XRW TEMP std::string StructProfileMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return StructProfileMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string StructProfileMover::mover_name(){
// XRW TEMP  return "StructProfileMover";
// XRW TEMP }


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StructProfileMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new StructProfileMover );
// XRW TEMP }



struct Hit{
	Real cend;
	Real cend_norm;
	Real rmsd;
	Real rmsd_norm;
	Real score;
	std::string aa;
	Hit(Real cend_i,Real rmsd_i,std::string aa_i){
		cend = cend_i;
		rmsd = rmsd_i;
		aa = aa_i;
	}
	void print(){
		std::cout << "cend " << cend << " cend_norm " << cend_norm << " rmsd " << rmsd << " rmsd_norm " << rmsd_norm << " score " << score << std::endl;
	}
};

struct less_then_match_rmsd
{
	inline bool operator() (const Hit& struct1, const Hit& struct2)
	{
		return (struct1.score < struct2.score);
	}
};

Size StructProfileMover::ss_type_convert(char ss_type){
	if ( ss_type == 'H' ) {
		return 1;
	} else if ( ss_type == 'L' ) {
		return 2;
	} else {
		return 3;
	}
}

void StructProfileMover::read_P_AA_SS_cen6(){
	utility::io::izstream stream;
	basic::database::open( stream,"scoring/score_functions/P_AA_SS_cen6/P_AA_SS_cen6.txt");
	Size SS_TYPES = 3;
	Size BURIAL_TYPES = 10;
	Size aa_types = aa_order_.size();
	P_AA_SS_burial_.resize(SS_TYPES);
	for ( Size ii = 1; ii <= SS_TYPES; ++ii ) {
		P_AA_SS_burial_[ii].resize(BURIAL_TYPES);
		for ( Size jj = 1; jj <= BURIAL_TYPES; ++jj ) {
			P_AA_SS_burial_[ii][jj].resize(aa_types);
		}
	}
	char ss_type;
	Size burial_type;
	Real tmp_probability;
	string line;
	stream >> skip;
	while ( getline(stream,line) ) {
		std::istringstream l(line);
		l >> ss_type >> burial_type;
		for ( Size ii=1; ii<=aa_types; ++ii ) {
			l >> tmp_probability >> skip(1);
			P_AA_SS_burial_[ss_type_convert(ss_type)][burial_type][ii]=tmp_probability;
			debug_assert( ( tmp_probability >= Probability( 0.0 ) ) && ( tmp_probability <= Probability( 1.0 ) ) );
			debug_assert( burial_type <= BURIAL_TYPES);
		}
	}
	stream.close();
}


vector1<std::string> StructProfileMover::get_closest_sequence_at_res(core::pose::Pose const pose, Size res,vector1<Real> cenList){
	vector1<string> top_hits_aa;
	//I want to normalize the rmsd & burial
	vector1<vector<Real> > hits_cen;
	vector1<Real> hits_rms;
	vector1<std::string> hits_aa;
	SSHashedFragmentStore_->get_hits_below_rms(pose,res,rmsThreshold_,hits_cen,hits_rms,hits_aa);
	if ( hits_cen.size() == 0 ) {
		return top_hits_aa;
	}
	vector1<Hit>hits;
	for ( Size ii=1; ii<=hits_cen.size(); ++ii ) {
		Real cen_deviation = get_cen_deviation(hits_cen[ii],cenList);
		struct Hit result_tmp(cen_deviation,hits_rms[ii],hits_aa[ii]);
		hits.push_back(result_tmp);
	}
	//step1 get max and min cen and rmsd
	Real maxRmsd = -9999;
	Real minRmsd = 9999;
	Real maxCend = -9999;
	Real minCend = 9999;
	for ( Size ii=1; ii<=hits.size(); ++ii ) {
		if ( hits[ii].cend>maxCend ) {
			maxCend = hits[ii].cend;
		}
		if ( hits[ii].cend<minCend ) {
			minCend = hits[ii].cend;
		}
		if ( hits[ii].rmsd>maxRmsd ) {
			maxRmsd = hits[ii].rmsd;
		}
		if ( hits[ii].rmsd<minRmsd ) {
			minRmsd = hits[ii].rmsd;
		}
	}
	//step2 set normatlized rmsd cend and set score
	for ( Size ii=1; ii<=hits.size(); ++ii ) {
		hits[ii].cend_norm = 1-(maxCend-hits[ii].cend)/(maxCend-minCend);
		hits[ii].rmsd_norm = 1-(maxRmsd-hits[ii].rmsd)/(maxRmsd-minRmsd);
		hits[ii].score = hits[ii].cend_norm*(burialWt_)+hits[ii].rmsd_norm*(1-burialWt_);
	}
	//step3 sort array based on score
	if ( consider_topN_frags_ < hits.size() ) {
		std::sort(hits.begin(), hits.end(), less_then_match_rmsd());
	}
	//for(Size ii=1; ii<=lookupResultsPlusV.size(); ++ii){
	// lookupResultsPlusV[ii].print();
	//}
	//step4 get AA from top hits
	for ( Size ii=1; ii<=hits.size()&&ii<consider_topN_frags_; ++ii ) {
		top_hits_aa.push_back(hits[ii].aa);
	}
	return(top_hits_aa);
}

vector1<vector1<std::string> > StructProfileMover::get_closest_sequences(core::pose::Pose const pose,vector1<Real> cenList){
	Size fragment_length = SSHashedFragmentStore_->get_fragment_length();
	vector1<vector1<std::string > > all_aa_hits;
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	for ( Size ii=1; ii<nres1-fragment_length+1; ++ii ) {
		vector1<Real>::const_iterator begin =cenList.begin();
		vector1<Real> shortCenList(begin+ii-1, begin+ii+fragment_length-1);
		vector1<std::string> aa_hits = get_closest_sequence_at_res(pose,ii,shortCenList);
		all_aa_hits.push_back(aa_hits);
	}
	return(all_aa_hits);
}

vector1<vector1<Size> > StructProfileMover::generate_counts(vector1<vector1<std::string> > top_frag_sequences,core::pose::Pose const pose){
	//step1---Initialize counts to zero
	vector1<vector1<Size> > counts;
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	for ( Size ii=1; ii<=nres1; ++ii ) {
		vector1<Size> counts_per_res;
		for ( Size jj=0; jj<aa_order_.size(); ++jj ) {
			counts_per_res.push_back(0);
		}
		counts.push_back(counts_per_res);
	}
	//step2--Gather counts for the variables
	for ( Size ii=1; ii<=top_frag_sequences.size(); ++ii ) {//@each residue in pose
		for ( Size jj=1; jj<=top_frag_sequences[ii].size(); ++jj ) {//@#hits per position
			for ( Size kk=0; kk<top_frag_sequences[ii][jj].size(); ++kk ) {//length of string
				char currentChar = top_frag_sequences[ii][jj].at(kk);
				Size charPosition = aa_order_.find(currentChar)+1; //+1 to convert to index 1 array
				Size residueIndex= ii+kk;
				counts[residueIndex][charPosition]++;
			}
		}
	}
	//step3--Print out first position
	/* for(Size ii=1; ii<=pose.size(); ++ii){
	for(Size kk=1; kk<=counts[1].size(); ++kk)
	std::cout << aa_order_.at(kk-1) << ":" << counts[ii][kk] << std::endl;
	std::cout << "____________________________________" << std::endl;
	}
	*/
	return(counts);
}

vector1<vector1<Real> > StructProfileMover::generate_profile_score(vector1<vector1<Size> > res_per_pos,Pose const pose){
	//step1---Get total counts for each position
	vector1<Size> total_cts;
	total_cts.resize(res_per_pos.size(),0);
	for ( Size ii=1; ii<= res_per_pos.size(); ++ii ) {
		for ( Size jj=1; jj<= res_per_pos[ii].size(); ++jj ) {
			total_cts[ii]+=res_per_pos[ii][jj];
		}
	}
	//step2--Generate score
	vector1<vector1<Real> > profile_score;
	for ( Size ii=1; ii<= res_per_pos.size(); ++ii ) {
		vector1<Real> pos_profile_score;
		for ( Size jj=1; jj<= res_per_pos[ii].size(); ++jj ) {
			Real tmp_score = -std::log((Real(res_per_pos[ii][jj]+1))/(Real(total_cts[ii]+20)));
			if ( pose.secstruct(ii) != 'L' && only_loops_ ) {
				tmp_score = 0.0;
			}
			pos_profile_score.push_back(tmp_score);
		}
		profile_score.push_back(pos_profile_score);
	}
	return(profile_score);
}

vector1<vector1<Real> > StructProfileMover::generate_profile_score_wo_background(vector1<vector1<Size> > res_per_pos, vector1<Real> cenList, Pose const pose){
	//step1---Get total counts for each position
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	vector1<Size> total_cts;
	total_cts.resize(res_per_pos.size(),0);
	for ( Size ii=1; ii<= res_per_pos.size(); ++ii ) {
		for ( Size jj=1; jj<= res_per_pos[ii].size(); ++jj ) {
			total_cts[ii]+=res_per_pos[ii][jj];
		}
	}
	//step2--Generate score
	vector1<vector1<Real> > profile_score;
	for ( Size ii=1; ii<= res_per_pos.size(); ++ii ) {
		vector1<Real> pos_profile_score;
		Size ssType = ss_type_convert(pose.secstruct(ii));
		Size burialType = round(cenList[ii]);
		if ( burialType>10 ) { //max for centype 6 was set to 10 atoms with 6ang because of low counts.
			burialType=10;
		}
		for ( Size jj=1; jj<= res_per_pos[ii].size(); ++jj ) {//represents the AA in the same order as in the file.
			Real rmsdProb = (Real(res_per_pos[ii][jj]+1))/(Real(total_cts[ii]+20));
			Size aaType = jj;
			Real backgroundProb = P_AA_SS_burial_[ssType][burialType][aaType];
			Real tmp_score = 0.0;
			if ( pose.secstruct(ii) == 'L' ) {
				if ( (rmsdProb-backgroundProb-allowed_deviation_loops_)>0.0 ) {
					tmp_score = -std::log(rmsdProb-backgroundProb);
				}
			} else {
				if ( (rmsdProb-backgroundProb-allowed_deviation_)>0.0 ) {
					tmp_score = -std::log(rmsdProb-backgroundProb);
				}
			}
			if ( ignore_terminal_res_ && (ii==1 || ii==nres1) ) {
				tmp_score = 0.0;  //phi is 0 in first position and psi and omega are 0 in last position. Can cause odd behavior to the profile so no weight is allowed
			}
			if ( pose.secstruct(ii) != 'L' && only_loops_ ) {
				tmp_score = 0.0;  //0's out when not in loops
			}
			pos_profile_score.push_back(tmp_score);
		}
		profile_score.push_back(pos_profile_score);
	}
	return(profile_score);
}


void StructProfileMover::save_MSAcst_file(vector1<vector1<Real> > profile_score,core::pose::Pose const pose){
	std::string profile_name( "profile" );
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	utility::io::ozstream profile_out(profile_name);
	profile_out << "aa     ";
	for ( char currentChar : aa_order_ ) {
		profile_out << currentChar << "       ";
	}
	profile_out << std::endl;
	for ( Size ii=1; ii<=profile_score.size(); ++ii ) {
		char aa_tmp = pose.residue(ii).name1();
		profile_out << aa_tmp;
		for ( Size jj=1; jj<=profile_score[ii].size(); ++jj ) {
			profile_out << F(8,2,profile_score[ii][jj]);
		}
		profile_out << std::endl;
	}
	profile_out.close();
	std::string msa_name( "MSAcst" );
	utility::io::ozstream msa_out(msa_name);
	for ( Size ii=1; ii<=nres1; ++ii ) {
		msa_out << "SequenceProfile " << ii << " profile" << std::endl;
	}
	msa_out.close();
}

void StructProfileMover::add_MSAcst_to_pose(vector1<vector1<Real> > profile_score,core::pose::Pose & pose){
	using namespace core::sequence;
	SequenceProfileOP profileOP = SequenceProfileOP(new SequenceProfile( profile_score, pose.sequence(), "structProfile" ));
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	for ( core::Size seqpos( 1 ), end( nres1 ); seqpos <= end; ++seqpos ) {
		pose.add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( pose, seqpos, profileOP ) ) ) );
	}
}

Real StructProfileMover::get_cen_deviation(vector<Real> cenListFrag,vector1<Real> cenListModel){
	Real total_deviation = 0;
	for ( Size ii=1; ii<=cenListFrag.size(); ++ii ) {
		total_deviation += (cenListModel[ii]-cenListFrag[ii-1])*(cenListModel[ii]-cenListFrag[ii-1]);
	}
	return(std::sqrt(total_deviation));
}

vector1< Real> StructProfileMover::calc_cenlist(Pose const pose){
	using namespace core::chemical;
	using namespace core::scoring;
	core::pose::PoseOP centroidPose = pose.clone();
	if ( centroidPose->is_fullatom() ) {
		core::util::switch_to_residue_type_set(*centroidPose, core::chemical::CENTROID_t );
	}
	ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		sfcen=core::scoring::symmetry::symmetrize_scorefunction(*sfcen);
	}
	sfcen->score(*centroidPose);
	vector1 <Real> cenlist;
	Size nres1 = centroidPose->size();
	if ( core::pose::symmetry::is_symmetric(*centroidPose) ) {
		nres1 = core::pose::symmetry::symmetry_info(*centroidPose)->num_independent_residues();
	}
	for ( Size ii = 1; ii <= nres1; ++ii ) {
		if ( cenType_ == 6 ) {
			Real fcen6( core::scoring::EnvPairPotential::cenlist_from_pose( *centroidPose ).fcen6(ii));
			cenlist.push_back(fcen6);
		}
		if ( cenType_ == 10 ) {
			Real const fcen10( core::scoring::EnvPairPotential::cenlist_from_pose( *centroidPose ).fcen10(ii));
			cenlist.push_back(fcen10);
		}
		if ( cenType_ == 12 ) {
			Real const fcen12( core::scoring::EnvPairPotential::cenlist_from_pose( *centroidPose ).fcen12(ii));
			cenlist.push_back(fcen12);
		}
	}
	return(cenlist);
}


void StructProfileMover::apply(core::pose::Pose & pose) {
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	vector1<Real> cenList = calc_cenlist(pose);
	vector1<vector1<std::string> > top_frag_sequences = get_closest_sequences(pose,cenList);
	vector1<vector1<Size> > res_per_pos = generate_counts(top_frag_sequences,pose);
	vector1<vector1<Real> > profile_score;
	if ( eliminate_background_ ) {
		profile_score = generate_profile_score_wo_background(res_per_pos,cenList,pose);
	} else {
		profile_score =generate_profile_score(res_per_pos,pose);
	}
	if ( outputProfile_ ) {
		save_MSAcst_file(profile_score,pose);
	}
	if ( add_csts_to_pose_ ) {
		add_MSAcst_to_pose(profile_score,pose);
	}
}


// XRW TEMP std::string StructProfileMover::get_name() const {
// XRW TEMP  return "StructProfileMover";
// XRW TEMP }

void
StructProfileMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	using namespace core::indexed_structure_store;
	rmsThreshold_ = tag->getOption< core::Real >( "RMSthreshold", 0.40 );
	burialWt_ =tag->getOption< Real > ("burialWt", 0.8); //other weight is toward RMSD
	consider_topN_frags_ =tag->getOption< Size > ("consider_topN_frags", 50);
	only_loops_=tag->getOption< bool > ("only_loops",false);
	allowed_deviation_=tag->getOption< Real >("allowed_deviation",0.10);
	allowed_deviation_loops_=tag->getOption< Real >("allowed_deviation_loops",0.10);
	eliminate_background_=tag->getOption< bool >("eliminate_background",true);
	SSHashedFragmentStore_ = SSHashedFragmentStore::get_instance();
	SSHashedFragmentStore_->set_threshold_distance(rmsThreshold_);
	cenType_ = tag->getOption<Size>("cenType",6); //Needs to match the datatabase. Likely I will find one I like and use that so this is an option that shouldn't be modified often
	outputProfile_ = tag->getOption<bool>("outputProfile",false);
	add_csts_to_pose_ = tag->getOption<bool>("add_csts_to_pose",true);
	ignore_terminal_res_ = tag->getOption<bool>("ignore_terminal_residue",true);
}

std::string StructProfileMover::get_name() const {
	return mover_name();
}

std::string StructProfileMover::mover_name() {
	return "StructProfileMover";
}

void StructProfileMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "RMSthreshold", xsct_real, "XRW TO DO", "0.40" )
		+ XMLSchemaAttribute::attribute_w_default( "burialWt", xsct_real, "XRW TO DO", "0.8" )
		+ XMLSchemaAttribute::attribute_w_default( "consider_topN_frags", xsct_non_negative_integer, "XRW TO DO", "50" )
		+ XMLSchemaAttribute::attribute_w_default( "only_loops", xsct_rosetta_bool, "XRW TO DO",  "false" )
		+ XMLSchemaAttribute::attribute_w_default( "allowed_deviation", xsct_real,  "XRW TO DO", "0.10" )
		+ XMLSchemaAttribute::attribute_w_default( "allowed_deviation_loops", xsct_real, "XRW TO DO", "0.10" )
		+ XMLSchemaAttribute::attribute_w_default( "eliminate_background", xsct_rosetta_bool, "XRW TO DO", "true")
		+ XMLSchemaAttribute::attribute_w_default( "cenType", xsct_non_negative_integer, "XRW TO DO",  "6")
		+ XMLSchemaAttribute::attribute_w_default( "outputProfile", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "add_csts_to_pose", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_terminal_residue", xsct_rosetta_bool, "XRW TO DO", "true" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Quickly generates a structure profile", attlist );
}

std::string StructProfileMoverCreator::keyname() const {
	return StructProfileMover::mover_name();
}

protocols::moves::MoverOP
StructProfileMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StructProfileMover );
}

void StructProfileMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StructProfileMover::provide_xml_schema( xsd );
}


}//simple_moves
}//protocols
