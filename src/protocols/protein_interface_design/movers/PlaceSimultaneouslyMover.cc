// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PlaceSimultaneouslyMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Lei Shi (shilei@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PlaceSimultaneouslyMover.hh>
#include <protocols/protein_interface_design/movers/PlaceSimultaneouslyMoverCreator.hh>

// Package headers
#include <protocols/protein_interface_design/movers/PlacementAuctionMover.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>
#include <basic/datacache/DataMap.hh>

#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>
//#include <core/pack/rotamer_trials.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilter.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>

//#include <protocols/docking/DockingProtocol.hh>
#include <protocols/moves/Mover.hh>
#include <core/chemical/ResidueType.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/moves/ResId.hh>
#include <protocols/rosetta_scripts/util.hh>


// Utility Headers
#include <utility/exit.hh>

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/design_utils.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <protocols/protein_interface_design/util.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/scoring/Interface.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <map>
#include <algorithm>
#include <utility>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace core::scoring;
using namespace protocols::protein_interface_design;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.PlaceSimultaneouslyMover" );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

// XRW TEMP std::string
// XRW TEMP PlaceSimultaneouslyMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PlaceSimultaneouslyMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PlaceSimultaneouslyMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PlaceSimultaneouslyMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PlaceSimultaneouslyMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "PlaceSimultaneously";
// XRW TEMP }

PlaceSimultaneouslyMover::~PlaceSimultaneouslyMover() {
}

protocols::moves::MoverOP
PlaceSimultaneouslyMover::clone() const {
	return( protocols::moves::MoverOP( new PlaceSimultaneouslyMover( *this ) ) );
}


/// @details preliminary rb minimization step in the presence of strong hotspot constraints. Reports failure
/// if pre-minimization bb_cst score is 0. Note that setting up constraints should be done outside
bool
PlaceSimultaneouslyMover::minimize_no_bb( core::pose::Pose & pose ) const {

	core::scoring::ScoreFunctionCOP stub_scorefxn( make_stub_scorefxn() );
	simple_filters::ScoreTypeFilter const stf( stub_scorefxn, backbone_stub_constraint, 1.0 );
	core::Real const before_min( stf.compute( pose ) );
	if ( before_min >= -0.0001 ) {
		TR<<"bb_cst evalutes to 0. Failing";
		return false;
	}
	//for minimization (rb and sc of previous placed stubs)
	utility::vector1< bool > const no_min( pose.size(), false );
	utility::vector1< core::Size > no_targets;
	MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, no_min/*sc_min*/, min_rb()/*rb*/, optimize_foldtree(), no_targets, true/*simultaneous optimization*/);
	return( true );
}

/// @details wraps around user-defined minimization movers. Note that setting up constraints should be done
/// outside this method.
void
PlaceSimultaneouslyMover::minimize_all( core::pose::Pose & pose, core::Size const minimization_steps ) const{
	using namespace protocols::hotspot_hashing;
	using namespace core::scoring::constraints;

	remove_hotspot_constraints_from_pose( pose );
	ConstraintCOPs const csts_before_min = pose.constraint_set()->get_all_constraints();
	core::Size const num_csts_before_min( csts_before_min.size() );

	//core::Size fixed_res(1);
	//if( host_chain_ == 1 ) fixed_res = pose.size();  // set but never used ~Labonte
	// core::id::AtomID const fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res ); // Unused variable causes a warning.

	core::scoring::ScoreFunctionCOP stub_scorefxn( make_stub_scorefxn() );

	for ( core::Size repeat( 1 ); repeat<=minimization_steps; ++repeat ) {
		for ( MoverRealPair const & curr : minimization_movers_ ) {
			using namespace core::scoring;

			simple_moves::DesignRepackMoverOP const curr_mover( curr.first );
			core::Real const bb_cst_weight( curr.second );
			TR<<"applying mover: "<<curr_mover->get_name()<<std::endl;
			//restricting movers for stub minimization
			curr_mover->prevent_repacking( prevent_repacking() );
			curr_mover->optimize_foldtree( false );
			curr_mover->design( false ); //we dont want any design to take place within any mover for stub minimization
			TR<<" design and repacking during stub minimization is prevented" << std::endl;
			TR<<"using weight: "<<bb_cst_weight<<" for the stub bb constraints" << std::endl;
			TR<<"and 1.0 for coordinate constraints"<<std::endl;
			ScoreFunctionOP minimize_mover_scorefxn_repack( curr_mover->scorefxn_repack() );
			if ( minimize_mover_scorefxn_repack ) {
				minimize_mover_scorefxn_repack->set_weight( backbone_stub_constraint, bb_cst_weight );
			}
			ScoreFunctionOP minimize_mover_scorefxn_minimize( curr_mover->scorefxn_minimize() );
			if ( minimize_mover_scorefxn_minimize ) {
				minimize_mover_scorefxn_minimize->set_weight( backbone_stub_constraint, bb_cst_weight );
			}
			curr_mover->apply( pose );
			utility::vector1< bool > sc_min( pose.size(), false );
			utility::vector1< bool > const no_min( pose.size(), false );
			utility::vector1< core::Size > targets;
			for ( StubSetStubPos const & stubset_pos_pair : stub_sets_ ) {
				core::Size const pos( stubset_pos_pair.second.second );
				if ( pos!=0 ) {
					targets.push_back( pos );
					sc_min[ pos ] = true;
				}
			}//foreach stubset_pos_pair
			using namespace core::scoring;
			ScoreFunctionCOP stub_scorefxn( make_stub_scorefxn() );
			ScoreFunctionOP scorefxn_mod( get_score_function());
			scorefxn_mod->set_weight( backbone_stub_constraint, 10.0 ); //This will not have any effect if the bbcsts are off
			scorefxn_mod->set_weight( coordinate_constraint, 1.0 );//similarly
			MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, optimize_foldtree(), targets, true/*simultaneous optimization*/ );
			utility::vector1< bool > min_host( pose.size(), false );
			core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ) );
			core::Size const host_chain_end  ( pose.conformation().chain_end( host_chain_ ) );
			for ( core::Size i( host_chain_begin ); i<=host_chain_end; ++i ) min_host[ i ] = true;
			utility::vector1< bool > const no_min_rb( pose.num_jump(), false );
			MinimizeInterface( pose, scorefxn_mod, min_host/*bb*/, sc_min, no_min_rb/*rb*/, optimize_foldtree(), targets, true/*simultaneous optimization*/ );
			TR<<"Doing rb minimization towards the stub and sc minimization of placed stubs" <<std::endl;
			MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, optimize_foldtree(), targets, true/*simultaneous optimization*/ );
		}//foreach minimization mover
	}//repeat

	ConstraintCOPs const csts_after_min = pose.constraint_set()->get_all_constraints();
	core::Size const num_csts_after_min( csts_after_min.size() );
	TR<<"Csts before min "<<num_csts_before_min<<" after min "<<num_csts_after_min<<std::endl;
	runtime_assert( num_csts_after_min == num_csts_before_min );//constraints shouldnt have changed
}

/// @details setup the residue level tasks for each of the paired positions in stub_sets_. This sets the explosion
/// level as well as the allowed identities at each position
/// Note that two things happen in this function: 1) the task is modified to reflect the rotamer explosion operations.
/// 2) residue_level_tasks_for_placed_hotspots_ is updated with these operations
core::pack::task::PackerTaskOP
PlaceSimultaneouslyMover::create_task_for_hotspot_packing( core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;
	using namespace core::pack::task;

	residue_level_tasks_for_placed_hotspots_->clear();
	if ( task_factory() ) {
		*residue_level_tasks_for_placed_hotspots_ = *(task_factory()); // this will allow PlaceStub's TaskAware paragraph to affect what is happening here, including trickling down to design movers
	}
	for ( StubSetStubPos const & stubset_pos : stub_sets_ ) {
		core::Size const pos( stubset_pos.second.second );
		runtime_assert( pos );
		HotspotStubSetCOP hs_set = stubset_pos.first;
		using namespace core::pack::task::operation;
		RotamerExplosionOP re_op( new RotamerExplosion( pos, EX_THREE_THIRD_STEP_STDDEVS, explosion_ ) );
		utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
		for ( ResidueAuctionItem const & item : auction_->auction_results() ) {
			HotspotStubSetCOP hs_set_curr( item.second.second.first );
			if ( hs_set_curr != hs_set ) continue;
			HotspotStubCOP hs_stub_curr( item.second.second.second );
			chemical::ResidueType const type( hs_stub_curr->residue()->type() );
			allowed_aas[ hs_stub_curr->residue()->type().aa() ] = true;
		}//foreach item in auction_->auction_results()
		RestrictAbsentCanonicalAASOP rac_op( new RestrictAbsentCanonicalAAS( pos, allowed_aas ) );

		residue_level_tasks_for_placed_hotspots_->push_back( rac_op );
		residue_level_tasks_for_placed_hotspots_->push_back( re_op );
	}//foreach stubset_pos
	PackerTaskOP task = residue_level_tasks_for_placed_hotspots_->create_task_and_apply_taskoperations( pose );
	return( task );
}


/// @details setup the residue level tasks for each of the paired positions in stub_sets_. Add all rotamers from the auction results
core::pack::task::PackerTaskOP
PlaceSimultaneouslyMover::create_task_for_allhotspot_packing( core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;
	using namespace core::pack::task;

	residue_level_tasks_for_placed_hotspots_->clear();
	if ( task_factory() ) {
		*residue_level_tasks_for_placed_hotspots_ = *(task_factory()); // this will allow PlaceStub's TaskAware paragraph to affect what is happening here, including trickling down to design movers
	}
	for ( StubSetStubPos const & stubset_pos : stub_sets_ ) {
		core::Size const pos( stubset_pos.second.second );
		runtime_assert( pos );
		HotspotStubSetCOP hs_set = stubset_pos.first;
		using namespace core::pack::task::operation;
		RotamerExplosionOP re_op( new RotamerExplosion( pos, EX_THREE_THIRD_STEP_STDDEVS, explosion_ ) );
		utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
		for ( ResidueAuctionItem const & item : auction_->auction_results() ) {
			HotspotStubSetCOP hs_set_curr( item.second.second.first );
			if ( pos==item.second.first ) {
				HotspotStubCOP hs_stub_curr( item.second.second.second );
				chemical::ResidueType const type( hs_stub_curr->residue()->type() );
				allowed_aas[ hs_stub_curr->residue()->type().aa() ] = true;
				TR << "RestrictAbsentCanonicalAAS: " << pos << " " << hs_stub_curr->residue()->type().aa() << " " << item.first <<  std::endl;
			}
		}//foreach item in auction_->auction_results()
		RestrictAbsentCanonicalAASOP rac_op( new RestrictAbsentCanonicalAAS( pos, allowed_aas ) );

		residue_level_tasks_for_placed_hotspots_->push_back( rac_op );
		residue_level_tasks_for_placed_hotspots_->push_back( re_op );
	}//foreach stubset_pos
	PackerTaskOP task = residue_level_tasks_for_placed_hotspots_->create_task_and_apply_taskoperations( pose );
	return( task );
}

/// @add coordinate constraints to pose
void PlaceSimultaneouslyMover::add_coordinatecst_for_hotspot_packing( core::pose::Pose & pose ) {
	using namespace protocols::hotspot_hashing;
	for ( StubSetStubPos const & stubset_pos : stub_sets_ ) {
		core::Size const pos( stubset_pos.second.second );
		runtime_assert( pos );
		HotspotStubSetCOP hs_set = stubset_pos.first;

		for ( ResidueAuctionItem const & item : auction_->auction_results() ) {
			HotspotStubSetCOP hs_set_curr( item.second.second.first );
			HotspotStubCOP hs_stub_curr( item.second.second.second );
			chemical::ResidueType const type( hs_stub_curr->residue()->type() );
			if ( pos==item.second.first ) {
				//loop through sidechain heavy atom of pos and the coordinates of hs_stub_curr->residue()->xyz(heavy_sidechain)
				//create coordinate constraints
				core::scoring::func::HarmonicFuncOP dummy_cst;
				add_coordinate_constraints( pose, *hs_stub_curr->residue(), host_chain_, pos, 0.5, dummy_cst );
				TR<<"applied coordinate constraints at position " << pos << std::endl;
			}
		}//foreach stubset_pos
	}
}

/// @details positions on the scaffold are auctioned to hotspot stub sets. The stubset that has a stub with the
/// lowest constraint energy with respect to that position wins the auction. If the number of positions that are
/// paired is less than the number of hotspot families, then failure is reported.
/// If pairing succeeds, each of the paired positions is allowed to adopt an identity from its parent stubset
/// that has below 0 constraint energy. Following which pack_rotamers is called with ex1, ex2, and rotamer explosion
/// if directed by the user.
/// Finally, each placed hotspot is tested for its total energy and for a user-defined filter. If one of these
/// filters fails, failure is reported.
/// Note that residues other than the lowest constraint energy identities might emerge from this process.
bool
PlaceSimultaneouslyMover::pair_sets_with_positions( core::pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;

	TR<<"Calling auction"<<std::endl;
	auction_->stub_sets( stub_sets_ );
	auction_->apply( pose );
	protocols::moves::MoverStatus const auction_stat( auction_->get_last_move_status() );
	if ( auction_stat != protocols::moves::MS_SUCCESS ) {
		TR<<"Auction failed => PlaceSimultaneously failing." << std::endl;
		return( false );
	}
	stub_sets_ = auction_->stub_sets();//copying the pairings information
	utility::vector1< bool > sc_min( pose.size(), false );
	utility::vector1< bool > const no_min( pose.size(), false );
	utility::vector1< core::Size > targets;
	for ( StubSetStubPos const & stubset_pos_pair : auction_->stub_sets() ) {
		core::Size const pos( stubset_pos_pair.second.second );
		TR << "pos: " << pos << std::endl;
		targets.push_back( pos );
		sc_min[ pos ] = true;
	}//foreach stubset_pos_pair

	//add to the prevent repacking
	utility::vector1< core::Size > prev_pack( prevent_repacking() );
	prev_pack.insert( prev_pack.begin(), targets.begin(), targets.end() );
	std::sort( prev_pack.begin(), prev_pack.end() );
	utility::vector1< core::Size >::iterator last = std::unique( prev_pack.begin(), prev_pack.end() );
	prev_pack.erase( last, prev_pack.end() );
	prevent_repacking( prev_pack );
	//end add to prevent repacking
	ScoreFunctionCOP scorefxn = get_score_function();

	if ( auction_->get_stub_scorefxn() == "backbone_stub_constraint" ) {
		using namespace core::pack;
		using namespace core::pack::task;
		PackerTaskOP task = create_task_for_hotspot_packing( pose );
		//  task->/*initialize_from_command_line().*/or_include_current( true ); // we don't want rotamer explosion with ex1 ex2
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( !pose.residue(i).is_protein() ) continue;
			if ( std::find( targets.begin(), targets.end(), i ) == targets.end() ) {
				task->nonconst_residue_task(i).prevent_repacking();
			}
		}//for residue i in pose
		using namespace core::scoring;
		pack_rotamers( pose, *scorefxn, task );
		using namespace core::scoring;
		ScoreFunctionCOP stub_scorefxn( make_stub_scorefxn() );
		MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, optimize_foldtree()/*optimize foldtree*/, targets, true/*simultaneous optimization*/ );
	} else if ( auction_->get_stub_scorefxn() == "backbone_stub_linear_constraint" ) {
		using namespace core::pack;
		using namespace core::pack::task;
		using core::pack::task::operation::TaskOperationCOP;
		protocols::scoring::Interface interface_obj(pose.num_jump());
		pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
		interface_obj.distance( 8 );
		interface_obj.calculate( pose );

		std::multimap< core::Real, std::pair< core::Size, StubsetStubPair > > saved_auction = auction_->auction_results();
		std::multimap< core::Real, std::pair< core::Size, StubsetStubPair > > new_auction;

		utility::vector1< Size > scanned_position;
		core::pose::Pose saved_pose=pose;

		for ( PlacementAuctionMover::ResidueAuction::iterator each_auction_result = saved_auction.begin(); each_auction_result != saved_auction.end(); ++each_auction_result ) {
			//old copy of pose and cleared auction
			auction_->clear();
			saved_pose=pose;
			auction_->insert(*each_auction_result);
			PackerTaskOP task = create_task_for_allhotspot_packing( saved_pose );
			//TR << "test loop through saved_auction: " << auction_->auction_results().size() << " , " << new_auction.size() << std::endl;

			//Size residue= each_auction_result->second.first;  // unused ~Labonte
			residue_level_tasks_for_placed_hotspots_->clear();
			core::pack::task::TaskFactoryOP pack_around_placed_hotspots_ = residue_level_tasks_for_placed_hotspots_->clone();
			pack_around_placed_hotspots_->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
			core::pack::task::operation::PreventRepackingOP prop( new core::pack::task::operation::PreventRepacking );

			for ( core::Size i=1; i<=saved_pose.size(); ++i ) {
				if ( !saved_pose.residue(i).is_protein() ) continue;
				if ( std::find( targets.begin(), targets.end(), i ) == targets.end() ) {
					//prevents everything else from repacking
					task->nonconst_residue_task(i).prevent_repacking();

					//if i in interface and in chain 1 and inpair with prev_pack
					if ( interface_obj.is_interface( i ) ) {
						bool contact_any=false;
						for ( core::Size j=1; j<=prev_pack.size(); ++j ) {
							if ( interface_obj.is_pair( saved_pose.residue(prev_pack[j]), saved_pose.residue(i) ) )  {
								contact_any=true;
								continue;
							}
						}
						if ( contact_any==true ) {
							;
						} else {
							prop->include_residue(i);
						}
					} else {
						prop->include_residue(i);
					} //non-interface residues

				}
			}//for residue i in saved_pose

			using namespace core::scoring;
			add_coordinatecst_for_hotspot_packing(saved_pose);
			pack_rotamers( saved_pose, *scorefxn, task );

			ScoreFunctionOP scorefxnc = get_score_function();
			scorefxnc->set_weight( coordinate_constraint, 1.0 );
			(*scorefxnc)(saved_pose);
			Real cst_score = saved_pose.energies().total_energies()[core::scoring::coordinate_constraint ];
			//TR << "residue: " << residue << " coordinate constraint energy: " << cst_score << std::endl;
			TR << std::endl;
			remove_coordinate_constraints_from_pose( saved_pose );

			//decide whether to keep it or not
			//complexity for the unknown number of hotspot
			if ( new_auction.empty() ) { //size()==0) {
				new_auction.insert( std::make_pair( cst_score, std::make_pair( each_auction_result->second.first, std::make_pair( each_auction_result->second.second.first, each_auction_result->second.second.second) ) ) );
			} else {
				Size status=0;
				// Modifying a container you're iterating over is perilous, as iterators may become invalid.
				// Luckily for a multimap only the iterator being deleted is invalidated. We can make a copy and increment the main iterator before we delete the copy.
				for ( PlacementAuctionMover::ResidueAuction::iterator selected_auction_iterator = new_auction.begin(); selected_auction_iterator != new_auction.end(); /*intentionally left blank*/ ) {
					PlacementAuctionMover::ResidueAuction::iterator selected_auction_result( selected_auction_iterator );
					++selected_auction_iterator; // Increment before delete, because it will be invalid after the deletion of the copy
					if ( each_auction_result->second.first == selected_auction_result->second.first ) {
						if ( cst_score < selected_auction_result->first ) {
							new_auction.erase( selected_auction_result );
							// selected_auction_result is now an invalid iterator
							// NOTE: Inserting into a multimap we're iterating over doesn't invalidate iterators, but it may mean we will encounter the inserted value during the current round of iteration.
							//       This may not be the behavior we're after here.
							new_auction.insert( std::make_pair( cst_score, std::make_pair( each_auction_result->second.first, std::make_pair( each_auction_result->second.second.first, each_auction_result->second.second.second) ) ) );
						}
						status=1;
					}
				}

				if ( status==0 ) {
					new_auction.insert( std::make_pair( cst_score, std::make_pair( each_auction_result->second.first, std::make_pair( each_auction_result->second.second.first, each_auction_result->second.second.second) ) ) );
				}

			} //insert when there is something

			remove_coordinate_constraints_from_pose( saved_pose );
		} //pack through each position

		//loop through each position, save best score for each position and insert to new_auction_results
		auction_->clear();
		for ( PlacementAuctionMover::ResidueAuction::iterator selected_auction_result = new_auction.begin(); selected_auction_result != new_auction.end(); ++selected_auction_result ) {
			auction_->insert(*selected_auction_result);
			TR << "selected coordinate constraint score: " << selected_auction_result->first << " residue: " << selected_auction_result->second.first << std::endl;
			if ( selected_auction_result->first >= coor_cst_cutoff_ ) {
				TR<<"coordinate constraint energy: " << selected_auction_result->first << " Failed cutoff "<<coor_cst_cutoff_<<", Can be bad rotamer" << std::endl;
				return( false );
			}
		}

		PackerTaskOP task = create_task_for_allhotspot_packing( pose );
		//PackerTaskOP taskc = task->clone();

		//create a new task that does not allow design
		core::pack::task::TaskFactoryOP pack_around_placed_hotspots_ = residue_level_tasks_for_placed_hotspots_->clone();
		//core::pack::task::TaskFactoryOP pack_around_placed_hotspots_ = new core::pack::task::TaskFactory;
		//if( task_factory() )
		//  *pack_around_placed_hotspots_ = *(task_factory());
		pack_around_placed_hotspots_->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
		core::pack::task::operation::PreventRepackingOP prop( new core::pack::task::operation::PreventRepacking );

		//  task->/*initialize_from_command_line().*/or_include_current( true ); // we don't want rotamer explosion with ex1 ex2

		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( !pose.residue(i).is_protein() ) continue;
			if ( std::find( targets.begin(), targets.end(), i ) == targets.end() ) {
				//prevents everything else from repacking
				task->nonconst_residue_task(i).prevent_repacking();

				//if i in interface and in chain 1 and inpair with prev_pack
				if ( interface_obj.is_interface( i ) ) {
					bool contact_any=false;
					for ( core::Size j=1; j<=prev_pack.size(); ++j ) {
						if ( interface_obj.is_pair( pose.residue(prev_pack[j]), pose.residue(i) ) )  {
							contact_any=true;
							continue;
						}
					}
					if ( contact_any==true ) {
						TR << "only repack " << i << std::endl;
					} else {
						prop->include_residue(i);
					}
					//taskc->nonconst_residue_task(i).prevent_repacking();
				} else {
					prop->include_residue(i);
					//taskc->nonconst_residue_task(i).prevent_repacking();
				} //non-interface residues

			}
		}//for residue i in pose

		using namespace core::scoring;
		TR << "Add hotspot" << std::endl;
		pack_rotamers( pose, *scorefxn, task );

	} // end of backbone_stub_linear_constraint

	//filter each placed hotspot
	for ( StubSetStubPos & stubset_pos_pair : stub_sets_ ) {
		using namespace core::scoring;
		core::Size const pos( stubset_pos_pair.second.second );
		HotspotStubSetOP stubset( stubset_pos_pair.first );
		HotspotStubCOP stub( stubset_pos_pair.second.first );
		using namespace protocols::filters;
		simple_filters::EnergyPerResidueFilter total_energy_filter( pos, scorefxn, total_score, stub_energy_threshold_ );
		bool const pass_tot_energy( total_energy_filter.apply( pose ) );
		protocols::filters::FilterOP modified_filter( stub_set_filters_[ stubset ]->clone() );
		protocols::moves::modify_ResId_based_object( modified_filter, pos );
		bool const pass_stub_set_filter( modified_filter->apply( pose ) );
		core::Real const distance( pose.residue( pos ).xyz( "CB" ).distance( stub->residue()->xyz( "CB" ) ) );
		if ( distance >= max_cb_cb_dist_ ) {
			TR<<"distance: " << distance << " Failed distance cutoff " << max_cb_cb_dist_ << std::endl;
			return( false );
		} else if ( !pass_tot_energy || !pass_stub_set_filter ) {
			TR<<"Failed stub filters " << std::endl;
			return( false );
		}
	}

	bool const after_placement_pass( after_placement_filter_->apply( pose ) );
	if ( !after_placement_pass ) {
		TR<<"Failed after_placement_filter" << std::endl;
		return( false );
	}

	for ( StubSetStubPos const & hs_set : stub_sets_ ) {
		core::Size const position( hs_set.second.second );
		TR<<"Paired position "<<position<<std::endl;
	}

	TR.flush();
	return( true );
}

/// @details this will be removed in the future
bool
PlaceSimultaneouslyMover::place_stubs( core::pose::Pose & pose ) const
{
	using namespace protocols::hotspot_hashing;
	core::Size const chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const chain_end( pose.conformation().chain_end( host_chain_ ) );

	for ( StubSetStubPos const & hs_set : stub_sets_ ) {
		HotspotStubCOP stub( hs_set.second.first );
		core::Size const res_num( hs_set.second.second );
		pose.replace_residue( res_num, *(stub->residue()), true );

		using namespace core::chemical;
		//removing variant types can only be done if the residue would still be connected to the chain
		if ( res_num < chain_end ) {
			core::pose::remove_upper_terminus_type_from_pose_residue( pose, res_num );
		}
		if ( res_num > chain_begin ) {
			core::pose::remove_lower_terminus_type_from_pose_residue( pose, res_num );
		}
		pose.conformation().update_polymeric_connection( res_num, true ); // o/w residues connections mess up
		if ( res_num > chain_begin ) {
			pose.conformation().update_polymeric_connection( res_num - 1, true );
		}
	}

	pose.update_residue_neighbors();
	bool const pass( after_placement_filter_->apply( pose ) );
	return( pass );
}

/// @details removes the coordinate constraints from the pose and reapplies them individually.
void
PlaceSimultaneouslyMover::refresh_coordinate_constraints( core::pose::Pose & pose, core::Real const coord_sdev )
{
	remove_coordinate_constraints_from_pose( pose );
	saved_coord_constraints_ = pose.add_constraints( saved_coord_constraints_ );
	if ( coord_sdev < 0.0001 ) {
		TR<<"no coordinate constraints applied" << std::endl;
		return;
	}
	for ( StubSetStubPos const & hs_set : stub_sets_ ) {
		using namespace protocols::hotspot_hashing;

		HotspotStubCOP stub( hs_set.second.first );
		core::Size const position( hs_set.second.second );
		core::scoring::func::HarmonicFuncOP dummy_cst;
		add_coordinate_constraints( pose, *stub->residue(), host_chain_, position, coord_sdev, dummy_cst );
	}
	TR<<"applied coordinate constraints" << std::endl;
}

/// @details wraps around the user-defined design movers. Applies coordinate constraints as specified by the user.
/// Prevents repacking of placed hotspots
void
PlaceSimultaneouslyMover::design( core::pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;
	TR<<"redesigning remainder of interface with user defined design movers" << std::endl;

	BuildAlaPose toAla( host_chain_ == 1/*partner1*/, host_chain_ == 2 /*partner2*/ );
	utility::vector1< core::Size > no_repack;
	if ( !prevent_repacking().empty() ) no_repack = prevent_repacking();
	if ( !no_repack.empty() ) {
		std::sort( no_repack.begin(), no_repack.end() );
		utility::vector1< core::Size >::iterator last = std::unique( no_repack.begin(), no_repack.end() );
		no_repack.erase( last, no_repack.end() );
		toAla.prevent_repacking( no_repack );
	}
	toAla.task_factory( residue_level_tasks_for_placed_hotspots_ );
	TR<<"switching interface to alanine" << std::endl;
	toAla.apply( pose );

	saved_coord_constraints_ = remove_coordinate_constraints_from_pose( pose );
	for ( MoverRealPair const & mover_coord_cst : design_movers_ ) {//design movers
		core::Real const sdev( mover_coord_cst.second );
		simple_moves::DesignRepackMoverOP mover( mover_coord_cst.first );
		TR<<"applying design mover "<<mover->get_name()<<std::endl;
		if ( sdev >= 0 ) { //use constraints
			core::Size const before_refresh( pose.constraint_set()->get_all_constraints().size() );
			refresh_coordinate_constraints( pose, sdev );
			core::Size const after_refresh( pose.constraint_set()->get_all_constraints().size() );
			TR.Debug<<"before refreshing coord cst "<<before_refresh<<" constraints. after: "<<after_refresh<<std::endl;
		} else { //use constraints
			TR<<"no constraints applied"<< std::endl;
		}

		core::Size const before_apply_design_mover( pose.constraint_set()->get_all_constraints().size() );
		mover->prevent_repacking( prevent_repacking() );
		using namespace core::pack::task;
		mover->task_factory( residue_level_tasks_for_placed_hotspots_ );
		mover->optimize_foldtree( false );
		TR<<"setting coordinate constraint weights to 1.0 in design movers" << std::endl;
		core::scoring::ScoreFunctionOP scorefxn_rep( mover->scorefxn_repack() );
		core::scoring::ScoreFunctionOP scorefxn_min( mover->scorefxn_minimize() );
		if ( scorefxn_rep ) scorefxn_rep->set_weight( coordinate_constraint, 1.0 );
		if ( scorefxn_min ) scorefxn_min->set_weight( coordinate_constraint, 1.0 );
		mover->apply( pose );
		core::Size const after_apply_design_mover( pose.constraint_set()->get_all_constraints().size() );
		TR.Debug<<"before applying design mover "<<before_apply_design_mover<<" constraints. After: "<<after_apply_design_mover<<std::endl;
		if ( before_apply_design_mover != after_apply_design_mover ) {
			TR<<"ERROR: This design mover changed the number of constraints on the pose. Before: "<< before_apply_design_mover<<" after: "<<after_apply_design_mover<<". This behaviour is unsupported."<<std::endl;
			runtime_assert( before_apply_design_mover != after_apply_design_mover );
		}
		TR<<"removing coordinate constraints" << std::endl;
		remove_coordinate_constraints_from_pose( pose );
	}//foreach Design movers
	remove_coordinate_constraints_from_pose( pose ); // just in case
	pose.add_constraints( saved_coord_constraints_ );
}

void
PlaceSimultaneouslyMover::apply( core::pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;

	bool const stub_score_filter_pass( stub_score_filter_->apply( pose ) );
	if ( !stub_score_filter_pass ) {
		TR<<"Stub score filter reported failure at the beginning of PlaceSimultaneously. Failing"<<std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	saved_pose_ = pose;
	saved_stub_sets_ = stub_sets_;
	// rescore the pose to ensure that backbone stub cst's are populated
	core::scoring::ScoreFunctionOP bbcst_scorefxn( new core::scoring::ScoreFunction );
	bbcst_scorefxn->reset();
	if ( auction_->get_stub_scorefxn() == "backbone_stub_linear_constraint" ) {
		bbcst_scorefxn->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
	}// else if ( auction_->get_stub_scorefxn() == "backbone_stub_linear_constraint" ) {
	//  bbcst_scorefxn->set_weight( core::scoring::backbone_stub_linear_constraint, 1.0 );
	//}
	(*bbcst_scorefxn)( pose );

	rbstub_minimization_->apply( pose );
	minimize_all( pose, minimization_repeats_before_placement_ );
	bool const pairing( pair_sets_with_positions( pose ) );
	if ( !pairing ) {
		TR<<"Failed to pair sets with individual scaffold positions"<<std::endl;
		final_cleanup( pose );
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}
	TR.flush();
	protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );
	minimize_all( pose, minimization_repeats_after_placement_ );
	TR.flush();
	design( pose );
	TR<<"Design done" << std::endl;
	TR.flush();
	set_last_move_status( protocols::moves::MS_SUCCESS );
}

// XRW TEMP std::string
// XRW TEMP PlaceSimultaneouslyMover::get_name() const {
// XRW TEMP  return PlaceSimultaneouslyMover::mover_name();
// XRW TEMP }

/// @details This should be called before failing placesimultaneously.
void
PlaceSimultaneouslyMover::final_cleanup( core::pose::Pose & pose )
{
	pose = saved_pose_;
	//don't accumulate state between ntrials
	prevent_repacking_.clear();
	restrict_to_repacking_.clear();
	saved_bb_constraints_.clear();
	stub_sets_ = saved_stub_sets_;
	auction_->clear();
	/// The statement below would clear task_factories of child design movers
	// if( residue_level_tasks_for_placed_hotspots_ )
	//  residue_level_tasks_for_placed_hotspots_->clear();
	TR.flush();
}

void
PlaceSimultaneouslyMover::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &filters,
	Movers_map const &movers,
	core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;
	using namespace protocols::filters;
	using namespace core::pack::task;

	TR<<"Parsing PlaceSimultaneouslyMover----"<<std::endl;
	/// auction, rbstub_minimization, and stub_score_filter are different, private
	/// instantiations and should parse the stubset, chain, etc. information
	auction_->parse_my_tag( tag, data, filters, movers, pose );
	rbstub_minimization_->parse_my_tag( tag, data, filters, movers, pose );
	stub_score_filter_->parse_my_tag( tag, data, filters, movers, pose );

	if ( tag->hasOption( "task_operations" ) ) {
		if ( task_factory() ) {
			TR<<"*****WARNING: ERASING existing task_factory, b/c of specifications for new task operations in " << std::endl<<tag<<std::endl;
		}
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}

	host_chain_ = tag->getOption<core::Size>( "chain_to_design", 2 );
	coor_cst_cutoff_ = tag->getOption< core::Real >( "coor_cst_cutoff", 100.0 );
	optimize_foldtree( tag->getOption< bool >( "optimize_fold_tree", true ) );
	TR<<"optimize_foldtree set to: "<<optimize_foldtree()<<std::endl;

	repack_non_ala_ = tag->getOption<bool>( "repack_non_ala", true );
	min_rb( true );

	std::string const after_placement_filter_name( tag->getOption<std::string>( "after_placement_filter", "true_filter" ) );
	Filters_map::const_iterator ap_filter( filters.find( after_placement_filter_name ) );
	if ( after_placement_filter_name == "true_filter" ) {
		after_placement_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
	} else {
		after_placement_filter_ = ap_filter->second->clone();
	}
	//parsing stub minimize movers and design movers for place stub
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );

	for ( TagCOP const btag : branch_tags ) {
		if ( btag->getName() == "StubMinimize" ) {
			minimization_repeats_before_placement_ = btag->getOption< core::Size >( "min_repeats_before_placement", 0 );
			minimization_repeats_after_placement_ = btag->getOption< core::Size >( "min_repeats_after_placement", 1 );
			utility::vector0< TagCOP > const & stub_min_tags( btag->getTags() );
			for ( TagCOP stub_m_tag : stub_min_tags ) {
				std::string const stub_mover_name( stub_m_tag->getOption<std::string>( "mover_name" ) );
				core::Real  const bb_stub_constraint_weight( stub_m_tag->getOption< core::Real > ( "bb_cst_weight", 10.0 ) );
				std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( stub_mover_name ));
				bool const stub_mover_found( find_mover != movers.end() );
				if ( stub_mover_found ) {
					simple_moves::DesignRepackMoverOP drSOP = utility::pointer::dynamic_pointer_cast< simple_moves::DesignRepackMover > ( find_mover->second->clone() );
					if ( !drSOP ) {
						TR<<"dynamic cast failed in tag "<<tag<<". Make sure that the mover is derived from DesignRepackMover"<<std::endl;
						runtime_assert( drSOP != 0 );
					}//done cast check
					minimization_movers_.push_back( std::make_pair( drSOP, bb_stub_constraint_weight) );
					TR<<"added stub minimize mover "<<stub_mover_name<<" to minimize towards the stub. Using this weight for the bb stub constraints: "<< bb_stub_constraint_weight<<std::endl;
				}
			}
		} else if ( btag->getName() == "DesignMovers" ) {
			utility::vector0< TagCOP > const & design_tags( btag->getTags() );
			for ( TagCOP const m_tag_ptr : design_tags ) {
				std::string const mover_name( m_tag_ptr->getOption< std::string >( "mover_name" ) );
				bool const apply_coord_constraints( m_tag_ptr->getOption< bool >( "use_constraints", true ) );
				core::Real const coord_cst_std( m_tag_ptr->getOption< core::Real >( "coord_cst_std", 0.5 ) );

				std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( mover_name ));
				bool const mover_found( find_mover != movers.end() );
				if ( mover_found ) {
					simple_moves::DesignRepackMoverOP drOP = utility::pointer::dynamic_pointer_cast< simple_moves::DesignRepackMover > ( find_mover->second );
					if ( !drOP ) {
						TR<<"dynamic cast failed in tag "<<tag<<". Make sure that the mover is derived from DesignRepackMover"<<std::endl;
						runtime_assert( drOP != 0 );
					}
					design_movers_.push_back( std::make_pair( drOP, ( apply_coord_constraints ? coord_cst_std : -1 ) ) );
					TR<<"added design mover "<<mover_name<<" to place simultaneously ";
					if ( apply_coord_constraints ) {
						TR<<"with with std "<< coord_cst_std<< std::endl;
					} else {
						TR<<"with no coord cst" << std::endl;
					}
				} else {
					TR<<"***WARNING WARNING! Mover defined for PlaceSimultaneouslyMoverMover not found in mover_list. EXITING ***"<<std::endl;
					runtime_assert( mover_found );
				}
			}
		} else if ( btag->getName() == "StubSets" ) {
			//PlaceSimultaneously doesn't use the utility function parse_stub_sets b/c it needs to connect filters with stubsets
			explosion_ = btag->getOption< core::Size >( "explosion", 0 );
			stub_energy_threshold_ = btag->getOption<core::Real>( "stub_energy_threshold", 1.0 );
			max_cb_cb_dist_ = btag->getOption< core::Real >( "max_cb_dist", 3.0 );

			utility::vector0< TagCOP > const & stubset_tags( btag->getTags() );
			for ( TagCOP const stubset_tag : stubset_tags ) {
				std::string const stub_fname = stubset_tag->getOption< std::string >( "stubfile" );
				HotspotStubSetOP stubset( new HotspotStubSet );
				if ( data.has( "hotspot_library", stub_fname ) ) {
					stubset = data.get_ptr<protocols::hotspot_hashing::HotspotStubSet>( "hotspot_library", stub_fname );
					TR<<"Associated mover with an already read stubset named "<<stub_fname<<std::endl;
				} else {
					stubset->read_data( stub_fname );
				}

				stub_sets_.push_back(std::make_pair(stubset, std::make_pair(HotspotStubOP( new HotspotStub() ), 0)));  // REQUIRED FOR WINDOWS
				//stub_sets_.push_back( StubSetStubPos( stubset, std::pair< HotspotStubOP, core::Size >( 0, 0 ) ) );
				core::pose::PoseOP ala_pose( new core::pose::Pose( pose ) );
				pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( *ala_pose ));
				task->initialize_from_command_line().or_include_current( true );

				utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
				allowed_aas[ chemical::aa_ala ] = true;

				core::Size const chain_begin( ala_pose->conformation().chain_begin( host_chain_ ) );
				core::Size const chain_end( ala_pose->conformation().chain_end( host_chain_ ) );

				for ( core::Size i = 1; i <= pose.size(); i++ ) {
					if ( !pose.residue(i).is_protein() ) continue;
					if ( i >= chain_begin && i <=chain_end ) {
						core::Size const restype( ala_pose->residue(i).aa() );
						if ( ( restype == chemical::aa_pro  && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) || restype == chemical::aa_gly ) {
							task->nonconst_residue_task(i).prevent_repacking();
						} else {
							task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
						}
					} else {
						task->nonconst_residue_task( i ).prevent_repacking();
					}
				}
				if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
					core::pack::task::parse_resfile(pose, *task);
				}

				core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
				pack::pack_rotamers( *ala_pose, *scorefxn, task);
				(*scorefxn)( *ala_pose );
				stubset->pair_with_scaffold( *ala_pose, host_chain_, protocols::filters::FilterCOP( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) ) );
				std::string const stub_set_filter_name( stubset_tag->getOption< std::string >( "filter_name", "true_filter" ) );
				Filters_map::const_iterator stub_set_filter( filters.find( stub_set_filter_name ) );
				runtime_assert( stub_set_filter != filters.end() );
				stub_set_filters_[ stubset ] = stub_set_filter->second->clone();
			}//foreach stubset_tag
		} else if ( btag->getName() != "NotifyMovers" ) { // fi stubsets
			utility_exit_with_message( "ERROR: tag in PlaceSimultaneouslyMover not defined\n" );
		}
		generate_taskfactory_and_add_task_awareness( btag, movers, data, task_factory() );//residue_level_tasks_for_placed_hotspots_ );
	}
	if ( minimization_movers_.size() == 0 ) {
		TR<<"No StubMinimize movers defined by user, defaulting to minimize_rb and _sc of stubs only" << std::endl;
	}

	auction_->task_factory( data.get_ptr< TaskFactory >( "TaskFactory", "placement" ) );
	rbstub_minimization_->task_factory( data.get_ptr< TaskFactory >( "TaskFactory", "placement" ) );
	// stub_score_filter_->parse_my_tag( tag, data, filters, movers, pose );
	// stub_score_filter_->stub_sets( stub_sets_ );
	TR<<"Using "<<minimization_repeats_before_placement_<<" minimization steps before placement (bbcst constraints on)" << std::endl;
	TR<<"Using "<<minimization_repeats_after_placement_<<" minimization steps after placement (no constraints on)" << std::endl;
	TR<<"max cb cb distance set to "<<max_cb_cb_dist_<<std::endl;
	TR<<"place simultaneously mover on chain "<<host_chain_<<" with repack_non_ala set to "<<repack_non_ala_<<std::endl;
	TR<<"Using auction energy function: "<<auction_->get_stub_scorefxn() << " with coordinate cst cutoff: " << coor_cst_cutoff_ <<std::endl;
}

PlaceSimultaneouslyMover::PlaceSimultaneouslyMover() :
	simple_moves::DesignRepackMover( PlaceSimultaneouslyMover::mover_name() )
{
	residue_level_tasks_for_placed_hotspots_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );//watch out! Never allocate a new task factory for this guy after parsing has been done, b/c in parsing all task aware movers will be watching it through their task_factory_
	// user_defined_auction_ = false;
	// user_defined_stub_score_filter_ = false;
	// user_defined_bbstub_minimization_ = false;
	auction_ = PlacementAuctionMoverOP( new PlacementAuctionMover );
	rbstub_minimization_ = PlacementMinimizationMoverOP( new PlacementMinimizationMover );
	stub_score_filter_ = protocols::protein_interface_design::filters::StubScoreFilterOP( new protocols::protein_interface_design::filters::StubScoreFilter );
	stub_sets_.clear();
	stub_set_filters_.clear();
	host_chain_=2;
	saved_stub_sets_.clear();
	stub_set_filters_.clear();
	max_cb_cb_dist_ = 3.0;
	stub_energy_threshold_=0.0;
	minimization_movers_.clear();
	minimization_repeats_before_placement_=0;
	minimization_repeats_after_placement_ = 1;
	saved_bb_constraints_.clear();
	design_movers_.clear();
	explosion_ = 0;
	saved_coord_constraints_.clear();
}

std::string PlaceSimultaneouslyMover::get_name() const {
	return mover_name();
}

std::string PlaceSimultaneouslyMover::mover_name() {
	return "PlaceSimultaneously";
}

//name mangling XML schema extravaganza
std::string PlaceSimultaneously_subelement_mangler( std::string const & element ) { return "PlaceSimultaneously_subelement" + element + "_type"; }
std::string PlaceSimultaneously_StubMinimize_subelement_mangler( std::string const & element ) { return "PlaceSimultaneously_StubMinimize_subelement" + element + "_type"; }
std::string PlaceSimultaneously_DesignMovers_subelement_mangler( std::string const & element ) { return "PlaceSimultaneously_DesignMovers_subelement" + element + "_type"; }
std::string PlaceSimultaneously_NotifyMovers_subelement_mangler( std::string const & element ) { return "PlaceSimultaneously_NotifyMovers_subelement" + element + "_type"; }
std::string PlaceSimultaneously_StubSets_subelement_mangler( std::string const & element ) { return "PlaceSimultaneously_StubSets_subelement" + element + "_type"; }

void PlaceSimultaneouslyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//import anything from these (their parse_my_tags are called)
	//PlacementAuctionMover
	//PlacementMinimizationMover
	//StubScoreFilter

	/*  This is PlacementAuctionMover
	attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_non_negative_integer, "Chain where design ought to be performed, numbered sequentially from 1", "2" )
	+ XMLSchemaAttribute::attribute_w_default( "max_cb_dist", xsct_real, "Maximum distance from ideal placement that is nonetheless considered a hit", "3.0" )
	+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Force to apply to CB atoms", "0.5" )
	+ XMLSchemaAttribute::attribute_w_default( "stubscorefxn", xs_string, "Scoring function to apply to the stubs being placed", "backbone_stub_constraint" );

	XMLSchemaSimpleSubelementList ssl;
	add_subelement_for_parse_stub_sets( ssl, xsd );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "A PlacementAuctionMover pays increasing prices to put hotspots in place on a chain", attlist, ssl );
	*/

	attlist + XMLSchemaAttribute( "name", xs_string, "XSD XRW: TO DO");
	attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_positive_integer, "Chain where design ought to be performed, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "max_cb_dist", xsct_real, "Maximum distance from ideal placement that is nonetheless considered a hit", "3.0" )
		+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Force to apply to CB atoms", "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "stubscorefxn", xs_string, "Scoring function to apply to the stubs being placed", "backbone_stub_constraint" );
	//will have to parse StubSets subelement later; fortunately this mover parses a superset of what PlacementAuctionMover sees

	/* PlacementMinimizationMover's schema

	attlist
	+ XMLSchemaAttribute::attribute_w_default( "host_chain", xsct_positive_integer, "Probably the chain where the stub goes", "2")
	+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Force to apply to CB atoms.  Must be positive.", "0.5" ) // XRW TODO, should be positive Real
	+ XMLSchemaAttribute::attribute_w_default( "optimize_foldtree", xsct_rosetta_bool, "setup new fold_tree for better numerical behaviour between the residue at the center of target_residues and the nearest residue on the partner", "false")
	+ XMLSchemaAttribute::attribute_w_default( "minimize_rb", xsct_rosetta_bool, "minimize the rigid body degree of freedom", "true");

	XMLSchemaSimpleSubelementList ssl;
	add_subelement_for_parse_stub_sets( ssl, xsd ); //XRW TODO the subelement is actually REQUIRED; might be for all uses of this function
	*/

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "host_chain", xsct_positive_integer, "Probably the chain where the stub goes", "2")
		+ XMLSchemaAttribute::attribute_w_default( "optimize_fold_tree", xsct_rosetta_bool, "setup new fold_tree for better numerical behaviour between the residue at the center of target_residues and the nearest residue on the partner", "false")
		+ XMLSchemaAttribute::attribute_w_default( "minimize_rb", xsct_rosetta_bool, "do we want to minimize the rb dof during stub placement? This will allow a previously placed stub to move a a little to accommodate the new stub. It's a good idea to use this with the previously placed stub adding its implied constraints.", "true");
	//will have to parse StubSets subelement later; fortunately this mover parses a superset of what PlacementAuctionMover sees

	//StubScoreFilter only has attributes we've already included from PlacementAuctionMover and PlacementMinimizationMover
	/*attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_non_negative_integer, "Chain that ought to be designed, numbered sequentially from 1", "2" )
	+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_non_negative_integer, "Chain that ought to be designed, numbered sequentially from 1", "2" );

	XMLSchemaSimpleSubelementList ssl;
	protein_interface_design::movers::add_subelement_for_parse_stub_sets( ssl, xsd );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "This filter checks whether in the current configuration the scaffold is 'feeling' any of the hotspot stub constraints. This is useful for quick triaging of hopeless configuration. If used in conjunction with the Mover PlaceSimultaneously, don't bother setting flags -- the Mover will take care of it.", attlist, ssl );*/


	//now we handle attributes actually in this class's parse_my_tag:
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist
		//+ XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_positive_integer, "Chain where design ought to be performed, numbered sequentially from 1", "2")
		+ XMLSchemaAttribute::attribute_w_default( "coor_cst_cutoff", xsct_real, "the threshold coordinate constraint energy between the added hotspot residues and the one in the stub library. Use with stubscorefxn=backbone_stub_linear_constraint. PlaceSimultaneously fails if placed residues deviates beyond this threshold.", "100.0")
		//+ XMLSchemaAttribute::attribute_w_default( "optimize_fold_tree", xsct_rosetta_bool, "setup new fold_tree for better numerical behaviour between the residue at the center of target_residues and the nearest residue on the partner", "true")
		+ XMLSchemaAttribute::attribute_w_default( "repack_non_ala", xsct_rosetta_bool, "clearly implies something about repacking, although the only comment I could find mentioned reDESIGN instead", "true")
		+ XMLSchemaAttribute::attribute_w_default( "after_placement_filter", xs_string, "The name of a filter to be applied immediately after stub placement and StubMinimize movers, but before the DesignMovers run. Useful for quick sanity check on the goodness of the stub.", "true_filter");
	if ( !attribute_w_name_in_attribute_list( "cb_force", attlist ) ) {
		attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	//see main/source/test/utility/tag/XMLSchemaGeneration.cxxtest.hh for example that clarifies
	//This is also copied from PlaceStubMover, with modifications
	//This is the ROOT NODE
	XMLSchemaRepeatableCTNodeOP PlaceSimultaneouslyMover_node( new XMLSchemaRepeatableCTNode );
	PlaceSimultaneouslyMover_node->set_element_w_attributes( mover_name(), attlist, "Places hotspot residues simultaneously on a scaffold, rather than iteratively as in PlaceStub. It is faster therefore allowing more backbone sampling, and should be useful in placing more than 2 hotspots." );
	PlaceSimultaneouslyMover_node->set_kids_naming_func( & PlaceSimultaneously_subelement_mangler );
	PlaceSimultaneouslyMover_node->set_root_node_naming_func( & complex_type_name_for_mover );

	//Subelement StubMinimize  //DIFFERENT from PlaceStubMover
	//attribute min_repeats_before_placement, 0 (nnint)
	//attribute min_repeats_after_placement, 1 (nnint)
	//Subelement Add (fake name)
	//attribute mover_name string
	//attribute decimal bb_cst_weight 10.0

	XMLSchemaRepeatableCTNodeOP StubMinimize_node( new XMLSchemaRepeatableCTNode );
	AttributeList StubMinimize_attlist;
	StubMinimize_attlist
		+ XMLSchemaAttribute::attribute_w_default( "min_repeats_before_placement", xsct_non_negative_integer, "minimization trials before placement", "0")
		+ XMLSchemaAttribute::attribute_w_default( "min_repeats_after_placement", xsct_non_negative_integer, "minimization trials before placement", "1");
	if ( !attribute_w_name_in_attribute_list( "cb_force", StubMinimize_attlist ) ) {
		StubMinimize_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}

	StubMinimize_node->set_element_w_attributes( "StubMinimize", StubMinimize_attlist, "Defines Movers used to minimize w/r/t stub placement" );
	StubMinimize_node->set_kids_naming_func( & PlaceSimultaneously_StubMinimize_subelement_mangler );

	XMLSchemaRepeatableCTNodeOP StubMinimize_Add_node( new XMLSchemaRepeatableCTNode );
	StubMinimize_node->add_child(StubMinimize_Add_node);

	AttributeList StubMinimize_Add_attlist;
	StubMinimize_Add_attlist
		+ XMLSchemaAttribute::required_attribute( "mover_name", xs_string, "Name of the Mover (defined elsewhere in the XML)")
		+ XMLSchemaAttribute::attribute_w_default( "bb_cst_weight", xsct_real, "determines the strength of the constraints derived from the stubs. This value is a weight on the cb_force, so larger values are stronger constraints.", "10.0");
	if ( !attribute_w_name_in_attribute_list( "cb_force", StubMinimize_Add_attlist ) ) {
		StubMinimize_Add_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	StubMinimize_Add_node->set_element_w_attributes( "Add", StubMinimize_Add_attlist, "Defines one Mover used to minimize w/r/t stub placement");

	//Subelement DesignMovers
	//Subelement Add
	//string mover_name
	//bool use_constraints, true
	//Real coord_cst_std, 0.5

	XMLSchemaRepeatableCTNodeOP DesignMovers_node( new XMLSchemaRepeatableCTNode );
	AttributeList empty_attlist;
	DesignMovers_node->set_element_w_attributes( "DesignMovers", empty_attlist, "Defines Movers used to do design after stub placement" );
	DesignMovers_node->set_kids_naming_func( & PlaceSimultaneously_DesignMovers_subelement_mangler );

	XMLSchemaRepeatableCTNodeOP DesignMovers_Add_node( new XMLSchemaRepeatableCTNode );
	DesignMovers_node->add_child(DesignMovers_Add_node);

	AttributeList DesignMovers_Add_attlist;
	DesignMovers_Add_attlist
		+ XMLSchemaAttribute::required_attribute( "mover_name", xs_string, "Name of the Mover (defined elsewhere in the XML)")
		+ XMLSchemaAttribute::attribute_w_default( "use_constraints", xsct_rosetta_bool, "obey constraints", "true")
		+ XMLSchemaAttribute::attribute_w_default( "coord_cst_std", xsct_real, "standard deviation (width) on coordinate constraint", "0.5");
	if ( !attribute_w_name_in_attribute_list( "cb_force", DesignMovers_Add_attlist ) ) {
		DesignMovers_Add_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	DesignMovers_Add_node->set_element_w_attributes( "Add", DesignMovers_Add_attlist, "Defines a Mover used to do design after stub placement");

	//This part ought to be in a shared function with PlaceSimultaneouslyMover, but fragility is the cost of these Movers doing things they ought not

	//Subelement NotifyMovers (this comes from a helper function)
	//Subelement Add
	//string mover_name

	XMLSchemaRepeatableCTNodeOP NotifyMovers_node( new XMLSchemaRepeatableCTNode );
	//AttributeList empty_attlist;
	NotifyMovers_node->set_element_w_attributes( "NotifyMovers", empty_attlist, "Defines Movers (named elsewhere in XML) that must be aware of the protected stub locations so they are not lost" );
	NotifyMovers_node->set_kids_naming_func( & PlaceSimultaneously_NotifyMovers_subelement_mangler );

	XMLSchemaRepeatableCTNodeOP NotifyMovers_Add_node( new XMLSchemaRepeatableCTNode );
	NotifyMovers_node->add_child(NotifyMovers_Add_node);

	AttributeList NotifyMovers_Add_attlist;
	NotifyMovers_Add_attlist
		+ XMLSchemaAttribute::required_attribute( "mover_name", xs_string, "Name of the Mover (defined elsewhere in the XML)");
	if ( !attribute_w_name_in_attribute_list( "cb_force", NotifyMovers_Add_attlist ) ) {
		NotifyMovers_Add_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	NotifyMovers_Add_node->set_element_w_attributes( "Add", NotifyMovers_Add_attlist, "Defines a Mover (named elsewhere in XML) that must be aware of the protected stub locations so they are not lost");

	//Subelement StubSets //MAY be similar to parse_stub_sets - actually a superset
	//attribute explosion, nnint, 0
	//attribute stub_energy_threshold, real, 1.0
	//attribute max_cb_dist, real, 3.0
	//Subelement Add
	//attribute stubfile, string, required
	//attribute filter_name, string, default true_filter

	XMLSchemaRepeatableCTNodeOP StubSets_node( new XMLSchemaRepeatableCTNode );
	AttributeList StubSets_attlist;
	StubSets_attlist
		+ XMLSchemaAttribute::attribute_w_default( "explosion", xsct_non_negative_integer, "which chis to explode, which probably means 'sample extensively'", "0") //XRW TODO: positive integer if chis?  is 0 a flag value?
		+ XMLSchemaAttribute::attribute_w_default( "stub_energy_threshold", xsct_real, "after placement and minimization, what energy cutoff to use for each of the hotspots", "1.0")
		+ XMLSchemaAttribute::attribute_w_default( "max_cb_dist", xsct_real, "Maximum distance from ideal placement that is nonetheless considered a hit", "3.0" );
	if ( !attribute_w_name_in_attribute_list( "cb_force", StubSets_attlist ) ) {
		StubSets_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	StubSets_node->set_element_w_attributes( "StubSets", StubSets_attlist, "What are the Stubs we are placing?" );
	StubSets_node->set_kids_naming_func( & PlaceSimultaneously_StubSets_subelement_mangler );

	XMLSchemaRepeatableCTNodeOP StubSets_Add_node( new XMLSchemaRepeatableCTNode );
	StubSets_node->add_child(StubSets_Add_node);

	AttributeList StubSets_Add_attlist;
	StubSets_Add_attlist
		+ XMLSchemaAttribute::required_attribute( "stubfile", xs_string, "File that has a stub in it")
		+ XMLSchemaAttribute( "filter_name", xs_string, "Stub placement tied to this Filter");

	if ( !attribute_w_name_in_attribute_list( "cb_force", StubSets_Add_attlist ) ) {
		StubSets_Add_attlist + XMLSchemaAttribute( "cb_force", xsct_real, "the force on CB atoms" );
	}
	StubSets_Add_node->set_element_w_attributes( "Add", StubSets_Add_attlist, "What is one Stub we are placing?");

	//now tie all the child nodes to the root node and report
	PlaceSimultaneouslyMover_node->add_child(StubMinimize_node);
	PlaceSimultaneouslyMover_node->add_child(DesignMovers_node);
	PlaceSimultaneouslyMover_node->add_child(NotifyMovers_node);
	PlaceSimultaneouslyMover_node->add_child(StubSets_node);

	PlaceSimultaneouslyMover_node->recursively_write_ct_to_schema( xsd );
}

std::string PlaceSimultaneouslyMoverCreator::keyname() const {
	return PlaceSimultaneouslyMover::mover_name();
}

protocols::moves::MoverOP
PlaceSimultaneouslyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PlaceSimultaneouslyMover );
}

void PlaceSimultaneouslyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PlaceSimultaneouslyMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
