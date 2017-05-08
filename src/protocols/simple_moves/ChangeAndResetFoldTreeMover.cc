// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ChangeAndResetFoldTreeMover.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/ChangeAndResetFoldTreeMover.hh>
#include <protocols/simple_moves/ChangeAndResetFoldTreeMoverCreator.hh>

#include <protocols/moves/ChangeFoldTreeMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/FoldTree.hh>

#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static THREAD_LOCAL basic::Tracer TR("ChangeAndResetFoldTreeMover");

namespace protocols {
namespace simple_moves {

using namespace protocols::moves;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;

ChangeAndResetFoldTreeMover::ChangeAndResetFoldTreeMover() :
	MoverApplyingMover("ChangeAndResetFoldTreeMover"),
	main_mover_(/* NULL */),
	ft_mover_(/* NULL */),
	scorefxn_(/* NULL */)
{

}

ChangeAndResetFoldTreeMover::ChangeAndResetFoldTreeMover(MoverOP main_mover):
	MoverApplyingMover("ChangeAndResetFoldTreeMover"),
	main_mover_(std::move(main_mover)),
	ft_mover_(/* NULL */),
	scorefxn_(/* NULL */)
{

}

ChangeAndResetFoldTreeMover::ChangeAndResetFoldTreeMover(MoverOP main_mover, ChangeFoldTreeMoverOP ft_mover):
	MoverApplyingMover("ChangeAndResetFoldTreeMover"),
	main_mover_(std::move(main_mover)),
	ft_mover_(std::move(ft_mover)),
	scorefxn_(/* NULL */)
{

}

ChangeAndResetFoldTreeMover::ChangeAndResetFoldTreeMover(
	MoverOP main_mover, ChangeFoldTreeMoverOP ft_mover, core::scoring::ScoreFunctionCOP scorefxn):
	MoverApplyingMover("ChangeAndResetFoldTreeMover"),
	main_mover_(std::move(main_mover)),
	ft_mover_(std::move(ft_mover)),
	scorefxn_(std::move(scorefxn))
{

}

ChangeAndResetFoldTreeMover::~ChangeAndResetFoldTreeMover() = default;

ChangeAndResetFoldTreeMover::ChangeAndResetFoldTreeMover(ChangeAndResetFoldTreeMover const & ) = default;

MoverOP
ChangeAndResetFoldTreeMover::clone() const {
	return MoverOP( new ChangeAndResetFoldTreeMover(*this) );
}

//ChangeAndResetFoldTreeMover & operator=(ChangeAndResetFoldTreeMover const & src) {
// return ChangeAndResetFoldTreeMover(src);
//}

MoverOP
ChangeAndResetFoldTreeMover::fresh_instance() const {
	return MoverOP( new ChangeAndResetFoldTreeMover );
}

// XRW TEMP std::string
// XRW TEMP ChangeAndResetFoldTreeMover::get_name() const {
// XRW TEMP  return "ChangeAndResetFoldTreeMover";
// XRW TEMP }

// XRW TEMP MoverOP
// XRW TEMP ChangeAndResetFoldTreeMoverCreator::create_mover() const {
// XRW TEMP  return MoverOP( new ChangeAndResetFoldTreeMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ChangeAndResetFoldTreeMoverCreator::keyname() const {
// XRW TEMP  return ChangeAndResetFoldTreeMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ChangeAndResetFoldTreeMover::mover_name() {
// XRW TEMP  return "ChangeAndResetFoldTreeMover";
// XRW TEMP }

//void
//ChangeAndResetFoldTreeMover::parse_my_tag(
//  TagCOP tag,
//  basic::datacache::DataMap&,
//  const Filters_map&,
//  const Movers_map&,
//  const Pose&) {
//
//}

void
ChangeAndResetFoldTreeMover::set_mover(protocols::moves::MoverOP main_mover){
	main_mover_ = main_mover;
}

protocols::moves::MoverOP
ChangeAndResetFoldTreeMover::mover() const {
	return main_mover_;
}

ChangeFoldTreeMoverOP
ChangeAndResetFoldTreeMover::ft_mover() const {
	return ft_mover_;
}

void
ChangeAndResetFoldTreeMover::set_ft_mover( protocols::moves::ChangeFoldTreeMoverOP ft_mover ) {
	ft_mover_ = ft_mover;
}

ScoreFunctionCOP
ChangeAndResetFoldTreeMover::scorefxn() const {
	return scorefxn_;
}

void
ChangeAndResetFoldTreeMover::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn;
}

void
ChangeAndResetFoldTreeMover::apply(core::pose::Pose& pose) {

	if ( (! main_mover_) || (! ft_mover_) ) {
		utility_exit_with_message("mover or FT mover not set in ChangeAndResetFoldTreeMover!");
	}


	TR << "Applying: " << main_mover_->get_name() << std::endl;

	//Print score info before
	if ( scorefxn_ ) {
		TR << "start: "<< scorefxn_->score(pose) << std::endl;
	}

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if ( ft_mover_ ) {
		ft_mover_->apply(pose);
	}
	main_mover_->apply(pose);
	pose.fold_tree(original_ft);

	//Print score info after
	if ( scorefxn_ ) {
		TR << "end: " << scorefxn_->score(pose) << std::endl;
	}


}

std::string ChangeAndResetFoldTreeMover::get_name() const {
	return mover_name();
}

std::string ChangeAndResetFoldTreeMover::mover_name() {
	return "ChangeAndResetFoldTreeMover";
}

void ChangeAndResetFoldTreeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Not valid for use with RosettaScripts; don't use it", attlist );
}

std::string ChangeAndResetFoldTreeMoverCreator::keyname() const {
	return ChangeAndResetFoldTreeMover::mover_name();
}

protocols::moves::MoverOP
ChangeAndResetFoldTreeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ChangeAndResetFoldTreeMover );
}

void ChangeAndResetFoldTreeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ChangeAndResetFoldTreeMover::provide_xml_schema( xsd );
}



} //simple_moves
} //Protocols
