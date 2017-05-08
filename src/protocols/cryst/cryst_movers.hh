// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_cryst_cryst_movers_hh
#define INCLUDED_protocols_cryst_cryst_movers_hh

#include <protocols/cryst/cryst_movers_creator.hh>

#include <utility>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/electron_density/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>


namespace protocols {
namespace cryst {


class ReportGradientsMover : public moves::Mover {
public:
	ReportGradientsMover() :
		Mover(), verbose_(false) {}

	ReportGradientsMover(core::scoring::ScoreFunctionOP sfin) :
		Mover(), verbose_(false), score_function_(std::move(sfin)) {}

	// XRW TEMP  std::string get_name() const override { return ReportGradientsMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new ReportGradientsMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;

	// compute gradients
	core::Real compute(core::pose::Pose & pose );

	// helper function normalizes gradients for a single atom
	core::Real normalization(core::pose::Pose & pose, core::id::AtomID atmid, core::scoring::ScoreFunctionOP sfxn );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool verbose_;
	core::scoring::ScoreFunctionOP score_function_;
	std::string outfile_;
};

class SetCrystWeightMover : public moves::Mover {
public:
	SetCrystWeightMover() :
		Mover(), autoset_wt_(true), cartesian_(false), weight_(0.0), weight_scale_(1.0), weight_min_(1.0) {}

	// XRW TEMP  std::string get_name() const override { return SetCrystWeightMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new SetCrystWeightMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool autoset_wt_;
	bool cartesian_;
	core::Real weight_, weight_scale_, weight_min_;

	core::scoring::ScoreFunctionOP score_function_;
	core::scoring::ScoreFunctionOP score_function_ref_;
	core::kinematics::MoveMapOP mm_;
};


class RecomputeDensityMapMover : public moves::Mover {
public:
	RecomputeDensityMapMover() : Mover(), keep_sidechains_(true) {}

	// XRW TEMP  std::string get_name() const override { return RecomputeDensityMapMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new RecomputeDensityMapMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	bool keep_sidechains_;
};


class LoadDensityMapMover : public moves::Mover {
public:
	LoadDensityMapMover() : Mover(), mapfile_("") {}

	// XRW TEMP  std::string get_name() const override { return LoadDensityMapMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new LoadDensityMapMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	std::string mapfile_;
	core::Real sc_scale_;
	core::Size window_;
};


class FitBfactorsMover : public moves::Mover {
public:
	FitBfactorsMover() : Mover(), adp_strategy_("individual"), b_min_(5.0), b_max_(5.0) {}

	// XRW TEMP  std::string get_name() const override { return FitBfactorsMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new FitBfactorsMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void randomize_bs( core::pose::Pose & pose );

	std::string adp_strategy_;
	core::Real b_min_, b_max_;
};

class UpdateSolventMover : public moves::Mover {
public:
	UpdateSolventMover() : Mover(), update_mask_(true), update_fcalc_(true), optimize_mask_(false), optimize_params_(false) {}

	// XRW TEMP  std::string get_name() const override { return UpdateSolventMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new UpdateSolventMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool update_mask_ , update_fcalc_, optimize_mask_, optimize_params_;
};


class TagPoseWithRefinementStatsMover : public moves::Mover {
public:
	TagPoseWithRefinementStatsMover() : Mover(), tag_(""), dump_pose_(false),report_grads_(false) {}

	// XRW TEMP  std::string get_name() const override { return TagPoseWithRefinementStatsMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new TagPoseWithRefinementStatsMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string tag_;
	bool dump_pose_,report_grads_;
};

class SetRefinementOptionsMover : public moves::Mover {
public:
	SetRefinementOptionsMover() :
		Mover(), res_high_(0.0), res_low_(0.0), sharpen_b_(0.0), twin_law_(""), algo_(""), target_(""), map_type_(""), setmap_type_(false)
	{}

	// XRW TEMP  std::string get_name() const override { return SetRefinementOptionsMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new SetRefinementOptionsMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real res_high_, res_low_;
	core::Real sharpen_b_;
	std::string twin_law_, algo_, target_, map_type_;
	utility::vector1<std::string> cif_files_;
	bool setmap_type_;
};


}
}

#endif
