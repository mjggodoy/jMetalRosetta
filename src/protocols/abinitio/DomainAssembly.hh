// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/abinitio/DomainAssembly.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_abinitio_DomainAssembly_hh
#define INCLUDED_protocols_abinitio_DomainAssembly_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

/// @brief insert fragments in a linker region. Very similar to what's in looprelax_main
class DomainAssembly : public protocols::moves::Mover
{
	typedef core::fragment::FragSetOP FragSetOP;
	typedef core::fragment::FragSetCOP FragSetCOP;
public:
	DomainAssembly();
	DomainAssembly(
		core::Size const linker_start,
		core::Size const linker_end,
		FragSetOP fragset_large,
		FragSetOP fragset_small
	);
	~DomainAssembly() override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void apply( core::pose::Pose & pose ) override;
	// XRW TEMP  std::string get_name() const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size linker_start_, linker_end_;
	core::fragment::FragSetOP fragset_large_, fragset_small_;
	bool fragments_set_;
};

} // abinitio
} // protocols


#endif /*INCLUDED_abinitio_movers_DomainAssembly_HH*/
