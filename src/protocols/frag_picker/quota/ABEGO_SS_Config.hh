// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Config_hh
#define INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Config_hh

#include <protocols/frag_picker/quota/QuotaConfig.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace frag_picker {
namespace quota {

/// @brief read a config file for quota selector
class ABEGO_SS_Config : public QuotaConfig {
public:

	/// @brief  Constructor reads a config file
	ABEGO_SS_Config(std::string & config_file_name);

	core::Size n_columns() { return bin_probs_[1].size(); }

	core::Size n_rows() { return bin_probs_.size(); }

	core::Size size() { return bin_probs_.size(); }

	std::string & source_file_name() { return source_file_name_; }

	core::Real probability(core::Size seq_pos,core::Size bin) { return bin_probs_[seq_pos][bin]; }

	core::Real highest_probability(core::Size);

	core::Size most_probable_bin(core::Size);

	utility::vector1< std::pair<core::Size,core::Size> > get_pool_bins(core::Size pool_id) { return pool_defs_[pool_id]; }
private:
	std::string source_file_name_;
	utility::vector1< utility::vector1< std::pair<core::Size,core::Size> > > pool_defs_;
	utility::vector1< utility::vector1<core::Real> > bin_probs_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Config_HH */
