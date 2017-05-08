// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/CacheableDataType.cc
/// @brief  enum indexing the data types stored in a Pose's internal DataCache

#include <core/pose/datacache/CacheableDataType.hh>

namespace core {
namespace pose {
namespace datacache {

utility::vector1< std::string > & CacheableDataType::enum2name_()
{
	static utility::vector1< std::string > enum_to_name;
	return enum_to_name;
}


std::map< std::string, CacheableDataType::Enum > & CacheableDataType::name2enum_()
{
	static std::map< std::string, CacheableDataType::Enum > name_to_enum;
	return name_to_enum;
}

std::string
CacheableDataType::get_name( CacheableDataType::Enum datatype) {
	initialize_name_map();
	return enum2name_()[ datatype ];
}


void
CacheableDataType::initialize_name_map() {
	static bool initialized_(false);

	if ( initialized_ ) return;

	name2enum_()["JOBDIST_OUTPUT_TAG"] = JOBDIST_OUTPUT_TAG;
	name2enum_()["ARBITRARY_FLOAT_DATA"] = ARBITRARY_FLOAT_DATA;
	name2enum_()["ARBITRARY_STRING_DATA"] = ARBITRARY_STRING_DATA;
	name2enum_()["STRING_MAP"] = STRING_MAP;
	name2enum_()["SCORE_LINE_STRINGS"] = SCORE_LINE_STRINGS;
	name2enum_()["SCORE_MAP"] = SCORE_MAP;
	name2enum_()["STM_STORED_TASKS"] = STM_STORED_TASKS;
	name2enum_()["STORED_RESIDUE_SUBSET"] = STORED_RESIDUE_SUBSET;
	name2enum_()["CONSTRAINT_GENERATOR"] = CONSTRAINT_GENERATOR;
	name2enum_()["POSE_BEFORE_CAVITIES_ADDED"] = POSE_BEFORE_CAVITIES_ADDED;
	name2enum_()["TEMPLATE_HYBRIDIZATION_HISTORY"] = TEMPLATE_HYBRIDIZATION_HISTORY;
	name2enum_()["NCS_RESIDUE_MAPPING"] = NCS_RESIDUE_MAPPING;
	name2enum_()["FULL_MODEL_INFO"] = FULL_MODEL_INFO;
	name2enum_()["VDW_REP_SCREEN_INFO"] = VDW_REP_SCREEN_INFO;
	name2enum_()["CDR_CLUSTER_INFO"] = CDR_CLUSTER_INFO;
	name2enum_()["NATIVE_ANTIBODY_SEQ"] = NATIVE_ANTIBODY_SEQ;
	name2enum_()["WRITEABLE_DATA"] = WRITEABLE_DATA;
	name2enum_()["POSITION_CONSERVED_RESIDUES"] = POSITION_CONSERVED_RESIDUES;
	name2enum_()["INTERFACE_DDG_MUTATION"] = INTERFACE_DDG_MUTATION;
	name2enum_()["BASE_PARTNER"] = BASE_PARTNER;
	name2enum_()["CEN_LIST_INFO"] = CEN_LIST_INFO;
	name2enum_()["SIGMOID_WEIGHTED_CEN_LIST"] = SIGMOID_WEIGHTED_CEN_LIST;
	name2enum_()["SIGMOID_WEIGHTED_D_CEN_LIST"] = SIGMOID_WEIGHTED_D_CEN_LIST;
	name2enum_()["RG_MINDATA"] = RG_MINDATA;
	name2enum_()["RG_LOCAL_MINDATA"] = RG_LOCAL_MINDATA;
	name2enum_()["MEMBRANE_TOPOLOGY"] = MEMBRANE_TOPOLOGY;
	name2enum_()["MEMBRANE_EMBED"] = MEMBRANE_EMBED;
	name2enum_()["INTERFACE_INFO"] = INTERFACE_INFO;
	name2enum_()["SS_INFO"] = SS_INFO;
	name2enum_()["SS_KILLHAIRPINS_INFO"] = SS_KILLHAIRPINS_INFO;
	name2enum_()["RNA_SCORING_INFO"] = RNA_SCORING_INFO;
	name2enum_()["RNA_SECSTRUCT_INFO"] = RNA_SECSTRUCT_INFO;
	name2enum_()["WATER_PACKING_INFO"] = WATER_PACKING_INFO;
	name2enum_()["HOLES_POSE_INFO"] = HOLES_POSE_INFO;
	name2enum_()["DAB_SASA_POSE_INFO"] = DAB_SASA_POSE_INFO;
	name2enum_()["DAB_SEV_POSE_INFO"] = DAB_SEV_POSE_INFO;
	name2enum_()["CHEMICAL_SHIFT_ANISOTROPY_DATA"] = CHEMICAL_SHIFT_ANISOTROPY_DATA;
	name2enum_()["RESIDUAL_DIPOLAR_COUPLING_DATA"] = RESIDUAL_DIPOLAR_COUPLING_DATA;
	name2enum_()["RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL"] = RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL;
	name2enum_()["RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA"] = RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA;
	name2enum_()["FIBER_DIFFRACTION_CEN_SCATTERING"] = FIBER_DIFFRACTION_CEN_SCATTERING;
	name2enum_()["FIBER_DIFFRACTION_FA_SCATTERING"] = FIBER_DIFFRACTION_FA_SCATTERING;
	name2enum_()["PSEUDOCONTACT_SHIFT_DATA"] = PSEUDOCONTACT_SHIFT_DATA;
	name2enum_()["TS1_PSEUDOCONTACT_SHIFT_DATA"] = TS1_PSEUDOCONTACT_SHIFT_DATA;
	name2enum_()["TS2_PSEUDOCONTACT_SHIFT_DATA"] = TS2_PSEUDOCONTACT_SHIFT_DATA;
	name2enum_()["TS3_PSEUDOCONTACT_SHIFT_DATA"] = TS3_PSEUDOCONTACT_SHIFT_DATA;
	name2enum_()["TS4_PSEUDOCONTACT_SHIFT_DATA"] = TS4_PSEUDOCONTACT_SHIFT_DATA;
	name2enum_()["PSEUDOCONTACT_SHIFT_MULTI_DATA"] = PSEUDOCONTACT_SHIFT_MULTI_DATA;
	name2enum_()["GEN_BORN_POSE_INFO"] = GEN_BORN_POSE_INFO;
	name2enum_()["MULTIPOLE_POSE_INFO"] = MULTIPOLE_POSE_INFO;
	name2enum_()["VDWTINKER_POSE_INFO"] = VDWTINKER_POSE_INFO;
	name2enum_()["FACTS_POSE_INFO"] = FACTS_POSE_INFO;
	name2enum_()["MEMBRANE_FAEMBED"] = MEMBRANE_FAEMBED;
	name2enum_()["LK_BALL_POSE_INFO"] = LK_BALL_POSE_INFO;
	name2enum_()["PB_LIFETIME_CACHE"] = PB_LIFETIME_CACHE;
	name2enum_()["FAELEC_CONTEXT_DATA"] = FAELEC_CONTEXT_DATA;
	name2enum_()["POSE_SEQUENCE"] = POSE_SEQUENCE;

	runtime_assert( name2enum_().size() == CacheableDataType::num_cacheable_data_types );

	enum2name_().resize( CacheableDataType::num_cacheable_data_types );
	for ( std::map< std::string, CacheableDataType::Enum >::const_iterator iter = name2enum_().begin(),
			iter_end = name2enum_().end(); iter != iter_end; ++iter ) {
		enum2name_()[ iter->second ] = iter->first;
	}

	initialized_ = true;
}

} // namespace datacache
} // namespace pose
} // namespace core
