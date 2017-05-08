// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/types.cc
/// @brief  Functions for manipulating hydrogen bond types
/// @author Matthew O'Meara

#include <core/chemical/types.hh>
#include <core/scoring/hbonds/types.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <sstream>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace hbonds {

static THREAD_LOCAL basic::Tracer TR("core.scoring.hbonds.types");

Size
hb_eval_type(
	HBDonChemType don_chem_type,
	HBAccChemType acc_chem_type,
	HBSeqSep seq_sep_type){

	return
		(don_chem_type-1)*(hbacc_MAX-2)*(seq_sep_MAX-2) +
		(acc_chem_type-1)*(seq_sep_MAX-2) +
		(seq_sep_type-1) + 1;
}


/// @brief converts a string into an HBEvalType
HBEvalType string_to_hb_eval_type(std::string const& hbe_str) {

	if ( hbe_str == "hbe_dH2OaHXL" ) {
		return hbe_dH2OaHXL;
	}
	if ( hbe_str == "hbe_GENERIC_SP3SCSC_LR" ) {
		return hbe_GENERIC_SP3SCSC_LR;
	}

	std::stringstream msg;
	msg << "Cannot convert " << hbe_str << " into an HBEvalType" << std::endl;
	throw( utility::excn::EXCN_Msg_Exception( msg.str()  ));
}


void
HBEval_lookup_initializer( ObjexxFCL::FArray3D<HBEvalType> & hbe )
{
	hbe = hbe_UNKNOWN; // General Default - If not otherwise redone below.

	hbe(hbdon_NONE, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_PBA, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_CXA, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_IMD, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_IME, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_IND, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_AMO, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_GDE, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_GDH, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_AHX, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_HXL, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_H2O, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_GENERIC_BB, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_GENERIC_SC, hbacc_NONE, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_PBA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_CXA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_CXL, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_IMD, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_IME, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_AHX, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_HXL, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_PCA_DNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_PES_DNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_RRI_DNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_PCA_RNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_PES_RNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_RRI_RNA, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_H2O, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_NONE;
	hbe(hbdon_NONE, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_NONE;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_M4) = hbe_dPBAaPBAsepM4helix;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_M3) = hbe_dPBAaPBAsepM3turn;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_M2) = hbe_dPBAaPBAsepM2turn;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_PM1) = hbe_dPBAaPBAsepPM1;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_P2) = hbe_dPBAaPBAsepP2turn;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_P3) = hbe_dPBAaPBAsepP3turn;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_P4) = hbe_dPBAaPBAsepP4helix;
	hbe(hbdon_PBA, hbacc_PBA, seq_sep_other) = hbe_dPBAaPBAsepother;
	hbe(hbdon_CXA, hbacc_PBA, seq_sep_PM1) = hbe_dCXAaPBAsepPM1;
	hbe(hbdon_IMD, hbacc_PBA, seq_sep_PM1) = hbe_dIMDaPBAsepPM1;
	hbe(hbdon_IME, hbacc_PBA, seq_sep_PM1) = hbe_dIMEaPBAsepPM1;
	hbe(hbdon_IND, hbacc_PBA, seq_sep_PM1) = hbe_dINDaPBAsepPM1;
	hbe(hbdon_AMO, hbacc_PBA, seq_sep_PM1) = hbe_dAMOaPBAsepPM1;
	hbe(hbdon_GDE, hbacc_PBA, seq_sep_PM1) = hbe_dGDEaPBAsepPM1;
	hbe(hbdon_GDH, hbacc_PBA, seq_sep_PM1) = hbe_dGDHaPBAsepPM1;
	hbe(hbdon_AHX, hbacc_PBA, seq_sep_PM1) = hbe_dAHXaPBAsepPM1;
	hbe(hbdon_HXL, hbacc_PBA, seq_sep_PM1) = hbe_dHXLaPBAsepPM1;
	hbe(hbdon_CXA, hbacc_PBA, seq_sep_other) = hbe_dCXAaPBAsepother;
	hbe(hbdon_IMD, hbacc_PBA, seq_sep_other) = hbe_dIMDaPBAsepother;
	hbe(hbdon_IME, hbacc_PBA, seq_sep_other) = hbe_dIMEaPBAsepother;
	hbe(hbdon_IND, hbacc_PBA, seq_sep_other) = hbe_dINDaPBAsepother;
	hbe(hbdon_AMO, hbacc_PBA, seq_sep_other) = hbe_dAMOaPBAsepother;
	hbe(hbdon_GDE, hbacc_PBA, seq_sep_other) = hbe_dGDEaPBAsepother;
	hbe(hbdon_GDH, hbacc_PBA, seq_sep_other) = hbe_dGDHaPBAsepother;
	hbe(hbdon_AHX, hbacc_PBA, seq_sep_other) = hbe_dAHXaPBAsepother;
	hbe(hbdon_HXL, hbacc_PBA, seq_sep_other) = hbe_dHXLaPBAsepother;
	hbe(hbdon_H2O, hbacc_PBA, seq_sep_other) = hbe_dH2OaPBA;
	hbe(hbdon_PBA, hbacc_CXA, seq_sep_PM1) = hbe_dPBAaCXAsepPM1;
	hbe(hbdon_PBA, hbacc_CXA, seq_sep_other) = hbe_dPBAaCXAsepother;
	hbe(hbdon_CXA, hbacc_CXA, seq_sep_other) = hbe_dCXAaCXA;
	hbe(hbdon_IMD, hbacc_CXA, seq_sep_other) = hbe_dIMDaCXA;
	hbe(hbdon_IME, hbacc_CXA, seq_sep_other) = hbe_dIMEaCXA;
	hbe(hbdon_IND, hbacc_CXA, seq_sep_other) = hbe_dINDaCXA;
	hbe(hbdon_AMO, hbacc_CXA, seq_sep_other) = hbe_dAMOaCXA;
	hbe(hbdon_GDE, hbacc_CXA, seq_sep_other) = hbe_dGDEaCXA;
	hbe(hbdon_GDH, hbacc_CXA, seq_sep_other) = hbe_dGDHaCXA;
	hbe(hbdon_AHX, hbacc_CXA, seq_sep_other) = hbe_dAHXaCXA;
	hbe(hbdon_HXL, hbacc_CXA, seq_sep_other) = hbe_dHXLaCXA;
	hbe(hbdon_H2O, hbacc_CXA, seq_sep_other) = hbe_dH2OaCXA;
	hbe(hbdon_PBA, hbacc_CXL, seq_sep_PM1) = hbe_dPBAaCXLsepPM1;
	hbe(hbdon_PBA, hbacc_CXL, seq_sep_other) = hbe_dPBAaCXLsepother;
	hbe(hbdon_CXA, hbacc_CXL, seq_sep_other) = hbe_dCXAaCXL;
	hbe(hbdon_IMD, hbacc_CXL, seq_sep_other) = hbe_dIMDaCXL;
	hbe(hbdon_IME, hbacc_CXL, seq_sep_other) = hbe_dIMEaCXL;
	hbe(hbdon_IND, hbacc_CXL, seq_sep_other) = hbe_dINDaCXL;
	hbe(hbdon_AMO, hbacc_CXL, seq_sep_other) = hbe_dAMOaCXL;
	hbe(hbdon_GDE, hbacc_CXL, seq_sep_other) = hbe_dGDEaCXL;
	hbe(hbdon_GDH, hbacc_CXL, seq_sep_other) = hbe_dGDHaCXL;
	hbe(hbdon_AHX, hbacc_CXL, seq_sep_other) = hbe_dAHXaCXL;
	hbe(hbdon_HXL, hbacc_CXL, seq_sep_other) = hbe_dHXLaCXL;
	hbe(hbdon_H2O, hbacc_CXL, seq_sep_other) = hbe_dH2OaCXL;
	hbe(hbdon_PBA, hbacc_IMD, seq_sep_PM1) = hbe_dPBAaIMDsepPM1;
	hbe(hbdon_PBA, hbacc_IMD, seq_sep_other) = hbe_dPBAaIMDsepother;
	hbe(hbdon_CXA, hbacc_IMD, seq_sep_other) = hbe_dCXAaIMD;
	hbe(hbdon_IMD, hbacc_IMD, seq_sep_other) = hbe_dIMDaIMD;
	hbe(hbdon_IME, hbacc_IMD, seq_sep_other) = hbe_dIMEaIMD;
	hbe(hbdon_IND, hbacc_IMD, seq_sep_other) = hbe_dINDaIMD;
	hbe(hbdon_AMO, hbacc_IMD, seq_sep_other) = hbe_dAMOaIMD;
	hbe(hbdon_GDE, hbacc_IMD, seq_sep_other) = hbe_dGDEaIMD;
	hbe(hbdon_GDH, hbacc_IMD, seq_sep_other) = hbe_dGDHaIMD;
	hbe(hbdon_AHX, hbacc_IMD, seq_sep_other) = hbe_dAHXaIMD;
	hbe(hbdon_HXL, hbacc_IMD, seq_sep_other) = hbe_dHXLaIMD;
	hbe(hbdon_H2O, hbacc_IMD, seq_sep_other) = hbe_dH2OaIMD;
	hbe(hbdon_PBA, hbacc_IME, seq_sep_PM1) = hbe_dPBAaIMEsepPM1;
	hbe(hbdon_PBA, hbacc_IME, seq_sep_other) = hbe_dPBAaIMEsepother;
	hbe(hbdon_CXA, hbacc_IME, seq_sep_other) = hbe_dCXAaIME;
	hbe(hbdon_IMD, hbacc_IME, seq_sep_other) = hbe_dIMDaIME;
	hbe(hbdon_IME, hbacc_IME, seq_sep_other) = hbe_dIMEaIME;
	hbe(hbdon_IND, hbacc_IME, seq_sep_other) = hbe_dINDaIME;
	hbe(hbdon_AMO, hbacc_IME, seq_sep_other) = hbe_dAMOaIME;
	hbe(hbdon_GDE, hbacc_IME, seq_sep_other) = hbe_dGDEaIME;
	hbe(hbdon_GDH, hbacc_IME, seq_sep_other) = hbe_dGDHaIME;
	hbe(hbdon_AHX, hbacc_IME, seq_sep_other) = hbe_dAHXaIME;
	hbe(hbdon_HXL, hbacc_IME, seq_sep_other) = hbe_dHXLaIME;
	hbe(hbdon_H2O, hbacc_IME, seq_sep_other) = hbe_dH2OaIME;
	hbe(hbdon_PBA, hbacc_AHX, seq_sep_PM1) = hbe_dPBAaAHXsepPM1;
	hbe(hbdon_PBA, hbacc_AHX, seq_sep_other) = hbe_dPBAaAHXsepother;
	hbe(hbdon_CXA, hbacc_AHX, seq_sep_other) = hbe_dCXAaAHX;
	hbe(hbdon_IMD, hbacc_AHX, seq_sep_other) = hbe_dIMDaAHX;
	hbe(hbdon_IME, hbacc_AHX, seq_sep_other) = hbe_dIMEaAHX;
	hbe(hbdon_IND, hbacc_AHX, seq_sep_other) = hbe_dINDaAHX;
	hbe(hbdon_AMO, hbacc_AHX, seq_sep_other) = hbe_dAMOaAHX;
	hbe(hbdon_GDE, hbacc_AHX, seq_sep_other) = hbe_dGDEaAHX;
	hbe(hbdon_GDH, hbacc_AHX, seq_sep_other) = hbe_dGDHaAHX;
	hbe(hbdon_AHX, hbacc_AHX, seq_sep_other) = hbe_dAHXaAHX;
	hbe(hbdon_HXL, hbacc_AHX, seq_sep_other) = hbe_dHXLaAHX;
	hbe(hbdon_H2O, hbacc_AHX, seq_sep_other) = hbe_dH2OaAHX;
	hbe(hbdon_PBA, hbacc_HXL, seq_sep_PM1) = hbe_dPBAaHXLsepPM1;
	hbe(hbdon_PBA, hbacc_HXL, seq_sep_other) = hbe_dPBAaHXLsepother;
	hbe(hbdon_CXA, hbacc_HXL, seq_sep_other) = hbe_dCXAaHXL;
	hbe(hbdon_IMD, hbacc_HXL, seq_sep_other) = hbe_dIMDaHXL;
	hbe(hbdon_IME, hbacc_HXL, seq_sep_other) = hbe_dIMEaHXL;
	hbe(hbdon_IND, hbacc_HXL, seq_sep_other) = hbe_dINDaHXL;
	hbe(hbdon_AMO, hbacc_HXL, seq_sep_other) = hbe_dAMOaHXL;
	hbe(hbdon_GDE, hbacc_HXL, seq_sep_other) = hbe_dGDEaHXL;
	hbe(hbdon_GDH, hbacc_HXL, seq_sep_other) = hbe_dGDHaHXL;
	hbe(hbdon_AHX, hbacc_HXL, seq_sep_other) = hbe_dAHXaHXL;
	hbe(hbdon_HXL, hbacc_HXL, seq_sep_other) = hbe_dHXLaHXL;
	hbe(hbdon_H2O, hbacc_HXL, seq_sep_other) = hbe_dH2OaHXL;
	hbe(hbdon_PBA, hbacc_PCA_DNA, seq_sep_PM1) = hbe_dPBAaPCA_DNAsepPM1;
	hbe(hbdon_PBA, hbacc_PCA_DNA, seq_sep_other) = hbe_dPBAaPCA_DNAsepother;
	hbe(hbdon_CXA, hbacc_PCA_DNA, seq_sep_other) = hbe_dCXAaPCA_DNA;
	hbe(hbdon_IMD, hbacc_PCA_DNA, seq_sep_other) = hbe_dIMDaPCA_DNA;
	hbe(hbdon_IME, hbacc_PCA_DNA, seq_sep_other) = hbe_dIMEaPCA_DNA;
	hbe(hbdon_IND, hbacc_PCA_DNA, seq_sep_other) = hbe_dINDaPCA_DNA;
	hbe(hbdon_AMO, hbacc_PCA_DNA, seq_sep_other) = hbe_dAMOaPCA_DNA;
	hbe(hbdon_GDE, hbacc_PCA_DNA, seq_sep_other) = hbe_dGDEaPCA_DNA;
	hbe(hbdon_GDH, hbacc_PCA_DNA, seq_sep_other) = hbe_dGDHaPCA_DNA;
	hbe(hbdon_AHX, hbacc_PCA_DNA, seq_sep_other) = hbe_dAHXaPCA_DNA;
	hbe(hbdon_HXL, hbacc_PCA_DNA, seq_sep_other) = hbe_dHXLaPCA_DNA;
	hbe(hbdon_H2O, hbacc_PCA_DNA, seq_sep_other) = hbe_dH2OaPCA_DNA;
	hbe(hbdon_PBA, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dPBAaPCA_RNAsepPM1;
	hbe(hbdon_PBA, hbacc_PCA_RNA, seq_sep_other) = hbe_dPBAaPCA_RNAsepother;
	hbe(hbdon_CXA, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dCXAaPCA_RNAsepPM1;
	hbe(hbdon_CXA, hbacc_PCA_RNA, seq_sep_other) = hbe_dCXAaPCA_RNAsepother;
	hbe(hbdon_IMD, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dIMDaPCA_RNAsepPM1;
	hbe(hbdon_IMD, hbacc_PCA_RNA, seq_sep_other) = hbe_dIMDaPCA_RNAsepother;
	hbe(hbdon_IME, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dIMEaPCA_RNAsepPM1;
	hbe(hbdon_IME, hbacc_PCA_RNA, seq_sep_other) = hbe_dIMEaPCA_RNAsepother;
	hbe(hbdon_IND, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dINDaPCA_RNAsepPM1;
	hbe(hbdon_IND, hbacc_PCA_RNA, seq_sep_other) = hbe_dINDaPCA_RNAsepother;
	hbe(hbdon_AMO, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dAMOaPCA_RNAsepPM1;
	hbe(hbdon_AMO, hbacc_PCA_RNA, seq_sep_other) = hbe_dAMOaPCA_RNAsepother;
	hbe(hbdon_GDE, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dGDEaPCA_RNAsepPM1;
	hbe(hbdon_GDE, hbacc_PCA_RNA, seq_sep_other) = hbe_dGDEaPCA_RNAsepother;
	hbe(hbdon_GDH, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dGDHaPCA_RNAsepPM1;
	hbe(hbdon_GDH, hbacc_PCA_RNA, seq_sep_other) = hbe_dGDHaPCA_RNAsepother;
	hbe(hbdon_AHX, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dAHXaPCA_RNAsepPM1;
	hbe(hbdon_AHX, hbacc_PCA_RNA, seq_sep_other) = hbe_dAHXaPCA_RNAsepother;
	hbe(hbdon_HXL, hbacc_PCA_RNA, seq_sep_PM1) = hbe_dHXLaPCA_RNAsepPM1;
	hbe(hbdon_HXL, hbacc_PCA_RNA, seq_sep_other) = hbe_dHXLaPCA_RNAsepother;
	hbe(hbdon_H2O, hbacc_PCA_RNA, seq_sep_other) = hbe_dH2OaPCA_RNA;
	hbe(hbdon_PBA, hbacc_PES_DNA, seq_sep_PM1) = hbe_dPBAaPES_DNAsepPM1;
	hbe(hbdon_PBA, hbacc_PES_DNA, seq_sep_other) = hbe_dPBAaPES_DNAsepother;
	hbe(hbdon_CXA, hbacc_PES_DNA, seq_sep_other) = hbe_dCXAaPES_DNA;
	hbe(hbdon_IMD, hbacc_PES_DNA, seq_sep_other) = hbe_dIMDaPES_DNA;
	hbe(hbdon_IME, hbacc_PES_DNA, seq_sep_other) = hbe_dIMEaPES_DNA;
	hbe(hbdon_IND, hbacc_PES_DNA, seq_sep_other) = hbe_dINDaPES_DNA;
	hbe(hbdon_AMO, hbacc_PES_DNA, seq_sep_other) = hbe_dAMOaPES_DNA;
	hbe(hbdon_GDE, hbacc_PES_DNA, seq_sep_other) = hbe_dGDEaPES_DNA;
	hbe(hbdon_GDH, hbacc_PES_DNA, seq_sep_other) = hbe_dGDHaPES_DNA;
	hbe(hbdon_AHX, hbacc_PES_DNA, seq_sep_other) = hbe_dAHXaPES_DNA;
	hbe(hbdon_HXL, hbacc_PES_DNA, seq_sep_other) = hbe_dHXLaPES_DNA;
	hbe(hbdon_H2O, hbacc_PES_DNA, seq_sep_other) = hbe_dH2OaPES_DNA;
	hbe(hbdon_PBA, hbacc_PES_RNA, seq_sep_PM1) = hbe_dPBAaPES_RNAsepPM1;
	hbe(hbdon_PBA, hbacc_PES_RNA, seq_sep_other) = hbe_dPBAaPES_RNAsepother;
	hbe(hbdon_CXA, hbacc_PES_RNA, seq_sep_PM1) = hbe_dCXAaPES_RNAsepPM1;
	hbe(hbdon_CXA, hbacc_PES_RNA, seq_sep_other) = hbe_dCXAaPES_RNAsepother;
	hbe(hbdon_IMD, hbacc_PES_RNA, seq_sep_PM1) = hbe_dIMDaPES_RNAsepPM1;
	hbe(hbdon_IMD, hbacc_PES_RNA, seq_sep_other) = hbe_dIMDaPES_RNAsepother;
	hbe(hbdon_IME, hbacc_PES_RNA, seq_sep_PM1) = hbe_dIMEaPES_RNAsepPM1;
	hbe(hbdon_IME, hbacc_PES_RNA, seq_sep_other) = hbe_dIMEaPES_RNAsepother;
	hbe(hbdon_IND, hbacc_PES_RNA, seq_sep_PM1) = hbe_dINDaPES_RNAsepPM1;
	hbe(hbdon_IND, hbacc_PES_RNA, seq_sep_other) = hbe_dINDaPES_RNAsepother;
	hbe(hbdon_AMO, hbacc_PES_RNA, seq_sep_PM1) = hbe_dAMOaPES_RNAsepPM1;
	hbe(hbdon_AMO, hbacc_PES_RNA, seq_sep_other) = hbe_dAMOaPES_RNAsepother;
	hbe(hbdon_GDE, hbacc_PES_RNA, seq_sep_PM1) = hbe_dGDEaPES_RNAsepPM1;
	hbe(hbdon_GDE, hbacc_PES_RNA, seq_sep_other) = hbe_dGDEaPES_RNAsepother;
	hbe(hbdon_GDH, hbacc_PES_RNA, seq_sep_PM1) = hbe_dGDHaPES_RNAsepPM1;
	hbe(hbdon_GDH, hbacc_PES_RNA, seq_sep_other) = hbe_dGDHaPES_RNAsepother;
	hbe(hbdon_AHX, hbacc_PES_RNA, seq_sep_PM1) = hbe_dAHXaPES_RNAsepPM1;
	hbe(hbdon_AHX, hbacc_PES_RNA, seq_sep_other) = hbe_dAHXaPES_RNAsepother;
	hbe(hbdon_HXL, hbacc_PES_RNA, seq_sep_PM1) = hbe_dHXLaPES_RNAsepPM1;
	hbe(hbdon_HXL, hbacc_PES_RNA, seq_sep_other) = hbe_dHXLaPES_RNAsepother;
	hbe(hbdon_H2O, hbacc_PES_RNA, seq_sep_other) = hbe_dH2OaPES_RNA;
	hbe(hbdon_PBA, hbacc_RRI_DNA, seq_sep_PM1) = hbe_dPBAaRRI_DNAsepPM1;
	hbe(hbdon_PBA, hbacc_RRI_DNA, seq_sep_other) = hbe_dPBAaRRI_DNAsepother;
	hbe(hbdon_CXA, hbacc_RRI_DNA, seq_sep_other) = hbe_dCXAaRRI_DNA;
	hbe(hbdon_IMD, hbacc_RRI_DNA, seq_sep_other) = hbe_dIMDaRRI_DNA;
	hbe(hbdon_IME, hbacc_RRI_DNA, seq_sep_other) = hbe_dIMEaRRI_DNA;
	hbe(hbdon_IND, hbacc_RRI_DNA, seq_sep_other) = hbe_dINDaRRI_DNA;
	hbe(hbdon_AMO, hbacc_RRI_DNA, seq_sep_other) = hbe_dAMOaRRI_DNA;
	hbe(hbdon_GDE, hbacc_RRI_DNA, seq_sep_other) = hbe_dGDEaRRI_DNA;
	hbe(hbdon_GDH, hbacc_RRI_DNA, seq_sep_other) = hbe_dGDHaRRI_DNA;
	hbe(hbdon_AHX, hbacc_RRI_DNA, seq_sep_other) = hbe_dAHXaRRI_DNA;
	hbe(hbdon_HXL, hbacc_RRI_DNA, seq_sep_other) = hbe_dHXLaRRI_DNA;
	hbe(hbdon_H2O, hbacc_RRI_DNA, seq_sep_other) = hbe_dH2OaRRI_DNA;
	hbe(hbdon_PBA, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dPBAaRRI_RNAsepPM1;
	hbe(hbdon_PBA, hbacc_RRI_RNA, seq_sep_other) = hbe_dPBAaRRI_RNAsepother;
	hbe(hbdon_CXA, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dCXAaRRI_RNAsepPM1;
	hbe(hbdon_CXA, hbacc_RRI_RNA, seq_sep_other) = hbe_dCXAaRRI_RNAsepother;
	hbe(hbdon_IMD, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dIMDaRRI_RNAsepPM1;
	hbe(hbdon_IMD, hbacc_RRI_RNA, seq_sep_other) = hbe_dIMDaRRI_RNAsepother;
	hbe(hbdon_IME, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dIMEaRRI_RNAsepPM1;
	hbe(hbdon_IME, hbacc_RRI_RNA, seq_sep_other) = hbe_dIMEaRRI_RNAsepother;
	hbe(hbdon_IND, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dINDaRRI_RNAsepPM1;
	hbe(hbdon_IND, hbacc_RRI_RNA, seq_sep_other) = hbe_dINDaRRI_RNAsepother;
	hbe(hbdon_AMO, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dAMOaRRI_RNAsepPM1;
	hbe(hbdon_AMO, hbacc_RRI_RNA, seq_sep_other) = hbe_dAMOaRRI_RNAsepother;
	hbe(hbdon_GDE, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dGDEaRRI_RNAsepPM1;
	hbe(hbdon_GDE, hbacc_RRI_RNA, seq_sep_other) = hbe_dGDEaRRI_RNAsepother;
	hbe(hbdon_GDH, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dGDHaRRI_RNAsepPM1;
	hbe(hbdon_GDH, hbacc_RRI_RNA, seq_sep_other) = hbe_dGDHaRRI_RNAsepother;
	hbe(hbdon_AHX, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dAHXaRRI_RNAsepPM1;
	hbe(hbdon_AHX, hbacc_RRI_RNA, seq_sep_other) = hbe_dAHXaRRI_RNAsepother;
	hbe(hbdon_HXL, hbacc_RRI_RNA, seq_sep_PM1) = hbe_dHXLaRRI_RNAsepPM1;
	hbe(hbdon_HXL, hbacc_RRI_RNA, seq_sep_other) = hbe_dHXLaRRI_RNAsepother;
	hbe(hbdon_H2O, hbacc_RRI_RNA, seq_sep_other) = hbe_dH2OaRRI_RNA;
	hbe(hbdon_PBA, hbacc_H2O, seq_sep_other) = hbe_dPBAaH2O;
	hbe(hbdon_CXA, hbacc_H2O, seq_sep_other) = hbe_dCXAaH2O;
	hbe(hbdon_IMD, hbacc_H2O, seq_sep_other) = hbe_dIMDaH2O;
	hbe(hbdon_IME, hbacc_H2O, seq_sep_other) = hbe_dIMEaH2O;
	hbe(hbdon_IND, hbacc_H2O, seq_sep_other) = hbe_dINDaH2O;
	hbe(hbdon_AMO, hbacc_H2O, seq_sep_other) = hbe_dAMOaH2O;
	hbe(hbdon_GDE, hbacc_H2O, seq_sep_other) = hbe_dGDEaH2O;
	hbe(hbdon_GDH, hbacc_H2O, seq_sep_other) = hbe_dGDHaH2O;
	hbe(hbdon_AHX, hbacc_H2O, seq_sep_other) = hbe_dAHXaH2O;
	hbe(hbdon_HXL, hbacc_H2O, seq_sep_other) = hbe_dHXLaH2O;
	hbe(hbdon_H2O, hbacc_H2O, seq_sep_other) = hbe_dH2OaH2O;
	hbe(hbdon_GENERIC_BB, hbacc_PBA, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_PBA, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BB_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BB_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBB_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBB_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_CXA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_CXA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_CXL, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_CXL, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_PCA_DNA, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_PCA_DNA, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_PCA_RNA, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_PCA_RNA, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_AHX, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_AHX, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_HXL, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_HXL, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_H2O, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_PES_DNA, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_PES_DNA, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_PES_RNA, seq_sep_PM1) = hbe_GENERIC_SP2BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_PES_RNA, seq_sep_other) = hbe_GENERIC_SP2BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_RRI_DNA, seq_sep_PM1) = hbe_GENERIC_SP3BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_RRI_DNA, seq_sep_other) = hbe_GENERIC_SP3BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_RRI_RNA, seq_sep_PM1) = hbe_GENERIC_SP3BB_SR;
	hbe(hbdon_GENERIC_BB, hbacc_RRI_RNA, seq_sep_other) = hbe_GENERIC_SP3BB_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_IMD, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_IMD, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_IME, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_IME, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GENERIC_BB, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_PBA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_PBA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_PES_DNA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_PES_DNA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_PES_RNA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_PES_RNA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP2BB, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP2BB, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_RRI_DNA, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_RRI_DNA, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_RRI_RNA, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_RRI_RNA, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP3BB, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP3BB, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_RINGBB, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_RINGBB, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_CXA, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_CXA, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_CXL, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_CXL, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_PCA_DNA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_PCA_DNA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_PCA_RNA, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_PCA_RNA, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2BSC_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2BSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP2SC, seq_sep_PM1) = hbe_GENERIC_SP2SCSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP2SC, seq_sep_other) = hbe_GENERIC_SP2SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_AHX, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_AHX, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_HXL, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_HXL, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_H2O, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_H2O, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3BSC_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3BSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_SP3SC, seq_sep_other) = hbe_GENERIC_SP3SCSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP3SC, seq_sep_PM1) = hbe_GENERIC_SP3SCSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_SP3SC, seq_sep_other) = string_to_hb_eval_type(
		basic::options::option[basic::options::OptionKeys::score::hbe_for_dH2O_aGEN_SP3SC_ssother]);
	hbe(hbdon_GENERIC_SC, hbacc_IMD, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_IMD, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_IME, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_IME, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_GENERIC_SC, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_PBA, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGBSC_SR;
	hbe(hbdon_PBA, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGBSC_LR;
	hbe(hbdon_CXA, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_CXA, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_IMD, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_IMD, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_IME, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_IME, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_IND, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_IND, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_AMO, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_AMO, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_GDE, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_GDE, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_GDH, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_GDH, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_AHX, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_AHX, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_HXL, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_HXL, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
	hbe(hbdon_H2O, hbacc_GENERIC_RINGSC, seq_sep_PM1) = hbe_GENERIC_RINGSCSC_SR;
	hbe(hbdon_H2O, hbacc_GENERIC_RINGSC, seq_sep_other) = hbe_GENERIC_RINGSCSC_LR;
}
utility::pointer::shared_ptr< ObjexxFCL::FArray3D<HBEvalType> const > HBEval_lookup;

/// @brief makes it explicit that the HBEval lookup table is initialized somewhere
void initialize_HBEval_lookup() {

	HBEval_lookup = utility::pointer::shared_ptr< ObjexxFCL::FArray3D< HBEvalType > const > (
		new ObjexxFCL::FArray3D<HBEvalType>(hbdon_MAX, hbacc_MAX, seq_sep_MAX, HBEval_lookup_initializer));
}

chemical::Hybridization
get_hbe_acc_hybrid( HBEvalType const & hbe )
{
	using namespace chemical;
	switch(hbe){
	case hbe_NONE : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepM4helix : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepM3turn : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepM2turn : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepP2turn : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepP3turn : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepP4helix : return SP2_HYBRID; break;
	case hbe_dPBAaPBAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMDaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMEaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dINDaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAMOaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDEaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDHaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAHXaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dHXLaPBAsepPM1 : return SP2_HYBRID; break;
	case hbe_dCXAaPBAsepother : return SP2_HYBRID; break;
	case hbe_dIMDaPBAsepother : return SP2_HYBRID; break;
	case hbe_dIMEaPBAsepother : return SP2_HYBRID; break;
	case hbe_dINDaPBAsepother : return SP2_HYBRID; break;
	case hbe_dAMOaPBAsepother : return SP2_HYBRID; break;
	case hbe_dGDEaPBAsepother : return SP2_HYBRID; break;
	case hbe_dGDHaPBAsepother : return SP2_HYBRID; break;
	case hbe_dAHXaPBAsepother : return SP2_HYBRID; break;
	case hbe_dHXLaPBAsepother : return SP2_HYBRID; break;
	case hbe_dH2OaPBA : return SP2_HYBRID; break;
	case hbe_dPBAaCXAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaCXAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaCXA : return SP2_HYBRID; break;
	case hbe_dIMDaCXA : return SP2_HYBRID; break;
	case hbe_dIMEaCXA : return SP2_HYBRID; break;
	case hbe_dINDaCXA : return SP2_HYBRID; break;
	case hbe_dAMOaCXA : return SP2_HYBRID; break;
	case hbe_dGDEaCXA : return SP2_HYBRID; break;
	case hbe_dGDHaCXA : return SP2_HYBRID; break;
	case hbe_dAHXaCXA : return SP2_HYBRID; break;
	case hbe_dHXLaCXA : return SP2_HYBRID; break;
	case hbe_dH2OaCXA : return SP2_HYBRID; break;
	case hbe_dPBAaCXLsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaCXLsepother : return SP2_HYBRID; break;
	case hbe_dCXAaCXL : return SP2_HYBRID; break;
	case hbe_dIMDaCXL : return SP2_HYBRID; break;
	case hbe_dIMEaCXL : return SP2_HYBRID; break;
	case hbe_dINDaCXL : return SP2_HYBRID; break;
	case hbe_dAMOaCXL : return SP2_HYBRID; break;
	case hbe_dGDEaCXL : return SP2_HYBRID; break;
	case hbe_dGDHaCXL : return SP2_HYBRID; break;
	case hbe_dAHXaCXL : return SP2_HYBRID; break;
	case hbe_dHXLaCXL : return SP2_HYBRID; break;
	case hbe_dH2OaCXL : return SP2_HYBRID; break;
	case hbe_dPBAaIMDsepPM1 : return RING_HYBRID; break;
	case hbe_dPBAaIMDsepother : return RING_HYBRID; break;
	case hbe_dCXAaIMD : return RING_HYBRID; break;
	case hbe_dIMDaIMD : return RING_HYBRID; break;
	case hbe_dIMEaIMD : return RING_HYBRID; break;
	case hbe_dINDaIMD : return RING_HYBRID; break;
	case hbe_dAMOaIMD : return RING_HYBRID; break;
	case hbe_dGDEaIMD : return RING_HYBRID; break;
	case hbe_dGDHaIMD : return RING_HYBRID; break;
	case hbe_dAHXaIMD : return RING_HYBRID; break;
	case hbe_dHXLaIMD : return RING_HYBRID; break;
	case hbe_dH2OaIMD : return RING_HYBRID; break;
	case hbe_dPBAaIMEsepPM1 : return RING_HYBRID; break;
	case hbe_dPBAaIMEsepother : return RING_HYBRID; break;
	case hbe_dCXAaIME : return RING_HYBRID; break;
	case hbe_dIMDaIME : return RING_HYBRID; break;
	case hbe_dIMEaIME : return RING_HYBRID; break;
	case hbe_dINDaIME : return RING_HYBRID; break;
	case hbe_dAMOaIME : return RING_HYBRID; break;
	case hbe_dGDEaIME : return RING_HYBRID; break;
	case hbe_dGDHaIME : return RING_HYBRID; break;
	case hbe_dAHXaIME : return RING_HYBRID; break;
	case hbe_dHXLaIME : return RING_HYBRID; break;
	case hbe_dH2OaIME : return RING_HYBRID; break;
	case hbe_dPBAaAHXsepPM1 : return SP3_HYBRID; break;
	case hbe_dPBAaAHXsepother : return SP3_HYBRID; break;
	case hbe_dCXAaAHX : return SP3_HYBRID; break;
	case hbe_dIMDaAHX : return SP3_HYBRID; break;
	case hbe_dIMEaAHX : return SP3_HYBRID; break;
	case hbe_dINDaAHX : return SP3_HYBRID; break;
	case hbe_dAMOaAHX : return SP3_HYBRID; break;
	case hbe_dGDEaAHX : return SP3_HYBRID; break;
	case hbe_dGDHaAHX : return SP3_HYBRID; break;
	case hbe_dAHXaAHX : return SP3_HYBRID; break;
	case hbe_dHXLaAHX : return SP3_HYBRID; break;
	case hbe_dH2OaAHX : return SP3_HYBRID; break;
	case hbe_dPBAaHXLsepPM1 : return SP3_HYBRID; break;
	case hbe_dPBAaHXLsepother : return SP3_HYBRID; break;
	case hbe_dCXAaHXL : return SP3_HYBRID; break;
	case hbe_dIMDaHXL : return SP3_HYBRID; break;
	case hbe_dIMEaHXL : return SP3_HYBRID; break;
	case hbe_dINDaHXL : return SP3_HYBRID; break;
	case hbe_dAMOaHXL : return SP3_HYBRID; break;
	case hbe_dGDEaHXL : return SP3_HYBRID; break;
	case hbe_dGDHaHXL : return SP3_HYBRID; break;
	case hbe_dAHXaHXL : return SP3_HYBRID; break;
	case hbe_dHXLaHXL : return SP3_HYBRID; break;
	case hbe_dH2OaHXL : return SP3_HYBRID; break;
	case hbe_dPBAaPCA_DNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaPCA_DNAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dIMDaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dIMEaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dINDaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dAMOaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dGDEaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dGDHaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dAHXaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dHXLaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dH2OaPCA_DNA : return SP2_HYBRID; break;
	case hbe_dPBAaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dCXAaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dIMDaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMDaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dIMEaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMEaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dINDaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dINDaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dAMOaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAMOaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dGDEaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDEaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dGDHaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDHaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dAHXaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAHXaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dHXLaPCA_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dHXLaPCA_RNAsepother : return SP2_HYBRID; break;
	case hbe_dH2OaPCA_RNA : return SP2_HYBRID; break;
	case hbe_dPBAaPES_DNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaPES_DNAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaPES_DNA : return SP2_HYBRID; break;
	case hbe_dIMDaPES_DNA : return SP2_HYBRID; break;
	case hbe_dIMEaPES_DNA : return SP2_HYBRID; break;
	case hbe_dINDaPES_DNA : return SP2_HYBRID; break;
	case hbe_dAMOaPES_DNA : return SP2_HYBRID; break;
	case hbe_dGDEaPES_DNA : return SP2_HYBRID; break;
	case hbe_dGDHaPES_DNA : return SP2_HYBRID; break;
	case hbe_dAHXaPES_DNA : return SP2_HYBRID; break;
	case hbe_dHXLaPES_DNA : return SP2_HYBRID; break;
	case hbe_dH2OaPES_DNA : return SP2_HYBRID; break;
	case hbe_dPBAaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dPBAaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dCXAaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dCXAaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dIMDaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMDaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dIMEaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dIMEaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dINDaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dINDaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dAMOaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAMOaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dGDEaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDEaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dGDHaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dGDHaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dAHXaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dAHXaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dHXLaPES_RNAsepPM1 : return SP2_HYBRID; break;
	case hbe_dHXLaPES_RNAsepother : return SP2_HYBRID; break;
	case hbe_dH2OaPES_RNA : return SP2_HYBRID; break;
	case hbe_dPBAaRRI_DNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dPBAaRRI_DNAsepother : return SP3_HYBRID; break;
	case hbe_dCXAaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dIMDaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dIMEaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dINDaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dAMOaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dGDEaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dGDHaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dAHXaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dHXLaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dH2OaRRI_DNA : return SP3_HYBRID; break;
	case hbe_dPBAaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dPBAaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dCXAaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dCXAaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dIMDaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dIMDaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dIMEaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dIMEaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dINDaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dINDaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dAMOaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dAMOaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dGDEaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dGDEaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dGDHaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dGDHaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dAHXaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dAHXaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dHXLaRRI_RNAsepPM1 : return SP3_HYBRID; break;
	case hbe_dHXLaRRI_RNAsepother : return SP3_HYBRID; break;
	case hbe_dH2OaRRI_RNA : return SP3_HYBRID; break;
	case hbe_dPBAaH2O : return SP3_HYBRID; break;
	case hbe_dCXAaH2O : return SP3_HYBRID; break;
	case hbe_dIMDaH2O : return SP3_HYBRID; break;
	case hbe_dIMEaH2O : return SP3_HYBRID; break;
	case hbe_dINDaH2O : return SP3_HYBRID; break;
	case hbe_dAMOaH2O : return SP3_HYBRID; break;
	case hbe_dGDEaH2O : return SP3_HYBRID; break;
	case hbe_dGDHaH2O : return SP3_HYBRID; break;
	case hbe_dAHXaH2O : return SP3_HYBRID; break;
	case hbe_dHXLaH2O : return SP3_HYBRID; break;
	case hbe_dH2OaH2O : return SP3_HYBRID; break;
	case hbe_GENERIC_SP2BB_SR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP2BB_LR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP3BB_SR : return SP3_HYBRID; break;
	case hbe_GENERIC_SP3BB_LR : return SP3_HYBRID; break;
	case hbe_GENERIC_RINGBB_SR : return RING_HYBRID; break;
	case hbe_GENERIC_RINGBB_LR : return RING_HYBRID; break;
	case hbe_GENERIC_SP2BSC_SR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP2BSC_LR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP3BSC_SR : return SP3_HYBRID; break;
	case hbe_GENERIC_SP3BSC_LR : return SP3_HYBRID; break;
	case hbe_GENERIC_RINGBSC_SR : return RING_HYBRID; break;
	case hbe_GENERIC_RINGBSC_LR : return RING_HYBRID; break;
	case hbe_GENERIC_SP2SCSC_SR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP2SCSC_LR : return SP2_HYBRID; break;
	case hbe_GENERIC_SP3SCSC_SR : return SP3_HYBRID; break;
	case hbe_GENERIC_SP3SCSC_LR : return SP3_HYBRID; break;
	case hbe_GENERIC_RINGSCSC_SR : return RING_HYBRID; break;
	case hbe_GENERIC_RINGSCSC_LR : return RING_HYBRID; break;
	case hbe_UNKNOWN :
		TR.Warning << "WARNING: Attempted to find hybrid status of Unknown HBEvalType" << std::endl;
		return UNKNOWN_HYBRID; break;
	}
	return UNKNOWN_HYBRID;
}


///////////////////////////////////////////////////////////
//:::DEPRICATION NOTICE::::
//
// Note backbone/sidechain identification is used
//1) packing: eg since the backbone is fixed, backbone-backbone, and
// backbone-sidechain interactions can be precomputed. This is found
// by looking at rsd.is_backbone()
//2) environmental dependence is only applied to side chains, (ask
//Tanja or Lin, because I don't know why)
//3) to take some account of variable correlation,
//sidechain-sidechain hbonds have two different angle potentials
//that depend on the AH distance so hbe_is_BB_type is used.
//
//
//Backbone/sidechain distinctions are tricky when one looks at
//modified-proteins or non-proteins.  Since these functions are only
//used in questional circumstances they are depricated and will be
//removed in subsiquent versions of the hbond-potential.
//

bool
hbe_is_BB_type( HBEvalType hbe )
{
	switch(hbe){
	case hbe_NONE : return false; break;
	case hbe_dPBAaPBAsepM4helix : return true; break;
	case hbe_dPBAaPBAsepM3turn : return true; break;
	case hbe_dPBAaPBAsepM2turn : return true; break;
	case hbe_dPBAaPBAsepPM1 : return true; break;
	case hbe_dPBAaPBAsepP2turn : return true; break;
	case hbe_dPBAaPBAsepP3turn : return true; break;
	case hbe_dPBAaPBAsepP4helix : return true; break;
	case hbe_dPBAaPBAsepother : return true; break;
	case hbe_dCXAaPBAsepPM1 : return false; break;
	case hbe_dIMDaPBAsepPM1 : return false; break;
	case hbe_dIMEaPBAsepPM1 : return false; break;
	case hbe_dINDaPBAsepPM1 : return false; break;
	case hbe_dAMOaPBAsepPM1 : return false; break;
	case hbe_dGDEaPBAsepPM1 : return false; break;
	case hbe_dGDHaPBAsepPM1 : return false; break;
	case hbe_dAHXaPBAsepPM1 : return false; break;
	case hbe_dHXLaPBAsepPM1 : return false; break;
	case hbe_dCXAaPBAsepother : return false; break;
	case hbe_dIMDaPBAsepother : return false; break;
	case hbe_dIMEaPBAsepother : return false; break;
	case hbe_dINDaPBAsepother : return false; break;
	case hbe_dAMOaPBAsepother : return false; break;
	case hbe_dGDEaPBAsepother : return false; break;
	case hbe_dGDHaPBAsepother : return false; break;
	case hbe_dAHXaPBAsepother : return false; break;
	case hbe_dHXLaPBAsepother : return false; break;
	case hbe_dH2OaPBA : return false; break;
	case hbe_dPBAaCXAsepPM1 : return false; break;
	case hbe_dPBAaCXAsepother : return false; break;
	case hbe_dCXAaCXA : return false; break;
	case hbe_dIMDaCXA : return false; break;
	case hbe_dIMEaCXA : return false; break;
	case hbe_dINDaCXA : return false; break;
	case hbe_dAMOaCXA : return false; break;
	case hbe_dGDEaCXA : return false; break;
	case hbe_dGDHaCXA : return false; break;
	case hbe_dAHXaCXA : return false; break;
	case hbe_dHXLaCXA : return false; break;
	case hbe_dH2OaCXA : return false; break;
	case hbe_dPBAaCXLsepPM1 : return false; break;
	case hbe_dPBAaCXLsepother : return false; break;
	case hbe_dCXAaCXL : return false; break;
	case hbe_dIMDaCXL : return false; break;
	case hbe_dIMEaCXL : return false; break;
	case hbe_dINDaCXL : return false; break;
	case hbe_dAMOaCXL : return false; break;
	case hbe_dGDEaCXL : return false; break;
	case hbe_dGDHaCXL : return false; break;
	case hbe_dAHXaCXL : return false; break;
	case hbe_dHXLaCXL : return false; break;
	case hbe_dH2OaCXL : return false; break;
	case hbe_dPBAaIMDsepPM1 : return false; break;
	case hbe_dPBAaIMDsepother : return false; break;
	case hbe_dCXAaIMD : return false; break;
	case hbe_dIMDaIMD : return false; break;
	case hbe_dIMEaIMD : return false; break;
	case hbe_dINDaIMD : return false; break;
	case hbe_dAMOaIMD : return false; break;
	case hbe_dGDEaIMD : return false; break;
	case hbe_dGDHaIMD : return false; break;
	case hbe_dAHXaIMD : return false; break;
	case hbe_dHXLaIMD : return false; break;
	case hbe_dH2OaIMD : return false; break;
	case hbe_dPBAaIMEsepPM1 : return false; break;
	case hbe_dPBAaIMEsepother : return false; break;
	case hbe_dCXAaIME : return false; break;
	case hbe_dIMDaIME : return false; break;
	case hbe_dIMEaIME : return false; break;
	case hbe_dINDaIME : return false; break;
	case hbe_dAMOaIME : return false; break;
	case hbe_dGDEaIME : return false; break;
	case hbe_dGDHaIME : return false; break;
	case hbe_dAHXaIME : return false; break;
	case hbe_dHXLaIME : return false; break;
	case hbe_dH2OaIME : return false; break;
	case hbe_dPBAaAHXsepPM1 : return false; break;
	case hbe_dPBAaAHXsepother : return false; break;
	case hbe_dCXAaAHX : return false; break;
	case hbe_dIMDaAHX : return false; break;
	case hbe_dIMEaAHX : return false; break;
	case hbe_dINDaAHX : return false; break;
	case hbe_dAMOaAHX : return false; break;
	case hbe_dGDEaAHX : return false; break;
	case hbe_dGDHaAHX : return false; break;
	case hbe_dAHXaAHX : return false; break;
	case hbe_dHXLaAHX : return false; break;
	case hbe_dH2OaAHX : return false; break;
	case hbe_dPBAaHXLsepPM1 : return false; break;
	case hbe_dPBAaHXLsepother : return false; break;
	case hbe_dCXAaHXL : return false; break;
	case hbe_dIMDaHXL : return false; break;
	case hbe_dIMEaHXL : return false; break;
	case hbe_dINDaHXL : return false; break;
	case hbe_dAMOaHXL : return false; break;
	case hbe_dGDEaHXL : return false; break;
	case hbe_dGDHaHXL : return false; break;
	case hbe_dAHXaHXL : return false; break;
	case hbe_dHXLaHXL : return false; break;
	case hbe_dH2OaHXL : return false; break;
	case hbe_dPBAaPCA_DNAsepPM1 : return false; break;
	case hbe_dPBAaPCA_DNAsepother : return false; break;
	case hbe_dCXAaPCA_DNA : return false; break;
	case hbe_dIMDaPCA_DNA : return false; break;
	case hbe_dIMEaPCA_DNA : return false; break;
	case hbe_dINDaPCA_DNA : return false; break;
	case hbe_dAMOaPCA_DNA : return false; break;
	case hbe_dGDEaPCA_DNA : return false; break;
	case hbe_dGDHaPCA_DNA : return false; break;
	case hbe_dAHXaPCA_DNA : return false; break;
	case hbe_dHXLaPCA_DNA : return false; break;
	case hbe_dH2OaPCA_DNA : return false; break;
	case hbe_dPBAaPCA_RNAsepPM1 : return true; break;
	case hbe_dPBAaPCA_RNAsepother : return true; break;
	case hbe_dCXAaPCA_RNAsepPM1 : return false; break;
	case hbe_dCXAaPCA_RNAsepother : return false; break;
	case hbe_dIMDaPCA_RNAsepPM1 : return false; break;
	case hbe_dIMDaPCA_RNAsepother : return false; break;
	case hbe_dIMEaPCA_RNAsepPM1 : return false; break;
	case hbe_dIMEaPCA_RNAsepother : return false; break;
	case hbe_dINDaPCA_RNAsepPM1 : return false; break;
	case hbe_dINDaPCA_RNAsepother : return false; break;
	case hbe_dAMOaPCA_RNAsepPM1 : return false; break;
	case hbe_dAMOaPCA_RNAsepother : return false; break;
	case hbe_dGDEaPCA_RNAsepPM1 : return false; break;
	case hbe_dGDEaPCA_RNAsepother : return false; break;
	case hbe_dGDHaPCA_RNAsepPM1 : return false; break;
	case hbe_dGDHaPCA_RNAsepother : return false; break;
	case hbe_dAHXaPCA_RNAsepPM1 : return false; break;
	case hbe_dAHXaPCA_RNAsepother : return false; break;
	case hbe_dHXLaPCA_RNAsepPM1 : return false; break;
	case hbe_dHXLaPCA_RNAsepother : return false; break;
	case hbe_dH2OaPCA_RNA : return false; break;
	case hbe_dPBAaPES_DNAsepPM1 : return false; break;
	case hbe_dPBAaPES_DNAsepother : return false; break;
	case hbe_dCXAaPES_DNA : return false; break;
	case hbe_dIMDaPES_DNA : return false; break;
	case hbe_dIMEaPES_DNA : return false; break;
	case hbe_dINDaPES_DNA : return false; break;
	case hbe_dAMOaPES_DNA : return false; break;
	case hbe_dGDEaPES_DNA : return false; break;
	case hbe_dGDHaPES_DNA : return false; break;
	case hbe_dAHXaPES_DNA : return false; break;
	case hbe_dHXLaPES_DNA : return false; break;
	case hbe_dH2OaPES_DNA : return false; break;
	case hbe_dPBAaPES_RNAsepPM1 : return true; break;
	case hbe_dPBAaPES_RNAsepother : return true; break;
	case hbe_dCXAaPES_RNAsepPM1 : return false; break;
	case hbe_dCXAaPES_RNAsepother : return false; break;
	case hbe_dIMDaPES_RNAsepPM1 : return false; break;
	case hbe_dIMDaPES_RNAsepother : return false; break;
	case hbe_dIMEaPES_RNAsepPM1 : return false; break;
	case hbe_dIMEaPES_RNAsepother : return false; break;
	case hbe_dINDaPES_RNAsepPM1 : return false; break;
	case hbe_dINDaPES_RNAsepother : return false; break;
	case hbe_dAMOaPES_RNAsepPM1 : return false; break;
	case hbe_dAMOaPES_RNAsepother : return false; break;
	case hbe_dGDEaPES_RNAsepPM1 : return false; break;
	case hbe_dGDEaPES_RNAsepother : return false; break;
	case hbe_dGDHaPES_RNAsepPM1 : return false; break;
	case hbe_dGDHaPES_RNAsepother : return false; break;
	case hbe_dAHXaPES_RNAsepPM1 : return false; break;
	case hbe_dAHXaPES_RNAsepother : return false; break;
	case hbe_dHXLaPES_RNAsepPM1 : return false; break;
	case hbe_dHXLaPES_RNAsepother : return false; break;
	case hbe_dH2OaPES_RNA : return false; break;
	case hbe_dPBAaRRI_DNAsepPM1 : return false; break;
	case hbe_dPBAaRRI_DNAsepother : return false; break;
	case hbe_dCXAaRRI_DNA : return false; break;
	case hbe_dIMDaRRI_DNA : return false; break;
	case hbe_dIMEaRRI_DNA : return false; break;
	case hbe_dINDaRRI_DNA : return false; break;
	case hbe_dAMOaRRI_DNA : return false; break;
	case hbe_dGDEaRRI_DNA : return false; break;
	case hbe_dGDHaRRI_DNA : return false; break;
	case hbe_dAHXaRRI_DNA : return false; break;
	case hbe_dHXLaRRI_DNA : return false; break;
	case hbe_dH2OaRRI_DNA : return false; break;
	case hbe_dPBAaRRI_RNAsepPM1 : return false; break;
	case hbe_dPBAaRRI_RNAsepother : return false; break;
	case hbe_dCXAaRRI_RNAsepPM1 : return false; break;
	case hbe_dCXAaRRI_RNAsepother : return false; break;
	case hbe_dIMDaRRI_RNAsepPM1 : return false; break;
	case hbe_dIMDaRRI_RNAsepother : return false; break;
	case hbe_dIMEaRRI_RNAsepPM1 : return false; break;
	case hbe_dIMEaRRI_RNAsepother : return false; break;
	case hbe_dINDaRRI_RNAsepPM1 : return false; break;
	case hbe_dINDaRRI_RNAsepother : return false; break;
	case hbe_dAMOaRRI_RNAsepPM1 : return false; break;
	case hbe_dAMOaRRI_RNAsepother : return false; break;
	case hbe_dGDEaRRI_RNAsepPM1 : return false; break;
	case hbe_dGDEaRRI_RNAsepother : return false; break;
	case hbe_dGDHaRRI_RNAsepPM1 : return false; break;
	case hbe_dGDHaRRI_RNAsepother : return false; break;
	case hbe_dAHXaRRI_RNAsepPM1 : return false; break;
	case hbe_dAHXaRRI_RNAsepother : return false; break;
	case hbe_dHXLaRRI_RNAsepPM1 : return false; break;
	case hbe_dHXLaRRI_RNAsepother : return false; break;
	case hbe_dH2OaRRI_RNA : return false; break;
	case hbe_dPBAaH2O : return false; break;
	case hbe_dCXAaH2O : return false; break;
	case hbe_dIMDaH2O : return false; break;
	case hbe_dIMEaH2O : return false; break;
	case hbe_dINDaH2O : return false; break;
	case hbe_dAMOaH2O : return false; break;
	case hbe_dGDEaH2O : return false; break;
	case hbe_dGDHaH2O : return false; break;
	case hbe_dAHXaH2O : return false; break;
	case hbe_dHXLaH2O : return false; break;
	case hbe_dH2OaH2O : return false; break;
	case hbe_GENERIC_SP2BB_SR : return true; break;
	case hbe_GENERIC_SP2BB_LR : return true; break;
	case hbe_GENERIC_SP3BB_SR : return true; break;
	case hbe_GENERIC_SP3BB_LR : return true; break;
	case hbe_GENERIC_RINGBB_SR : return true; break;
	case hbe_GENERIC_RINGBB_LR : return true; break;
	case hbe_GENERIC_SP2BSC_SR : return false; break;
	case hbe_GENERIC_SP2BSC_LR : return false; break;
	case hbe_GENERIC_SP3BSC_SR : return false; break;
	case hbe_GENERIC_SP3BSC_LR : return false; break;
	case hbe_GENERIC_RINGBSC_SR : return false; break;
	case hbe_GENERIC_RINGBSC_LR : return false; break;
	case hbe_GENERIC_SP2SCSC_SR : return false; break;
	case hbe_GENERIC_SP2SCSC_LR : return false; break;
	case hbe_GENERIC_SP3SCSC_SR : return false; break;
	case hbe_GENERIC_SP3SCSC_LR : return false; break;
	case hbe_GENERIC_RINGSCSC_SR : return false; break;
	case hbe_GENERIC_RINGSCSC_LR : return false; break;
	case hbe_UNKNOWN :
		TR.Error << "ERROR: Attempted to find backbone status of Unknown HBEvalType" << std::endl;
		break;
	}
	utility_exit_with_message("Unhandled HBEvalType");
	return false; // <- to assure gcc a bool will be returned
}
bool hbe_is_SC_type( HBEvalType hbe ){return !hbe_is_BB_type(hbe);}


HBondWeightType
get_hbond_weight_type( HBEvalType const & hbe_type )
{
	switch(hbe_type){
	case hbe_NONE : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepM4helix : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepM3turn : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepM2turn : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepPM1 : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepP2turn : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepP3turn : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepP4helix : return hbw_SR_BB; break;
	case hbe_dPBAaPBAsepother : return hbw_LR_BB; break;
	case hbe_dCXAaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMDaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMEaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dINDaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAMOaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDEaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDHaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAHXaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dHXLaPBAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dCXAaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMDaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMEaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dINDaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAMOaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDEaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDHaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAHXaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dHXLaPBAsepother : return hbw_LR_BB_SC; break;
	case hbe_dH2OaPBA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaCXAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaCXAsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaCXA : return hbw_SC; break;
	case hbe_dIMDaCXA : return hbw_SC; break;
	case hbe_dIMEaCXA : return hbw_SC; break;
	case hbe_dINDaCXA : return hbw_SC; break;
	case hbe_dAMOaCXA : return hbw_SC; break;
	case hbe_dGDEaCXA : return hbw_SC; break;
	case hbe_dGDHaCXA : return hbw_SC; break;
	case hbe_dAHXaCXA : return hbw_SC; break;
	case hbe_dHXLaCXA : return hbw_SC; break;
	case hbe_dH2OaCXA : return hbw_SC; break;
	case hbe_dPBAaCXLsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaCXLsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaCXL : return hbw_SC; break;
	case hbe_dIMDaCXL : return hbw_SC; break;
	case hbe_dIMEaCXL : return hbw_SC; break;
	case hbe_dINDaCXL : return hbw_SC; break;
	case hbe_dAMOaCXL : return hbw_SC; break;
	case hbe_dGDEaCXL : return hbw_SC; break;
	case hbe_dGDHaCXL : return hbw_SC; break;
	case hbe_dAHXaCXL : return hbw_SC; break;
	case hbe_dHXLaCXL : return hbw_SC; break;
	case hbe_dH2OaCXL : return hbw_SC; break;
	case hbe_dPBAaIMDsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaIMDsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaIMD : return hbw_SC; break;
	case hbe_dIMDaIMD : return hbw_SC; break;
	case hbe_dIMEaIMD : return hbw_SC; break;
	case hbe_dINDaIMD : return hbw_SC; break;
	case hbe_dAMOaIMD : return hbw_SC; break;
	case hbe_dGDEaIMD : return hbw_SC; break;
	case hbe_dGDHaIMD : return hbw_SC; break;
	case hbe_dAHXaIMD : return hbw_SC; break;
	case hbe_dHXLaIMD : return hbw_SC; break;
	case hbe_dH2OaIMD : return hbw_SC; break;
	case hbe_dPBAaIMEsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaIMEsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaIME : return hbw_SC; break;
	case hbe_dIMDaIME : return hbw_SC; break;
	case hbe_dIMEaIME : return hbw_SC; break;
	case hbe_dINDaIME : return hbw_SC; break;
	case hbe_dAMOaIME : return hbw_SC; break;
	case hbe_dGDEaIME : return hbw_SC; break;
	case hbe_dGDHaIME : return hbw_SC; break;
	case hbe_dAHXaIME : return hbw_SC; break;
	case hbe_dHXLaIME : return hbw_SC; break;
	case hbe_dH2OaIME : return hbw_SC; break;
	case hbe_dPBAaAHXsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaAHXsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaAHX : return hbw_SC; break;
	case hbe_dIMDaAHX : return hbw_SC; break;
	case hbe_dIMEaAHX : return hbw_SC; break;
	case hbe_dINDaAHX : return hbw_SC; break;
	case hbe_dAMOaAHX : return hbw_SC; break;
	case hbe_dGDEaAHX : return hbw_SC; break;
	case hbe_dGDHaAHX : return hbw_SC; break;
	case hbe_dAHXaAHX : return hbw_SC; break;
	case hbe_dHXLaAHX : return hbw_SC; break;
	case hbe_dH2OaAHX : return hbw_SC; break;
	case hbe_dPBAaHXLsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaHXLsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaHXL : return hbw_SC; break;
	case hbe_dIMDaHXL : return hbw_SC; break;
	case hbe_dIMEaHXL : return hbw_SC; break;
	case hbe_dINDaHXL : return hbw_SC; break;
	case hbe_dAMOaHXL : return hbw_SC; break;
	case hbe_dGDEaHXL : return hbw_SC; break;
	case hbe_dGDHaHXL : return hbw_SC; break;
	case hbe_dAHXaHXL : return hbw_SC; break;
	case hbe_dHXLaHXL : return hbw_SC; break;
	case hbe_dH2OaHXL : return hbw_SC; break;
	case hbe_dPBAaPCA_DNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaPCA_DNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaPCA_DNA : return hbw_SC; break;
	case hbe_dIMDaPCA_DNA : return hbw_SC; break;
	case hbe_dIMEaPCA_DNA : return hbw_SC; break;
	case hbe_dINDaPCA_DNA : return hbw_SC; break;
	case hbe_dAMOaPCA_DNA : return hbw_SC; break;
	case hbe_dGDEaPCA_DNA : return hbw_SC; break;
	case hbe_dGDHaPCA_DNA : return hbw_SC; break;
	case hbe_dAHXaPCA_DNA : return hbw_SC; break;
	case hbe_dHXLaPCA_DNA : return hbw_SC; break;
	case hbe_dH2OaPCA_DNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaPCA_RNAsepPM1 : return hbw_SR_BB; break;
	case hbe_dPBAaPCA_RNAsepother : return hbw_LR_BB; break;
	case hbe_dCXAaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dCXAaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMDaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMDaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMEaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMEaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dINDaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dINDaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAMOaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAMOaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDEaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDEaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDHaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDHaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAHXaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAHXaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dHXLaPCA_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dHXLaPCA_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dH2OaPCA_RNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaPES_DNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaPES_DNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaPES_DNA : return hbw_SC; break;
	case hbe_dIMDaPES_DNA : return hbw_SC; break;
	case hbe_dIMEaPES_DNA : return hbw_SC; break;
	case hbe_dINDaPES_DNA : return hbw_SC; break;
	case hbe_dAMOaPES_DNA : return hbw_SC; break;
	case hbe_dGDEaPES_DNA : return hbw_SC; break;
	case hbe_dGDHaPES_DNA : return hbw_SC; break;
	case hbe_dAHXaPES_DNA : return hbw_SC; break;
	case hbe_dHXLaPES_DNA : return hbw_SC; break;
	case hbe_dH2OaPES_DNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaPES_RNAsepPM1 : return hbw_SR_BB; break;
	case hbe_dPBAaPES_RNAsepother : return hbw_LR_BB; break;
	case hbe_dCXAaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dCXAaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMDaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMDaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMEaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMEaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dINDaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dINDaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAMOaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAMOaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDEaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDEaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDHaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDHaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAHXaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAHXaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dHXLaPES_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dHXLaPES_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dH2OaPES_RNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaRRI_DNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dPBAaRRI_DNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dCXAaRRI_DNA : return hbw_SC; break;
	case hbe_dIMDaRRI_DNA : return hbw_SC; break;
	case hbe_dIMEaRRI_DNA : return hbw_SC; break;
	case hbe_dINDaRRI_DNA : return hbw_SC; break;
	case hbe_dAMOaRRI_DNA : return hbw_SC; break;
	case hbe_dGDEaRRI_DNA : return hbw_SC; break;
	case hbe_dGDHaRRI_DNA : return hbw_SC; break;
	case hbe_dAHXaRRI_DNA : return hbw_SC; break;
	case hbe_dHXLaRRI_DNA : return hbw_SC; break;
	case hbe_dH2OaRRI_DNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaRRI_RNAsepPM1 : return hbw_SR_BB; break;
	case hbe_dPBAaRRI_RNAsepother : return hbw_LR_BB; break;
	case hbe_dCXAaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dCXAaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMDaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMDaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dIMEaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dIMEaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dINDaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dINDaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAMOaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAMOaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDEaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDEaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dGDHaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dGDHaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dAHXaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dAHXaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dHXLaRRI_RNAsepPM1 : return hbw_SR_BB_SC; break;
	case hbe_dHXLaRRI_RNAsepother : return hbw_LR_BB_SC; break;
	case hbe_dH2OaRRI_RNA : return hbw_LR_BB_SC; break;
	case hbe_dPBAaH2O : return hbw_LR_BB_SC; break;
	case hbe_dCXAaH2O : return hbw_SC; break;
	case hbe_dIMDaH2O : return hbw_SC; break;
	case hbe_dIMEaH2O : return hbw_SC; break;
	case hbe_dINDaH2O : return hbw_SC; break;
	case hbe_dAMOaH2O : return hbw_SC; break;
	case hbe_dGDEaH2O : return hbw_SC; break;
	case hbe_dGDHaH2O : return hbw_SC; break;
	case hbe_dAHXaH2O : return hbw_SC; break;
	case hbe_dHXLaH2O : return hbw_SC; break;
	case hbe_dH2OaH2O : return hbw_SC; break;
	case hbe_GENERIC_SP2BB_SR : return hbw_SR_BB; break;
	case hbe_GENERIC_SP2BB_LR : return hbw_LR_BB; break;
	case hbe_GENERIC_SP3BB_SR : return hbw_SR_BB; break;
	case hbe_GENERIC_SP3BB_LR : return hbw_LR_BB; break;
	case hbe_GENERIC_RINGBB_SR : return hbw_SR_BB; break;
	case hbe_GENERIC_RINGBB_LR : return hbw_LR_BB; break;
	case hbe_GENERIC_SP2BSC_SR : return hbw_SR_BB_SC; break;
	case hbe_GENERIC_SP2BSC_LR : return hbw_LR_BB_SC; break;
	case hbe_GENERIC_SP3BSC_SR : return hbw_SR_BB_SC; break;
	case hbe_GENERIC_SP3BSC_LR : return hbw_LR_BB_SC; break;
	case hbe_GENERIC_RINGBSC_SR : return hbw_SR_BB_SC; break;
	case hbe_GENERIC_RINGBSC_LR : return hbw_LR_BB_SC; break;
	case hbe_GENERIC_SP2SCSC_SR : return hbw_SC; break;
	case hbe_GENERIC_SP2SCSC_LR : return hbw_SC; break;
	case hbe_GENERIC_SP3SCSC_SR : return hbw_SC; break;
	case hbe_GENERIC_SP3SCSC_LR : return hbw_SC; break;
	case hbe_GENERIC_RINGSCSC_SR : return hbw_SC; break;
	case hbe_GENERIC_RINGSCSC_LR : return hbw_SC; break;
	case hbe_UNKNOWN :
		TR.Error << "ERROR: Attempted to find weight type Unknown HBEvalType" << std::endl;
		break;
	}
	utility_exit_with_message("Unhandled HBEvalType");
	return hbw_NONE; // <- to assure gcc an HBondWeightType is returned
}

} // hbonds
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::HBondDerivs::save( Archive & arc ) const {
	arc( CEREAL_NVP( don_deriv ) ); // class core::scoring::DerivVectorPair
	arc( CEREAL_NVP( h_deriv ) ); // class core::scoring::DerivVectorPair
	arc( CEREAL_NVP( acc_deriv ) ); // class core::scoring::DerivVectorPair
	arc( CEREAL_NVP( abase_deriv ) ); // class core::scoring::DerivVectorPair
	arc( CEREAL_NVP( abase2_deriv ) ); // class core::scoring::DerivVectorPair
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::HBondDerivs::load( Archive & arc ) {
	arc( don_deriv ); // class core::scoring::DerivVectorPair
	arc( h_deriv ); // class core::scoring::DerivVectorPair
	arc( acc_deriv ); // class core::scoring::DerivVectorPair
	arc( abase_deriv ); // class core::scoring::DerivVectorPair
	arc( abase2_deriv ); // class core::scoring::DerivVectorPair
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::HBondDerivs );
#endif // SERIALIZATION
