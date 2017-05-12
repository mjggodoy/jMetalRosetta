// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file    core/chemical/VariantType.hh
/// @brief   Enumeration definition for VariantType.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    DO NOT EDIT THIS FILE DIRECTLY!  It is auto-generated.
/// If you wish to edit it, modify the update_ResidueType_enum_files.py script.

#ifndef INCLUDED_core_chemical_VariantType_HH
#define INCLUDED_core_chemical_VariantType_HH

namespace core {
namespace chemical {

// If adding a new variant, DO NOT MANUALLY EDIT THIS ENUMERATOR DEFINITION.
// Instead, add it to variant_types.list, and run the update_ResidueType_enum_files.py script.
/// @brief   Enumerators for all the ResidueType variants.
/// @details VariantTypes are primarily utilized by the patch system.  All the type does is add an identifier
/// that can be used later on in different protocols.  It also helps the patch system keep track of which residues are
/// patched with which patches.
enum VariantType {
	NO_VARIANT = 0,
	FIRST_VARIANT = 1,
	UPPER_TERMINUS_VARIANT = 1,  // C-terminus cap
	LOWER_TERMINUS_VARIANT,  // N-terminus cap
	UPPERTERM_TRUNC_VARIANT,  // C-terminus truncation
	LOWERTERM_TRUNC_VARIANT,  // N-terminus truncation
	NTERM_CONNECT,  // for cyclic peptides
	CTERM_CONNECT,  // for cyclic peptides
	BRANCH_LOWER_TERMINUS_VARIANT,  // used for branched polymers and glycosylations
	CUTPOINT_LOWER,  // for use during loop modeling, at positions before a cutpoint
	CUTPOINT_UPPER,  // for use during loop modeling, at positions after a cutpoint
	SC_BRANCH_POINT,  // for branched polymers and glycosylations from amino acid side chains
	C1_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C2_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C3_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C4_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C5_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C6_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C7_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C8_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	C9_BRANCH_POINT,  // for branched polymers and glycosylations from position 1 of saccharide residues
	SIDECHAIN_CONJUGATION,  // for chemically conjugatable residues and side-chain conjugation (like ubiquitination)
	PROTONATED,
	DEPROTONATED,
	ALTERNATIVE_PROTONATION,
	SG_CONNECT,
	NE2_CONNECT,
	ZN_CONNECT,
	METHYLATION,
	METHYLATED_NTERM_VARIANT,  // N-terminal methylation
	PHOSPHORYLATION,  // used for PTM patches
	ACETYLATION,  // used for PTM patches
	SULFATION,  // used for PTM patches
	CARBOXYLATION,  // used for PTM patches
	HYDROXYLATION,  // used for PTM patches
	DIMETHYLATION,  // used for PTM patches
	TRIMETHYLATION,  // used for PTM patches
	DIIODINATION,  // used for PTM patches
	N_METHYLATION,  // Methylation on the peptide bond nitrogen
	ACETYLATED_NTERMINUS_VARIANT,  // for creating amino acid dipeptides for NCAA rotamer libraries
	METHYLATED_CTERMINUS_VARIANT,  // for creating amino acid dipeptides for NCAA rotamer libraries
	PHOSPHONATE_UPPER_VARIANT,  // for terminal amino phosphonic acids
	VIRTUAL_DNA_PHOSPHATE,
	VIRTUAL_PHOSPHATE,  // in wide use for RNA -- 'default' 5' ending for fragment assembly & stepwise
	REPL_PHOSPHATE,  // in wide use for RNA to represent steric but 'unstructured' atoms
	VIRTUAL_RIBOSE,  // in use for RNA during stepwise assembly
	VIRTUAL_BB,  // useful for docking nucleobases
	VIRTUAL_BACKBONE_EXCEPT_C1PRIME,  // not in use, may deprecated along with turner_rules_test
	VIRTUAL_BASE,  // not in wide use, may deprecated in future
	VIRTUAL_BASE_HEAVY_ATOM,  // not in use, may deprecated along with old swa_rna_main in 2015
	VIRTUAL_RNA_RESIDUE,  // not in use, may deprecated along with old swa_rna_main in 2015
	VIRTUAL_RNA_RESIDUE_EXCLUDE_PHOSPHATE,  // not in use, may deprecated along with old swa_rna_main in 2015
	BULGE,  // not in use, may deprecated along with old swa_rna_main in 2015
	VIRTUAL_O2PRIME_HYDROGEN,  // important for stepwise RNA code
	VIRTUAL_SIDE_CHAIN,  // important for stepwise protein code.
	THREE_PRIME_END_OH,  // alternative terminal variant in real RNAs, keep this.
	THREE_PRIME_PHOSPHATE,  // alternative terminal variant in real RNAs, keep this.
	FIVE_PRIME_END_OH,  // alternative terminal variant in real RNAs, keep this.
	FIVE_PRIME_END_PHOSPHATE,  // alternative terminal variant in real RNAs, keep this.
	FIVE_PRIME_PHOSPHATE,  // alternative terminal variant in real RNAs, keep this.
	DEOXY_O2PRIME,  // Make RNA into DNA
	UPPER_CONNECTION_RNA,
	THREE_PRIME_PACKABLE_PHOSPHATE,  // experimental terminal variant -- deprecate in 2015 if not in use by then.
	FIVE_PRIME_PACKABLE_PHOSPHATE,  // experimental terminal variant -- deprecate in 2015 if not in use by then.
	PROTONATED_N1_ADENOSINE,  // in use in RNA -- functionally important protonation state
	PROTONATED_N3_ADENOSINE,  // in use in RNA -- functionally important protonation state
	THREE_PRIME_FIVE_PRIME_METHYL_PHOSPHATE,  // capping for optimization
	BLOCK_STACK_ABOVE,  // Puts a steric block above a base to disallow stacking; for solving RNA motifs
	BLOCK_STACK_BELOW,  // Puts a steric block above a base to disallow stacking; for solving RNA motifs
	N_ACETYLATION,  // for stepwise assembly (SWA) code; different geometry/atoms then ACETYLATED_NTERMINUS above
	N_FORMYLATION,
	C_METHYLAMIDATION,  // for stepwise assembly (SWA) code; distinct from METHYLATED_CTERMINUS above
	CTERM_AMIDATION,
	ALDONIC_ACID_VARIANT,  // C1 by definition
	C2_KETOALDONIC_ACID,
	C3_KETOALDONIC_ACID,
	C4_KETOALDONIC_ACID,
	C5_KETOALDONIC_ACID,
	C6_KETOALDONIC_ACID,
	C7_KETOALDONIC_ACID,
	C8_KETOALDONIC_ACID,
	URONIC_ACID_VARIANT,  // Cn by definition
	C1_DEOXY_SUGAR,
	C2_DEOXY_SUGAR,
	C3_DEOXY_SUGAR,
	C4_DEOXY_SUGAR,
	C5_DEOXY_SUGAR,
	C6_DEOXY_SUGAR,
	C7_DEOXY_SUGAR,
	C8_DEOXY_SUGAR,
	C9_DEOXY_SUGAR,
	C1_AMINO_SUGAR,
	C2_AMINO_SUGAR,
	C3_AMINO_SUGAR,
	C4_AMINO_SUGAR,
	C5_AMINO_SUGAR,
	C6_AMINO_SUGAR,
	C7_AMINO_SUGAR,
	C8_AMINO_SUGAR,
	C9_AMINO_SUGAR,
	C1_ACETYLAMINO_SUGAR,
	C2_ACETYLAMINO_SUGAR,
	C3_ACETYLAMINO_SUGAR,
	C4_ACETYLAMINO_SUGAR,
	C5_ACETYLAMINO_SUGAR,
	C6_ACETYLAMINO_SUGAR,
	C7_ACETYLAMINO_SUGAR,
	C8_ACETYLAMINO_SUGAR,
	C9_ACETYLAMINO_SUGAR,
	C1_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C2_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C3_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C4_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C5_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C6_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C7_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C8_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C9_R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR,
	C1_SULFATED_SUGAR,
	C2_SULFATED_SUGAR,
	C3_SULFATED_SUGAR,
	C4_SULFATED_SUGAR,
	C5_SULFATED_SUGAR,
	C6_SULFATED_SUGAR,
	C7_SULFATED_SUGAR,
	C8_SULFATED_SUGAR,
	C9_SULFATED_SUGAR,
	C1_SULFOAMINO_SUGAR,
	C2_SULFOAMINO_SUGAR,
	C3_SULFOAMINO_SUGAR,
	C4_SULFOAMINO_SUGAR,
	C5_SULFOAMINO_SUGAR,
	C6_SULFOAMINO_SUGAR,
	C7_SULFOAMINO_SUGAR,
	C8_SULFOAMINO_SUGAR,
	C9_SULFOAMINO_SUGAR,
	C1_METHYLATED_SUGAR,
	C2_METHYLATED_SUGAR,
	C3_METHYLATED_SUGAR,
	C4_METHYLATED_SUGAR,
	C5_METHYLATED_SUGAR,
	C6_METHYLATED_SUGAR,
	C7_METHYLATED_SUGAR,
	C8_METHYLATED_SUGAR,
	C9_METHYLATED_SUGAR,
	C1_PHOSPHATE,
	C2_PHOSPHATE,
	C3_PHOSPHATE,
	C4_PHOSPHATE,
	C5_PHOSPHATE,
	C6_PHOSPHATE,
	C7_PHOSPHATE,
	C8_PHOSPHATE,
	C9_PHOSPHATE,
	METHYL_GLYCOSIDE,  // for methylated saccharide lower termini
	OOP_PRE,  // used for oligooxopiperazines (OOPs)
	OOP_POST,  // used for oligooxopiperazines (OOPs)
	HBS_PRE,  // used for hydrogen bond surrogates
	HBS_POST,  // used for hydrogen bond surrogates
	A3B_HBS_PRE,  // used for a3b hydrogen bond surrogates
	A3B_HBS_POST,  // used for a3b hydrogen bond surrogates
	TRIAZOLAMERN,
	TRIAZOLAMERC,
	DISULFIDE,
	ADDUCT_VARIANT,
	CENTROID_WITH_HA,  // used in NOE apps.
	SPECIAL_ROT,  // generic VariantType that allows for differential scoring of a set of residues/rotamers
	HYDROXYLATION1,  // used specifically for hydroxylated Pro
	HYDROXYLATION2,  // used specifically for hydroxylated Pro
	REPLONLY,  // Only the repulsive energy will be considered during structure calculations on this residue.
	REPLS_BB,
	SC_FRAGMENT,  // used by PlaceProbeMover.cc and hotspot hashing
	SHOVE_BB,  // used by MapHotSpot.cc, HotspotStubSet.cc, ShoveResidueMover.cc, and TryRotamers.cc
	VIRTUAL_RESIDUE_VARIANT,
	VIRTUAL_NTERM,
	N_VARIANTS = VIRTUAL_NTERM
};

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_VariantType_HH