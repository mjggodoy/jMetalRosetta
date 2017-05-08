// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_ScoringManager_hh
#define INCLUDED_core_scoring_ScoringManager_hh

// Unit headers
#include <core/scoring/ScoringManager.fwd.hh>

// Package headers
#include <core/scoring/AtomVDW.fwd.hh>
#include <core/scoring/CenRotEnvPairPotential.fwd.hh>
#include <core/scoring/CenHBPotential.fwd.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/GenBornPotential.fwd.hh>
#include <core/scoring/HydroxylTorsionPotential.fwd.hh>
#include <core/scoring/MultipoleElecPotential.fwd.hh>
#include <core/scoring/SASAPotential.fwd.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh> //pba
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/P_AA_ss.fwd.hh>
#include <core/scoring/PairEPotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <core/scoring/ProQPotential.fwd.hh>
#include <core/scoring/RamaPrePro.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/SecondaryStructurePotential.fwd.hh>
#include <core/scoring/SmoothEnvPairPotential.fwd.hh>
#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/VdWTinkerPotential.fwd.hh>
#include <core/scoring/WaterAdductHBondPotential.fwd.hh>

#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.fwd.hh>

#include <core/scoring/carbohydrates/CHIEnergyFunction.fwd.hh>
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.fwd.hh>

#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>

#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>

#include <core/scoring/elec/CPRepMapType.fwd.hh>

#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.fwd.hh>

#include <core/scoring/dna/DNA_BasePotential.fwd.hh>
#include <core/scoring/dna/DNA_DihedralPotential.fwd.hh>
#include <core/scoring/dna/DNA_EnvPairPotential.fwd.hh>
#include <core/scoring/dna/DNABFormPotential.fwd.hh>
#include <core/scoring/dna/DNATorsionPotential.fwd.hh>
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>

#include <core/scoring/facts/FACTSPotential.fwd.hh>

#include <core/scoring/geometric_solvation/DatabaseOccSolEne.fwd.hh>

#include <core/scoring/interface_/DDPlookup.fwd.hh>

#include <core/scoring/membrane/MembraneData.hh>

#include <core/scoring/mm/MMLJLibrary.fwd.hh>
#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondLengthLibrary.fwd.hh>

#include <core/scoring/nv/NVlookup.fwd.hh>

#include <core/scoring/orbitals/OrbitalsLookup.fwd.hh>

#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <core/scoring/rna/RNA_SuitePotential.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/rna/RNP_LowResPotential.fwd.hh>
#include <core/scoring/rna/RNP_LowResStackData.fwd.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.fwd.hh>

#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.fwd.hh>

#if defined(WIN32) || defined(WIN_PYROSETTA)
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh> // WIN32 INCLUDE
#include <core/scoring/etable/Etable.hh>
#endif

#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh> //pba
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>

#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>


// C++ headers
#include <map>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

//singelton class
class ScoringManager : public utility::SingletonBase< ScoringManager >
{
public:
	friend class utility::SingletonBase< ScoringManager >;

public:

	void factory_register( methods::EnergyMethodCreatorOP creator );

	//P_AA const & get_P_AA() const;
	//ReferenceEnergyPotential const & get_ReferenceEnergyPotnential() const;
	PairEPotential const & get_PairEPotential() const;

	GenBornPotential const & get_GenBornPotential() const;

	HydroxylTorsionPotential const & get_HydroxylTorsionPotential() const;

	VdWTinkerPotential const & get_VdWTinkerPotential() const;

	MultipoleElecPotential const & get_MultipoleElecPotential( methods::EnergyMethodOptions const & options ) const;

	SASAPotential const & get_SASAPotential() const;

	FACTSPotential const & get_FACTSPotential() const;

	dna::DNA_BasePotential const & get_DNA_BasePotential() const;

	///RotamerLibrary & get_RotamerLibrary() const;

	RamaPrePro const & get_RamaPrePro() const;

	Ramachandran2BCOP get_Ramachandran2B_ptr() const;
	Ramachandran2B const & get_Ramachandran2B() const;

	RamachandranCOP get_Ramachandran_ptr() const;

	/// @brief Get a non-const owning pointer to the Ramachandran object.
	/// @details Warning!  Probably not threadsafe!  Needed, though, for
	/// lazy loading.
	RamachandranOP get_Ramachandran_ptr_nonconst();

	Ramachandran const & get_Ramachandran() const;

	/// @brief Get a non-const reference to the Ramachandran object.
	/// @details Warning!  Probably not threadsafe!  Needed, though, for
	/// lazy loading.
	Ramachandran & get_Ramachandran_nonconst();

	dna::DNABFormPotential const & get_DNABFormPotential() const;

	dna::DNATorsionPotential const & get_DNATorsionPotential() const;

	OmegaTether const & get_OmegaTether() const;

	SmoothEnvPairPotential const & get_SmoothEnvPairPotential() const;

	CenRotEnvPairPotential const & get_CenRotEnvPairPotential() const;

	CenHBPotential const & get_CenHBPotential() const;

	EnvPairPotential const & get_EnvPairPotential() const;

	dna::DNA_EnvPairPotential const & get_DNA_EnvPairPotential() const;
	dna::DNA_DihedralPotential const & get_DNA_DihedralPotential() const;

	SecondaryStructurePotential const & get_SecondaryStructurePotential() const;

	AtomVDW const & get_AtomVDW( std::string const & atom_type_set_name ) const;

	rna::RNA_AtomVDW const & get_RNA_AtomVDW() const;

	geometric_solvation::DatabaseOccSolEne const &
	get_DatabaseOccSolEne(
		std::string const & atom_type_set_name,
		Real const & min_occ_energy
	) const;

	carbon_hbonds::CarbonHBondPotential const & get_CarbonHBondPotential() const;

	loop_graph::evaluator::SixDTransRotPotential const & get_LoopCloseSixDPotential( std::string const & database_file ) const;

	rna::RNA_LowResolutionPotential const & get_RNA_LowResolutionPotential() const;

	rna::RNP_LowResPotential const & get_RNP_LowResPotential() const;

	rna::RNP_LowResStackData const & get_RNP_LowResStackData() const;

	rna::chemical_shift::RNA_ChemicalShiftPotential const & get_RNA_ChemicalShiftPotential() const;

	rna::data::RNA_DMS_Potential & get_RNA_DMS_Potential() const;

	rna::data::RNA_DMS_LowResolutionPotential & get_RNA_DMS_LowResolutionPotential() const;

	dna::DirectReadoutPotential const & get_DirectReadoutPotential() const;

	mm::MMLJLibrary const & get_MMLJLibrary() const;

	mm::MMLJEnergyTable const & get_MMLJEnergyTable() const;

	mm::MMTorsionLibrary const & get_MMTorsionLibrary() const;

	mm::MMBondAngleLibrary const & get_MMBondAngleLibrary() const;

	mm::MMBondLengthLibrary const & get_MMBondLengthLibrary() const;

	nv::NVlookup const & get_NVLookupTable() const;
	core::scoring::orbitals::OrbitalsLookup const & get_OrbitalsLookupTable() const;

	interface_::DDPlookup const & get_DDPLookupTable() const;

	P_AA const & get_P_AA() const;

	P_AA_ss const & get_P_AA_ss() const;

	UnfoldedStatePotential const & get_UnfoldedStatePotential( std::string const & type ) const;

	WaterAdductHBondPotential const & get_WaterAdductHBondPotential() const;

	MembranePotential const & get_MembranePotential() const;

	membrane::MembraneData const & get_MembraneData() const;

	Membrane_FAPotential const & get_Membrane_FAPotential() const; //pba

	ProQPotential const & get_ProQPotential() const;

	PoissonBoltzmannPotential const & get_PoissonBoltzmannPotential() const;

	SplitUnfoldedTwoBodyPotential const & get_SplitUnfoldedTwoBodyPotential(std::string const & label_type,std::string const & value_type, std::string const & score_func_type ) const;

	disulfides::FullatomDisulfidePotential &
	get_FullatomDisulfidePotential() const;

	disulfides::CentroidDisulfidePotential &
	get_CentroidDisulfidePotential() const;

	disulfides::DisulfideMatchingPotential &
	get_DisulfideMatchingPotential() const;

	carbohydrates::CHIEnergyFunction const & get_CHIEnergyFunction(
		bool setup_sampling_data = false,
		Real sampling_step_size = 0.1 ) const;

	carbohydrates::OmegaPreferencesFunction const & get_OmegaPreferencesFunction(
		bool setup_sampling_data = false,
		Real sampling_step_size = 0.1 ) const;


	bool
	has_energy_method( ScoreType t ) const;


	methods::EnergyMethodOP
	energy_method( ScoreType const & t, methods::EnergyMethodOptions const & options ) const;

	/// @brief make etable for extra partially softies
	etable::EtableOP
	make_partially_soft_etable( std::string const & name, etable::EtableOptions etable_options ) const;

	///pba
	void
	add_memb_etable( std::string const & name, etable::MembEtableOP etable );

	///pba
	etable::MembEtableCAP
	memb_etable( std::string const & table_id ) const;

	/// @brief Request an etable specifying an EnergyMethodOptions object; internally
	/// this will retrieve the EtableOptions object, and invoke the EtableOptions
	/// version of this function.
	etable::EtableCAP
	etable( methods::EnergyMethodOptions const & options_in ) const;

	/// @brief Request an etable specifying an EtableOptions; internally this will
	/// query the ScoringManager's map from EtableOptions to Etables for the desired
	/// Etable, and construct a new one if needed.
	etable::EtableCAP
	etable( etable::EtableOptions const & options_in ) const;

	/// @brief Create and return an etable specified only by the etable_type of the
	/// etable::EtableOptions class.  This, internally, will create an EtableOptions object,
	/// initialized from the command line, set the etable_type of this object, and then
	/// invoke the etable( EtableOptions ) method.
	etable::EtableCAP
	etable( std::string const & etable_id ) const;

	/// @brief Get an owning pointer to data used by the FA_ElecEnergy in beta_nov15 mode.
	/// @details If the data have not been loaded, this loads the data (lazy loading).  NOT THREADSAFE.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::scoring::elec::CPRepMapTypeCOP get_cp_rep_map_byname() const;

	/// @brief Get a vector of owning pointers to data used by the AACompositionEnergy score term.
	/// @details If this vector has not yet been populated, this loads the data from disk (lazy loading).  NOT THREADSAFE.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP > get_cloned_aa_comp_setup_helpers( core::scoring::methods::EnergyMethodOptions const &options ) const;

	/// @brief Get a particular MainchainScoreTable for the rama_prepro score term, for a particular residue type.
	/// @details If this has not yet been populated, loads the data from disk (lazy loading).  NOT THREADSAFE.
	/// @note Each restype stores separate tables for general and pre-proline scoring.  The prepro_table parameter determines
	/// whether we're loading the default scoring table or the version for residues occurring before a proline.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::chemical::mainchain_potential::MainchainScoreTableCOP
	get_rama_prepro_mainchain_torsion_potential(
		core::chemical::ResidueTypeCOP restype,
		bool const use_polycubic_interpolation,
		bool const prepro_table
	) const;

private:

	//private constructor
	ScoringManager();
	~ScoringManager();

private:

	// WARNING -- if you add something here don't forget to initialize to 0 in the constructor
	mutable VdWTinkerPotentialOP vdw_tinker_potential_;
	mutable PairEPotentialOP pairE_potential_;
	mutable RamachandranOP rama_;
	mutable Ramachandran2BOP rama2b_;
	mutable RamaPreProOP rama_pp_;
	mutable OmegaTetherOP omega_;
	mutable EnvPairPotentialOP env_pair_potential_;
	mutable SmoothEnvPairPotentialOP smooth_env_pair_potential_;
	mutable CenRotEnvPairPotentialOP cen_rot_pair_potential_;
	mutable CenHBPotentialOP cen_hb_potential_;
	mutable SecondaryStructurePotentialOP secondary_structure_potential_;
	mutable std::map< std::string, AtomVDWOP > atom_vdw_;
	mutable rna::RNA_AtomVDWOP rna_atom_vdw_;
	mutable geometric_solvation::DatabaseOccSolEneOP occ_hbond_sol_database_;
	mutable dna::DirectReadoutPotentialOP dna_dr_potential_;
	mutable mm::MMLJLibraryOP mm_lj_library_;
	mutable mm::MMLJEnergyTableOP mm_lj_energy_table_;
	mutable mm::MMTorsionLibraryOP mm_torsion_library_;
	mutable mm::MMBondAngleLibraryOP mm_bondangle_library_;
	mutable mm::MMBondLengthLibraryOP mm_bondlength_library_;
	mutable dna::DNA_EnvPairPotential * dna_env_pair_potential_;
	mutable dna::DNA_DihedralPotential * dna_dihedral_potential_;
	mutable dna::DNABFormPotentialOP dnabform_;
	mutable dna::DNATorsionPotentialOP dna_torsion_potential_;
	mutable dna::DNA_BasePotentialOP DNA_base_potential_;
	mutable carbon_hbonds::CarbonHBondPotentialOP carbon_hbond_potential_;
	mutable rna::RNA_LowResolutionPotentialOP rna_low_resolution_potential_;
	mutable rna::RNP_LowResPotentialOP rnp_low_res_potential_;
	mutable rna::RNP_LowResStackDataOP rnp_low_res_stack_data_;
	mutable rna::chemical_shift::RNA_ChemicalShiftPotentialOP rna_chemical_shift_potential_;
	mutable rna::data::RNA_DMS_PotentialOP rna_dms_potential_;
	mutable rna::data::RNA_DMS_LowResolutionPotentialOP rna_dms_low_resolution_potential_;
	mutable std::map< std::string, loop_graph::evaluator::SixDTransRotPotentialCOP > loop_close_six_d_potential_;
	mutable P_AAOP p_aa_;
	mutable P_AA_ssOP p_aa_ss_;
	mutable WaterAdductHBondPotentialOP water_adduct_hbond_potential_;
	mutable GenBornPotentialOP gen_born_potential_;
	mutable HydroxylTorsionPotentialOP hxl_tors_potential_;
	mutable MultipoleElecPotentialOP multipole_elec_potential_;
	mutable SASAPotentialOP sasa_potential_;
	mutable FACTSPotentialOP facts_potential_;
	mutable disulfides::FullatomDisulfidePotentialOP fa_disulfide_potential_;
	mutable disulfides::CentroidDisulfidePotentialOP cen_disulfide_potential_;
	mutable disulfides::DisulfideMatchingPotentialOP disulfide_matching_potential_;
	mutable MembranePotentialOP membrane_potential_;
	mutable membrane::MembraneDataOP mp_base_potential_;
	mutable Membrane_FAPotentialOP membrane_fapotential_; //pba
	mutable ProQPotentialOP ProQ_potential_;
	mutable PoissonBoltzmannPotentialOP PB_potential_;
	mutable SplitUnfoldedTwoBodyPotentialOP sutbp_;
	//ReferenceEnergyPotential referenceEnergyPotential_;
	mutable UnfoldedStatePotentialOP unf_state_;
	mutable carbohydrates::CHIEnergyFunctionOP CHI_energy_function_;
	mutable carbohydrates::OmegaPreferencesFunctionOP carbohydrate_omega_preferences_function_;

	mutable nv::NVlookupOP NV_lookup_table_;
	mutable orbitals::OrbitalsLookupOP orbitals_lookup_table_;
	mutable interface_::DDPlookupOP DDP_lookup_table_;
	// data

	// original map using string as key, let's keep it
	mutable std::map< std::string, etable::EtableOP > etables_by_string_;

	// new map for etables using EtableOptions as key
	mutable std::map< etable::EtableOptions, etable::EtableOP > etables_by_options_;

	mutable std::map< std::string, etable::MembEtableOP > memb_etables_; //pba

	/// @brief Cached data used by FA_ElecEnergy with beta_nov15.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable core::scoring::elec::CPRepMapTypeOP cp_rep_map_byname_;

	/// @brief Cached data used by the AACompositionEnergy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP > aa_composition_setup_helpers_;

	/// @brief Cached mainchain torsion potentials, used by rama_prepro.
	/// @details This one is for potentials for residues NOT occurring before proline.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::chemical::mainchain_potential::MainchainScoreTableOP > rama_prepro_mainchain_potentials_;

	/// @brief Cached mainchain torsion potentials, used by rama_prepro.
	/// @details This one is for potentials for residues occurring before proline.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::chemical::mainchain_potential::MainchainScoreTableOP > rama_prepro_mainchain_potentials_beforeproline_;


	utility::vector1< methods::EnergyMethodCreatorOP > method_creator_map_;


};
} // namespace core
} // namespace scoring


#endif
