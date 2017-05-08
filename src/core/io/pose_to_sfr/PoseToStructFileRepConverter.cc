// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/pose_to_sfr/PoseToStructFileRepConverter.cc
/// @brief A class to convert a pose to a StructFileRep.
/// @details This conversion is a first step in PDB or mmCIF output.  It could be useful for other
/// input/output, too.
/// @author Vikram K. Mulligan (vmullig@uw.edu), XRW 2016 Team
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

// Unit headers
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/Remarks.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/io/HeaderInformation.hh>
#include <core/io/ResidueInformation.hh>

// Core headers
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <sstream>
#include <vector>


namespace core {
namespace io {
namespace pose_to_sfr {
using utility::to_string;

static THREAD_LOCAL basic::Tracer TR( "core.io.pose_to_sfr.PoseToStructFileRepConverter" );

// Pose to StructFileRep methods ///////////////////////////////////////////////////

/// @brief Constructor.
/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
/// (using the sfr() object).  The options_ object is created and initialized from the options system.
PoseToStructFileRepConverter::PoseToStructFileRepConverter() :
	utility::pointer::ReferenceCount(),
	options_() //Initialize from options system.
{
	new_sfr();
}

/// @brief Constructor with options input.
/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
/// (using the sfr() object).  The options_ object is duplicated from the input options object.
PoseToStructFileRepConverter::PoseToStructFileRepConverter( StructFileRepOptions const & options_in ) :
	utility::pointer::ReferenceCount(),
	options_( options_in )
{
	new_sfr();
}


/// @brief Resets the PoseToStructFileRepConverter object, and reinitializes
/// it with a fresh StruftFileRep object, returning an owning pointer to the
/// new object.
core::io::StructFileRepOP
PoseToStructFileRepConverter::new_sfr() {
	sfr_ = StructFileRepOP( new core::io::StructFileRep() );
	return sfr_;
}


/// @details Non-const access to the StructFileRep object.
StructFileRepOP
PoseToStructFileRepConverter::sfr() {
	return sfr_;
}


/// @details Read atoms/residue information from Pose object and put it in StructFileRep object.
void
PoseToStructFileRepConverter::init_from_pose( core::pose::Pose const & pose )
{
	init_from_pose( pose, options_ );
}

/// @details Read atoms/residue information from Pose object and put it in StructFileRep object using options defined in
/// StructFileRepOptions.
void
PoseToStructFileRepConverter::init_from_pose( core::pose::Pose const & pose, StructFileRepOptions const & options )
{
	id::AtomID_Mask mask = id::AtomID_Mask( pose.size() );
	for ( Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		for ( Size atom_index = 1; atom_index <= pose.residue_type( resnum ).natoms(); ++atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			mask.set( atm, true );
		}
	}

	init_from_pose( pose, mask, options );
}

/// @details A lightweight, direct way of limiting pose pdb output to a subset of residues.
/// @note The alternative of constructing new subposes for output only would be unnecessary/less efficient (?)
void
PoseToStructFileRepConverter::init_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residue_indices
)
{
	using namespace core;

	id::AtomID_Mask mask = id::AtomID_Mask( pose.size() );
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {

		for ( Size atom_index = 1; atom_index <= pose.residue_type( resnum ).natoms(); ++ atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			if ( residue_indices.has_value( resnum ) ) {
				mask.set( atm, true );
			} else {
				mask.set( atm, false );
			}
		}

	}

	init_from_pose( pose, mask );
}

void
PoseToStructFileRepConverter::init_from_pose( core::pose::Pose const & pose, id::AtomID_Mask const & mask )
{
	init_from_pose( pose, mask, options_ );
}


// Main Init Function /////////////////////////////////////////////////////////
void
PoseToStructFileRepConverter::init_from_pose(
	core::pose::Pose const & pose,
	id::AtomID_Mask const & mask,
	StructFileRepOptions const & options )
{
	using namespace core;
	using core::pose::PDBInfo;

	new_sfr();
	if ( &options_ != &options ) { // Yes, address-of comparison
		// In an if block to avoid self (re)assignment, say from the init_from_pose(pose) call.
		options_ = options;
	}

	// Get Title Section information.
	if ( ( options.preserve_header() == true || options.preserve_crystinfo() == true ) && pose.pdb_info() ) {
		*(sfr_->remarks()) = pose.pdb_info()->remarks();  // Get OP to PDBInfo object for remarks.
		if ( pose.pdb_info()->header_information() ) {
			sfr_->header() =
				io::HeaderInformationOP( new io::HeaderInformation( *( pose.pdb_info()->header_information() ) ) );
		} else {
			sfr_->header() = io::HeaderInformationOP( new io::HeaderInformation() );
		}
	} else {
		sfr_->header() = io::HeaderInformationOP( new io::HeaderInformation() );
	}

	// Get parametric information
	if ( options.write_pdb_parametric_info() ) {
		get_parametric_info( sfr_->remarks(), pose );
	}

	// Get Connectivity Annotation Section information.
	if ( options.write_pdb_link_records() ) {
		get_connectivity_annotation_info( pose );
	}

	// Get Crystallographic and Coordinate Transformation Section information.
	if ( options.preserve_crystinfo() && pose.pdb_info() ) {
		sfr_->crystinfo() = pose.pdb_info()->crystinfo();
	}

	if ( options.output_secondary_structure() ) {
		generate_secondary_structure_informations( pose );
	}

	// Get membrane, other data:
	grab_additional_pose_data( pose );

	// Setup options.
	bool renumber_chains( false );
	if ( options.per_chain_renumbering() ) {
		renumber_chains = true;
	}

	//TR << " WE " << ( renumber_chains ? " are " : " aren't " ) << " renumebring chains" << std::endl;
	//utility::vector1< bool > res_info_added( pose.size(), false );

	core::Size new_atom_num( 1 );
	core::Size new_tercount( 0 );
	for ( Size resnum = 1, resnum_max = pose.size(); resnum <= resnum_max; ++resnum ) {
		conformation::Residue const & rsd( pose.residue( resnum ) );
		bool const use_pdb_info( use_pdb_info_for_num( pose, resnum ) );
		//TR << "use pdb info for num ? " << ( use_pdb_info ?  "yes" : "no" ) << std::endl;
		if ( ( resnum > 1 ) && ( pose.chain( resnum ) != pose.chain( resnum - 1 ) ) ) {
			++new_tercount;
		}

		ResidueInformation res_info;
		get_residue_information( pose, resnum, use_pdb_info, renumber_chains, new_tercount, res_info );
		//TR << "In init_from_pose " << resnum << " " << pose.residue_type( resnum ).name() << " " << pose.pdb_info()->chain( resnum ) << std::endl;
		//TR << "In init_from_pose " << resnum << " " << res_info.chainID() << std::endl;

		// Add res information only if we havn't done so already.
		//if ( ! res_info_added[ resnum ] ) {
		append_residue_info_to_sfr( res_info, rsd );
		//res_info_added[ resnum ] = true;
		//}

		for ( Size atom_index = 1; atom_index <= rsd.natoms(); ++atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			if ( ( ! mask[ atm ] ) || ( ! mask.has( atm ) ) ) { continue; }

			bool const success( append_atom_info_to_sfr(
				pose, res_info, rsd, atom_index, use_pdb_info, new_atom_num, new_tercount ) );
			if ( success ) { ++new_atom_num; }
		}
	}
}


// Append pdb information to StructFileRep for a single residue.
void
PoseToStructFileRepConverter::append_residue_to_sfr( core::pose::Pose const & pose, core::Size const resnum )
{
	core::Size new_atom_num_start = get_new_atom_serial_num();
	core::Size const new_tercount( pose.chain( resnum ) - 1 );
	append_residue_to_sfr( pose, resnum, new_atom_num_start, new_tercount );
}

/// @brief Append pdb information to StructFileRep for a single residue.
/// @details Start atom numbering from given n
void
PoseToStructFileRepConverter::append_residue_to_sfr(
	pose::Pose const & pose,
	core::Size const resnum,
	core::Size & new_atom_num,
	core::Size const new_tercount )
{
	using namespace core;
	using namespace utility;

	conformation::Residue const & rsd = pose.residue( resnum );

	bool use_pdb_info = use_pdb_info_for_num( pose, resnum );
	bool renumber_chains(false);
	if ( options_.per_chain_renumbering() ) {
		renumber_chains = true;
	}
	//TR << " append_residue " << resnum << " " << pose.pdb_info()->segmentID( resnum ) << std::endl;
	ResidueInformation res_info;
	get_residue_information( pose, resnum, use_pdb_info, renumber_chains, new_tercount, res_info );
	//TR << " append_residue " << resnum << " " << res_info.segmentID() << std::endl;
	append_residue_info_to_sfr( res_info, rsd );

	// Loop through each atom in the residue and generate ATOM or HETATM data.
	for ( Size j = 1; j <= rsd.natoms(); ++j ) {
		bool const success( append_atom_info_to_sfr(
			pose, res_info, rsd, j, use_pdb_info, new_atom_num, new_tercount ) );
		if ( success ) { ++new_atom_num; }
	}
}


/// @brief Append just residue-based info to StructFileRep
void PoseToStructFileRepConverter::append_residue_info_to_sfr(
	ResidueInformation const & res_info,
	conformation::Residue const & rsd )
{
	// Determine residue identifier information.

	// Generate HETNAM data, if applicable.
	if ( ! ( rsd.is_protein() ||
			rsd.is_NA() ||
			rsd.is_virtual_residue() ||
			rsd.type().is_membrane() ) ) {
		String const & hetID( rsd.name3() );
		String const & base_name( rsd.type().base_name() );

		String resSeq( utility::pad_left( res_info.resSeq(), 4 ) ); //("%4d", res_info.resSeq);
		String const resID = String( 1, res_info.chainID() ) + resSeq + String( 1, res_info.iCode() );

		String const text = resID + " " + base_name;

		if ( ! sfr_->heterogen_names().count( hetID ) ) {
			sfr_->heterogen_names()[ hetID ] = base_name;
		}
		sfr_->residue_type_base_names()[ res_info.resid() ] = std::make_pair( hetID, base_name );
	}
}


bool
PoseToStructFileRepConverter::append_atom_info_to_sfr(
	core::pose::Pose const & pose,
	ResidueInformation const & res_info,
	core::conformation::Residue const & rsd,
	core::Size const atom_index,
	bool const use_pdb_info)
{
	return append_atom_info_to_sfr( pose, res_info, rsd, atom_index, use_pdb_info, get_new_atom_serial_num() /*atom index*/, pose.chain( rsd.seqpos() ) - 1 /*Number of termini before this atom*/);
}


bool
PoseToStructFileRepConverter::append_atom_info_to_sfr(
	core::pose::Pose const & pose,
	ResidueInformation const & res_info,
	core::conformation::Residue const & rsd,
	core::Size const atom_index,
	bool const use_pdb_info,
	core::Size const new_atom_num,
	core::Size const new_tercount )
{
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	conformation::Atom const & atom = rsd.atom( atom_index ) ;

	if ( ! add_atom_to_sfr( pose, rsd, atom_index, use_pdb_info ) ) { return false; }

	AtomInformation ai;
	AtomInformation orb;  //have to initialize this out here.

	ai.isHet = (!rsd.is_polymer() || rsd.is_ligand());
	ai.chainID = res_info.chainID();
	ai.resSeq = res_info.resSeq();
	ai.iCode = res_info.iCode();
	ai.serial = new_atom_num;
	ai.name = rsd.atom_name( atom_index );
	ai.resName = rsd.name3();
	ai.x = atom.xyz()(1);
	ai.y = atom.xyz()(2);
	ai.z = atom.xyz()(3);
	ai.occupancy = 1.0; // dummy occupancy, can be overridden by PDBInfo
	ai.segmentID = res_info.segmentID();
	ai.terCount = new_tercount;

	// Output with pdb-specific info if possible.
	if ( use_pdb_info ) {
		if ( pdb_info->is_het( rsd.seqpos(), atom_index ) ) { // override standard het only if .is_het() is true
			ai.isHet = true;
		}
		ai.altLoc = pdb_info->alt_loc( rsd.seqpos(), atom_index );
		ai.occupancy = pdb_info->occupancy( rsd.seqpos(), atom_index );
		ai.temperature = pdb_info->temperature( rsd.seqpos(), atom_index );
		ai.segmentID = pdb_info->segmentID( rsd.seqpos() );
	}

	if ( options_.output_virtual_zero_occ() ) {
		if ( rsd.atom_type( atom_index ).is_virtual() ) {
			ai.occupancy = 0.0;
		}
	}

	grab_conect_records_for_atom( pose, rsd.seqpos(), atom_index, ai );

	// Element
	// (written by fpd; moved here by Labonte)
	core::chemical::AtomTypeSet const &ats = rsd.type().atom_type_set();
	ai.element = ats[atom.type()].element();
	if ( ai.element.length() == 1 ) ai.element = " "+ai.element;

	// 'chains' is member data
	if ( sfr_->chains().size() < static_cast <core::Size> (rsd.chain() + 1) ) sfr_->chains().resize( rsd.chain() + 1 );
	sfr_->chains()[rsd.chain()].push_back(ai);

	return true;
}

core::Size
PoseToStructFileRepConverter::get_new_atom_serial_num() const {
	core::Size new_atom_num;
	if ( total_sfr_atoms( *sfr_ ) == 0 ) {
		new_atom_num = 1;
	} else {
		core::Size total_chains = sfr_->chains().size();
		core::Size total_entries_last_chain = total_sfr_atoms( *sfr_, total_chains - 1);
		new_atom_num = sfr_->chains()[ total_chains - 1][ total_entries_last_chain - 1].serial;
	}

	return new_atom_num;
}

bool
PoseToStructFileRepConverter::add_atom_to_sfr(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	core::Size atom_index,
	bool use_pdb_info ) const
{
	//skip outputting virtual atom unless specified
	if ( !options_.output_virtual() &&
			rsd.atom_type( atom_index ).is_virtual() ) return false;

	//fpd optionally don't output centroids
	if ( options_.no_output_cen() &&
			rsd.atom_name( atom_index ) == " CEN" ) return false;

	// skip outputting zero occupancy atoms if specified
	if ( use_pdb_info && options_.suppress_zero_occ_pdb_output() &&
			( rsd.seqpos() <= pose.pdb_info()->nres() ) ) {
		if ( pose.pdb_info()->occupancy( rsd.seqpos(), atom_index ) < 0.0001 ) return false;
	}

	return true;
}

bool
PoseToStructFileRepConverter::use_pdb_info_for_num( pose::Pose const & pose, Size resnum )
{
	//TR << "use pdb info? " << resnum << " " << pose.pdb_info()->nres()  << " " << pose.pdb_info()->obsolete() << " " << options_.renumber_pdb() << std::endl;
	// Setup options.
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	if ( pdb_info
			&& !( pdb_info->obsolete() )
			&& resnum <= pdb_info->nres()
			&& !( options_.renumber_pdb() ) ) {
		return true;
	} else {
		return false;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


LinkInformation
PoseToStructFileRepConverter::get_link_record( core::pose::Pose const & pose, core::Size ii, core::Size conn ) {

	using namespace id;
	LinkInformation link;
	Size jj = pose.residue( ii ).connected_residue_at_resconn( conn );
	if ( jj == 0 ) {
		link.name1 = "ABORT";
		return link;
	}
	Size jj_conn = pose.residue( ii ).residue_connection_conn_id( conn );

	link.name1 = pose.residue_type( ii ).atom_name( pose.residue_type( ii ).residue_connect_atom_index( conn ) );
	link.resName1 = pose.residue_type( ii ).name3();
	link.chainID1 = pose.pdb_info()->chain( ii );
	link.resSeq1 = pose.pdb_info()->number( ii );
	link.iCode1 = pose.pdb_info()->icode( ii );
	std::stringstream ss;
	ss.width(6);
	ss << std::right << pose.pdb_info()->number( ii );
	link.resID1 = ss.str() + link.iCode1 + link.chainID1;

	link.name2 =  pose.residue_type( jj ).atom_name( pose.residue_type( jj ).residue_connect_atom_index( jj_conn ) );
	link.resName2 = pose.residue_type( jj ).name3();
	link.chainID2 = pose.pdb_info()->chain( jj );
	link.resSeq2 = pose.pdb_info()->number( jj );
	link.iCode2 = pose.pdb_info()->icode( jj );
	std::stringstream ss2;
	ss2.width(6);
	ss2 << std::right << pose.pdb_info()->number( jj );
	link.resID2 = ss2.str() + link.iCode2 + link.chainID2;

	// Calculate bond distance.
	uint start_atom_index = pose.residue_type( ii ).atom_index( link.name1 );
	uint stop_atom_index = pose.residue_type( jj ).atom_index( link.name2 );
	link.length = pose.conformation().bond_length(
		AtomID( start_atom_index, ii ),
		AtomID( stop_atom_index, jj ) );
	return link;
}

SSBondInformation
PoseToStructFileRepConverter::get_ssbond_record( core::pose::Pose const & pose, core::Size ii, core::Size conn ) {

	using namespace id;
	SSBondInformation ssbond;

	Size jj = pose.residue( ii ).connected_residue_at_resconn( conn );

	ssbond.resName1 = pose.residue_type( ii ).name3();
	ssbond.chainID1 = pose.pdb_info()->chain( ii );
	ssbond.resSeq1 = pose.pdb_info()->number( ii );
	ssbond.iCode1 = pose.pdb_info()->icode( ii );
	std::stringstream ss;
	ss.width(6);
	ss << std::right << pose.pdb_info()->number( ii );
	ssbond.resID1 = ss.str() + ssbond.iCode1 + ssbond.chainID1;

	ssbond.resName2 = pose.residue_type( jj ).name3();
	ssbond.chainID2 = pose.pdb_info()->chain( jj );
	ssbond.resSeq2 = pose.pdb_info()->number( jj );
	ssbond.iCode2 = pose.pdb_info()->icode( jj );
	std::stringstream ss2;
	ss2.width(6);
	ss2 << std::right << pose.pdb_info()->number( jj );
	ssbond.resID2 = ss2.str() + ssbond.iCode2 + ssbond.chainID2;

	// Calculate bond distance.
	uint start_atom_index = pose.residue_type( ii ).atom_index( pose.residue_type( ii ).get_disulfide_atom_name() );
	uint stop_atom_index = pose.residue_type( jj ).atom_index( pose.residue_type( ii ).get_disulfide_atom_name() );
	ssbond.length = pose.conformation().bond_length(
		AtomID( start_atom_index, ii ),
		AtomID( stop_atom_index, jj ) );

	return ssbond;
}

/// @brief Get connectivity annotation information from the Pose object and create LinkInformation and
/// SSBondInformation data as appropriate.
/// @author Watkins
void PoseToStructFileRepConverter::get_connectivity_annotation_info( core::pose::Pose const & pose ) {

	using namespace utility;
	using namespace id;
	using namespace kinematics;
	using namespace conformation;

	// In the past, we walked through the fold_tree and found the termini of all
	// the chemical edges.
	// This technique is not very general because many complicated branched
	// topologies would end up creating cycles and thus aren't included.
	// Similarly, it doesn't permit the use of SSBOND records, as the
	// FoldTree doesn't go through those (except in exotic and generally
	// temporary circumstances.

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		conformation::Residue const & ii_res = pose.residue( ii );

		if ( ii_res.has_lower_connect() ) {
			Size lower = ii_res.lower_connect().index();

			// if bonded to not ii - 1 or bonded to not ii - 1's upper
			if ( ii == 1 || ii_res.connected_residue_at_resconn( lower ) != ii - 1 ||
					ii_res.residue_connection_conn_id( lower ) != static_cast<Size>(pose.residue( ii - 1 ).upper_connect().index()) ) {

				LinkInformation link = get_link_record( pose, ii, lower );
				if ( link.name1 != "ABORT" ) {
					vector1<LinkInformation> links;
					// If key is found in the links map, add this new linkage information to the links already keyed to
					// this residue.
					if ( sfr_->link_map().count(link.resID1) ) {
						links = sfr_->link_map()[link.resID1];
					}
					links.push_back(link);

					sfr_->link_map()[link.resID1] = links;
				}
			}
		}

		if ( ii_res.has_upper_connect() ) {
			Size upper = ii_res.upper_connect().index();

			// if bonded to not ii + 1 or bonded to not ii + 1's lower
			if ( ii == pose.total_residue() || ii_res.connected_residue_at_resconn( upper ) != ii + 1 ||
					ii_res.residue_connection_conn_id( upper ) != static_cast<Size>( pose.residue( ii + 1 ).lower_connect().index()) ) {

				// Escape if it's bonded to residue 1's lower--we don't want to double-count cyclization here.
				// If jj < ii, we already counted it
				if ( ii_res.connected_residue_at_resconn( upper ) < ii ) continue;

				LinkInformation link = get_link_record( pose, ii, upper );
				if ( link.name1 != "ABORT" ) {
					vector1<LinkInformation> links;

					// If key is found in the links map, add this new linkage information to the links already keyed to
					// this residue.
					if ( sfr_->link_map().count(link.resID1) ) {
						links = sfr_->link_map()[link.resID1];
					}
					links.push_back(link);

					sfr_->link_map()[link.resID1] = links;
				}
			}
		}

		if ( ii_res.n_non_polymeric_residue_connections() == 0 ) continue;

		for ( Size conn = ii_res.n_polymeric_residue_connections()+1; conn <= ii_res.n_possible_residue_connections(); ++conn ) {

			Size jj = ii_res.connected_residue_at_resconn( conn );
			Size jj_conn = ii_res.residue_connection_conn_id( conn );

			// Either LINK or SSBOND
			// Note that it's not an SSBOND if you are a cysteine bonded
			// to the UPPER of another cysteine!
			// Are the two atoms both the get_disulfide_atom_name() of the
			// connected residue types?
			if ( ( pose.residue_type( ii ).has_property( chemical::DISULFIDE_BONDED ) || pose.residue_type( ii ).has_property( chemical::SIDECHAIN_THIOL ) ) &&
					( pose.residue_type( jj ).has_property( chemical::DISULFIDE_BONDED ) || pose.residue_type( jj ).has_property( chemical::SIDECHAIN_THIOL ) ) &&
					pose.residue_type( ii ).residue_connect_atom_index( conn ) ==
					pose.residue_type( ii ).atom_index( pose.residue_type( ii ).get_disulfide_atom_name() ) &&
					pose.residue_type( jj ).residue_connect_atom_index( jj_conn ) ==
					pose.residue_type( jj ).atom_index( pose.residue_type( jj ).get_disulfide_atom_name() ) ) {
				// Disulfide.

				// If jj < ii, we already counted it
				if ( jj < ii ) continue;

				SSBondInformation ssbond = get_ssbond_record( pose, ii, conn );
				vector1<SSBondInformation> ssbonds;

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( sfr_->ssbond_map().count(ssbond.resID1) ) {
					ssbonds = sfr_->ssbond_map()[ssbond.resID1];
				}
				ssbonds.push_back(ssbond);

				sfr_->ssbond_map()[ssbond.resID1] = ssbonds;

			} else {

				LinkInformation link = get_link_record( pose, ii, conn );
				if ( link.name1 != "ABORT" ) {
					vector1<LinkInformation> links;
					// If key is found in the links map, add this new linkage information to the links already keyed to
					// this residue.
					if ( sfr_->link_map().count(link.resID1) ) {
						links = sfr_->link_map()[link.resID1];
					}
					// If this link is found under the OTHER record...
					bool skip = false;
					if ( sfr_->link_map().count( link.resID2 ) ) {
						for ( Size i = 1; i <= sfr_->link_map()[link.resID2].size(); ++i ) {
							if ( link.resID1 == sfr_->link_map()[link.resID2][i].resID2 ) {
								skip = true;
								break;
							}
						}
					}
					if ( skip )  continue;
					// Make sure it isn't a dupe--for example, make sure that
					// we didn't already push this back as a lower to upper thing.
					// Right now, we assume you can't have two noncanonical connections to the same residue.
					// This is not strictly necessarily true, but until we implement
					// LinkInformation ==, it's good enough.
					bool push_it = true;
					for ( Size i = 1; i <= links.size(); ++i ) {
						if ( ( link.resID1 == links[i].resID1 && link.resID2 == links[i].resID2 )
								|| ( link.resID2 == links[i].resID1 && link.resID1 == links[i].resID2 ) ) {
							push_it = false;
						}
					}
					if ( push_it ) {
						links.push_back(link);
					}

					sfr_->link_map()[link.resID1] = links;
				}
			}
		}
	}
}

/// @brief Get parametric information from the Pose object and add it to the PDB remarks.
void PoseToStructFileRepConverter::get_parametric_info(
	core::io::RemarksOP remarks,
	core::pose::Pose const & pose
) {

	using namespace core::conformation::parametric;

	core::Size const nsets(pose.conformation().n_parameters_sets()); //How many ParametersSet objects are there in the pose?
	if ( nsets == 0 ) { return; }  //No need to proceed if this isn't a parametric conformation.

	for ( core::Size iset=1; iset<=nsets; ++iset ) { //Loop through all of the ParametersSet objects.
		ParametersSetCOP curset( pose.conformation().parameters_set(iset) );
		std::stringstream curset_summary;
		curset->get_pdb_remark( curset_summary );

		int cur_remark_number(1); //int instead of core::Size, to match the Remarks class.
		if ( remarks->size() > 0 ) cur_remark_number = (*sfr_->remarks())[remarks->size()-1].num + 1;

		std::string cur_remark_str;
		while ( std::getline(curset_summary,cur_remark_str) ) {
			core::io::RemarkInfo cur_remark;
			cur_remark.value = cur_remark_str; // Ugh.  The RemarkInfo class provides no getters or setters -- only public access to its members.
			cur_remark.num=cur_remark_number;
			remarks->push_back( cur_remark );
			cur_remark_str.clear();
		}
	}

	return;
}

/// @brief Get additional pose data.
/// @details This is rewritten from the old core/io/pdb/file_data.cc:write_additional_pdb_data() function.
/// This was a catch-all for dumping out all sorts of protocol-specific stuff (e.g. membrane info).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
PoseToStructFileRepConverter::grab_additional_pose_data( core::pose::Pose const &pose )
{
	grab_membrane_info( pose, options_.normalize_to_thk() );
	grab_foldtree( pose, options_.fold_tree_io() );
	grab_pdb_parents( pose, options_.pdb_parents() );
	grab_pdb_comments( pose, options_.pdb_comments() );
	grab_torsion_records( pose, options_.output_torsions() );
	grab_pdbinfo_labels( pose );

	//Here due to packstat.
	if ( options_.output_pose_energies_table() ) {
		grab_pose_energies_table( pose);
	}
	if ( options_.output_pose_cache() ) {
		grab_pose_cache_data( pose);
	}
}


/// @brief Set whether to write the fold tree, in the
/// StructFileRepOptions object (options_).
void
PoseToStructFileRepConverter::set_fold_tree_io(
	bool const setting
) {
	options_.set_fold_tree_io( setting );
}

/********* PRIVATE FUNCTIONS ************/

/// @brief Get the membrane information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @param [in] pose The pose.
/// @param [in] normalize_to_thk Normalized MEM lines, useful for visualizing the boundaries
/// of the membrane by coupling the NORM and THK coordinates.  Added by Rebecca.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
PoseToStructFileRepConverter::grab_membrane_info(
	core::pose::Pose const &pose,
	bool const normalize_to_thk
) {

	if ( pose.conformation().is_membrane() && normalize_to_thk == true ) {

		// Grab membrane residue & current data
		core::Real const thkn( pose.conformation().membrane_info()->membrane_thickness() );
		core::Vector const cntr( pose.conformation().membrane_info()->membrane_center(pose.conformation()) );
		core::Vector norm( pose.conformation().membrane_info()->membrane_normal(pose.conformation()) );

		// Actually normalize the membrane residue to thk
		norm.normalize();
		norm.x( norm.x() * thkn);
		norm.y( norm.y() * thkn);
		norm.z( norm.z() * thkn);

		// Get rescount, chain count
		core::Size const chaincount( sfr_->chains().size() - 1 );
		core::Size const rescount( sfr_->chains()[chaincount][ sfr_->chains()[chaincount].size() - 1 ].serial );
		char const newchain( sfr_->chains()[chaincount][ sfr_->chains()[chaincount].size() - 1 ].chainID + 1 );

		AtomInformation ai1, ai2, ai3;
		ai1.isHet = true;
		ai2.isHet = true;
		ai3.isHet = true;
		ai1.serial = rescount + 1;
		ai2.serial = rescount + 2;
		ai3.serial = rescount + 3;
		ai1.name = "THKN";
		ai2.name = "CNTR";
		ai3.name = "NORM";
		ai1.resName = "MEM";
		ai2.resName = "MEM";
		ai3.resName = "MEM";
		ai1.chainID = newchain;
		ai2.chainID = newchain;
		ai3.chainID = newchain;
		ai1.resSeq = rescount + 1;
		ai2.resSeq = rescount + 2;
		ai3.resSeq = rescount + 3;
		ai1.x = thkn; ai1.y = 0; ai1.z = 0;
		ai2.x = cntr.x(); ai2.y = cntr.y(); ai2.z = cntr.z();
		ai3.x = norm.x(); ai3.y = norm.y(); ai3.z = norm.z();
		ai1.occupancy=1.0;
		ai2.occupancy=1.0;
		ai3.occupancy=1.0;
		ai1.element=" H";
		ai2.element=" H";
		ai3.element=" H";

		sfr_->chains().push_back( utility::vector0< AtomInformation >() );
		sfr_->chains()[ chaincount + 1 ].push_back( ai1 );
		sfr_->chains()[ chaincount + 1 ].push_back( ai2 );
		sfr_->chains()[ chaincount + 1 ].push_back( ai3 );

		//out << "HETATM XXXX THKN MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, thkn) << "   0.000   0.000 \n";
		//out << "HETATM XXXX CNTR MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, cntr.x()) << F(8, 3, cntr.y()) << F(8, 3, cntr.z()) << "\n";
		//out << "HETATM XXXX NORM MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, norm.x()) << F(8, 3, norm.y()) << F(8, 3, norm.z()) << "\n";

		//              HETATM XXXX THKN MEM X  81      15.000   0.000   0.000
		//    HETATM XXXX CNTR MEM X  81       0.000   0.000   0.000
		//              HETATM XXXX NORM MEM X  81       0.000   0.000  15.000

	}

}

/// @brief Get the conect record information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @param [in] pose The pose, for determining what's bonded to what.
/// @param [in] atom_index_in_pose The index (in the whole pose) of the atom that we're setting up.
/// @param [in] res_index The current residue index.
/// @param [in] atom_index_rsd The index (in the residue) of the atom that we're setting up.
/// @param [in] dist_cutoff The atom separation, above which a CONECT record is written.
/// @param [in] write_virtuals Are virtual atoms being written out?
/// @param [in/out] ai The AtomInformation object that's being set up.  Modified by this function
void
PoseToStructFileRepConverter::grab_conect_records_for_atom(
	core::pose::Pose const &pose,
	core::Size const res_index,
	core::Size const atom_index_in_rsd,
	core::io::AtomInformation &ai
) const {

	if ( options_.skip_connect_info() ) return;
	if ( pose.size() == 0 ) return; //Probably unnecesary, but why not?

	debug_assert( res_index <= pose.size() );
	debug_assert( atom_index_in_rsd <= pose.residue( res_index ).natoms() );

	bool const write_virtuals( options_.output_virtual() );

	//Return if this is virtual and we're not writing virtuals:
	if ( ( pose.residue_type(res_index).is_virtual_residue() || pose.residue_type(res_index).atom_type(atom_index_in_rsd).is_virtual() ) && !write_virtuals ) return;

	bool const writeall( options_.write_all_connect_info() );
	bool const this_res_is_canonical_or_solvent( pose.residue_type(res_index).is_canonical() || pose.residue_type(res_index).is_solvent() );
	core::Real const dist_cutoff_sq( options_.connect_info_cutoff()*options_.connect_info_cutoff() );
	core::Size count(0);
	core::id::AtomID const this_atom_id( atom_index_in_rsd, res_index ); //The AtomID of this atom.
	utility::vector1<core::id::AtomID> const & bonded_ids(  pose.conformation().bonded_neighbor_all_res( this_atom_id, write_virtuals ) ); //List of AtomIDs of atoms bound to this atom.

	for ( core::Size i=1; i<=res_index; ++i ) { //Loop through all residues
		for ( core::Size j=1, jmax=pose.residue_type(i).natoms(); j<=jmax; ++j ) { //Loop through all atoms in this residue
			if ( i == res_index && j == atom_index_in_rsd ) {
				break;
			}

			if ( !write_virtuals && pose.residue_type(i).atom_type(j).is_virtual() ) continue; //Skip writing virtuals if we should do so.
			++count;
			if ( !writeall ) {
				if ( i == res_index ) { //If this is an intra-residue bond:
					if ( pose.residue_type(i).is_canonical() || pose.residue_type(i).is_solvent() ) continue; //Skip canonical and solvent intra-res bonds.
				} else { //If this is an inter-residue bond:
					if ( this_res_is_canonical_or_solvent && ( pose.residue_type(i).is_canonical() || pose.residue_type(i).is_solvent() ) ) continue; //Skip canonical or solvent inter-res bonds.
				}
			}

			core::id::AtomID const other_atom_id( j, i ); //Candidate other atom to which this one might be bonded.
			for ( core::Size n=1, nmax=bonded_ids.size(); n<=nmax; ++n ) {
				if ( bonded_ids[n] == other_atom_id ) { //If the candidate atom is in the list of bonded atoms for this atom, add CONECT data.
					if ( pose.xyz( this_atom_id ).distance_squared( pose.xyz( other_atom_id ) ) >= dist_cutoff_sq ) { //Are we over the distance cutoff?
						ai.connected_indices.push_back( count );
					}
					break;
				}
			} //Loop through bonded atoms

		} //Loop through all atoms
	} //Loop through all residues

}

/// @brief Get the foldtree from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_foldtree(
	core::pose::Pose const &pose,
	bool const output_foldtree
) {
	if ( !output_foldtree ) return;
	std::stringstream sstr;
	sstr << pose.fold_tree();
	sfr_->foldtree_string() = sstr.str();
}

/// @brief Get the pdb parents from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_pdb_parents(
	core::pose::Pose const &pose,
	bool const output_parents
) {
	if ( !output_parents ) return;
	std::string value;
	bool has_parents = core::pose::get_comment( pose, "parents", value );
	if ( has_parents ) {
		sfr_->additional_string_output() = sfr_->additional_string_output() + "REMARK PARENT    " + value.substr(0,5) + "\n";
	}
}

/// @brief Get the pdb comments from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_pdb_comments(
	core::pose::Pose const &pose,
	bool const output_comments
) {
	if ( !output_comments ) return;
	sfr_->pdb_comments() = core::pose::get_all_comments( pose );
}

/// @brief Get the torsion information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @details This is an incredibly ugly way to handle the output of these
/// records.  This abuses the REMARK lines in the PDB format, and should be
/// refactored.  Keeping as-is for now to preserve legacy behaviour (VKM
/// 29 January 2016, Chemical XRW).
void
PoseToStructFileRepConverter::grab_torsion_records(
	core::pose::Pose const &pose,
	bool const output_torsions
) {
	using namespace ObjexxFCL::format;

	if ( !output_torsions ) return;

	std::stringstream out("");
	core::pose::PDBInfoCOP pdb_info(pose.pdb_info());

	if ( !core::pose::is_ideal_pose(pose) ) {
		TR << "Ignoring out::file::output_torsions option because pose is non-ideal!" << std::endl;
	} else {
		ObjexxFCL::FArray1D_char dssp_reduced_secstruct(pose.size());
		scoring::dssp::Dssp(pose).dssp_reduced(dssp_reduced_secstruct);
		if ( pdb_info ) {
			out << "REMARK torsions: res pdbres pdbchain seq dssp phi psi omega" << std::endl;
		} else {
			out << "REMARK torsions: res    res    chain seq dssp phi psi omega" << std::endl;
		}
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( pdb_info ) {
				out << "REMARK " << I( 4, i ) << " " << I( 4, pose.pdb_info()->number(i)) << " " << pose.pdb_info()->chain(i) << " " << pose.residue_type( i ).name1() << " " <<
					dssp_reduced_secstruct(i) << " " << F( 9, 3, pose.phi(i)) << " " << F( 9, 3, pose.psi(i)) << " " << F( 9, 3, pose.omega(i)) << std::endl;
			} else {
				out << "REMARK " << I( 4, i ) << " " << I( 4, i) << " " << pose.chain(i) << " " << pose.residue_type( i ).name1() << " " <<
					dssp_reduced_secstruct(i) << " " << F( 9, 3, pose.phi(i)) << " " << F( 9, 3, pose.psi(i)) << " " << F( 9, 3, pose.omega(i)) << std::endl;
			}
		}
	}

	//Append the results of all of this to the additional_string_output in the StructFileRep:
	sfr_->additional_string_output() = sfr_->additional_string_output() + out.str();
}


/// @brief Get the pdbinfo labels from the pose and store them in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @details This also abuses REMARK lines, and should be rewritten.
void
PoseToStructFileRepConverter::grab_pdbinfo_labels(
	core::pose::Pose const &pose
) {
	using namespace ObjexxFCL::format;

	// Added by Daniel-Adriano Silva, used to write the PDBInfoLabels to the REMARK
	// First test that the pdb_info() is not empty
	if ( pose.pdb_info() ) {
		// Then output the labels
		std::stringstream out;
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			utility::vector1 < std::string > const tmp_v_reslabels( pose.pdb_info()->get_reslabels(i) ); //Ugh.  Passing a vector of strings by return value.
			core::Size const numLables( tmp_v_reslabels.size() );
			//Only write if the residue has any label (keep the file as small as possible)
			if ( numLables > 0 ) {
				out << "REMARK PDBinfo-LABEL: " << I( 4, i );
				for ( core::Size lndx=1; lndx <= numLables; ++lndx ) {
					out << " " << tmp_v_reslabels[lndx];
				}
				out << std::endl;
			}
		}
		sfr_->additional_string_output() = sfr_->additional_string_output() + out.str();
	}
}

void
PoseToStructFileRepConverter::grab_pose_energies_table(
	core::pose::Pose const & pose
){

	utility::vector1< std::string > labels;
	utility::vector1< core::Real >  weights;
	utility::vector1< std::vector< std::string > > table;

	core::scoring::EnergyMap emap_weights = pose.energies().weights();

	typedef utility::vector1<core::scoring::ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if ( emap_weights[ii] != 0 ) score_types.push_back(ii);
	}

	for ( core::scoring::ScoreType const & score_type : score_types ) {
		labels.push_back( name_from_score_type(score_type) );
	}
	labels.push_back( "total" );

	for ( core::scoring::ScoreType const & score_type : score_types ) {
		weights.push_back( emap_weights[score_type] );
	}

	using utility::to_string;
	core::Real pose_total = 0.0;
	if ( pose.energies().energies_updated() ) {
		std::vector<std::string> line;
		line.push_back( "pose" );
		for ( core::scoring::ScoreType const & score_type : score_types ) {
			core::Real score = (emap_weights[score_type] * pose.energies().total_energies()[ score_type ]);
			line.push_back( restrict_prec(score) );
			pose_total += score;
		}
		line.push_back( restrict_prec(pose_total)); //end first for overall pose energy;
		table.push_back(line);

		for ( core::Size j = 1, end_j = pose.size(); j <= end_j; ++j ) {
			line.clear();
			core::Real rsd_total = 0.0;
			line.push_back( pose.residue(j).name() + "_" + to_string( j ) );
			for ( core::scoring::ScoreType const & score_type : score_types ) {
				core::Real score = (emap_weights[score_type] * pose.energies().residue_total_energies(j)[ score_type ]);

				line.push_back( restrict_prec(score));
				rsd_total += score;
			}
			line.push_back( restrict_prec(rsd_total) ); // end line;
			table.push_back( line );
		}
	}


	sfr_->score_table_labels() =   labels;
	sfr_->score_table_weights() =  weights;
	sfr_->score_table_lines() =    table;
}

void
PoseToStructFileRepConverter::grab_pose_cache_data(const core::pose::Pose &pose){

	std::map< std::string, std::string > string_data;
	std::map< std::string, float >  float_data;

	// ARBITRARY_STRING_DATA
	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		basic::datacache::CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		assert( data.get() != NULL );

		for ( std::map< std::string, std::string >::const_iterator it( data->map().begin() ), end( data->map().end() );
				it != end;
				++it ) {
			//TR << it->first << " " << it->second << std::endl;
			string_data[ it->first ] = it->second;
		}
	}

	// ARBITRARY_FLOAT_DATA
	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		basic::datacache::CacheableStringFloatMapCOP data
			= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringFloatMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		assert( data.get() != NULL );

		for ( std::map< std::string, float >::const_iterator it( data->map().begin() ), end( data->map().end() );
				it != end;
				++it ) {
			//TR << it->first << " " << it->second << std::endl;
			float_data [ it->first ] = it->second;
		}
	}

	sfr_->pose_cache_string_data() = string_data;
	sfr_->pose_cache_float_data() =  float_data;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Get the total number of atoms in the SFR.
core::Size
PoseToStructFileRepConverter::total_sfr_atoms(
	StructFileRep const & sfr
) const {
	core::Size total = 0;
	for ( core::Size chain = 0, chain_max=sfr.chains().size(); chain < chain_max ; ++chain ) {
		total += total_sfr_atoms( sfr, chain );
	}
	return total;
}

/// @brief Get the total number of atoms in a chain in the SFR.
core::Size
PoseToStructFileRepConverter::total_sfr_atoms(
	StructFileRep const & sfr,
	core::Size const chain_num
) const {
	runtime_assert_string_msg( chain_num < sfr.chains().size(), "Error in core::io::pose_to_sfr::PoseToStructFileRepConverter::total_sfr_atoms(): The chain index is out of range.");
	return sfr.chains()[ chain_num ].size();
}

/// @brief Return the PDB resName, chainID, resSeq, and iCode for the given Rosetta sequence position.
/// @details Output is res_info.
void
PoseToStructFileRepConverter::get_residue_information(
	core::pose::Pose const & pose,
	core::uint const seqpos,
	bool const use_PDB,
	bool const renumber_chains,
	core::Size const new_tercount,
	ResidueInformation & res_info ) const
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::pose;
	using core::chemical::chr_chains;

	res_info.resName( pose.residue_type(seqpos).name3() );

	// Use PDB-specific information?
	if ( use_PDB ) {
		PDBInfoCOP pdb_info = pose.pdb_info();

		res_info.chainID( pdb_info->chain( seqpos ) );
		if ( res_info.chainID() == PDBInfo::empty_record() ) {  // safety
			TR.Warning << "PDBInfo chain ID was left as character '" << PDBInfo::empty_record()
				<< "', denoting an empty record; for convenience, replacing with space." << endl;
			res_info.chainID( ' ' );
		}
		res_info.resSeq(    pdb_info->number( seqpos ) );
		res_info.iCode(     pdb_info->icode( seqpos ) );
		res_info.segmentID( pdb_info->segmentID( seqpos ) );
		res_info.terCount( new_tercount );
		// ...or not?
	} else {
		uint const chain_num = pose.chain( seqpos );
		runtime_assert( chain_num > 0 );

		res_info.chainID( chr_chains[ ( chain_num - 1 ) % chr_chains.size() ] );
		res_info.resSeq( seqpos );
		res_info.iCode( ' ' );
		res_info.terCount( new_tercount );

		// If option is specified, renumber per-chain.
		if ( renumber_chains ) {
			vector1< uint > const & chn_ends = pose.conformation().chain_endings();
			for ( uint i = 1; i <= chn_ends.size(); ++i ) {
				if ( chn_ends[ i ] < seqpos ) {
					res_info.resSeq( seqpos - chn_ends[ i ] );
				}
			}
		}

		// Fix for >10k residues.
		res_info.resSeq( res_info.resSeq() % 10000 );
	}
	//TR << " res_info chain is " << res_info.chainID() << std::endl;
}

/* OLDer - DELETE - left for reference
core::Real restrict_prec( core::Real inval )
{
if ( inval >= 1 || inval <= -1 ) { // Don't alter value, as the default precision of 6 works fine, and we avoid rounding artifacts
return inval;
}
core::Real outval;
std::stringstream temp;
temp << std::fixed << std::setprecision(5) << inval;
temp >> outval;
return outval;
}
*/


// Left over from pose energies table.  Remove if possible! JAB
std::string restrict_prec( core::Real inval )
{
	if ( inval >= 1 || inval <= -1 ) { // Don't alter value, as the default precision of 6 works fine, and we avoid rounding artifacts
		return utility::to_string(inval);
	}
	core::Real outval;
	std::stringstream temp;
	temp << std::fixed << std::setprecision(5) << inval;
	temp >> outval;
	return utility::to_string(outval);

}

/// @brief fills HELIXInformation and SHEETInformation for SFR - wrappers for individual functions
/// @author Steven Lewis smlewi@gmail.com, but cribbed from code of Yifan Song as part of Cyrus Biotechnology
/// option dependencies in this function: options_.output_secondary_structure() controls if it runs at all, options_.per_chain_renumbering() is used with ResidueInformation, and it also calls get_residue_information()
void PoseToStructFileRepConverter::generate_secondary_structure_informations(
	core::pose::Pose const & pose
) {

	//skip if option not requested
	if ( !options_.output_secondary_structure() ) return;

	//generate ss string ourselves unless requested to use the one in the pose
	std::string secstructs(pose.secstruct());
	if ( !options_.do_not_autoassign_SS() ) { //this defaults FALSE, so usually this first if DOES run
		secstructs = core::scoring::dssp::Dssp(pose).get_dssp_secstruct();
	} else if ( secstructs == "" ) {  //I think this never works, because pose seems to set all L secstruct by default
		TR.Error << "PoseToStructFileRepConverter::generate_secondary_structure_informations:: you have requested secondary structure output (-out:file:output_secondary_structure true) without automatically generating the secstruct string (-out:file:do_not_autoassign_SS true), but your pose does not have its secstruct string set; skipping secondary structure output" << std::endl;
		return;
	} else if ( secstructs.size() != pose.size() ) {
		TR.Error << "PoseToStructFileRepConverter::generate_secondary_structure_informations:: you have requested secondary structure output (-out:file:output_secondary_structure true) without automatically generating the secstruct string (-out:file:do_not_autoassign_SS true), but your pose's secstruct string is a different length from your pose; skipping secondary structure output" << std::endl;
		return;
	} else { //we are checking if the string consists of only LLLLLLLLL... or if H and E occur
		std::size_t const found_H = secstructs.find("H");
		std::size_t const found_E = secstructs.find("E");
		if ( found_H == std::string::npos && found_E == std::string::npos ) {
			TR.Error << "PoseToStructFileRepConverter::generate_secondary_structure_informations:: you have requested secondary structure output (-out:file:output_secondary_structure true) without automatically generating the secstruct string (-out:file:do_not_autoassign_SS true), but your pose appears to have an all-L secondary structure.  It is likely that your pose never had its secondary structure set by DSSP.  No secondary structure will be output." << std::endl;
			return; //not clear that this return is relevant since the rest of this code is no-op in this circumstance; OK to change, later code-reviewer!
		}
	}

	//count numbers of SS elements for their ID field
	core::Size n_helix = 0, n_sheet = 0;//, n_loop = 0;
	core::Size new_tercount( 0 ); //we have to track this for ResidueInformation

	//Now we are going to iterate through the pose, identifying secondary structure elements
	for ( Size ires=1; ires<pose.size(); ++ires ) {
		char secstruct = secstructs[ires-1]; //H, E, or L; note indexing fix
		Size chain = pose.residue(ires).chain();
		Size jres = ires;

		//this is tracked for certain option combinations in ResidueInformation
		//this was copied from elsewhere in the file and NOT TESTED CAREFULLY
		if ( ( ires > 1 ) && ( pose.chain( ires ) != pose.chain( ires - 1 ) ) ) {
			++new_tercount;
		}

		// iterate jres until jres+1's secstruct mismatches, so that ires and jres mark the beginning and end of an ss segment
		while (
				jres < pose.size() //not past end of pose
				&& secstructs[jres] == secstruct //ss matches
				&& pose.residue(jres+1).chain() == chain /*still on same chain*/ ) {
			jres += 1;
		}

		ResidueInformation ires_info, jres_info;
		get_residue_information(pose, ires, use_pdb_info_for_num(pose, ires), options_.per_chain_renumbering(), new_tercount, ires_info);
		get_residue_information(pose, jres, use_pdb_info_for_num(pose, jres), options_.per_chain_renumbering(), new_tercount, jres_info);

		if ( secstruct == 'H' ) {
			n_helix += 1;
			core::Size const helix_length(jres-ires+1);
			generate_HELIXInformation(ires_info, jres_info, n_helix, helix_length);
		} else if ( secstruct == 'E' ) {
			n_sheet += 1;
			generate_SHEETInformation(ires_info, jres_info, n_sheet);
		}

		ires = jres; //increment ires to end of this ss element; loop iterator will move us to the start of the next ss element
	}
	return;
}

/// @brief fills one HELIXInformation into SFR
/// @author Steven Lewis smlewi@gmail.com
void PoseToStructFileRepConverter::generate_HELIXInformation(
	ResidueInformation const & start_info,
	ResidueInformation const & stop_info,
	core::Size const index,
	core::Size const length
){

	HELIXInformation helix;
	helix.helixID = index;
	helix.helix_name = std::string(ObjexxFCL::format::I(3, index)); //3-width string
	helix.name3_1 = start_info.resName();
	helix.chainID1 = start_info.chainID();
	helix.seqNum1 = start_info.resSeq();
	helix.icode1 = start_info.iCode();
	helix.name3_2 = stop_info.resName();
	helix.chainID2 = stop_info.chainID();
	helix.seqNum2 = stop_info.resSeq();
	helix.icode2 = stop_info.iCode();
	//IGNORING helixClass
	//IGNORING comment
	helix.length = length;

	sfr_->HELIXInformations().push_back(helix);

	return;
}

/// @brief fills one SHEETInformation into SFR
/// @author Steven Lewis smlewi@gmail.com
void PoseToStructFileRepConverter::generate_SHEETInformation(
	ResidueInformation const & start_info,
	ResidueInformation const & stop_info,
	core::Size const index
){

	SHEETInformation sheet;
	//IGNORING strand_num
	sheet.sheetID = std::string(ObjexxFCL::format::I(3, index)); //3-width string
	//IGNORING num_strands
	sheet.name3_1 = start_info.resName();
	sheet.chainID1 = start_info.chainID();
	sheet.seqNum1 = start_info.resSeq();
	sheet.icode1 = start_info.iCode();
	sheet.name3_2 = stop_info.resName();
	sheet.chainID2 = stop_info.chainID();
	sheet.seqNum2 = stop_info.resSeq();
	sheet.icode2 = stop_info.iCode();
	//IGNORING strandClass
	//IGNORING the remainder of the record, this is enough to get it to show nicely in pymol and PV

	sfr_->SHEETInformations().push_back(sheet);


	return;
}

} // namespace pose_to_sfr
} // namespace io
} // namespace core
