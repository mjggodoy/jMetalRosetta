// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Jared Adolf-Bryfogle <jadolfbr@gmail.com>

// Unit Headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package Headers
#include <core/io/carbohydrates/pose_io.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/carbohydrates/GlycanTreeSetObserver.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>

// Utility Headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

// External Headers
#include <boost/lexical_cast.hpp>
#include <basic/basic.hh>

// C++ Header
#include <list>
#include <utility>

// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.pose.carbohydrates.util" );


namespace core {
namespace pose {
namespace carbohydrates {
	using std::map;
	using namespace conformation::carbohydrates;
	
// Glycosylate the Pose at the given sequence position and atom using an IUPAC sequence.
/// @details   Format for <iupac_sequence>:\n
/// Prefixes apply to the residue to which they are attached, below indicated by residue n.\n
/// Residues are listed from N to 1, where N is the total number of residues in the saccharide.\n
/// The sequence is parsed by reading to the next hyphen, so hyphens are crucial.\n
/// Linkage indication: "(a->x)-" specifies the linkage of residue n, where a is the anomeric carbon number of residue
/// (n+1) and x is the oxygen number of residue n.  The first residue listed in the annotated sequence (residue N)
/// need not have the linkage prefix.  A ->4) ResidueType will automatically be assigned by default if not specified.\n
/// Anomer indication: The strings "alpha-" or "beta-" are supplied next, which determines the stereochemistry of the
/// anomeric carbon of the residue to which it is prefixed.  An alpha ResidueType will automatically be assigned by
/// default.\n
/// Stereochemical indication: "L-" or "D-" specifies whether residue n is an L- or D-sugar.  The default is "D-".\n
/// 3-Letter code: A three letter code (in sentence case) MUST be supplied next.  This specifies the "base sugar name",
/// e.g., Glc is for glucose.  (A list of all recognized 3-letter codes for sugars can be found in
/// database/chemical/carbohydrates/codes_to_roots.map.)\n
/// 1-Letter suffix: If no suffix follows, residue n will be linear.  If a letter is present, it indicates the ring
/// size, where "f" is furanose, "p" is puranose, and "s" is septanose.\n
/// Branches are indicated using nested brackets and are best explained by example:\n
/// beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc- is:\n
/// beta-D-Galp-(1->4)-D-GlcpNAc-\n
///                       |\n
///     alpha-L-Fucp-(1->3)
/// @remarks   The order of residues in the created pose is in the opposite direction as the annotated sequence of
/// saccharide residues, as sugars are named with the 1st residue as the "main chain", with all other residues named as
/// substituents and written as prefixes.  In other words, sugars are usually drawn and named with residue 1 to the
/// right.
/// At present time, param files only exist for a limited number of sugars! ~ Labonte
/// Also, I've not written this to handle creation of glycolipids yet.
void
glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & iupac_sequence,
	bool const idealize_linkages /*true*/ )
{
	using namespace std;
	using namespace utility;
	using namespace chemical;
	using namespace conformation;

	conformation::Residue const & residue( pose.residue( sequence_position ) );
	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose( residue.type().mode() ) );

	// Get list of carbohydrate ResidueTypes from which to construct the Pose.
	ResidueTypeCOPs residue_types( residue_types_from_saccharide_sequence( iupac_sequence, *residue_set ) );

	Size const n_types( residue_types.size() );
	if ( ! n_types ) {
		TR.Warning << "No saccharide residues in sequence to append to pose!" << endl;
		return;
	}

	// We need to reorder the ResidueTypes, such that each chain is complete before a new branch is added.
	reorder_saccharide_residue_types( residue_types );

	// Prepare the glycosylation site for branching.
	// TODO: Add code to attach to lipids.
	if ( residue.is_carbohydrate() ) {
		chemical::VariantType const variant(
			chemical::carbohydrates::CarbohydrateInfoManager::branch_variant_type_from_atom_name( atom_name ) );
		add_variant_type_to_pose_residue( pose, variant, sequence_position );
	} else {  // It's a typical branch point, in that it has a single type of branch point variant.
		add_variant_type_to_pose_residue( pose, SC_BRANCH_POINT, sequence_position );
	}


	// Now we can extend the Pose.
	// Keep track of branch points as we go.
	list< pair< uint, string > > branch_points;

	//JAB - reverting this till we can fix it
	// TODO: Fix it.
	/*
	Size const initial_sizes( pose.size() );
	PDBInfoOP info( new PDBInfo( *pose.pdb_info() ) );//Copy because as we add residues it will change PDB Info - and not the way we want.
	if ( ! info ) {
	info = PDBInfoOP( new PDBInfo( pose ) );
	info->name( pose.sequence() );  // Use the sequence as the default name.
	//pose.pdb_info( info ); //This copies PDBInfo into the pose, not the actuall OP.  Can't set it here.
	}

	char const last_chain_id( info->chain( initial_sizes ) );  // Get the chain ID of the last residue.
	char new_chain_id;
	uint last_chain_index( utility::ALPHANUMERICS.find( last_chain_id ) );
	if ( ( last_chain_index == string::npos ) ||  // not a standard chain ID
	( last_chain_index == utility::ALPHANUMERICS.size() - 1 ) ) {  // out of chain IDs
	new_chain_id = 'A';  // Wrap around to the beginning.
	} else {
	new_chain_id = utility::ALPHANUMERICS[ last_chain_index + 1 ];
	}
	*/

	// Begin with the first sugar.
	ResidueType const & first_sugar_type( *residue_types.front() );
	ResidueOP first_sugar( ResidueFactory::create_residue( first_sugar_type ) );
	string const upper_atom( first_sugar->carbohydrate_info()->anomeric_carbon_name() );
	pose.append_residue_by_atoms( *first_sugar, true, upper_atom, sequence_position, atom_name, true );

	// Build any other sugars.
	append_pose_with_glycan_residues( pose, ResidueTypeCOPs( residue_types.begin() + 1, residue_types.end() ) );

	// Let the Conformation know that it (now) contains sugars.
	pose.conformation().contains_carbohydrate_residues( true );




	/*
	// Reverting PDBInfo changes as they don't quite work.
	Size const sizes( pose.size() );
	uint new_seqpos( 0 );
	for ( uint i( initial_sizes + 1 ); i <= sizes; ++i ) {
	++new_seqpos;
	info->append_res( i - 1, pose.residue( i ).natoms() );
	info->chain( i, new_chain_id );
	info->number( i, new_seqpos );
	}
	info->name( info->name() + "_glycosylated" );
	info->obsolete( false );
	pose.pdb_info(info);
	*/

	//JAB - this leaves an intact PDBInfo, which we absolutely need for Link Records to be written out properly.
	// However, the PDBInfo records we start with are completely wiped out.
	// We tried a fix above, but this leaves the PDB unreadable by Rosetta...
	// Figure out a way to preserve some of the original info that got us here.

	PDBInfoOP info( new PDBInfo( pose ) );
	info->name( pose.sequence() );  // Use the sequence as the default name.
	pose.pdb_info( info );

	//TR << "InitialNRES: " << initial_sizes << " NRES: "<< pose.size() << " PDBINFO: "<< pose.pdb_info()->nres() << std::endl;
	TR << "Glycosylated pose with " << iupac_sequence << '-' << atom_name <<
		pose.residue( sequence_position ).name3() << sequence_position << endl;
	
	//Make sure that the new pose has a GlycanTreeSetObserver.
	if ( ! pose.glycan_tree_set() ){
		GlycanTreeSetObserverOP observer = GlycanTreeSetObserverOP(new GlycanTreeSetObserver(pose.conformation()));
		pose.observer_cache().set( datacache::GLYCAN_TREE_OBSERVER, observer, true /*auto_attach */ );
	}

	if ( idealize_linkages ) {
		TR.Debug << "Idealizing glycosidic torsions." << endl;
		idealize_last_n_glycans_in_pose( pose, n_types );
	}
	
}

// Glycosylate the Pose at the given sequence position using an IUPAC sequence.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, or Ser.
void
glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & iupac_sequence,
	bool const idealize_linkages /*true*/ )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose( pose, sequence_position, "ND2", iupac_sequence, idealize_linkages );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose( pose, sequence_position, "OG", iupac_sequence, idealize_linkages );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose( pose, sequence_position, "OG1", iupac_sequence, idealize_linkages );
	} else if ( glycosylation_site == "TRP" ) {
		glycosylate_pose( pose, sequence_position, "CD1", iupac_sequence, idealize_linkages );
	} else {
		utility_exit_with_message( glycosylation_site + " is not a common site of glycosylation or else it is "
			"ambiguous; Rosetta cannot determine attachment atom.  Use glycosylate_pose( Pose & pose, uint const "
			"sequence_position, std::string const & atom_name, std::string const & iupac_sequence ) instead." );
	}
}


// Glycosylate the Pose at the given sequence position and atom using a .GWS or IUPAC sequence file.
void
glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & filename,
	bool const idealize_linkages /*true*/ )
{
	using namespace std;
	using namespace io::carbohydrates;

	string const & sequence_file( core::chemical::carbohydrates::find_glycan_sequence_file( filename ) );
	if ( sequence_file.find( ".gws" ) != string::npos ) {
		utility_exit_with_message( "Rosetta cannot yet read .gws format sequence files." );  // TODO: Change this.
		//string const & gws_sequence( core::chemical::carbohydrates::read_glycan_sequence_file( sequence_file ) );
	} else {  // Assume any other file type contains an IUPAC sequence.
		string const & iupac_sequence( core::chemical::carbohydrates::read_glycan_sequence_file( sequence_file ) );
		glycosylate_pose( pose, sequence_position, atom_name, iupac_sequence, idealize_linkages );
	}
}

// Glycosylate the Pose at the given sequence position using a .GWS or IUPAC sequence file.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, or Ser.
void
glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & filename,
	bool const idealize_linkages /*true*/ )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose_by_file( pose, sequence_position, "ND2", filename, idealize_linkages );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG", filename, idealize_linkages );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG1", filename, idealize_linkages );
	} else if ( glycosylation_site == "TRP" ) {
		glycosylate_pose_by_file( pose, sequence_position, "CD1", filename, idealize_linkages );
	} else {
		utility_exit_with_message( glycosylation_site + " is not a common site of glycosylation or else it is "
			"ambiguous; Rosetta cannot determine attachment atom.  Use glycosylate_pose_by_file( Pose & pose, uint "
			"const sequence_position, std::string const & atom_name, std::string const & filename ) instead." );
	}
}
	


// Return pointers to the two residues of the glycosidic bond.
/// @return  Pointers to the residue at <sequence_position> and its parent or else the same pointer twice if undefined.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP >
get_glycosidic_bond_residues( Pose const & pose, uint const sequence_position )
{
	using namespace std;
	using namespace conformation;

	// Get the 1st residue of interest.
	ResidueCOP res_n( pose.residue( sequence_position ).get_self_ptr() );

	if ( res_n->is_lower_terminus() ) {
		if ( TR.Info.visible() ) {
			TR.Info << "Glycosidic torsions are undefined for the first polysaccharide residue of a chain unless part "
				"of a branch." << endl;
		}
		return make_pair( res_n, res_n );
	}

	// Get the 2nd residue of interest.
	// (res_n_minus_1 is a misnomer for the lower termini of branches.)
	ResidueCOP res_n_minus_1;
	if ( pose.glycan_tree_set() ){
		 res_n_minus_1 = pose.residue( pose.glycan_tree_set()->get_parent( sequence_position ) ).get_self_ptr();
	}
	else {
		res_n_minus_1 = pose.residue( find_seqpos_of_saccharides_parent_residue( *res_n ) ).get_self_ptr();
	}

	return make_pair( res_n, res_n_minus_1 );
}



// Return the AtomIDs of the four phi torsion reference atoms.
/// @details For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...\n
utility::vector1< id::AtomID >
get_reference_atoms_for_phi( Pose const & pose, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}
	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is O(cyclic) for cyclic saccharides, C? for linear saccharides.
	// Because the cyclic O is not connected by the atom tree,
	// we actually will return the AtomID for the virtual atom that superimposes with the real atom.
	AtomID ref1;
	if ( residues.first->carbohydrate_info()->is_cyclic() ) {
		ref1 = AtomID( residues.first->carbohydrate_info()->virtual_cyclic_oxygen_index(), residues.first->seqpos() );
	} else /* is linear */ {
		utility_exit_with_message("Linear Polysacharides are not yet handled!");  // TODO: Figure out how linear polysaccharides are handled by IUPAC.
	}
	ids[ 1 ] = ref1;

	// Reference 2 is always the anomeric carbon.
	AtomID const ref2( residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is OX(n-1) for polysaccharides.
	AtomID const ref3( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for phi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}

// Return the AtomIDs of the four psi torsion reference atoms.
/// @details For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_psi( Pose const & pose, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is always the anomeric carbon.
	AtomID const ref1( residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos() );
	ids[ 1 ] = ref1;

	// Reference 2 is OX(n-1) for polysaccharides.
	AtomID const ref2( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is CX(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-1(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for psi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}

// Return the AtomIDs of the four omega torsion reference atoms.
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-CX-2(n-1),
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_1st_omega( Pose const & pose, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}
	
	
	//if ( residues.second->is_carbohydrate() && ( ! has_exocyclic_glycosidic_linkage( *residues.first, *residues.second ) ) ) {
	bool exocylic_linkage;
	if ( pose.glycan_tree_set() ){
		exocylic_linkage = pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage( residues.first->seqpos() );
	}
	else {
		exocylic_linkage = has_exocyclic_glycosidic_linkage(pose.conformation(), residues.first->seqpos() );
	}

	if (residues.second->is_carbohydrate() && ( ! exocylic_linkage)){
		TR.Debug << "Omega is undefined for residue " << sequence_position <<
			" because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}
	
	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is OX(n-1) for polysaccharides.
	AtomID const ref1( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 1 ] = ref1;

	// Reference 2 is CX(n-1) for polysaccharides.
	AtomID const ref2( residues.second->type().atom_base( ref1.atomno() ), residues.second->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is CX-1(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-2(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for omega: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}

// Return the AtomIDs of the four omega2 torsion reference atoms.
/// For carbohydrates glycosylated at positions with more than 1 exocyclic
/// position, or for linkages to longer amino acid residues,
/// omega2 of residue n is defined as CX(n-1)-CX-1(n-1)-CX-2(n-1)-CX-3(n-1),
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_2nd_omega( Pose const & pose, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}

	bool exocyclic_linkage;
	if ( pose.glycan_tree_set() ){
		exocyclic_linkage = pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage( residues.first->seqpos() );
	}
	else {
		exocyclic_linkage = has_exocyclic_glycosidic_linkage(pose.conformation(), residues.first->seqpos() );
	}
	
	if ( residues.second->is_carbohydrate() && ( ! exocyclic_linkage ) ) {
		TR.Debug << "Omega2 is undefined for residue " << sequence_position <<
			" because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}

	// Set the atom names of the four reference atoms.
	// Reference 0 is OX(n-1) for polysaccharides.
	AtomID const ref0( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );

	// Reference 1 is CX(n-1) for polysaccharides.
	AtomID const ref1( residues.second->type().atom_base( ref0.atomno() ), residues.second->seqpos() );


	// Reference 2 is CX-1(n-1) for polysaccharides.
	AtomID const ref2( residues.second->type().atom_base( ref1.atomno() ), residues.second->seqpos() );

	//If Ref2 is a ring atom and we do not have an omega2 angle
	if ( residues.second->is_carbohydrate() && residues.second->type().is_ring_atom( 1, ref2.atomno() ) ) {
		return ids;
	} else if ( (! residues.second->is_carbohydrate() ) && (residues.second->aa() == core::chemical::aa_ser || residues.second->aa() == core::chemical::aa_thr ) ) {
		//O-linked Glycosylation has 3 dihedrals.  We should address this better.
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.
	ids[ 1 ] = ref1;
	ids[ 2 ] = ref2;

	// Reference 3 is CX-2(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-3(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for omega2: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}


// Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1< id::AtomID >
get_reference_atoms( uint const torsion_id, Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace utility;

	vector1< AtomID > ref_atoms;
	switch ( torsion_id ) {
	case phi_torsion :
		ref_atoms = get_reference_atoms_for_phi( pose, sequence_position );
		break;
	case psi_torsion :
		ref_atoms = get_reference_atoms_for_psi( pose, sequence_position );
		break;
	case omega_torsion :
		ref_atoms = get_reference_atoms_for_1st_omega( pose, sequence_position );
		break;
	case 4 :
		ref_atoms = get_reference_atoms_for_2nd_omega( pose, sequence_position );
		break;
	default :
		utility_exit_with_message( "An invalid torsion angle was requested." );
	}
	return ref_atoms;
}


// Virtual Atom Alignment /////////////////////////////////////////////////////
// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
// conformation.
/// @details  This method aligns virtual atom VOX, where X is the position of the cyclic oxygen, OY and HOY, where Y is
/// the position of the anomeric carbon, (provided the residue is not the reducing end, where OY and HOY would be real
/// atoms), and HOZ, where Z is the mainchain glycosidic bond location.  OY and HOY are aligned with the last two "main
/// chain" atoms of the parent residue.  This ensures that torsion angles with duplicate names, e.g., chi1 and phi for
/// internal linked aldoses, will always return the same values.  The same concept applies for HOZ, which aligns with
/// the anomeric carbon of the downstream residue.  This method should be called after any coordinate change for a sac-
/// charide residue and after loading a saccharide residue from a file or sequence for the first time.
/// @note     Do I need to worry about aligning the virtual atoms left over from modified sugar patches?  Can such
/// virtual atoms be deleted?
void
align_virtual_atoms_in_carbohydrate_residue( conformation::Conformation & conf, uint const sequence_position ) {
	using namespace std;
	using namespace id;
	using namespace conformation;

	TR.Debug << " Aligning virtual atoms on residue " << sequence_position << "..." << endl;

	ResidueCOP res( conf.residue( sequence_position ).get_self_ptr() );

	// Find and align VOX, if applicable.
	if ( res->carbohydrate_info()->is_cyclic() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning VOX..." << endl;
		}
		uint const x( res->carbohydrate_info()->cyclic_oxygen() );
		uint const OX( res->atom_index( res->carbohydrate_info()->cyclic_oxygen_name() ) );
		uint const VOX( res->atom_index( "VO" + string( 1, x + '0' ) ) );

		conf.set_xyz( AtomID( VOX, sequence_position ), conf.xyz( AtomID( OX, sequence_position ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  VOX aligned." << endl;
		}
	}

	// Find and align OY and HOY, if applicable.
	if ( ! res->is_lower_terminus() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning OY and HOY..." << endl;
		}
		uint const y( res->carbohydrate_info()->anomeric_carbon() );
		uint const OY( res->atom_index( "O" + string( 1, y + '0') ) );
		uint const HOY( res->atom_index( "HO" + string( 1, y + '0' ) ) );

		uint const parent_res_seqpos( find_seqpos_of_saccharides_parent_residue( *res ) );
		ResidueCOP parent_res( conf.residue( parent_res_seqpos ).get_self_ptr() );
		uint const OY_ref( parent_res->connect_atom( *res ) );
		uint const HOY_ref( parent_res->first_adjacent_heavy_atom( OY_ref ) );

		conf.set_xyz( AtomID( HOY, sequence_position ), conf.xyz( AtomID( HOY_ref, parent_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  HOY aligned with atom " << parent_res->atom_name( HOY_ref ) <<
				" of residue " << parent_res_seqpos << endl;

			TR.Debug << "   Updating torsions..." << endl;
		}
		ResidueCOP dummy( conf.residue( sequence_position ).get_self_ptr() );  // to trigger private method commented below
		//conf.update_residue_torsions( res->seqpos(), false );
		if ( TR.Debug.visible() ) {
			TR.Debug << "   Torsions updated." << endl;
		}

		conf.set_xyz( AtomID( OY, sequence_position ), conf.xyz( AtomID( OY_ref, parent_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  OY aligned with atom " << parent_res->atom_name( OY_ref ) <<
				" of residue " << parent_res_seqpos << endl;
		}
	}

	// Find and align HOZ(s), if applicable.
	if ( ! res->is_upper_terminus() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning HOZ..." << endl;
		}
		uint const z( res->carbohydrate_info()->mainchain_glycosidic_bond_acceptor() );
		uint const HOZ( res->atom_index( "HO" + string( 1, z + '0' ) ) );

		uint const downstream_res_seqpos( sequence_position + 1 );
		ResidueCOP downstream_res( conf.residue( downstream_res_seqpos ).get_self_ptr() );
		uint const HOZ_ref( downstream_res->atom_index( downstream_res->carbohydrate_info()->anomeric_carbon_name() ) );

		conf.set_xyz( AtomID( HOZ, sequence_position ), conf.xyz( AtomID( HOZ_ref, downstream_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  HOZ aligned." << endl;
		}
	}
	Size const n_branches( res->carbohydrate_info()->n_branches() );
	for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
		uint const z( res->carbohydrate_info()->branch_point( branch_num ) );
		uint const OZ( res->atom_index( "O" + string( 1, z + '0' ) ) );
		uint const HOZ( res->atom_index( "HO" + string( 1, z + '0' ) ) );

		uint const branch_connection_id( res->type().residue_connection_id_for_atom( OZ ) );
		uint const branch_res_seqpos( res->residue_connection_partner( branch_connection_id ) );
		ResidueCOP branch_res( conf.residue( branch_res_seqpos ).get_self_ptr() );
		uint const HOZ_ref( branch_res->atom_index( branch_res->carbohydrate_info()->anomeric_carbon_name() ) );

		conf.set_xyz( AtomID( HOZ, sequence_position ), conf.xyz( AtomID( HOZ_ref, branch_res_seqpos ) ) );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << " All virtual atoms aligned." << endl;
	}
}


// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given pose.
void
align_virtual_atoms_in_carbohydrate_residue( Pose & pose, uint const sequence_position ) {
	align_virtual_atoms_in_carbohydrate_residue( pose.conformation(), sequence_position );
}


// TorsionID Queries //////////////////////////////////////////////////////////
// Is this is the phi torsion angle of a glycosidic linkage?
/// @details  Carbohydrate linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 1 ) is phi.  For a sugar, phi (of the next
/// residue) is the last torsion, and the number of main-chain torsions varies per saccharide residue.
bool
is_glycosidic_phi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == BRANCH ) {  // Phi for a branched saccharide has a BRANCH type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case  BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() );  // The last BB is phi.
				}
			}
			break;
		case BRANCH :
			if ( torsion_id.torsion() <= residue.n_non_polymeric_residue_connections() ) {
				Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					return true;  // If it's a branch to a sugar, it must be the phi from the branching residue.
				}
			}
			break;
		default :
			break;
		}
	}
	return false;
}

// Is this is the psi torsion angle of a glycosidic linkage?
/// @details  Carbohydrates linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 2 ) is psi.  For a sugar, psi (of the next
/// residue) is the penultimate torsion, and the number of main-chain torsions varies per saccharide residue.
/// @note  This function will return false if the TorsionID is a CHI and the torsion is already covered by an equivalent
/// BB.  In other words, TorsionIDs for CHIs only return true for branch connections, where the ONLY way to access that
/// torsion is through CHI.
bool
is_glycosidic_psi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Psi for a branched saccharide has a CHI type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() - 1 );
				}
			}
			break;
		case CHI :
			{
			Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
			// A psi angle will always have the third atom of its definition be a connect atom.
			uint const third_atom( residue.chi_atoms( torsion_id.torsion() )[ 3 ] );
			Size const n_branches( residue.n_non_polymeric_residue_connections() );
			for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + branch_num );
				conformation::Residue const & next_rsd( pose.residue( next_rsd_num ) );
				if ( ( next_rsd.is_carbohydrate() ) && ( residue.connect_atom( next_rsd ) == third_atom ) ) {
					return true;
				}
			}
		}
			break;
		default :
			break;
		}
	}
	return false;
}

// Is this is an omega torsion angle of a glycosidic linkage?
/// @details  Carbohydrates linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 3 ) is omega.  For a sugar, omega (of the next
/// residue) is the 3rd-to-last torsion, and the number of main-chain torsions varies per saccharide residue.
/// @note     This function currently will not recognize additional omegas beyond omega1.\n
/// Also, this function will return false if the TorsionID is a CHI and the torsion is already covered by an equivalent
/// BB.  In other words, TorsionIDs for CHIs only return true for branch connections, where the ONLY way to access that
/// torsion is through CHI.
bool
is_glycosidic_omega_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Omega for a branched saccharide has a CHI type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					chemical::carbohydrates::CarbohydrateInfoCOP info( residue.carbohydrate_info() );
					if ( info->has_exocyclic_linkage_to_child_mainchain() ) {
						return ( torsion_id.torsion() == residue.n_mainchain_atoms() - 2 );
					}
				}
			}
			break;
		case CHI :
			{
			Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
			// An omega angle will always have the fourth atom of its definition be a connect atom.
			uint const fourth_atom( residue.chi_atoms( torsion_id.torsion() )[ 4 ] );
			Size const n_branches( residue.n_non_polymeric_residue_connections() );
			for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + branch_num );
				conformation::Residue const & next_rsd( pose.residue( next_rsd_num ) );
				if ( ( next_rsd.is_carbohydrate() ) && ( residue.connect_atom( next_rsd ) == fourth_atom ) ) {
					return true;
				}
			}
		}
			break;
		default :
			break;
		}
	}
	return false;
}

// Return the sequence position of the immediate downstream (child) residue affected by this torsion.
/// @note  This method is designed for saccharide residues, because it makes assumptions about CHI and BRANCH numbering.
core::uint
get_downstream_residue_that_this_torsion_moves( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	conformation::Residue const & rsd( pose.residue( torsion_id.rsd() ) );
	uint next_rsd( 0 );
	if ( torsion_id.type() == BB ) {
		debug_assert( torsion_id.torsion() <= rsd.n_mainchain_atoms() );
		// If this is a main-chain torsion, we can be confident that the next residue on this chain MUST be n+1.
		next_rsd = torsion_id.rsd() + 1;
	} else if ( torsion_id.type() == BRANCH ) {
		debug_assert( torsion_id.torsion() <= rsd.n_non_polymeric_residue_connections() );
		Size const n_mainchain_connections( rsd.n_polymeric_residue_connections() );
		next_rsd = rsd.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
	} else if ( torsion_id.type() == CHI ) {
		next_rsd = find_seqpos_of_saccharides_child_residue_at( rsd, torsion_id.torsion() );
	}
	if ( ! next_rsd ) {
		TR.Error << "Torsion " << torsion_id << " does not move any downstream residues!  Returning 0." << std::endl;
	}
	return next_rsd;
}


// Torsion Access /////////////////////////////////////////////////////////////
// Getters ////////////////////////////////////////////////////////////////////

///@brief Get the number of glycosidic torsions for this residue.  Up to 4 (omega2).
Size
get_n_glycosidic_torsions_in_res( Pose & pose, uint const sequence_position )
{
	core::Size n_torsions( 0 );
	for ( core::Size torsion_id = 1; torsion_id <= 4; ++torsion_id ) {
		utility::vector1< id::AtomID > const ref_atoms = get_reference_atoms( torsion_id, pose, sequence_position );
		if ( ref_atoms.size() != 0 ) {
			n_torsions+=1;
		}
	}
	return n_torsions;
}


// Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Residue::mainchain_torsion() since the main-chain torsions of saccharide
/// residues only make sense in the context of two residues.  Moreover, the reference atoms as defined by the IUPAC are
/// different from the ones that Rosetta uses by default for mainchain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
core::Angle
get_glycosidic_torsion( uint const torsion_id, Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	utility::vector1< id::AtomID > const ref_atoms = get_reference_atoms( torsion_id, pose, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		TR.Debug << "Returning zero." << std::endl;
		return 0.0;
	}

	Angle const angle_in_radians(
		pose.conformation().torsion_angle( ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ]) );
	return principal_angle_degrees( conversions::degrees( angle_in_radians ) );
}


// Setters ////////////////////////////////////////////////////////////////////
// Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Conformation::set_torsion() since the reference atoms as defined by the
/// IUPAC are different from the ones that Rosetta uses by default for main-chain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3\n
/// <setting> is in degrees.
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
void
set_glycosidic_torsion( uint const torsion_id, Pose & pose, uint const sequence_position, core::Angle const setting )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1< AtomID > const ref_atoms = get_reference_atoms( torsion_id, pose, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		return;
	}

	Angle const setting_in_radians( conversions::radians( setting ) );
	pose.conformation().set_torsion_angle(
		ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ], setting_in_radians );
	align_virtual_atoms_in_carbohydrate_residue( pose, sequence_position );
}


// Glycosylation //////////////////////////////////////////////////////////////
// Idealize the glycosidic torsion angles for the last n glycan residues added or built.
void
idealize_last_n_glycans_in_pose( Pose & pose, Size const n_glycans_added )
{
	using namespace std;
	using namespace conformation;

	// Work backward through the list of glycans and fix their glycosidic bonds to be ideal, if such data is known.
	Size const sizes( pose.size() );
	uint const glycosylation_site( sizes - n_glycans_added );  // 0, if a free oligo-/polysaccharide was built.
	for ( uint i( sizes ); i > glycosylation_site; --i ) {
		Residue const & non_red_end_res( pose.residue( i ) );
		if ( ! non_red_end_res.is_carbohydrate() ) {
			TR.Warning << "Cannot idealize a non-saccharide residue!" << endl;
			continue;
		}
		uint parent_seqpos( pose.glycan_tree_set()->get_parent( i ) );
		if ( parent_seqpos == 0 ) { continue; }
		Residue const & red_end_res( pose.residue( parent_seqpos ) );
		string const & non_red_end_short_name( non_red_end_res.carbohydrate_info()->short_name() );
		string red_end_short_name;
		if ( red_end_res.is_carbohydrate() ) {
			red_end_short_name = red_end_res.carbohydrate_info()->short_name();  // 3-letter code not enough
			uint const link_pos( pose.glycan_tree_set()->get_linkage_position( i ) );
			red_end_short_name[ 2 ] = '0' + link_pos;  // Set the correct connectivity.
		} else {
			red_end_short_name = red_end_res.name3();
		}
		TR.Trace << non_red_end_short_name << "(?-" << red_end_short_name << " linkage: ";
		if ( chemical::carbohydrates::CarbohydrateInfoManager::pair_has_linkage_statistics(
				red_end_short_name, non_red_end_short_name ) ) {
			TR.Trace << "found" << endl;
			// We will assume that the 1st conformer has the highest population and use that.
			chemical::carbohydrates::LinkageConformerData const conformer(
				chemical::carbohydrates::CarbohydrateInfoManager::linkages_from_pair(
				red_end_short_name, non_red_end_short_name )[ 1 ] );
			carbohydrates::set_dihedrals_from_linkage_conformer_data(
				pose, i, conformer );
		} else {
			TR.Trace << "not found" << endl;
			// TODO: Idealize from the sugar_bb scoring term.
		}
	}
}


// Set the dihedral angles involved in a glycosidic linkage based on statistical data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
void set_dihedrals_from_linkage_conformer_data( Pose & pose,
	uint const upper_residue,
	core::chemical::carbohydrates::LinkageConformerData const & conformer,
	bool idealize /* true */,
	bool use_prob_for_sd /* false */ )
{
	using namespace core::id;

	core::Size const n_torsions = conformer.n_torsions();

	if ( idealize ) {
		for ( core::Size i = 1; i <= n_torsions; ++i ) {
			TR.Debug << "torsion "<< i << std::endl;
			Angle torsion_mean;
			if ( i <= 3 ) {  // phi, psi, or omega1
				MainchainTorsionType torsion_type = static_cast< MainchainTorsionType >( i );
				torsion_mean = basic::periodic_range( conformer.get_torsion_mean( torsion_type ), 360);
			} else {
				torsion_mean = basic::periodic_range( conformer.get_torsion_mean( omega_dihedral, i - 2 ), 360);
			}
			TR.Debug << "Setting from " << get_glycosidic_torsion( i, pose, upper_residue ) <<
				" to " << torsion_mean << std::endl;  // TEMP
			set_glycosidic_torsion( i, pose, upper_residue, torsion_mean );
		}
	} else if ( use_prob_for_sd ) {
		TR << "This is not yet implemented." << std::endl;

	} else {
		Real mean;
		Real sd;
		Real conformer_angle;
		for ( core::Size i = 1; i <= n_torsions; ++i ) {
			TR.Debug << "torsion: " << i << std::endl;


			if ( i <= 3 ) { // phi, psi, or omega1
				MainchainTorsionType torsion_type = static_cast< MainchainTorsionType >( i );
				mean = conformer.get_torsion_mean( torsion_type );
				sd = conformer.get_torsion_sd( torsion_type );
				conformer_angle = basic::periodic_range( mean - sd + numeric::random::rg().uniform() * sd * 2, 360.0 );
			} else {
				mean = conformer.get_torsion_mean( omega_dihedral, i - 2 );
				sd = conformer.get_torsion_sd( omega_dihedral, i - 2);
				conformer_angle = basic::periodic_range( mean - sd + numeric::random::rg().uniform() * sd * 2, 360.0 );
			}

			//Debugging for now.
			TR.Debug << "Current: " << get_glycosidic_torsion(i, pose, upper_residue) << std::endl;
			TR.Debug << "Mean: " << mean << std::endl;
			TR.Debug << "New conformer: " << conformer_angle << std::endl;

			set_glycosidic_torsion( i, pose, upper_residue, conformer_angle );
		}
	}
}


/// @brief Remove branch points from a carbohydrate or aa residue.
void
remove_carbohydrate_branch_point_variants( Pose & pose, core::Size const seqpos)
{
	using namespace core::chemical;

	if ( pose.residue(seqpos).is_branch_point() ) {
		if ( pose.residue(seqpos).is_carbohydrate() ) {
			utility::vector1< VariantType > branch_types;

			branch_types.push_back( C1_BRANCH_POINT );
			branch_types.push_back( C2_BRANCH_POINT );
			branch_types.push_back( C3_BRANCH_POINT );
			branch_types.push_back( C4_BRANCH_POINT );
			branch_types.push_back( C5_BRANCH_POINT );
			branch_types.push_back( C6_BRANCH_POINT );
			branch_types.push_back( C7_BRANCH_POINT );
			branch_types.push_back( C8_BRANCH_POINT );
			branch_types.push_back( C9_BRANCH_POINT );


			for ( core::Size i = 1; i <= branch_types.size(); ++i ) {

				//Must check current residue since we are editing - something gets strange with valgrind
				//   since we replace the residue on variant type removal
				if ( pose.residue(seqpos).has_variant_type( branch_types[ i ] ) ) {
					remove_variant_type_from_pose_residue( pose, branch_types[ i ], seqpos);
				}
			}

		} else {
			remove_variant_type_from_pose_residue( pose, SC_BRANCH_POINT, seqpos);
		}
	}
}

/// @brief find connected atom in neighboring residue
// pair: <residue number,main chain atom number
void 
find_neighbor( Pose &pose, core::Size res_pos, core::Size atom_nr, std::pair<core::Size, core::Size> &neighbor )
{
	// TODO: set proper bond length cut offs
	core::Real max_cutoff = 3.0;
	core::Real min_cutoff = 1.0;
	for ( core::Size i = res_pos + 1; i <= res_pos + 20; ++res_pos ) {
		conformation::ResidueOP res( new conformation::Residue( pose.residue( i ) ) );
		chemical::ResidueTypeOP RT( new chemical::ResidueType( pose.residue_type( i ) ) );
		// Loop through main chain atoms of potential neighbor residue
		for ( utility::vector1< core::Size >::const_iterator mc_atom = RT->mainchain_atoms().begin(); mc_atom != RT->mainchain_atoms().end(); ++mc_atom ){
			core::Real distance = pose.residue(res_pos).xyz( atom_nr ).distance(pose.residue( i ).xyz( *mc_atom ));
			if ( (distance < max_cutoff) | (distance > min_cutoff) ){
				neighbor = std::make_pair(i,*mc_atom);
			}
		}
	}
	// Apparently no neighbor was found
	
}

/// @brief Go through pose from non-Rosetta PDB file and infer and properly set up branching
void
setup_existing_glycans(Pose &pose)
{
	
	//chain_map = chainID, residue number //
	map< core::Size/*chainID*/,utility::vector1<core::Size>  > chain_map;
	for ( core::Size res_pos = 1; res_pos <= pose.total_residue(); ++res_pos ){
		conformation::ResidueOP res( new conformation::Residue( pose.residue( res_pos ) ) );
		if ( res->is_carbohydrate() ) {
			//core::Size chain = res->chain();
			chemical::ResidueTypeOP RT( new chemical::ResidueType( pose.residue_type(res_pos) ) );
			for ( utility::vector1< Size >::const_iterator atom = RT->mainchain_atoms().begin(); atom != RT->mainchain_atoms().end(); ++atom){
				// Find the fisrt exocyclic oxigen
				std::string atom_name = RT->atom_name(*atom);
				if ( atom_name.find('O') != std::string::npos ){  // is oxigen
					std::pair< core::Size,core::Size  > neighbor = std::make_pair(0,0);
					find_neighbor(pose,res_pos,*atom,neighbor);
					if ( neighbor.first == 0 || neighbor.second == 0 ) {
						// No neighbor neighbor found
						// Must be a lower terminus
						// TODO: add terminus variant type here
					} else {
						// Main chain continues at residue neighbor.first
						// Add branch variant type
						chemical::VariantType const variant(
							chemical::carbohydrates::CarbohydrateInfoManager::branch_variant_type_from_atom_name( atom_name ) );
						add_variant_type_to_pose_residue( pose, variant, res_pos );
					}
						
				}		
			}
		}
	}
	PDBInfoOP info( new PDBInfo( pose ) );
	info->name( pose.sequence() );  // Use the sequence as the default name.
	pose.pdb_info( info );
}
		




////////////////////////////////////////    Branch Information    /////////////////////////////////////////
///
///
///








/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_residues_of_branch(
	Pose const & pose,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(pose, starting_position).first;
}
	


/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_tips_of_branch(
	Pose const & pose,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(pose, starting_position).second;
}
	


/// @brief Get residues further down the branch from this residue.  starting_position ->
///  May require a better name.
///  Returns pair of all_upstream_residues, tips.
///  Tips are the ends of linear glycan branches.

std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	Pose const & pose,
	uint const starting_position)

{
	return conformation::carbohydrates::get_carbohydrate_residues_and_tips_of_branch(pose.conformation(), starting_position);
}





///@brief Get the carbohydrate residue connecting the protein branch point.
core::Size
get_glycan_connecting_protein_branch_point(Pose const & pose, core::Size const protein_branch_point_resnum){

	return conformation::carbohydrates::get_glycan_connecting_protein_branch_point( pose.conformation(), protein_branch_point_resnum );
	
}


///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
///
core::Size
get_resnum_from_glycan_position(Pose const & pose, core::Size const glycan_one, core::Size const glycan_position){
	return conformation::carbohydrates::get_resnum_from_glycan_position( pose.conformation(), glycan_one, glycan_position);
	
}

core::Size
get_glycan_position_from_resnum(Pose const & pose, core::Size const glycan_one, core::Size const glycan_residue ){
	
	return conformation::carbohydrates::get_glycan_position_from_resnum( pose.conformation(), glycan_one, glycan_residue );
	
}



select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(Pose const & pose, core::Size const glycan_one, utility::vector1< core::Size > const & glycan_positions){
	
	utility::vector1< bool > subset(pose.total_residue(), false);
	
	core::Size glycan_length = pose.glycan_tree_set()->get_tree(glycan_one)->get_size();
	for (core::Size index = 1; index <= glycan_positions.size(); ++index){
		
		Size glycan_position = glycan_positions[ index ];
		if (glycan_position <= glycan_length){
	
			Size glycan_residue = get_resnum_from_glycan_position(pose, glycan_one, glycan_position);
			
			if (glycan_residue <= pose.total_residue()){
				subset[ glycan_residue ] = true;
			}
			else {
				utility_exit_with_message("Glycan residue not in pose: "+utility::to_string( glycan_residue) + " corresponding to glycan position "+utility::to_string( glycan_position));
			}

		}
	}
	return subset;
}

select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(Pose const & pose, utility::vector1< core::Size > const & glycan_positions){
	select::residue_selector::OrResidueSelector combine_subsets = select::residue_selector::OrResidueSelector();
	
	utility::vector1< bool > subset(pose.total_residue(), false);
	utility::vector1< bool > glycan_start_points =  pose.glycan_tree_set()->get_start_points();
	
	for (Size i = 1; i <= glycan_start_points.size(); ++i){
		if (! glycan_start_points[ i ]) continue;
		
		select::residue_selector::ResidueSubset single_glycan = get_resnums_from_glycan_positions(pose, i, glycan_positions);
		combine_subsets.apply_or_to_subset(single_glycan, subset);
	}
	return subset;
}


////////////////////////////////////////    Branch Deletion    /////////////////////////////////////////
///
///
///


/// @brief Delete the glycan from this residue onward toward the end of the branch.  Like chopping off a tree trunk at position resnum (not including resnum). Also known as defoliating.
///  If resnum is the protein branch point, will change variant.
//   If no more carbohydrates exist in the pose, will change the pose status.
void
delete_carbohydrate_branch( Pose & pose, uint const delete_to )
{
	using namespace utility;
	using namespace chemical;
	using namespace conformation;

	if ( ! pose.residue( delete_to ).is_carbohydrate() && ! pose.residue( delete_to ).is_branch_point() ) {
		TR << "Delete to residue is not carbohydrate and not a branch point.  Nothing to be done." << std::endl;
		return;
	}
	
	TR << "Delete Carb Branch Starting Points " << pose.glycan_tree_set()->get_start_points() << std::endl;
	
	std::pair< vector1< Size >, vector1< Size > > resnums_and_tips;
	resnums_and_tips = get_carbohydrate_residues_and_tips_of_branch( pose, delete_to );

	vector1< Size > resnums_to_delete = resnums_and_tips.first;
	vector1< Size > tips              = resnums_and_tips.second;

	Size i = 0;
	std::string ref_pose_name = "temp_ref_pose" ;

	Size current_delete_to = delete_to;
	
	while ( tips.size() > 0 ) {

		pose.reference_pose_from_current( ref_pose_name, true );

		vector1< vector1< Size > > leafs;

		//Get linear carbohydrate leafs
		for ( core::Size x = 1; x <= tips.size(); ++x ) {
			Size tip = tips[ x ];
			vector1< Size > leaf = get_resnums_in_leaf_on_the_fly(pose, tip, delete_to);
			leafs.push_back(leaf);
		}

		//Delete the leafs
		for ( core::Size x = 1; x <= leafs.size(); ++x ) {
			vector1< Size > leaf = leafs[ x ];
			delete_leaf( pose, leaf, ref_pose_name);
		}

		//Get new tips to delete.
		current_delete_to = pose.corresponding_residue_in_current( current_delete_to, ref_pose_name);
		resnums_and_tips = get_carbohydrate_residues_and_tips_of_branch( pose, current_delete_to );

		resnums_to_delete = resnums_and_tips.first;
		tips              = resnums_and_tips.second;


		//TR << pose << std::endl; //Debugging for now.

		i+=1;
	}

	// Update branching of starting_position residue.  After all this deletion, it will not have any branches and should be a tip.
	//   -> Deal with branch of aa side chain or carbohydrate branch properly!

	remove_carbohydrate_branch_point_variants( pose, current_delete_to );


	// Update contents of Pose - does it have any carbohydrates left?
	bool found_carbohydrate = false;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).is_carbohydrate() ) {
			found_carbohydrate = true;
			break;
		}
	}
	pose.conformation().contains_carbohydrate_residues( found_carbohydrate );
}

/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
get_branching_residues( Pose const & pose,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	return conformation::carbohydrates::get_branching_residues( pose.conformation(), parent_residue, children_residues, list_of_residues, tips );
}


/// @brief  Find all children residues, list of residues, and any found tips from a given residue not including parent
///
///  Children Residues:  Filled in list of children residues found if not tips.
///  List Of  Residues:  All the residue nums found.
///  Tips:  All 'ends' of of children found.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
fill_upstream_children_res_and_tips( Pose const & pose,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	return conformation::carbohydrates::fill_upstream_children_res_and_tips(pose.conformation(), res, parent_residue, children_residues, list_of_residues, tips);
}


/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
utility::vector1< core::Size >
get_resnums_in_leaf( Pose const & pose, Size tip_residue, Size stop_at_residue )
{
	// This logic can be simplified even further, but this works.
	
	utility::vector1< Size > resnums;

	//If residue is not carbohydrate, we can't find the parent residue.
	/// We really shouldn't be here in the first place!
	if ( ! pose.residue( tip_residue ).is_carbohydrate() ) {
		return resnums;
	}

	resnums.push_back( tip_residue );
	Size res = tip_residue;

	while ( pose.residue( res ).n_current_residue_connections() < 3 ) {
		
		Size parent = pose.glycan_tree_set()->get_parent( res );

		if ( res != tip_residue ) { resnums.push_back( res ); }
		if ( parent == 0 || parent == stop_at_residue ) { break; }

		res = parent;
	}
	return resnums;
}

/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
utility::vector1< core::Size >
get_resnums_in_leaf_on_the_fly( Pose const & pose, Size tip_residue, Size stop_at_residue )
{
	// This logic can be simplified even further, but this works.

	utility::vector1< Size > resnums;

	//If residue is not carbohydrate, we can't find the parent residue.
	/// We really shouldn't be here in the first place!
	if ( ! pose.residue( tip_residue ).is_carbohydrate() ) {
		return resnums;
	}

	resnums.push_back( tip_residue );
	Size res = tip_residue;

	while ( pose.residue( res ).n_current_residue_connections() < 3 ) {

		Size parent = find_seqpos_of_saccharides_parent_residue( pose.residue( res ) );

		if ( res != tip_residue ) { resnums.push_back( res ); }
		if ( parent == 0 || parent == stop_at_residue ) { break; }

		res = parent;
	}
	return resnums;
}

/// @brief Delete a leaf of glycan residues.  Use the ReferencePose to associate residue numbers.
/// @details
///  Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
///  This is a bug in the FoldTree where chemical edges are being treated as Jumps.
///  This is being addressed in a branch and a FoldTree.
void
delete_leaf( Pose & pose, utility::vector1< Size > leaf, std::string ref_pose_name )
{
	for ( core::Size i = 1; i <= leaf.size(); ++i ) {
		Size resnum = leaf[ i ];
		Size new_resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name);

		// Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
		// This exit checks for foldtree jump, and it looks like it counts chemical edges as a foldtree jump.
		// UNTIL WE CAN RESOLVE THIS, use delete_residue_slow here.

		//pose.delete_polymer_residue( new_resnum );
		pose.delete_residue_slow( new_resnum );
	}
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
