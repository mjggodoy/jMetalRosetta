// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/bond_support.hh
/// @brief support functions for class Bond; functions that
/// should not be included as part of the class.
/// @author Steven Combs
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/bond_support.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

#include <utility/graph/RingDetection.hh>
#include <utility/graph/ring_detection.hh>

#include <basic/Tracer.hh>

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.bond_support" );

/// @brief convert bond order or aromatic into the corresponding radius
/// @param BOND_ORDER_OR_AROMATIC bond type in notation: 1=single, 2=double, 3=triple, 4=aromatic
gasteiger::GasteigerAtomTypeData::Properties bond_order_to_property( const core::Size &BOND_ORDER_OR_AROMATIC)
{
	static gasteiger::GasteigerAtomTypeData::Properties properties[ 5] =
		{
		gasteiger::GasteigerAtomTypeData::VdWaalsRadiusCSD,
		gasteiger::GasteigerAtomTypeData::CovalentRadiusSingleBond,
		gasteiger::GasteigerAtomTypeData::CovalentRadiusDoubleBond,
		gasteiger::GasteigerAtomTypeData::CovalentRadiusTripleBond,
		gasteiger::GasteigerAtomTypeData::CovalentRadiusAromaticBond
		};
	if ( BOND_ORDER_OR_AROMATIC > 5 ) {
		utility_exit_with_message("Bond order must be between 1 and 4 (4 = aromatic)");
	}
	return properties[ BOND_ORDER_OR_AROMATIC];
}

void find_bonds_in_rings(ResidueType & res, bool const complex_ring){
	//first, we assign all the bonds in the residue to having no rings
	EIter edge_begin, edge_end;
	boost::tie(edge_begin, edge_end) = boost::edges( res.graph() );
	for ( EIter edge_iter = edge_begin; edge_iter != edge_end; ++edge_iter ) {
		Bond & bond = res.bond(*edge_iter);
		//for now, nothing is in the ring. Defaulted to not known when constructed, after ring detection, either
		//in a ring or not in a ring
		bond.ringness( BondNotInRing );
	}

	complex_ring ? complex_ring_detection( res) : quick_ring_detection( res);

}
void
complex_ring_detection( ResidueType & res){

	//first we get the light weight residue graph
	LightWeightResidueGraph lwrg = convert_residuetype_to_light_graph(res);
	//now get the property maps
	boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	//boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);

	//for ring systems that are complex, you only want to annotate edges that are in the given system.
	//if you do not only annotate the edges, there is a combinatorial explosion, which results in
	//lots of times enumerating all rings in the system.
	utility::graph::RingDetection<LightWeightResidueGraph> ring_detect( lwrg); //initialize the ring detector. Automatically assigns rings


	utility::vector1<utility::vector1<lwrg_VD> > rings = ring_detect.GetRings(); //these are the path of the rings

	//iterate through the rings, then assign the bonds for ringness
	for ( core::Size i=1; i<= rings.size(); ++i ) {
		utility::vector1< ED > just_the_edges;
		for ( core::Size j=1; j<  rings[i].size(); ++j ) { //not less than, see explanation below
			VD source =  lwrg_vd_to_VD[ rings[i][j] ];
			//next set of code is a little convulted. The ring code returns all the vertex, but we need
			//the edges to assign (bond), not the vertex. (atom, at least this point in time).
			//Therefore, we have to move one past the vector to get the edge.
			VD target = lwrg_vd_to_VD[ rings[i][ j+1] ];
			ED edge;
			bool edge_exists;
			boost::tie(edge, edge_exists) = boost::edge(source, target, res.graph());
			if ( edge_exists ) { //if there is an edge, mark it as being a ring
				just_the_edges.push_back(edge); //get the edge
				Bond & bond = res.bond(edge);
				bond.ringness(BondInRing);
			} else {
				utility_exit_with_message("In ring detection, cannot find bond for " + res.atom_name( source ) + " to " + res.atom_name( target ) );
			}
		}
	}
}

void quick_ring_detection( ResidueType & res){
	std::map< VD, std::map<VD, bool > >  ring_edges( utility::graph::annotate_ring_edges( res.graph()) );
	//std::map< VD,std::map<VD, bool > >::const_iterator it_start, it_end;
	//it_start = ring_edges.begin();
	//it_end = ring_edges.end();
	for ( auto const & vd_ring_edge_elem : ring_edges ) {
		for ( auto const & vd_bool_elem : vd_ring_edge_elem.second ) {
			if ( ! vd_bool_elem.second ) continue;

			ED bond_edge;
			bool edge_exists;
			boost::tie( bond_edge, edge_exists) = boost::edge( vd_ring_edge_elem.first, vd_bool_elem.first, res.graph());
			if ( edge_exists ) {
				Bond & bond = res.bond( bond_edge);
				bond.ringness( BondInRing);
			} else {
				utility_exit_with_message("In quick ring detection, cannot find bond for " + res.atom_name( vd_ring_edge_elem.first ) + " to " + res.atom_name( vd_bool_elem.first ) );
			}
		}
	}
}

utility::vector1<VD> get_connecting_atoms(ResidueType const & res, ED const & edge) {
	return get_connecting_atoms(res.graph(), edge);
}

utility::vector1<VD> get_connecting_atoms(ResidueGraph const & graph, ED const & edge){
	utility::vector1<VD> connecting_atoms;
	connecting_atoms.push_back( boost::source( edge, graph ) );
	connecting_atoms.push_back( boost::target( edge, graph ) );
	return connecting_atoms;
}


ED get_bond(ResidueType const & res, VD const & source, VD const & target){
	bool bond_there(false);
	ED edge;
	boost::tie(edge,bond_there) = boost::edge(source, target, res.graph() );
	debug_assert(bond_there);
	return edge;
}

Real create_bond_length(
	gasteiger::GasteigerAtomTypeData const & atom1,
	gasteiger::GasteigerAtomTypeData const & atom2,
	BondName bond_type)
{
	return atom1.get_atom_type_property( bond_order_to_property(bond_type) ) + atom2.get_atom_type_property( bond_order_to_property(bond_type));
}

/// @brief Find which bonds are rotatatable (chi) bonds
/// Returns a list of four vds representing the chi
/// @details Assumes:
///
/// * Complete Atom/bond graph
/// * All element types have been set
/// * All bond_names have been set.
/// * Ringness has been set for all bonds in rings
/// * The Icoor graph has been set.

utility::vector1<VDs> find_chi_bonds( ResidueType const & restype ) {
	using namespace core::chemical;
	utility::vector1<VDs> found_chis;

	// Step 1) Find all rotatable bonds
	core::chemical::EIter eiter, eiter_end;
	for ( boost::tie(eiter, eiter_end) = restype.bond_iterators(); eiter != eiter_end; ++eiter ) {
		Bond const & bond( restype.bond( *eiter ) );
		VD source(boost::source(*eiter,restype.graph()));
		VD target(boost::target(*eiter,restype.graph()));
		// Step 1.1) Bonds which aren't single, which are in rings,
		//     or which are terminal/to hydrogen aren't rotatable.
		if ( bond.bond_name() != SingleBond || // Should this be bond order instead?
				bond.ringness() == BondInRing ||
				restype.atom(source).element_type()->element() == element::H  ||
				restype.atom(target).element_type()->element() == element::H ||
				boost::out_degree(source,restype.graph()) == 1 || boost::out_degree(target,restype.graph()) == 1 ) {
			continue;
		}
		// Step 1.2) Chi bonds need to be oriented along the ICOOR tree. Swap the order if they're reversed.
		if ( restype.atom_base(source) == target &&
				source != restype.root_atom() ) { // Root atom has atom_base as it's child atom - don't swap there.
			// Swap target and source such that source is nearer root.
			std::swap( source, target );
		} else if ( restype.atom_base(target) != source ) {
			TR << "Found non-tree bond " << restype.atom_name(source) << " --- " << restype.atom_name(target) << std::endl;
			TR << "\t   Expected tree bond " << restype.atom_name( restype.atom_base(target) ) << " --- " << restype.atom_name(target) << std::endl;
			utility_exit_with_message("Error: Non-ring bond not found in ResidueType atom tree.");
		}

		// Step 1.3a) Validity of rotatable bonds also depends on what's attached to them.
		//    Specifically, the number of heavy atoms and hydrogens
		//    For efficiency, we also notate selected heavy atoms and hydrogens for later use.
		core::Size targ_heavy(0), targ_hydro(0);
		VD selected_targ_heavy( boost::graph_traits<ResidueGraph>::null_vertex() );
		VD last_targ_hydro( boost::graph_traits<ResidueGraph>::null_vertex() );
		ResidueGraph::adjacency_iterator aiter, aiter_end;
		for ( boost::tie(aiter, aiter_end) = boost::adjacent_vertices(target, restype.graph()); aiter != aiter_end; ++aiter ) {
			if ( restype.atom(*aiter).element_type()->element() == element::H ) {
				++targ_hydro;
				last_targ_hydro = *aiter;
			} else {
				++targ_heavy;
				// Residue::set_chi() requires that chi angles have the center and "target" bonds be along the along the atom base graph.
				// We also don't want to pick the source atom, or the root atom (again, because the atom_base of the root atom is it's child).
				if ( selected_targ_heavy == boost::graph_traits<ResidueGraph>::null_vertex() &&
						restype.atom_base(*aiter) == target &&
						*aiter != source && *aiter != restype.root_atom() ) {
					selected_targ_heavy = *aiter;
				}
			}
		}
		// Step 1.3b) For rotatable bonds, we want at least one other (not the source) heavy atom on the target, or we want a single, non-carbon-attached hydrogen
		element::Elements const & target_element( restype.atom(target).element_type()->element() );
		if ( targ_heavy < 2 && (targ_hydro != 1 || target_element == element::C) ) {
			continue;
		}

		// Step 2) Pick the other two atoms which will make up the chi
		VD d, c(source), b(target), a;
		// Step 2.1) Pick the terminal (built/distal) atom for the chi
		if ( targ_heavy >= 2 ) {
			// Step 2.1a) A regular all-heavy atom chi: use the selected heavy atom.
			if ( selected_targ_heavy == boost::graph_traits<ResidueGraph>::null_vertex() ) {
				TR.Error << "Issues with finding chi angles with the ICOOR graph!" << std::endl;
				TR.Error << "    Chi bond for " << restype.atom_name( source ) << " -- " << restype.atom_name( target ) << std::endl;
				utility_exit_with_message("Can't autodetermine chi angles.");
			}
			a = selected_targ_heavy;
		} else {
			// Step 2.1b) Proton chi
			debug_assert( targ_hydro == 1 && target_element != element::C );
			debug_assert( last_targ_hydro != boost::graph_traits<ResidueGraph>::null_vertex() );
			a = last_targ_hydro;
		}

		// Step 1.3c) If the terminal bond is co-linear with the central bond
		//    (e.g. it's a triple bond) the chi isn't rotatable.
		if ( restype.bond( b, a ).order() ==  TripleBondOrder ) {
			continue;
		}

		// Step 2.2) Pick the dihedral reference atom.
		if ( source != restype.root_atom() ) {
			// Step 2.2a) Source isn't the root atom - we can (and should) use the atom_base as reference
			d = restype.atom_base(source);
		} else if ( restype.atom_base(source) != target ) {
			// Step 2.2b) Source is root atom, but the atom_base isn't part of the bond
			d = restype.atom_base(source);
		} else {
			// Step 2.2b) Source is root atom, and atom_base is part of the bond.
			// We need to find another heavy atom to base things off of.
			d = boost::graph_traits<ResidueGraph>::null_vertex();
			ResidueGraph::adjacency_iterator aiter2, aiter2_end;
			for ( boost::tie(aiter2, aiter2_end) = boost::adjacent_vertices(source, restype.graph()); aiter2 != aiter2_end; ++aiter2 ) {
				if ( *aiter2 != target && restype.atom(*aiter2).element_type()->element() != element::H ) {
					d = *aiter2;
					break;
				}
			}
			if ( d == boost::graph_traits<ResidueGraph>::null_vertex() ) {
				// Can't find a suitable reference atom - this probably isn't a rotatable chi after all.
				continue;
			}
		}
		// Step 2 + 1i) Chis propagate out from the root, and the distal two bonds follow the ICOOR tree.
		// Double check for sanity.
		debug_assert( a != restype.root_atom() );
		debug_assert( b != restype.root_atom() );
		debug_assert( restype.atom_base( a ) == b );
		debug_assert( restype.atom_base( b ) == c );
		// Step 3) Valid chi bond
		VDs chi;
		chi.push_back(d);
		chi.push_back(c);
		chi.push_back(b);
		chi.push_back(a);
		TR.Debug << "Found chi: " << restype.atom_name(d) << " --- " << restype.atom_name(c) << " --- " <<  restype.atom_name(b) << " --- " << restype.atom_name(a) << std::endl;
		found_chis.push_back( chi );
	} // For all edges
	return found_chis;
}

/// @brief Is the given chi a proton chi with the proton attached to an atom attached to an non-sp3 atom?
/// @details The use case is to see if the proton chi should flat or staggered with rotamers
bool is_sp2_proton_chi( core::Size chi, ResidueType const & restype ) {
	VDs atoms( restype.chi_atom_vds(chi) );
	debug_assert( atoms.size() == 4 );
	// Note this is used in setting up proton chis, so we can't assume is_proton_chi() and associated are valid yet.
	VD gp( atoms[2] );
	OutEdgeIter iter, iter_end;
	for ( boost::tie(iter, iter_end) = restype.bond_iterators(gp); iter != iter_end; ++iter ) {
		core::chemical::BondName bt( restype.bond(*iter).bond_name() );
		if ( bt == DoubleBond || bt == SingleBond || bt == AromaticBond ) {
			return true;
		}
	}
	return false;
}

}
}
