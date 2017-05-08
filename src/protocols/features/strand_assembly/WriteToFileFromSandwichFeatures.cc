// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/WriteToFileFromSandwichFeatures.cc
/// @brief Write to a file after SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)
/// @overview
///  Write AA distribution


//Devel
#include <protocols/features/strand_assembly/WriteToFileFromSandwichFeatures.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

static THREAD_LOCAL basic::Tracer TR( "protocols.features.strand_assembly.WriteToFileFromSandwichFeatures" );

//get_vector_of_loop_AA_distribution
utility::vector1<Size>
get_vector_of_loop_AA_distribution (
	StructureID struct_id, // needed argument
	sessionOP db_session, // needed argument
	string loop_kind // like 'hairpin' or 'inter-sheet-loop'
)
{
	string sum_string =
		"SELECT\n"
		"\tsum(A), sum(C), sum(D), sum(E), sum(F), \n"
		"\tsum(G), sum(H), sum(I), sum(K), sum(L), \n"
		"\tsum(M), sum(N), sum(P), sum(Q), sum(R), \n"
		"\tsum(S), sum(T), sum(V), sum(W), sum(Y) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tloop_kind = ? \n"
		"\tAND struct_id = ? ;";

	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1, loop_kind);
	sum_statement.bind(2, struct_id);
	result res(basic::database::safely_read_from_database(sum_statement));

	Size num_A, num_C, num_D, num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y;

	utility::vector1<Size> vector_loop_AA_distribution;

	while ( res.next() )
			{
		res >> num_A >>  num_C >>  num_D >>  num_E >>  num_F >>  num_G >>  num_H >>  num_I >>  num_K >>  num_L >>  num_M >>  num_N >>  num_P >>  num_Q >>  num_R >>  num_S >>  num_T >>  num_V >>  num_W >>  num_Y;
		vector_loop_AA_distribution.push_back(num_A);
		vector_loop_AA_distribution.push_back(num_C);
		vector_loop_AA_distribution.push_back(num_D);
		vector_loop_AA_distribution.push_back(num_E);
		vector_loop_AA_distribution.push_back(num_F);

		vector_loop_AA_distribution.push_back(num_G);
		vector_loop_AA_distribution.push_back(num_H);
		vector_loop_AA_distribution.push_back(num_I);
		vector_loop_AA_distribution.push_back(num_K);
		vector_loop_AA_distribution.push_back(num_L);

		vector_loop_AA_distribution.push_back(num_M);
		vector_loop_AA_distribution.push_back(num_N);
		vector_loop_AA_distribution.push_back(num_P);
		vector_loop_AA_distribution.push_back(num_Q);
		vector_loop_AA_distribution.push_back(num_R);

		vector_loop_AA_distribution.push_back(num_S);
		vector_loop_AA_distribution.push_back(num_T);
		vector_loop_AA_distribution.push_back(num_V);
		vector_loop_AA_distribution.push_back(num_W);
		vector_loop_AA_distribution.push_back(num_Y);
	}

	return vector_loop_AA_distribution;
} // get_vector_of_loop_AA_distribution


//prepare_to_write_number_of_electrostatic_interactions_of_residues_to_files
core::Size
prepare_and_write_number_of_electrostatic_interactions_of_residues_to_files(
	string tag,
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	Real distance_cutoff_for_electrostatic_interactions_,
	Real CB_b_factor_cutoff_for_electrostatic_interactions_,
	Size min_primary_seq_distance_diff_for_electrostatic_interactions_,
	bool write_electrostatic_interactions_of_surface_residues_in_a_strand_,
	bool write_electrostatic_interactions_of_all_residues_in_a_strand_,
	bool write_electrostatic_interactions_of_all_residues_,
	Size rkde_PK_id_counter,
	Size rkde_in_strands_PK_id_counter)
{
	if ( write_electrostatic_interactions_of_surface_residues_in_a_strand_ || write_electrostatic_interactions_of_all_residues_in_a_strand_ ) {
		// <begin> store RKDE in strands to a database table
		for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
			Size residue_begin = bs_of_sw_can_by_sh[ii].get_start();
			Size residue_end = bs_of_sw_can_by_sh[ii].get_end();
			for ( Size residue_num = residue_begin; residue_num <= residue_end; residue_num++ ) {
				if (
						(pose.residue_type(residue_num).name3() == "ARG")
						|| (pose.residue_type(residue_num).name3() == "LYS")
						|| (pose.residue_type(residue_num).name3() == "ASP")
						|| (pose.residue_type(residue_num).name3() == "GLU")
						) {
					string heading =
						determine_heading_direction_by_vector (
						struct_id,
						db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), residue_begin, residue_end,
						residue_num);
					WriteToDB_rkde_in_strands(
						struct_id,
						db_session,
						rkde_in_strands_PK_id_counter,
						tag,
						bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), //sw_can_by_sh_id
						residue_num, //residue_number,
						pose.residue_type(residue_num).name3(), //residue_type,
						heading // "surface" or "core"
					);
					rkde_in_strands_PK_id_counter++;
				}
			}
		}
		// <end> store RKDE in strands to a database table
	}

	if ( write_electrostatic_interactions_of_surface_residues_in_a_strand_ ) {
		// <begin> report number of electrostatic_interactions_of_surface_residues_in_a_strand
		write_number_of_electrostatic_interactions_of_residues_to_files (
			tag,
			struct_id,
			db_session,
			pose,
			"E",
			"surface",
			distance_cutoff_for_electrostatic_interactions_,
			CB_b_factor_cutoff_for_electrostatic_interactions_,
			min_primary_seq_distance_diff_for_electrostatic_interactions_);
		// <end> report number of electrostatic_interactions_of_surface_residues_in_a_strand
	}

	if ( write_electrostatic_interactions_of_all_residues_in_a_strand_ ) {
		// <begin> report number of electrostatic_interactions_of_all_residues_in_a_strand
		write_number_of_electrostatic_interactions_of_residues_to_files (
			tag,
			struct_id,
			db_session,
			pose,
			"E",
			"all_direction",
			distance_cutoff_for_electrostatic_interactions_,
			CB_b_factor_cutoff_for_electrostatic_interactions_,
			min_primary_seq_distance_diff_for_electrostatic_interactions_);
		// <end> report number of electrostatic_interactions_of_all_residues_in_a_strand
	}

	if ( write_electrostatic_interactions_of_all_residues_ ) {
		// <begin> store RKDE to a database table
		for ( Size ii=1; ii<=pose.size(); ii++ ) {
			if (
					(pose.residue_type(ii).name3() == "ARG")
					|| (pose.residue_type(ii).name3() == "LYS")
					|| (pose.residue_type(ii).name3() == "ASP")
					|| (pose.residue_type(ii).name3() == "GLU")
					) {
				WriteToDB_rkde(
					struct_id,
					db_session,
					rkde_PK_id_counter,
					tag,
					ii, //residue_number,
					pose.residue_type(ii).name3() //residue_type,
				);
				rkde_PK_id_counter++;
			}
		}
		// <end> store RKDE to a database table

		// <begin> report number of electrostatic_interactions_of_all_residues
		write_number_of_electrostatic_interactions_of_residues_to_files (
			tag,
			struct_id,
			db_session,
			pose,
			"all_dssp",
			"all_direction",
			distance_cutoff_for_electrostatic_interactions_,
			CB_b_factor_cutoff_for_electrostatic_interactions_,
			min_primary_seq_distance_diff_for_electrostatic_interactions_);
		// <end> report number of electrostatic_interactions_of_all_residues

	}

	return 0;
} //prepare_to_write_number_of_electrostatic_interactions_of_residues_to_files

//write_AA_kind_to_a_file
core::Size
write_AA_kind_to_a_file(
	string tag,
	StructureID struct_id, // needed argument
	utility::sql_database::sessionOP db_session, // needed argument
	Size sw_can_by_sh_id,
	Size sw_res_size)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	string AA_kind_file_name = pdb_file_name + "_AA_kind.txt";
	ofstream AA_kind_file;

	AA_kind_file.open(AA_kind_file_name.c_str());

	utility::vector1<Size> vec_AA_kind = get_vec_AA_kind(struct_id, db_session,
		sw_can_by_sh_id // used to be vec_sw_can_by_sh_id[ii]
	);

	// positive, negative, polar, aromatic, hydrophobic
	AA_kind_file << "pdb_name\tPos_percent\tNeg_percent\tPol_percent\tAro_percent\tHydropho_percent\t" ; // as Jenny's thesis
	AA_kind_file << "Pos_raw_count\tNeg_raw_count\tPol_raw_count\tAro_raw_count\tHydropho_raw_count" << endl; // to calculate net_charge_of_sw

	AA_kind_file << pdb_file_name << "\t" ;

	for ( Size i =1; i<=(vec_AA_kind.size()); i++ ) {
		Real percent = vec_AA_kind[i]*100/static_cast<Real>(sw_res_size);
		Size rounded_percent = round_to_Size(percent);
		AA_kind_file << rounded_percent << "\t" ;
	}

	for ( Size i =1; i<=(vec_AA_kind.size()); i++ ) {
		AA_kind_file << vec_AA_kind[i] << "\t" ;
	}
	AA_kind_file << endl; // to marshal
	AA_kind_file.close();

	return 0;
} //write_AA_kind_to_a_file

//write_AA_distribution_with_direction_to_a_file
core::Size
write_AA_distribution_with_direction_to_a_file(
	string tag,
	StructureID struct_id, // needed argument
	utility::sql_database::sessionOP db_session // needed argument
)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string AA_dis_file_name = pdb_file_name + "_AA_distribution_of_strands_sorted_alphabetically.txt";
	ofstream AA_dis_file;

	AA_dis_file.open(AA_dis_file_name.c_str());
	utility::vector1<Size> vec_core_heading_at_core_strand = get_vector_of_strand_AA_distribution (struct_id, db_session, "core_heading", "core");
	// struct_id, db_session, heading_direction, strand_location
	utility::vector1<Size> vec_surface_heading_at_core_strand = get_vector_of_strand_AA_distribution (struct_id, db_session, "surface_heading", "core");
	// struct_id, db_session, heading_direction, strand_location
	utility::vector1<Size> vec_core_heading_at_edge_strand = get_vector_of_strand_AA_distribution (struct_id, db_session, "core_heading", "edge");
	// struct_id, db_session, heading_direction, strand_location
	utility::vector1<Size> vec_surface_heading_at_edge_strand = get_vector_of_strand_AA_distribution (struct_id, db_session, "surface_heading", "edge");
	// struct_id, db_session, heading_direction, strand_location
	AA_dis_file << "core_heading_at_core_strand\tsurface_heading_at_core_strand\tcore_heading_at_edge_strand\tsurface_heading_at_edge_strand" << endl;
	for ( Size i =1; i<=(vec_core_heading_at_core_strand.size()); i++ ) {
		AA_dis_file << vec_core_heading_at_core_strand[i] << "\t" << vec_surface_heading_at_core_strand[i] << "\t" << vec_core_heading_at_edge_strand[i] << "\t" << vec_surface_heading_at_edge_strand[i] << endl;
	}
	AA_dis_file.close();

	return 0;
} //write_AA_distribution_with_direction_to_a_file


//write_AA_distribution_without_direction_to_a_file
core::Size
write_AA_distribution_without_direction_to_a_file(
	string tag,
	StructureID struct_id, // needed argument
	utility::sql_database::sessionOP db_session // needed argument
)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string AA_dis_file_name = pdb_file_name + "_AA_distribution_of_loops_sorted_alphabetically.txt";
	ofstream AA_dis_file;

	AA_dis_file.open(AA_dis_file_name.c_str());
	utility::vector1<Size> vector_of_hairpin_AA = get_vector_of_loop_AA_distribution (struct_id, db_session, "hairpin");
	utility::vector1<Size> vector_of_inter_sheet_loop_AA = get_vector_of_loop_AA_distribution (struct_id, db_session, "loop_connecting_two_sheets");

	AA_dis_file << "hairpin_AA\tinter_sheet_loop_AA" << endl;
	for ( Size i =1; i<=(vector_of_hairpin_AA.size()); i++ ) {
		AA_dis_file << vector_of_hairpin_AA[i] << "\t" << vector_of_inter_sheet_loop_AA[i] << endl;
	}
	AA_dis_file.close();

	return 0;
} //write_AA_distribution_without_direction_to_a_file


//write_beta_sheet_capping_info_to_a_file
core::Size
write_beta_sheet_capping_info_to_a_file(
	string tag,
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	int primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_,
	int primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_,
	int primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_,
	int primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_)
{
	int number_of_strands_where_capping_is_checked = 0;
	int number_of_strands_with_2_cappings = 0;
	int number_of_strands_with_1_capping = 0;

	utility::vector1<Size> residue_begin_of_strands_with_1_capping;
	utility::vector1<Size> residue_begin_of_strands_with_0_capping;

	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
		int residue_begin = bs_of_sw_can_by_sh[ii].get_start();
		int residue_end = bs_of_sw_can_by_sh[ii].get_end();
		if ( residue_end - residue_begin < 2 ) {
			continue; // this strand should be too short like "short_edge"
		}
		number_of_strands_where_capping_is_checked++;

		bool capping_near_begin = false;
		bool capping_near_end = false;
		for ( int jj = residue_begin - primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_;
				jj <= residue_begin + primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_;
				jj++ ) { // per each beta-strand
			if ( jj <= 0 ) {
				capping_near_begin = true; // I assume that capping_near_begin is not necessary for the first strand
				break;
			}
			if (
					(pose.residue_type(jj).name3() == "GLY")
					|| (pose.residue_type(jj).name3() == "ASP")
					|| (pose.residue_type(jj).name3() == "ASN")
					|| (pose.residue_type(jj).name3() == "PRO")
					) {
				capping_near_begin = true;
				break;
			}
		}
		for ( int jj = residue_end - primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_;
				jj <= residue_end + primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_;
				jj++ ) { // per each beta-strand
			if ( jj > static_cast<int>(pose.size()) ) {
				capping_near_end = true; // I assume that capping_near_end is not necessary for the last strand
				break;
			}
			if (
					(pose.residue_type(jj).name3() == "GLY")
					|| (pose.residue_type(jj).name3() == "ASP")
					|| (pose.residue_type(jj).name3() == "ASN")
					|| (pose.residue_type(jj).name3() == "PRO")
					) {
				capping_near_end = true;
				break;
			}
		}
		if ( capping_near_begin && capping_near_end ) {
			number_of_strands_with_2_cappings++;
		} else if ( capping_near_begin || capping_near_end ) {
			number_of_strands_with_1_capping++;
			residue_begin_of_strands_with_1_capping.push_back(residue_begin);
		} else {
			residue_begin_of_strands_with_0_capping.push_back(residue_begin);
		}
	}
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string info_file_name = pdb_file_name + "_beta_sheet_capping_info.txt";
	ofstream info_file;

	info_file.open(info_file_name.c_str());
	info_file << "number_of_strands_where_capping_is_checked\tnumber_of_strands_with_2_cappings\tnumber_of_strands_with_1_capping\t\tnumber_of_strands_with_0_capping\tresidue_begin_of_strands_with_1_capping\tresidue_begin_of_strands_with_0_capping" << endl;

	info_file
		<< number_of_strands_where_capping_is_checked << "\t"
		<< number_of_strands_with_2_cappings << "\t"
		<< number_of_strands_with_1_capping << "\t"
		<< number_of_strands_where_capping_is_checked - number_of_strands_with_2_cappings - number_of_strands_with_1_capping << "\t";

	for ( Size kk = 1;  kk <= residue_begin_of_strands_with_1_capping.size(); kk++ ) {
		if ( kk == residue_begin_of_strands_with_1_capping.size() ) {
			info_file << residue_begin_of_strands_with_1_capping[kk];
		} else {
			info_file << residue_begin_of_strands_with_1_capping[kk] << ",";
		}
	}

	info_file << "\t";

	for ( Size kk = 1;  kk <= residue_begin_of_strands_with_0_capping.size(); kk++ ) {
		if ( kk == residue_begin_of_strands_with_0_capping.size() ) {
			info_file << residue_begin_of_strands_with_0_capping[kk];
		} else {
			info_file << residue_begin_of_strands_with_0_capping[kk] << ",";
		}
	}

	info_file << endl;
	info_file.close();

	return 0;
} //write_beta_sheet_capping_info_to_a_file


//write_chain_B_resNum_to_a_file
core::Size
write_chain_B_resNum_to_a_file(
	string tag,
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string report_file_name = pdb_file_name + "_chain_B_resNum.txt";
	ofstream report_file;
	report_file.open(report_file_name.c_str());
	utility::vector1<SandwichFragment> chain_B_resNum =
		get_chain_B_resNum(
		struct_id,
		db_session,
		sw_can_by_sh_id // was 'vec_sw_can_by_sh_id[ii]'
	);
	for ( Size i=1; i<=chain_B_resNum.size(); ++i ) {
		report_file << chain_B_resNum[i].get_resNum() << endl;
	}
	report_file.close();
	return 0;
}//write_chain_B_resNum_to_a_file

//write_heading_direction_of_all_AA_in_a_strand_to_a_file
core::Size
write_heading_direction_of_all_AA_in_a_strand_to_a_file(
	string tag,
	StructureID struct_id, // needed argument
	sessionOP db_session,
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string heading_file_name = pdb_file_name + "_heading_direction_of_all_AA_in_a_strand.txt";
	ofstream heading_file;

	heading_file.open(heading_file_name.c_str());
	heading_file << "residue_begin\tresidue_end\theading_directions" << endl;

	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
		string heading_directions =
			report_heading_directions_of_all_AA_in_a_strand (
			struct_id,
			db_session,
			pose,
			bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),
			bs_of_sw_can_by_sh[ii].get_sheet_id(),
			bs_of_sw_can_by_sh[ii].get_start(),
			bs_of_sw_can_by_sh[ii].get_end());
		heading_file << bs_of_sw_can_by_sh[ii].get_start() << "\t" << bs_of_sw_can_by_sh[ii].get_end() << "\t" << heading_directions << endl;
	}

	heading_file.close();

	return 0;
} //write_heading_direction_of_all_AA_in_a_strand_to_a_file






// write_number_of_electrostatic_interactions_of_residues_to_files
//
// detail: according to
// write_electrostatic_interactions_of_surface_residues_in_a_strand_,
// write_electrostatic_interactions_of_all_residues_in_a_strand_,
// write_electrostatic_interactions_of_all_residues_,
// count number of electrostatic_interactions among selected residues,

// for example, count number of electrostatic_interactions among surface_residues_in_strands if
// write_electrostatic_interactions_of_surface_residues_in_a_strand_ is selected (2014/09/30)

core::Size
write_number_of_electrostatic_interactions_of_residues_to_files(
	string tag,
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	string dssp_code,
	string heading_direction,
	Real distance_cutoff_for_electrostatic_interactions_,
	Real CB_b_factor_cutoff_for_electrostatic_interactions_,
	Size min_primary_seq_distance_diff_for_electrostatic_interactions_)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	string ElectroStatics_txt_file_name;
	string SaltBridges_pymol_ready_file_name;
	if ( dssp_code == "all_dssp" ) {
		ElectroStatics_txt_file_name = pdb_file_name + "_electrostatic_interactions_among_all_residues.txt";
		SaltBridges_pymol_ready_file_name = pdb_file_name + "_salt_bridges_among_all_residues.pymolready";
	} else { // dssp_code = "E"
		if ( heading_direction == "surface" ) {
			ElectroStatics_txt_file_name = pdb_file_name + "_electrostatic_interactions_among_surface_residues_in_a_strand.txt";
			SaltBridges_pymol_ready_file_name = pdb_file_name + "_salt_bridges_among_surface_residues_in_a_strand.pymolready";
		} else { // (heading_direction == "all_direction")
			ElectroStatics_txt_file_name = pdb_file_name + "_electrostatic_interactions_among_all_residues_in_a_strand.txt";
			SaltBridges_pymol_ready_file_name = pdb_file_name + "_salt_bridges_among_all_residues_in_a_strand.pymolready";
		}
	}

	ofstream ElectroStatics_txt_file;
	ofstream SaltBridges_pymol_ready_file;

	ElectroStatics_txt_file.open(ElectroStatics_txt_file_name.c_str());
	SaltBridges_pymol_ready_file.open(SaltBridges_pymol_ready_file_name.c_str());

	string number_of_attractions_title;
	string head_attrac = "attrac_by_centroid_w_";
	number_of_attractions_title.append(head_attrac);

	string number_of_repulsions_title;
	string head_repul = "repul_by_centroid_w_";
	number_of_repulsions_title.append(head_repul);

	std::ostringstream cutoff;
	cutoff << distance_cutoff_for_electrostatic_interactions_;
	std::string cutoff_str = cutoff.str();
	number_of_attractions_title.append(cutoff_str);
	number_of_repulsions_title.append(cutoff_str);

	string tail = "_A";
	number_of_attractions_title.append(tail);
	number_of_repulsions_title.append(tail);

	ElectroStatics_txt_file << "resNum\ttype\t" << number_of_attractions_title << "\t" << number_of_repulsions_title << "\tattrac_minus_repul\tsalt_bridge\tCC_bridge\tNO_bridge\tsum_of_salt_CC_NO_bridges\tlonger_range_attrac_ion_pair" << endl;

	SaltBridges_pymol_ready_file << "select salt_bridges, resi ";
	bool first_addition_for_SaltBridges_pymol_ready_file = true;

	utility::vector1<Size> vector_of_unique_distinct_sw_ids = get_distinct_sw_id_from_sandwich_table (struct_id, db_session);

	std::vector<int> vec_number_of_attractions_by_centroid;
	std::vector<int> vec_number_of_repulsions_by_centroid;
	std::vector<int> vec_net_attrac_by_centroid;
	std::vector<int> vec_number_of_salt_bridges;
	std::vector<int> vec_number_of_CC_bridges;
	std::vector<int> vec_number_of_NO_bridges;
	std::vector<int> vec_sum_of_salt_CC_NO_bridges;
	std::vector<int> vec_number_of_longer_range_attrac_ion_pair;

	for ( Size sw_ii=1; sw_ii<=vector_of_unique_distinct_sw_ids.size(); ++sw_ii ) { // per each beta-sandwich
		utility::vector1<int> vector_of_residue_num_of_rkde =
			retrieve_residue_num_of_rkde(
			struct_id,
			db_session,
			vector_of_unique_distinct_sw_ids[sw_ii], //sw_can_by_sh_id
			dssp_code,
			heading_direction
		);
		for ( Size residue_i=1; residue_i<=vector_of_residue_num_of_rkde.size(); residue_i++ ) {
			int residue_num = vector_of_residue_num_of_rkde[residue_i];

			// <begin> check whether "current" residue has low atom position uncertainty
			pose::PDBInfoCOP info = pose.pdb_info();
			Real B_factor_of_CB = info->temperature( residue_num, 5 ); // '5' atom will be 'H' for Gly, otherwise typycal CB
			if ( B_factor_of_CB > CB_b_factor_cutoff_for_electrostatic_interactions_ ) {
				TR << "residue_num:\t" << residue_num << " is dropped by too high CB_b_factor" << endl;
				// this "current" residue has too high atom position uncertainty,
				// so let's not use this residue when counting number of electrostatic interactions
				continue;
			}
			// <end> check whether "current" residue has low atom position uncertainty

			ElectroStatics_txt_file << residue_num << "\t" << pose.residue_type(residue_num).name3();

			numeric::xyzVector< core::Real > xyz_of_centroid_of_RKDE;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_R;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_R;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_of_K;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_DE;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_DE;

			if ( pose.residue_type(residue_num).name3() == "ARG" ) {
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x() =
					(pose.residue(residue_num).atom(" NE ").xyz().x()
					+ pose.residue(residue_num).atom(" CZ ").xyz().x()
					+ pose.residue(residue_num).atom(" NH1").xyz().x()
					+ pose.residue(residue_num).atom(" NH2").xyz().x())/4;

				xyz_of_centroid_of_RKDE.y() =
					(pose.residue(residue_num).atom(" NE ").xyz().y()
					+ pose.residue(residue_num).atom(" CZ ").xyz().y()
					+ pose.residue(residue_num).atom(" NH1").xyz().y()
					+ pose.residue(residue_num).atom(" NH2").xyz().y())/4;

				xyz_of_centroid_of_RKDE.z() =
					(pose.residue(residue_num).atom(" NE ").xyz().z()
					+ pose.residue(residue_num).atom(" CZ ").xyz().z()
					+ pose.residue(residue_num).atom(" NH1").xyz().z()
					+ pose.residue(residue_num).atom(" NH2").xyz().z())/4;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_R = pose.residue(residue_num).atom(" NH1").xyz();
				xyz_of_terminal_atom_2_of_R = pose.residue(residue_num).atom(" NH2").xyz();
			} else if ( pose.residue_type(residue_num).name3() == "LYS" ) {
				xyz_of_centroid_of_RKDE = pose.residue(residue_num).atom(" NZ ").xyz();
				xyz_of_terminal_atom_of_K = pose.residue(residue_num).atom(" NZ ").xyz();
			} else if ( pose.residue_type(residue_num).name3() == "ASP" ) {
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x() =
					(pose.residue(residue_num).atom(" CG ").xyz().x()
					+ pose.residue(residue_num).atom(" OD1").xyz().x()
					+ pose.residue(residue_num).atom(" OD2").xyz().x())/3;

				xyz_of_centroid_of_RKDE.y() =
					(pose.residue(residue_num).atom(" CG ").xyz().y()
					+ pose.residue(residue_num).atom(" OD1").xyz().y()
					+ pose.residue(residue_num).atom(" OD2").xyz().y())/3;

				xyz_of_centroid_of_RKDE.z() =
					(pose.residue(residue_num).atom(" CG ").xyz().z()
					+ pose.residue(residue_num).atom(" OD1").xyz().z()
					+ pose.residue(residue_num).atom(" OD2").xyz().z())/3;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_DE = pose.residue(residue_num).atom(" OD1").xyz();
				xyz_of_terminal_atom_2_of_DE = pose.residue(residue_num).atom(" OD2").xyz();
			} else { //if (pose.residue_type(residue_num).name3() == "GLU")
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x() =
					(pose.residue(residue_num).atom(" CD ").xyz().x()
					+ pose.residue(residue_num).atom(" OE1").xyz().x()
					+ pose.residue(residue_num).atom(" OE2").xyz().x())/3;

				xyz_of_centroid_of_RKDE.y() =
					(pose.residue(residue_num).atom(" CD ").xyz().y()
					+ pose.residue(residue_num).atom(" OE1").xyz().y()
					+ pose.residue(residue_num).atom(" OE2").xyz().y())/3;

				xyz_of_centroid_of_RKDE.z() =
					(pose.residue(residue_num).atom(" CD ").xyz().z()
					+ pose.residue(residue_num).atom(" OE1").xyz().z()
					+ pose.residue(residue_num).atom(" OE2").xyz().z())/3;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_DE = pose.residue(residue_num).atom(" OE1").xyz();
				xyz_of_terminal_atom_2_of_DE = pose.residue(residue_num).atom(" OE2").xyz();
			}

			int number_of_attractions_by_centroid = 0;
			int number_of_repulsions_by_centroid = 0;

			Size number_of_salt_bridge = 0;
			Size number_of_C_C_bridge = 0;
			Size number_of_N_O_bridge = 0;
			Size number_of_longer_range_attrac_ion_pair = 0;

			for ( Size other_residue_i=1; other_residue_i<=vector_of_residue_num_of_rkde.size(); ++other_residue_i ) {
				if ( other_residue_i == residue_i ) {
					continue;
				}

				int other_residue_num = vector_of_residue_num_of_rkde[other_residue_i];

				//     TR << "residue_num: " << residue_num << endl;
				//     TR << "other_residue_num: " << other_residue_num << endl;

				if ( core::Size(std::abs(other_residue_num - residue_num)) <
						min_primary_seq_distance_diff_for_electrostatic_interactions_ ) {
					continue;
				}

				// <begin> check whether "other_residue" has low atom position uncertainty
				B_factor_of_CB = info->temperature( other_residue_num, 5 ); // '5' atom will be 'H' for Gly
				if ( B_factor_of_CB > CB_b_factor_cutoff_for_electrostatic_interactions_ ) {
					TR << "other_residue_num:\t" << other_residue_num << " is dropped by too high CB_b_factor" << endl;
					continue; // this "other_residue" has too high atom position uncertainty, so let's not use this residue when counting number of electrostatic interactions
				}
				// <end> check whether "other_residue" has low atom position uncertainty

				numeric::xyzVector< core::Real > xyz_of_centroid_of_other_RKDE;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_other_R;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_other_R;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_of_other_K;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_other_DE;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_other_DE;

				if ( pose.residue_type(other_residue_num).name3() == "ARG" ) {
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x() =
						(pose.residue(other_residue_num).atom(" NE ").xyz().x()
						+ pose.residue(other_residue_num).atom(" CZ ").xyz().x()
						+ pose.residue(other_residue_num).atom(" NH1").xyz().x()
						+ pose.residue(other_residue_num).atom(" NH2").xyz().x())/4;

					xyz_of_centroid_of_other_RKDE.y() =
						(pose.residue(other_residue_num).atom(" NE ").xyz().y()
						+ pose.residue(other_residue_num).atom(" CZ ").xyz().y()
						+ pose.residue(other_residue_num).atom(" NH1").xyz().y()
						+ pose.residue(other_residue_num).atom(" NH2").xyz().y())/4;

					xyz_of_centroid_of_other_RKDE.z() =
						(pose.residue(other_residue_num).atom(" NE ").xyz().z()
						+ pose.residue(other_residue_num).atom(" CZ ").xyz().z()
						+ pose.residue(other_residue_num).atom(" NH1").xyz().z()
						+ pose.residue(other_residue_num).atom(" NH2").xyz().z())/4;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_R = pose.residue(other_residue_num).atom(" NH1").xyz();
					xyz_of_terminal_atom_2_of_other_R = pose.residue(other_residue_num).atom(" NH2").xyz();
				} else if ( pose.residue_type(other_residue_num).name3() == "LYS" ) {
					xyz_of_centroid_of_other_RKDE = pose.residue(other_residue_num).atom(" NZ ").xyz();
					xyz_of_terminal_atom_of_other_K = pose.residue(other_residue_num).atom(" NZ ").xyz();
				} else if ( pose.residue_type(other_residue_num).name3() == "ASP" ) {
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x() =
						(pose.residue(other_residue_num).atom(" CG ").xyz().x()
						+ pose.residue(other_residue_num).atom(" OD1").xyz().x()
						+ pose.residue(other_residue_num).atom(" OD2").xyz().x())/3;

					xyz_of_centroid_of_other_RKDE.y() =
						(pose.residue(other_residue_num).atom(" CG ").xyz().y()
						+ pose.residue(other_residue_num).atom(" OD1").xyz().y()
						+ pose.residue(other_residue_num).atom(" OD2").xyz().y())/3;

					xyz_of_centroid_of_other_RKDE.z() =
						(pose.residue(other_residue_num).atom(" CG ").xyz().z()
						+ pose.residue(other_residue_num).atom(" OD1").xyz().z()
						+ pose.residue(other_residue_num).atom(" OD2").xyz().z())/3;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_DE = pose.residue(other_residue_num).atom(" OD1").xyz();
					xyz_of_terminal_atom_2_of_other_DE = pose.residue(other_residue_num).atom(" OD2").xyz();
				} else { //if (pose.residue_type(other_residue_num).name3() == "GLU")
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x() =
						(pose.residue(other_residue_num).atom(" CD ").xyz().x()
						+ pose.residue(other_residue_num).atom(" OE1").xyz().x()
						+ pose.residue(other_residue_num).atom(" OE2").xyz().x())/3;

					xyz_of_centroid_of_other_RKDE.y() =
						(pose.residue(other_residue_num).atom(" CD ").xyz().y()
						+ pose.residue(other_residue_num).atom(" OE1").xyz().y()
						+ pose.residue(other_residue_num).atom(" OE2").xyz().y())/3;

					xyz_of_centroid_of_other_RKDE.z() =
						(pose.residue(other_residue_num).atom(" CD ").xyz().z()
						+ pose.residue(other_residue_num).atom(" OE1").xyz().z()
						+ pose.residue(other_residue_num).atom(" OE2").xyz().z())/3;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_DE = pose.residue(other_residue_num).atom(" OE1").xyz();
					xyz_of_terminal_atom_2_of_other_DE = pose.residue(other_residue_num).atom(" OE2").xyz();
				}

				Real distance_between_centroid = xyz_of_centroid_of_RKDE.distance(xyz_of_centroid_of_other_RKDE);
				// TR << "distance_between_centroid: " << distance_between_centroid << endl;

				if ( pose.residue_type(residue_num).name3() == "ARG" ) {
					if ( pose.residue_type(other_residue_num).name3() == "ASP" || pose.residue_type(other_residue_num).name3() == "GLU" ) {
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_attractions_by_centroid++;
						}
						Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_1_of_R.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_1_of_R.distance(xyz_of_terminal_atom_2_of_other_DE);
						Real distance_3_between_terminal_atoms = xyz_of_terminal_atom_2_of_R.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_4_between_terminal_atoms = xyz_of_terminal_atom_2_of_R.distance(xyz_of_terminal_atom_2_of_other_DE);
						if ( distance_between_centroid < 7.0 ) {
							if ( distance_between_centroid < 4.0 ) {
								if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0
										|| distance_3_between_terminal_atoms < 4.0 || distance_4_between_terminal_atoms < 4.0 ) {
									number_of_salt_bridge++;
								} else {
									number_of_C_C_bridge++;
								}
							} else {
								if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0
										|| distance_3_between_terminal_atoms < 4.0 || distance_4_between_terminal_atoms < 4.0 ) {
									number_of_N_O_bridge++;
								} else {
									number_of_longer_range_attrac_ion_pair++;
								}
							}
						}
					} else {
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_repulsions_by_centroid++;
						}
					}
				} else if ( pose.residue_type(residue_num).name3() == "LYS" ) {
					if ( pose.residue_type(other_residue_num).name3() == "ASP" || pose.residue_type(other_residue_num).name3() == "GLU" ) {
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_attractions_by_centroid++;
						}
						Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_of_K.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_of_K.distance(xyz_of_terminal_atom_2_of_other_DE);
						if ( distance_between_centroid < 7.0 ) {
							if ( distance_between_centroid < 4.0 ) {
								if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 ) {
									number_of_salt_bridge++;
								} else {
									number_of_C_C_bridge++;
								}
							} else {
								if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 ) {
									number_of_N_O_bridge++;
								} else {
									number_of_longer_range_attrac_ion_pair++;
								}
							}
						}
					} else {
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_repulsions_by_centroid++;
						}
					}
				} else { // (pose.residue_type(residue_num).name3() == "ASP") || (pose.residue_type(residue_num).name3() == "GLU")
					if ( (pose.residue_type(other_residue_num).name3() == "ASP") || (pose.residue_type(other_residue_num).name3() == "GLU") ) {
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_repulsions_by_centroid++;
						}
					} else { // (pose.residue_type(other_residue_num).name3() == "ARG") || (pose.residue_type(other_residue_num).name3() == "LYS")
						if ( distance_between_centroid < distance_cutoff_for_electrostatic_interactions_ ) {
							number_of_attractions_by_centroid++;
						}
						if ( pose.residue_type(other_residue_num).name3() == "ARG" ) {
							Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_1_of_other_R.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_1_of_other_R.distance(xyz_of_terminal_atom_2_of_DE);
							Real distance_3_between_terminal_atoms = xyz_of_terminal_atom_2_of_other_R.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_4_between_terminal_atoms = xyz_of_terminal_atom_2_of_other_R.distance(xyz_of_terminal_atom_2_of_DE);

							if ( distance_between_centroid < 7.0 ) {
								if ( distance_between_centroid < 4.0 ) {
									if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 || distance_3_between_terminal_atoms < 4.0  || distance_4_between_terminal_atoms < 4.0 ) {
										number_of_salt_bridge++;
									} else {
										number_of_C_C_bridge++;
									}
								} else {
									if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 || distance_3_between_terminal_atoms < 4.0  || distance_4_between_terminal_atoms < 4.0 ) {
										number_of_N_O_bridge++;
									} else {
										number_of_longer_range_attrac_ion_pair++;
									}
								}
							}
						} else { // (pose.residue_type(other_residue_num).name3() == "LYS")
							Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_of_other_K.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_of_other_K.distance(xyz_of_terminal_atom_2_of_DE);
							if ( distance_between_centroid < 7.0 ) {
								if ( distance_between_centroid < 4.0 ) {
									if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 ) {
										number_of_salt_bridge++;
									} else {
										number_of_C_C_bridge++;
									}
								} else {
									if ( distance_1_between_terminal_atoms < 4.0 || distance_2_between_terminal_atoms < 4.0 ) {
										number_of_N_O_bridge++;
									} else {
										number_of_longer_range_attrac_ion_pair++;
									}
								}
							}
						}
					}
				}
			}

			vec_number_of_attractions_by_centroid.push_back(number_of_attractions_by_centroid);
			vec_number_of_repulsions_by_centroid.push_back(number_of_repulsions_by_centroid);
			vec_net_attrac_by_centroid.push_back(number_of_attractions_by_centroid-number_of_repulsions_by_centroid);
			vec_number_of_salt_bridges.push_back(number_of_salt_bridge);
			vec_number_of_CC_bridges.push_back(number_of_C_C_bridge);
			vec_number_of_NO_bridges.push_back(number_of_N_O_bridge);
			vec_sum_of_salt_CC_NO_bridges.push_back(number_of_salt_bridge+number_of_C_C_bridge+number_of_N_O_bridge);
			vec_number_of_longer_range_attrac_ion_pair.push_back(number_of_longer_range_attrac_ion_pair);

			ElectroStatics_txt_file << "\t" << number_of_attractions_by_centroid << "\t" << number_of_repulsions_by_centroid << "\t" << number_of_attractions_by_centroid-number_of_repulsions_by_centroid << "\t" << number_of_salt_bridge << "\t" << number_of_C_C_bridge << "\t" << number_of_N_O_bridge << "\t" << number_of_salt_bridge+number_of_C_C_bridge+number_of_N_O_bridge << "\t" <<  number_of_longer_range_attrac_ion_pair << endl;

			if ( number_of_salt_bridge + number_of_C_C_bridge + number_of_N_O_bridge > 0 ) {

				if ( first_addition_for_SaltBridges_pymol_ready_file ) {
					SaltBridges_pymol_ready_file << residue_num;
				} else {
					SaltBridges_pymol_ready_file << "+" << residue_num;
				}
				first_addition_for_SaltBridges_pymol_ready_file = false;
			}

		} // per each residue
	}

	SaltBridges_pymol_ready_file << endl;

	// report avg (min~max)
	// warning: this avg (min~max) assumes that I deal with 1 sandwich per 1 pdb file
	float avg_attrac = std::accumulate(vec_number_of_attractions_by_centroid.begin(),  vec_number_of_attractions_by_centroid.end(), 0) / static_cast<float>(vec_number_of_attractions_by_centroid.size());

	float avg_repul = std::accumulate(vec_number_of_repulsions_by_centroid.begin(),  vec_number_of_repulsions_by_centroid.end(), 0) / static_cast<float>(vec_number_of_repulsions_by_centroid.size());

	float avg_net_attrac = std::accumulate(vec_net_attrac_by_centroid.begin(),  vec_net_attrac_by_centroid.end(), 0) / static_cast<float>(vec_net_attrac_by_centroid.size());

	float avg_salt_bridges = std::accumulate(vec_number_of_salt_bridges.begin(),  vec_number_of_salt_bridges.end(), 0) / static_cast<float>(vec_number_of_salt_bridges.size());

	float avg_CC_bridges = std::accumulate(vec_number_of_CC_bridges.begin(),  vec_number_of_CC_bridges.end(), 0) / static_cast<float>(vec_number_of_CC_bridges.size());

	float avg_NO_bridges = std::accumulate(vec_number_of_NO_bridges.begin(),  vec_number_of_NO_bridges.end(), 0) / static_cast<float>(vec_number_of_NO_bridges.size());

	float avg_sum_of_salt_CC_NO_bridges = std::accumulate(vec_sum_of_salt_CC_NO_bridges.begin(),  vec_sum_of_salt_CC_NO_bridges.end(), 0) / static_cast<float>(vec_sum_of_salt_CC_NO_bridges.size());

	float avg_number_of_longer_range_attrac_ion_pair = std::accumulate(vec_number_of_longer_range_attrac_ion_pair.begin(),  vec_number_of_longer_range_attrac_ion_pair.end(), 0) / static_cast<float>(vec_number_of_longer_range_attrac_ion_pair.size());


	int min_attrac;
	int max_attrac;
	if ( vec_number_of_attractions_by_centroid.size() == 0 ) {
		min_attrac = 0;
		max_attrac = 0;
	} else {
		min_attrac = *std::min_element(vec_number_of_attractions_by_centroid.begin(), vec_number_of_attractions_by_centroid.end());
		max_attrac = *std::max_element(vec_number_of_attractions_by_centroid.begin(), vec_number_of_attractions_by_centroid.end());
	}


	int min_repul;
	int max_repul;
	if ( vec_number_of_repulsions_by_centroid.size() == 0 ) {
		min_repul = 0;
		max_repul = 0;
	} else {
		min_repul = *std::min_element(vec_number_of_repulsions_by_centroid.begin(), vec_number_of_repulsions_by_centroid.end());
		max_repul = *std::max_element(vec_number_of_repulsions_by_centroid.begin(), vec_number_of_repulsions_by_centroid.end());
	}

	int min_net_attrac;
	int max_net_attrac;
	if ( vec_net_attrac_by_centroid.size() == 0 ) {
		min_net_attrac = 0;
		max_net_attrac = 0;
	} else {
		min_net_attrac = *std::min_element(vec_net_attrac_by_centroid.begin(), vec_net_attrac_by_centroid.end());
		max_net_attrac = *std::max_element(vec_net_attrac_by_centroid.begin(), vec_net_attrac_by_centroid.end());
	}

	int min_salt_bridge;
	int max_salt_bridge;
	if ( vec_number_of_salt_bridges.size() == 0 ) {
		min_salt_bridge = 0;
		max_salt_bridge = 0;
	} else {
		min_salt_bridge = *std::min_element(vec_number_of_salt_bridges.begin(), vec_number_of_salt_bridges.end());
		max_salt_bridge = *std::max_element(vec_number_of_salt_bridges.begin(), vec_number_of_salt_bridges.end());
	}

	int min_CC_bridge;
	int max_CC_bridge;
	if ( vec_number_of_CC_bridges.size() == 0 ) {
		min_CC_bridge = 0;
		max_CC_bridge = 0;
	} else {
		min_CC_bridge = *std::min_element(vec_number_of_CC_bridges.begin(), vec_number_of_CC_bridges.end());
		max_CC_bridge = *std::max_element(vec_number_of_CC_bridges.begin(), vec_number_of_CC_bridges.end());
	}

	int min_NO_bridge;
	int max_NO_bridge;
	if ( vec_number_of_NO_bridges.size() == 0 ) {
		min_NO_bridge = 0;
		max_NO_bridge = 0;
	} else {
		min_NO_bridge = *std::min_element(vec_number_of_NO_bridges.begin(), vec_number_of_NO_bridges.end());
		max_NO_bridge = *std::max_element(vec_number_of_NO_bridges.begin(), vec_number_of_NO_bridges.end());
	}

	int min_sum_of_salt_CC_NO_bridges;
	int max_sum_of_salt_CC_NO_bridges;
	if ( vec_sum_of_salt_CC_NO_bridges.size() == 0 ) {
		min_sum_of_salt_CC_NO_bridges = 0;
		max_sum_of_salt_CC_NO_bridges = 0;
	} else {
		min_sum_of_salt_CC_NO_bridges = *std::min_element(vec_sum_of_salt_CC_NO_bridges.begin(), vec_sum_of_salt_CC_NO_bridges.end());
		max_sum_of_salt_CC_NO_bridges = *std::max_element(vec_sum_of_salt_CC_NO_bridges.begin(), vec_sum_of_salt_CC_NO_bridges.end());
	}

	int min_number_of_longer_range_attrac_ion_pair;
	int max_number_of_longer_range_attrac_ion_pair;
	if ( vec_number_of_longer_range_attrac_ion_pair.size() == 0 ) {
		min_number_of_longer_range_attrac_ion_pair = 0;
		max_number_of_longer_range_attrac_ion_pair = 0;
	} else {
		min_number_of_longer_range_attrac_ion_pair = *std::min_element(vec_number_of_longer_range_attrac_ion_pair.begin(), vec_number_of_longer_range_attrac_ion_pair.end());
		max_number_of_longer_range_attrac_ion_pair = *std::max_element(vec_number_of_longer_range_attrac_ion_pair.begin(), vec_number_of_longer_range_attrac_ion_pair.end());
	}

	ElectroStatics_txt_file << "avg\t(min~max)\t";
	ElectroStatics_txt_file << round_to_float(avg_attrac) << " (" << min_attrac << "~" << max_attrac << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_repul) << " (" << min_repul << "~" << max_repul << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_net_attrac) << " (" << min_net_attrac << "~" << max_net_attrac << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_salt_bridges) << " (" << min_salt_bridge << "~" << max_salt_bridge << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_CC_bridges) << " (" << min_CC_bridge << "~" << max_CC_bridge << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_NO_bridges) << " (" << min_NO_bridge << "~" << max_NO_bridge << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_sum_of_salt_CC_NO_bridges) << " (" << min_sum_of_salt_CC_NO_bridges << "~" << max_sum_of_salt_CC_NO_bridges << ")\t";
	ElectroStatics_txt_file << round_to_float(avg_number_of_longer_range_attrac_ion_pair) << " (" << min_number_of_longer_range_attrac_ion_pair << "~" << max_number_of_longer_range_attrac_ion_pair << ")\t" << endl;

	ElectroStatics_txt_file.close();

	return 0;
} // write_number_of_electrostatic_interactions_of_residues_to_files


// write_p_aa_pp_of_AAs_to_a_file
// (Probability of amino acid at phipsi)
// (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )
core::Size
write_p_aa_pp_of_AAs_to_a_file(
	string tag,
	Pose & dssp_pose)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string p_aa_pp_file_name = pdb_file_name + "_p_aa_pp_at_each_AA.txt";
	ofstream p_aa_pp_file;

	p_aa_pp_file.open(p_aa_pp_file_name.c_str());
	p_aa_pp_file << "residue_number\tres_type\tp_aa_pp" << endl;

	for ( Size ii=1; ii<=dssp_pose.size(); ii++ ) {
		core::scoring::EnergyMap em1 = dssp_pose.energies().residue_total_energies(ii);
		Real resi_p_aa_pp = em1[core::scoring::p_aa_pp];

		p_aa_pp_file << ii << "\t" << dssp_pose.residue_type(ii).name3() << "\t" << round_to_Real(resi_p_aa_pp) << endl;
	}
	p_aa_pp_file.close();
	return 0;
} //write_p_aa_pp_of_AAs_to_a_file


// (begin) write_phi_psi_of_each_residue
core::Size
write_phi_psi_of_each_residue_to_a_file(
	string tag,
	Pose & dssp_pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	bool write_phi_psi_of_E_,
	bool write_phi_psi_of_all_,
	Size max_num_sw_per_pdb_,
	StructureID struct_id, // needed argument
	sessionOP db_session,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_
)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	if ( write_phi_psi_of_E_ ) {
		string phi_psi_file_name = pdb_file_name + "_phi_psi_of_strand_res.txt";
		ofstream phi_psi_file;
		phi_psi_file.open(phi_psi_file_name.c_str());
		phi_psi_file << "tag\tres_num\tres_type\tres_at_terminal\tsheet_is_antiparallel\tstrand_is_at_edge\tphi\tpsi" << endl;

		for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii ) {
			if ( bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_ ) {
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());

			string strand_is_at_edge = is_this_strand_at_edge (
				dssp_pose,
				struct_id,
				db_session,
				bs_of_sw_can_by_sh[ii].get_sheet_id(),
				bs_of_sw_can_by_sh[ii].get_start(),
				bs_of_sw_can_by_sh[ii].get_end(),
				min_CA_CA_dis_,
				max_CA_CA_dis_);

			Size res_at_terminal;
			for ( Size res_num = bs_of_sw_can_by_sh[ii].get_start(); res_num <= bs_of_sw_can_by_sh[ii].get_end(); res_num++ ) {
				if ( res_num == bs_of_sw_can_by_sh[ii].get_start() || res_num == bs_of_sw_can_by_sh[ii].get_end() ) {
					res_at_terminal = 1;
				} else {
					res_at_terminal = 0;
				}
				Real phi = dssp_pose.phi(res_num);
				Real psi = dssp_pose.psi(res_num);
				phi_psi_file << tag << "\t" << res_num << "\t" << dssp_pose.residue_type(res_num).name3() << "\t" << res_at_terminal << "\t" << sheet_antiparallel << "\t" << strand_is_at_edge << "\t" << round_to_Real(phi) << "\t" << round_to_Real(psi) << endl;
			}
		}
		phi_psi_file.close();
	} //write_phi_psi_of_E_


	if ( !write_phi_psi_of_E_ && write_phi_psi_of_all_ ) { // 2014/10/11 it should be like this to avoid "ERROR: seqpos <= size() ERROR:: Exit from: src/core/conformation/Conformation.hh line: 361", when I have time, it should be debugged
		//if (write_phi_psi_of_all_)
		string phi_psi_file_name = pdb_file_name + "_phi_psi_of_all_res.txt";
		ofstream phi_psi_file;
		phi_psi_file.open(phi_psi_file_name.c_str());
		phi_psi_file << "tag\tres_num\tres_type\tdssp\tphi\tpsi" << endl;

		for ( core::Size ii=1; ii<=dssp_pose.size(); ii++ ) {
			char res_ss( dssp_pose.secstruct( ii ) ) ;
			Real phi = dssp_pose.phi(ii);
			Real psi = dssp_pose.psi(ii);
			phi_psi_file << tag << "\t" << ii << "\t" << dssp_pose.residue_type(ii).name3() << "\t" << res_ss << "\t" << round_to_Real(phi) << "\t" << round_to_Real(psi) << endl;
		}
	}//write_phi_psi_of_all_

	return 0;
}
// (end) write_phi_psi_of_each_residue



//write_rama_of_AAs_to_a_file
// (ramachandran preferences)
// (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )
core::Size
write_rama_of_AAs_to_a_file(
	string tag,
	Pose & dssp_pose)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string rama_file_name = pdb_file_name + "_rama_at_each_AA.txt";
	ofstream rama_file;

	rama_file.open(rama_file_name.c_str());
	rama_file << "residue_number\t\tres_type\trama" << endl;

	for ( Size ii=1; ii<=dssp_pose.size(); ii++ ) {
		core::scoring::EnergyMap em1 = dssp_pose.energies().residue_total_energies(ii);
		Real rama_at_this_AA = em1[core::scoring::rama];

		rama_file << ii << "\t" << dssp_pose.residue_type(ii).name3() << "\t" << round_to_Real(rama_at_this_AA) << endl;
	}
	rama_file.close();
	return 0;
} //write_rama_of_AAs_to_a_file


// write_resfile_to_a_file_when_seq_rec_is_bad
core::Size
write_resfile_to_a_file_when_seq_rec_is_bad(
	string tag,
	StructureID struct_id, // needed argument
	utility::sql_database::sessionOP db_session, // needed argument
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	bool write_resfile_to_minimize_too_much_hydrophobic_surface_,
	bool write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_,
	bool write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_
)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	// <begin> make a resfile to design all residues
	string resfile_name = pdb_file_name + "_resfile_to_design_all_residues_when_seq_rec_is_bad.txt";
	ofstream resfile_stream;

	resfile_stream.open(resfile_name.c_str());

	resfile_stream << "EX 1 NOTAA C" << endl;
	resfile_stream << "USE_INPUT_SC" << endl;
	resfile_stream << "start" << endl;

	resfile_stream << "# final resfile rule update: 05/08/2014" << endl;
	resfile_stream << "#\tbased on 203 native sandwiches and designed His at core-heading core-strands" << endl;
	resfile_stream << "# NOTAA\tCFPWY for surface_heading residues at core strands" << endl;
	resfile_stream << "# NOTAA\tCFWY for surface_heading residues at edge strands" << endl;

	if ( write_resfile_to_minimize_too_much_hydrophobic_surface_ ) {
		resfile_stream << "# PIKAA\tDENQST\tfor every 2nd surface-heading\tresidues " << endl;
	}

	resfile_stream << "# NOTAA\tCDEHKNPQR for core_heading residues at core strands" << endl;
	if ( write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_ ) {
		resfile_stream << "# NOTAA\tCDEFHKNPQRWY for every 3rd core-heading\tresidues\ton\tcore\tstrands" << endl;
	}

	resfile_stream << "# NOTAA\tCDNP for core_heading residues at edge strands" << endl;
	if ( write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_ ) {
		resfile_stream << "# NOTAA\tCDFNPWY\tfor every 2nd core-heading\tresidues\ton\tedge\tstrands" << endl;
	}

	int count_to_minimize_too_much_hydrophobic_surface = 0;
	int count_to_minimize_too_many_core_heading_FWY_on_core_strands = 0;
	int count_to_minimize_too_many_core_heading_FWY_on_edge_strands = 0;

	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
		Size residue_begin = bs_of_sw_can_by_sh[ii].get_start();
		Size residue_end = bs_of_sw_can_by_sh[ii].get_end();
		for ( Size residue_num = residue_begin; residue_num <= residue_end; residue_num++ ) {
			{
				string heading = determine_heading_direction_by_vector (struct_id, db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), residue_begin, residue_end, residue_num);
				if ( heading == "surface" ) {
					if ( write_resfile_to_minimize_too_much_hydrophobic_surface_ ) {
						// see_edge_or_core_or_loop_or_short_edge is refactored
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "loop" || edge_or_core == "short_edge" ) {
							continue; // I will write resfile for this residue later
						}

						count_to_minimize_too_much_hydrophobic_surface++;
						if ( count_to_minimize_too_much_hydrophobic_surface == 2 ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tPIKAA\tDENQST" << endl;
							count_to_minimize_too_much_hydrophobic_surface = 0;
						} else {
							if ( edge_or_core == "core" ) {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFPWY" << endl;
							} else if ( edge_or_core == "edge" ) {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFWY" << endl;
							}
						}
					} else { //(!write_resfile_to_minimize_too_much_hydrophobic_surface_)
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "core" ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFPWY" << endl;
						} else if ( edge_or_core == "edge" ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFWY" << endl;
						}
					}
				} else if ( heading == "core" ) {
					if ( (write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_) && (write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_) ) {
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "core" ) {
							count_to_minimize_too_many_core_heading_FWY_on_core_strands++;
							if ( count_to_minimize_too_many_core_heading_FWY_on_core_strands == 2 ) {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDEFHKNPQRWY" << endl;
								count_to_minimize_too_many_core_heading_FWY_on_core_strands = 0;
							} else {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDEHKNPQR" << endl;
							}
						} else if ( edge_or_core == "edge" ) {
							count_to_minimize_too_many_core_heading_FWY_on_edge_strands++;
							if ( count_to_minimize_too_many_core_heading_FWY_on_edge_strands == 2 ) {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDFNPWY" << endl;
								count_to_minimize_too_many_core_heading_FWY_on_edge_strands = 0;
							} else {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDNP" << endl;
							}
						}
					} else if ( (!write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_) && (write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_) ) {
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "core" ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDEHKNPQR" << endl;
						} else if ( edge_or_core == "edge" ) {
							count_to_minimize_too_many_core_heading_FWY_on_edge_strands++;
							if ( count_to_minimize_too_many_core_heading_FWY_on_edge_strands == 2 ) {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDFNPWY" << endl;
								count_to_minimize_too_many_core_heading_FWY_on_edge_strands = 0;
							} else {
								resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDNP" << endl;
							}
						}
					} else { //(!write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_)
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "core" ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDEHKNPQR" << endl;
						} else if ( edge_or_core == "edge" ) {
							resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCDNP" << endl;
						}
					}
				}
			}
		}
	}

	resfile_stream << "# NOTAA\tCFMWY for loop residues" << endl;
	for ( Size i =1; i<=(pose.size()); i++ ) {
		string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, i);
		if ( edge_or_core == "loop" || edge_or_core == "short_edge" ) {
			resfile_stream << i << "\tA\tEX\t1\tNOTAA\tCFMWY" << endl; // I think that both hairpin-loop and inter-sheet-loop can be treated with 'NOTAA CFMWY'

		}
	}

	resfile_stream.close();
	// <end> make a resfile to design all residues


	// <begin> make a resfile to design surface heading or loop residues
	string surface_loop_resfile_name = pdb_file_name + "_resfile_to_design_surface_heading_or_loop_residues_when_seq_rec_is_bad.txt";
	ofstream surface_loop_resfile_stream;

	surface_loop_resfile_stream.open(surface_loop_resfile_name.c_str());

	surface_loop_resfile_stream << "EX 1 NOTAA C" << endl;
	surface_loop_resfile_stream << "USE_INPUT_SC" << endl;
	surface_loop_resfile_stream << "start" << endl;

	surface_loop_resfile_stream << "# final resfile rule update: 05/06/2014" << endl;
	surface_loop_resfile_stream << "#\tbased on 203 native sandwiches" << endl;

	surface_loop_resfile_stream << "# NOTAA\tCFPWY for surface_heading residues at core strands" << endl;
	surface_loop_resfile_stream << "# NOTAA\tCFWY for surface_heading residues at edge strands" << endl;

	if ( write_resfile_to_minimize_too_much_hydrophobic_surface_ ) {
		surface_loop_resfile_stream << "# PIKAA\tDENQST\tfor every 2nd surface-heading\tresidues " << endl;
	}

	surface_loop_resfile_stream << "# NATAA\tfor core_heading residues " << endl;


	count_to_minimize_too_much_hydrophobic_surface = 0;
	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
		Size residue_begin = bs_of_sw_can_by_sh[ii].get_start();
		Size residue_end = bs_of_sw_can_by_sh[ii].get_end();
		for ( Size residue_num = residue_begin; residue_num <= residue_end; residue_num++ ) {
			{
				string heading = determine_heading_direction_by_vector (struct_id, db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), residue_begin, residue_end, residue_num);
				if ( heading == "surface" ) {
					if ( write_resfile_to_minimize_too_much_hydrophobic_surface_ ) {
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "loop" || edge_or_core == "short_edge" ) {
							continue; // I will write resfile for this residue later
						}
						count_to_minimize_too_much_hydrophobic_surface++;
						if ( count_to_minimize_too_much_hydrophobic_surface == 2 ) {
							surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tPIKAA\tDENQST" << endl;
							count_to_minimize_too_much_hydrophobic_surface = 0;
						} else {
							//        string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
							if ( edge_or_core == "core" ) {
								surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFPWY" << endl;
							} else if ( edge_or_core == "edge" ) {
								surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFWY" << endl;
							}
						}
					} else { //(!write_resfile_to_minimize_too_much_hydrophobic_surface_)
						string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, residue_num);
						if ( edge_or_core == "core" ) {
							surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFPWY" << endl;
						} else if ( edge_or_core == "edge" ) {
							surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tNOTAA\tCFWY" << endl;
						}
					}
				} else if ( heading == "core" ) {
					surface_loop_resfile_stream << residue_num << "\tA\tEX\t1\tNATAA" << endl;
				}
			}
		}
	}

	surface_loop_resfile_stream << "# NOTAA\tCFMWY for loop residues" << endl;
	for ( Size i =1; i<=(pose.size()); i++ ) {
		string edge_or_core = see_edge_or_core_or_loop_or_short_edge (struct_id, db_session, i);
		if ( edge_or_core == "loop" || edge_or_core == "short_edge" ) {
			surface_loop_resfile_stream << i << "\tA\tEX\t1\tNOTAA\tCFMWY" << endl; // I think that both hairpin-loop and inter-sheet-loop can be treated with 'NOTAA CFMWY'
		}
	}

	surface_loop_resfile_stream.close();
	// <end> make a resfile to design surface heading or loop residues

	return 0;
}//write_resfile_to_a_file_when_seq_rec_is_bad


// write_resfile_to_a_file
core::Size
write_resfile_to_a_file(
	string tag,
	StructureID struct_id, // needed argument
	utility::sql_database::sessionOP db_session, // needed argument
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	bool write_resfile_NOT_FWY_on_surface_
)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	// <begin> make a resfile to design all residues
	string resfile_name = pdb_file_name + "_resfile_to_design_all_residues.txt";
	ofstream resfile_stream;

	resfile_stream.open(resfile_name.c_str());

	resfile_stream << "NOTAA C" << endl;
	resfile_stream << "start" << endl;
	resfile_stream << "# final resfile rule update: 07/29/2014" << endl;
	resfile_stream << "# NOTAA CFWY for loops and surface_heading residues" << endl;

	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++ ) { // per each beta-strand
		Size residue_begin = bs_of_sw_can_by_sh[ii].get_start();
		Size residue_end = bs_of_sw_can_by_sh[ii].get_end();
		for ( Size residue_num = residue_begin; residue_num <= residue_end; residue_num++ ) {
			if ( write_resfile_NOT_FWY_on_surface_ ) {
				string  heading = determine_heading_direction_by_vector (struct_id, db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(),  residue_begin,  residue_end,  residue_num);
				if ( heading == "surface" ) {
					resfile_stream << residue_num << " A NOTAA CFWY #surface heading" << endl;
				}
			}
		}
	}

	for ( Size residue_num = 1;  residue_num <=  pose.size(); residue_num++ ) {
		string id = see_edge_or_core_or_loop_or_short_edge (struct_id,  db_session, residue_num);
		if ( id == "loop" ) {
			resfile_stream << residue_num << " A NOTAA CFWY #loop" << endl;
		}
	}

	resfile_stream.close();
	// <end> make a resfile to design all residues

	return 0;
}//write_resfile_to_a_file

} //namespace strand_assembly
} //namespace features
} //namespace protocols
