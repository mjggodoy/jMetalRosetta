// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/graph/Graph.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/Stub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/interpolate.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

static THREAD_LOCAL basic::Tracer TR( "core.pack.dunbrack.cenrot" );

Size const CentroidRotamerSampleData::NUMBER_OF_PARAMS = 7;

Size const SingleResidueCenrotLibrary::N_PHIPSI_BINS = 36;
Real const SingleResidueCenrotLibrary::PHIPSI_BINRANGE = 10.0;
Real const SingleResidueCenrotLibrary::NEUTRAL_PHI = -90; //dun -60 //-60
Real const SingleResidueCenrotLibrary::NEUTRAL_PSI = 130; //dun 60? //-60
Size const SingleResidueCenrotLibrary::RSD_PHI_INDEX = 1;
Size const SingleResidueCenrotLibrary::RSD_PSI_INDEX = 2;
Real const SingleResidueCenrotLibrary::MAX_ROT_ENERGY = 16;
Real const SingleResidueCenrotLibrary::MIN_ROT_PROB = 1e-6;

using namespace basic::options;

SingleResidueCenrotLibrary::SingleResidueCenrotLibrary(AA const aa)
:aa_(aa),max_rot_num(0),ref_energy_(0.0){
	all_rots_bb_.dimension(N_PHIPSI_BINS, N_PHIPSI_BINS);
	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		entropy_.dimension(N_PHIPSI_BINS, N_PHIPSI_BINS);
	}
}

SingleResidueCenrotLibrary::~SingleResidueCenrotLibrary(){}

std::string SingleResidueCenrotLibrary::read_from_file(
	utility::io::izstream & infile,
	bool first_line_three_letter_code_already_read
) {
	std::string const my_name (chemical::name_from_aa( aa() ));
	std::string next_name; //empty string to start
	std::string three_letter_code;
	Real phi, psi;
	Real prob, dis, ang, dih, vdis, vang, vdih;
	Size count;

	//make up a dummy sample for GLY/ALA
	CentroidRotamerSampleData crsd(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
	dummy_sample_.push_back(crsd);

	while ( infile ) {
		/// 1. peek at the line; if it starts with #, skip to the next line.
		char first_char = infile.peek();
		if ( first_char == '#' ) {
			std::string line;
			infile.getline( line );
			continue;
		}

		/// currently no bb-dependent info, do we need that?
		/// if so, how to cal, grab from dunlib directly? average for each chi1?
		/// 2. Read the line.  Format is:
		/// a. three-letter-code
		/// b. count
		/// c. probability
		/// d. distance
		/// e. angle
		/// f. dihedral
		/// g. std

		if ( first_line_three_letter_code_already_read ) {
			/// The last library to read from this file already ate my three letter code;
			/// skip directly to reading phi/psi values...
			first_line_three_letter_code_already_read = false;
		} else {
			infile >> three_letter_code;
			if ( three_letter_code != my_name ) {
				next_name = three_letter_code;
				break;
			}
		}

		infile >> phi >> psi >> count;
		infile >> prob >> dis >> ang >> dih >> vdis >> vang >> vdih;

		if ( count>max_rot_num ) max_rot_num=count;

		Size phi_bin, psi_bin;
		get_phipsi_bins(phi, psi, phi_bin, psi_bin);

		if ( count<=all_rots_bb_(phi_bin,psi_bin).size() ) break;

		// ang,  dih  -- degree
		// vang, vdih -- radian^2
		//convert deg to rad
		CentroidRotamerSampleData crsd(
			prob, dis,
			ang*numeric::constants::r::pi_over_180,
			dih*numeric::constants::r::pi_over_180,
			vdis, vang, vdih);
		all_rots_bb_(phi_bin,psi_bin).push_back(crsd);

		if ( all_rots_bb_(phi_bin,psi_bin).size()!=count ) {
			TR.Debug << "Warning: index (col-2) in cenrotlib file is wrong!" << std::endl;
			TR.Debug << "Warning: " << my_name << ": ndx="
				<< count << " real=" << all_rots_bb_(phi_bin,psi_bin).size()
				<< " phi=" << phi << " psi=" << psi << std::endl;

		}
	}

	//reference energy:
	ref_energy_ = log(float(max_rot_num));

	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		setup_entropy_correction();
	}

	TR.Debug << "Cenrotlib for AA " << my_name << " has been loaded!" << std::endl;
	return next_name;
}

void SingleResidueCenrotLibrary::setup_entropy_correction()
{
	//calculate the entropy for each bb bin
	for ( Size nphi=1; nphi<=N_PHIPSI_BINS; nphi++ ) {
		for ( Size npsi=1; npsi<=N_PHIPSI_BINS; npsi++ ) {
			//for each bin
			Real sum_plogp = 0.0;
			Real sum_p = 0.0;
			for ( Size nrot=1; nrot<=max_rot_num; nrot++ ) {
				sum_p += all_rots_bb_(nphi, npsi)[nrot].prob();
				//energy = -logP
				sum_plogp -= all_rots_bb_(nphi, npsi)[nrot].prob()*all_rots_bb_(nphi, npsi)[nrot].energy();
			}
			entropy_(nphi, npsi) = sum_plogp / sum_p; //in case P is not normalized
			//std::cout << entropy_(nphi, npsi) << std::endl;
		}
	}
}

Real SingleResidueCenrotLibrary::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const {
	return eval_rotameric_energy_deriv(rsd, scratch, true);
}

Real SingleResidueCenrotLibrary::rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const {
	return eval_rotameric_energy_deriv(rsd, scratch, false);
}

Real SingleResidueCenrotLibrary::eval_rotameric_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const {
	Real p(0.0);
	Real factori[10]; //for each rot (max=9)
	const utility::vector1< CentroidRotamerSampleData > rotamer_sample_data( get_rotamer_samples( rsd ) );

	////////////////////////////////////////////////////////
	//entropy correction
	////////////////////////////////////////////////////////
	Real entropy(0.0);
	Real dS_dphi(0.0);
	Real dS_dpsi(0.0);
	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		Size phibin, psibin;
		Size phibin_next, psibin_next;
		Real phi_alpha, psi_alpha;

		Real phi(get_phi_from_rsd(rsd));
		Real psi(get_psi_from_rsd(rsd));

		get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

		//calculate entropy for phi, psi
		Real s00(entropy_(phibin, psibin));
		Real s01(entropy_(phibin, psibin_next));
		Real s10(entropy_(phibin_next, psibin));
		Real s11(entropy_(phibin_next, psibin_next));

		basic::interpolate_bilinear_by_value(
			s00, s10, s01, s11,
			phi_alpha, psi_alpha, PHIPSI_BINRANGE, true,
			entropy, dS_dphi, dS_dpsi
		);
	}
	////////////////////////////////////////////////////////

	for ( Size nr=1; nr<=max_rot_num; nr++ ) {
		Real d_sq(0.0), a_sq(0.0), w_sq(0.0);

		Real cur_ang = rotamer_sample_data[nr].cal_delta_internal_coordinates_squared(rsd, d_sq, a_sq, w_sq);
		Real sin_ang = sin( cur_ang );
		Real weight = sin_ang * sin_ang;

		factori[nr] = rotamer_sample_data[nr].prob() \
			* rotamer_sample_data[nr].norm_factor()  \
			* exp(
			-d_sq / ( 2.0 * rotamer_sample_data[nr].sd_dis() )
			-a_sq / ( 2.0 * rotamer_sample_data[nr].sd_ang() )
			-w_sq * weight / ( 2.0 * rotamer_sample_data[nr].sd_dih() )
		);
		p += factori[nr];
	}

	if ( p<MIN_ROT_PROB ) return MAX_ROT_ENERGY; //too far away

	Real e = -log(p);

	//check meaningless score
	if ( e!=e || e>MAX_ROT_ENERGY+4.0 ) {
		TR.Error << "Dunbrack term calculation fail!" << std::endl;
		TR.Error << rsd.seqpos() << " res: " << rsd.name() << std::endl;
		for ( Size nr=1; nr<=max_rot_num; nr++ ) {
			TR.Error << "nr " << nr << ": " << factori[nr] << " p: " << rotamer_sample_data[nr].prob() << std::endl;
			Real d_sq(0.0), a_sq(0.0), w_sq(0.0);
			Real cur_ang = rotamer_sample_data[nr].cal_delta_internal_coordinates_squared(rsd, d_sq, a_sq, w_sq);
			Real sin_ang = sin( cur_ang );
			TR.Error << "delta: " << d_sq << " " << a_sq << " " << w_sq << std::endl;
			TR.Error << "sigma: " << rotamer_sample_data[nr].sd_dis() << " "
				<< rotamer_sample_data[nr].sd_ang() << " "
				<< rotamer_sample_data[nr].sd_dih() << std::endl;
			TR.Error << "Sin(ang): " << sin_ang << std::endl;
			utility_exit();
		}
		return MAX_ROT_ENERGY;
	}

	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		e += entropy;
	} else {
		e -= ref_energy_;
	}

	if ( !eval_deriv ) return e ;

	//cal devriv
	Real dE_ddis(0.0), dE_dang(0.0), dE_ddih(0.0);
	for ( Size nr=1; nr<=max_rot_num; nr++ ) {

		CentroidRotamerSampleData const &sample(rotamer_sample_data[nr]);

		Real delta_dis, delta_ang, delta_dih;
		Real cur_ang = rotamer_sample_data[nr].cal_delta_internal_coordinates(rsd, delta_dis, delta_ang, delta_dih);
		Real sin_ang = sin(cur_ang);
		Real cos_ang = cos(cur_ang);

		dE_ddis += factori[nr] * delta_dis / sample.sd_dis();
		dE_dang -= factori[nr] * (
			delta_ang / sample.sd_ang() +
			delta_dih * delta_dih * sin_ang * cos_ang / sample.sd_dih()
		);
		dE_ddih += factori[nr] * delta_dih * sin_ang * sin_ang / sample.sd_dih();
	}

	//hack save derivs (dis, ang, dih) in scratch.dE_dchi
	Real4 & dE_dchi(scratch.dE_dchi());
	dE_dchi[1] = dE_ddis / p;
	dE_dchi[2] = dE_dang / p;
	dE_dchi[3] = dE_ddih / p;

	return e ;
}

Real SingleResidueCenrotLibrary::eval_rotameric_energy_bb_dof_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const {
	Real p(0.0);
	Real factori[10]; //for each rot (max=9)
	const utility::vector1< CentroidRotamerSampleData > rotamer_sample_data( get_rotamer_samples( rsd ) );

	////////////////////////////////////////////////////////
	//entropy correction
	////////////////////////////////////////////////////////
	Real entropy(0.0);
	Real dS_dphi(0.0);
	Real dS_dpsi(0.0);
	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		Size phibin, psibin;
		Size phibin_next, psibin_next;
		Real phi_alpha, psi_alpha;

		Real phi(get_phi_from_rsd(rsd));
		Real psi(get_psi_from_rsd(rsd));

		get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

		//calculate entropy for phi, psi
		Real s00(entropy_(phibin, psibin));
		Real s01(entropy_(phibin, psibin_next));
		Real s10(entropy_(phibin_next, psibin));
		Real s11(entropy_(phibin_next, psibin_next));

		basic::interpolate_bilinear_by_value(
			s00, s10, s01, s11,
			phi_alpha, psi_alpha, PHIPSI_BINRANGE, true,
			entropy, dS_dphi, dS_dpsi
		);
	}
	////////////////////////////////////////////////////////

	for ( Size nr=1; nr<=max_rot_num; nr++ ) {
		Real d_sq(0.0), a_sq(0.0), w_sq(0.0);
		Real cur_ang = rotamer_sample_data[nr].cal_delta_internal_coordinates_squared(rsd, d_sq, a_sq, w_sq);
		Real sin_ang = sin( cur_ang );
		Real weight = sin_ang * sin_ang;

		//debug
		//std::cout << "sq: " << nr << "- " << d_sq << " " << a_sq << " " << w_sq << std::endl;
		factori[nr] = rotamer_sample_data[nr].norm_factor() * exp(
			-d_sq / ( 2.0 * rotamer_sample_data[nr].sd_dis() )
			-a_sq / ( 2.0 * rotamer_sample_data[nr].sd_ang() )
			-w_sq * weight / ( 2.0 * rotamer_sample_data[nr].sd_dih() )
		);

		p +=  rotamer_sample_data[nr].prob() * factori[nr];
	}

	if ( p<MIN_ROT_PROB ) return MAX_ROT_ENERGY; //too far away

	Real e = -log(p);

	//check meaningless score
	if ( e!=e || e>MAX_ROT_ENERGY+4.0 ) {
		TR.Error << "Dunbrack term calculation fail!" << std::endl;
		utility_exit();
	}

	utility::fixedsizearray1<Real, FIVE> & dE_dbb(scratch.dE_dbb());

	for ( Size nr=1; nr<=max_rot_num; nr++ ) {
		CentroidRotamerSampleData const &sample(rotamer_sample_data[nr]);
		const Size NDOF = 3;
		Real delta[NDOF];
		Real cur_ang = sample.cal_delta_internal_coordinates(rsd, delta[0], delta[1], delta[2]);
		Real sin_ang = sin(cur_ang);

		Real tmp_phi(0.0), tmp_psi(0.0);
		Size dat_shift = 1; // skip prob
		Size dev_shift = 1 + NDOF; // skip prob + dat
		for ( Size i=0; i<NDOF; i++ ) {

			Real weight = 1.0;
			//for dih(2) only, dis->0, ang->1
			if ( i==2 ) weight = sin_ang * sin_ang;

			Real f = delta[i] * delta[i] * weight;
			Real df_dphi = - 2.0 * delta[i] * sample.deriv_phi_[ i + dat_shift ] * weight;
			Real df_dpsi = - 2.0 * delta[i] * sample.deriv_psi_[ i + dat_shift ] * weight;
			Real g = sample.data_[ i + dev_shift ];
			Real dg_dphi = sample.deriv_phi_[ i + dev_shift ];
			Real dg_dpsi = sample.deriv_psi_[ i + dev_shift ];

			tmp_phi += ( df_dphi * g - dg_dphi * f ) / ( g * g ) / 2.0;
			tmp_psi += ( df_dpsi * g - dg_dpsi * f ) / ( g * g ) / 2.0;
		}

		dE_dbb[RSD_PHI_INDEX] += factori[nr] * ( tmp_phi * sample.prob() - sample.deriv_phi_[0] );
		dE_dbb[RSD_PSI_INDEX] += factori[nr] * ( tmp_psi * sample.prob() - sample.deriv_psi_[0] );
	}

	dE_dbb[RSD_PHI_INDEX] /= p;
	dE_dbb[RSD_PSI_INDEX] /= p;

	if ( option[ OptionKeys::corrections::score::dun_entropy_correction ] ) {
		e += entropy;
		dE_dbb[RSD_PHI_INDEX] += dS_dphi;
		dE_dbb[RSD_PSI_INDEX] += dS_dpsi;
	} else {
		e -= ref_energy_;
	}

	return e ;
}

CentroidRotamerSampleData const & SingleResidueCenrotLibrary::get_closest_rotamer(
	conformation::Residue const & rsd,
	Size &nrot, Real &dis_sq
) const {
	if ( rsd.aa()==core::chemical::aa_gly || rsd.aa()==core::chemical::aa_ala ) {
		nrot = 1;
		dis_sq = 0.0;
		return dummy_sample_[nrot];
	}

	//get the rotamer int
	const utility::vector1< CentroidRotamerSampleData > rotamer_sample_data( get_rotamer_samples( rsd ) );

	/// find the closest rotamer
	Real closest_dis = 9999;
	Size closest_rot = 0;
	for ( Size nr=1; nr<=max_rot_num; nr++ ) {
		Real d2 = rotamer_sample_data[nr].cal_distance_squared(rsd);
		if ( closest_dis>d2 ) {
			closest_dis = d2;
			closest_rot = nr;
		}
	}

	dis_sq = closest_dis;
	nrot = closest_rot;

	return rotamer_sample_data[nrot];
}

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
Real SingleResidueCenrotLibrary::best_rotamer_energy(
	conformation::Residue const &, //rsd,
	bool , //curr_rotamer_only,
	RotamerLibraryScratchSpace & //scratch
) const {
	return 0.0;
}

/// @brief Pick a rotamer for the input residue according to the rotamer probability
/// distribution and assign chi angles to the input rsd.
/// -- currently no perturbation allowed
void SingleResidueCenrotLibrary::assign_random_rotamer_with_bias(
	conformation::Residue const &, //rsd,
	pose::Pose const &, //pose,
	RotamerLibraryScratchSpace &, //scratch,
	numeric::random::RandomGenerator &, //RG,
	ChiVector &, //new_chi_angles,
	bool //perturb_from_rotamer_center
) const {
	TR.Debug << "do not use this function anymore, CHIs are not the only dof" << std::endl;
}

void SingleResidueCenrotLibrary::fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const &, //scorefxn,
	pack::task::PackerTask const & task,
	utility::graph::GraphCOP, //packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const &, //extra_chi_steps,
	bool, //buried,
	rotamers::RotamerVector & rotamers
) const {
	Size const max_rots_that_can_be_built = max_rot_num;
	Size count_rotamers_built = 0;
	Real accumulated_probability(0.0);
	Real const requisit_probability(1.0);

	// ** new logic: since we use interperlate value, should call get_rotamer_samples here
	const utility::vector1< CentroidRotamerSampleData > samples( get_rotamer_samples( existing_residue ) );

	while ( accumulated_probability<requisit_probability-MIN_ROT_PROB ) {
		++count_rotamers_built;

		//build it
		pack::task::ResidueLevelTask const & rtask( task.residue_task( existing_residue.seqpos() ) );
		conformation::ResidueOP rotamer = conformation::ResidueFactory::create_residue(
			*concrete_residue, existing_residue, pose.conformation(), rtask.preserve_c_beta() );

		if ( concrete_residue->aa()==core::chemical::aa_gly
				|| concrete_residue->aa()==core::chemical::aa_ala ) {
			// no other rotamer, no suppose to be here
			rotamers.push_back( rotamer );
			break;
		}

		//get sample from CentroidRotamerSampleData
		utility::fixedsizearray1< Real, 3 > sample;
		samples[count_rotamers_built].assign_best_rotamer(sample);

		core::kinematics::Stub::Vector a, b, c;
		a = rotamer->atom("N").xyz();
		b = rotamer->atom("CA").xyz();
		c = rotamer->atom("CB").xyz();


		core::kinematics::Stub stub(c,b,a);
		core::kinematics::Stub::Matrix M(stub.M * numeric::x_rotation_matrix_radians( sample[3] ));
		M *= numeric::z_rotation_matrix_radians( sample[2] );
		core::kinematics::Stub::Vector new_v(stub.v + sample[1] * M.col_x());
		rotamer->set_xyz("CEN", new_v);

		accumulated_probability += samples[count_rotamers_built].prob();

		// ** discard low prob rot
		// ** 0.16 optimized by sc recovery benchmark
		Real rot_cut =  option[ OptionKeys::packing::cenrot_cutoff ]();

		if ( samples[count_rotamers_built].prob()> rot_cut/max_rots_that_can_be_built ) {
			rotamers.push_back( rotamer );
		}

		if ( count_rotamers_built == max_rots_that_can_be_built ) {
			break;
		}
	}
}

void SingleResidueCenrotLibrary::write_to_file( utility::io::ozstream & ) const
{
	TR.Debug << "cenrot write func" << std::endl;
}

const utility::vector1< CentroidRotamerSampleData >
SingleResidueCenrotLibrary::get_rotamer_samples(
	conformation::Residue const & rsd
) const {

	Size phibin, psibin;
	Size phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;

	Real phi(get_phi_from_rsd(rsd));
	Real psi(get_psi_from_rsd(rsd));

	//get_phipsi_bins(phi, psi, phi_bin, psi_bin);
	get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	//prepare to call interpolate_bilinear_by_value
	utility::vector1< CentroidRotamerSampleData > r00( all_rots_bb_(phibin, psibin) );
	utility::vector1< CentroidRotamerSampleData > r01( all_rots_bb_(phibin, psibin_next) );
	utility::vector1< CentroidRotamerSampleData > r10( all_rots_bb_(phibin_next, psibin) );
	utility::vector1< CentroidRotamerSampleData > r11( all_rots_bb_(phibin_next, psibin_next) );

	utility::vector1< CentroidRotamerSampleData > samples(max_rot_num);

	for ( Size nr=1; nr<=max_rot_num; nr++ ) {
		//load
		r00[nr].private_data_to_public_array();
		r01[nr].private_data_to_public_array();
		r10[nr].private_data_to_public_array();
		r11[nr].private_data_to_public_array();

		//interp
		//CentroidRotamerSampleData sample;
		for ( Size i=0; i<7; i++ ) {
			if ( i==3 ) {
				//only for i==3, dih, treat as pbc 360!! unit!!
				basic::interpolate_bilinear_by_value(
					r00[nr].data_[i]*numeric::constants::r::radians_to_degrees,
					r10[nr].data_[i]*numeric::constants::r::radians_to_degrees,
					r01[nr].data_[i]*numeric::constants::r::radians_to_degrees,
					r11[nr].data_[i]*numeric::constants::r::radians_to_degrees,
					phi_alpha, psi_alpha, PHIPSI_BINRANGE, true,
					samples[nr].data_[i], samples[nr].deriv_phi_[i], samples[nr].deriv_psi_[i]
				);
				samples[nr].data_[i]*=numeric::constants::r::pi_over_180;
				samples[nr].deriv_phi_[i]*=numeric::constants::r::pi_over_180;
				samples[nr].deriv_psi_[i]*=numeric::constants::r::pi_over_180;
			} else {
				basic::interpolate_bilinear_by_value(
					r00[nr].data_[i], r10[nr].data_[i],
					r01[nr].data_[i], r11[nr].data_[i],
					phi_alpha, psi_alpha, PHIPSI_BINRANGE, false,
					samples[nr].data_[i], samples[nr].deriv_phi_[i], samples[nr].deriv_psi_[i]
				);
			}
		}

		//save
		samples[nr].public_array_to_private_data();
	}

	return samples;
}

void SingleResidueCenrotLibrary::get_phipsi_bins(
	Real phi,
	Real psi,
	Size & phibin,
	Size & psibin,
	Size & phibin_next,
	Size & psibin_next,
	Real & phi_alpha,
	Real & psi_alpha
) const {
	bin_angle( -180.0, PHIPSI_BINRANGE, 360.0, N_PHIPSI_BINS, basic::periodic_range( phi, 360 ), phibin, phibin_next, phi_alpha );
	bin_angle( -180.0, PHIPSI_BINRANGE, 360.0, N_PHIPSI_BINS, basic::periodic_range( psi, 360 ), psibin, psibin_next, psi_alpha );

	verify_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next );
}

void SingleResidueCenrotLibrary::get_phipsi_bins(
	Real phi,
	Real psi,
	Size & phibin,
	Size & psibin
) const {
	Size phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );
}

Real SingleResidueCenrotLibrary::get_phi_from_rsd(
	conformation::Residue const & rsd
) const {
	debug_assert( rsd.is_protein() );
	if ( rsd.is_lower_terminus() ) return NEUTRAL_PHI;
	else return rsd.mainchain_torsion( RSD_PHI_INDEX );
}

Real SingleResidueCenrotLibrary::get_psi_from_rsd(
	conformation::Residue const & rsd
) const {
	debug_assert( rsd.is_protein() );
	if ( rsd.is_upper_terminus() ) return NEUTRAL_PSI;
	else return rsd.mainchain_torsion( RSD_PSI_INDEX );
}

void SingleResidueCenrotLibrary::verify_phipsi_bins(
	Real phi,
	Real psi,
	Size const phibin,
	Size const psibin,
	Size const phibin_next,
	Size const psibin_next
) const {
	if ( ( phibin < 1 || phibin > 36 ) || (psibin < 1 || psibin > 36 ) ||
			( phibin_next < 1 || phibin_next > 36 ) || (psibin_next < 1 || psibin_next > 36 ) ) {
		std::cerr << "ERROR: phi/psi bin out of range: " <<
			aa() << " " << phi << " " << psi;
		std::cerr << phibin << " " << phibin_next << " " << psibin << " " << psibin_next <<  std::endl;
		utility_exit();
	}
}


/// CentroidRotamerSampleData
Real CentroidRotamerSampleData::cal_delta_internal_coordinates_squared(
	const conformation::Residue & rsd,
	Real & d_sq, Real & a_sq, Real & w_sq ) const
{
	Real ddis, dang, ddih;
	Real cur_ang = cal_delta_internal_coordinates(rsd, ddis, dang, ddih);

	d_sq = ddis * ddis;
	a_sq = dang * dang;
	w_sq = ddih * ddih;

	return cur_ang;
}

Real CentroidRotamerSampleData::cal_delta_internal_coordinates(
	const conformation::Residue & rsd,
	Real & ddis, Real & dang, Real & ddih ) const
{
	using namespace numeric::constants::f;

	//get rsd.sidechain coordinates
	core::kinematics::Stub::Vector const a(rsd.atom("N").xyz());
	core::kinematics::Stub::Vector const b(rsd.atom("CA").xyz());
	core::kinematics::Stub::Vector const c(rsd.atom("CB").xyz());
	core::kinematics::Stub::Vector const d(rsd.atom("CEN").xyz());

	//get centroid_rot int
	Real cur_ang = numeric::constants::r::pi-numeric::angle_radians(b,c,d);
	ddis = (d-c).length() - distance_;
	dang = cur_ang - angle_;
	ddih = basic::periodic_range( numeric::dihedral_radians(a,b,c,d) - dihedral_, pi_2 );

	return cur_ang;
}

Real CentroidRotamerSampleData::cal_distance_squared( const conformation::Residue & rsd ) const
{
	//get rsd.sidechain coordinates
	core::kinematics::Stub::Vector const a(rsd.atom("N").xyz());
	core::kinematics::Stub::Vector const b(rsd.atom("CA").xyz());
	core::kinematics::Stub::Vector const c(rsd.atom("CB").xyz());
	core::kinematics::Stub::Vector const d(rsd.atom("CEN").xyz());

	//get centroid_rot int
	DOF3 sample;
	sample[1] = (d-c).length();
	sample[2] = numeric::constants::r::pi-numeric::angle_radians(b,c,d);
	sample[3] = numeric::dihedral_radians(a,b,c,d);

	//std::cout << sample[1] << " " << sample[2] << " " << sample[3] << std::endl;
	// DOF3 sample;
	// core::kinematics::Stub::Vector const cenxyz(rsd.atom("CEN").xyz());
	// sample[1] = cenxyz[1];
	// sample[2] = cenxyz[2];
	// sample[3] = cenxyz[3];
	return cal_distance_squared(sample, false);
}

Real CentroidRotamerSampleData::cal_distance( const DOF3 & sample, bool use_xyz) const
{
	return sqrt(cal_distance_squared(sample, use_xyz));
}

Real CentroidRotamerSampleData::cal_distance_squared( const DOF3 & sample, bool use_xyz) const
{
	Real x1, y1, z1;

	if ( use_xyz ) {
		x1 = sample[1];
		y1 = sample[2];
		z1 = sample[3];
	} else {
		Real ang1 = numeric::constants::r::pi - sample[2];
		x1 = sample[1]*sin(ang1)*cos(sample[3]);
		y1 = sample[1]*sin(ang1)*sin(sample[3]);
		z1 = sample[1]*cos(ang1);
	}

	Real ang2 = numeric::constants::r::pi - angle_;
	Real x2 = distance_*sin(ang2)*cos(dihedral_);
	Real y2 = distance_*sin(ang2)*sin(dihedral_);
	Real z2 = distance_*cos(ang2);

	Real dx = x1-x2;
	Real dy = y1-y2;
	Real dz = z1-z2;
	return dx*dx+dy*dy+dz*dz;
}

void CentroidRotamerSampleData::assign_random_rotamer(
	DOF3 &sample, numeric::random::RandomGenerator & ) const
{
	//temp just return best
	sample[1] = distance_;
	sample[2] = angle_;
	sample[3] = dihedral_;
}

void CentroidRotamerSampleData::assign_best_rotamer( DOF3 & sample ) const
{
	sample[1] = distance_;
	sample[2] = angle_;
	sample[3] = dihedral_;
}

void CentroidRotamerSampleData::private_data_to_public_array() {
	data_[0] = prob_;
	data_[1] = distance_;
	data_[2] = angle_;
	data_[3] = dihedral_;
	data_[4] = sd_dis_;
	data_[5] = sd_ang_;
	data_[6] = sd_dih_;
}
void CentroidRotamerSampleData::public_array_to_private_data() {
	prob_   = data_[0];
	distance_  = data_[1];
	angle_   = data_[2];
	dihedral_  = data_[3];
	sd_dis_  = data_[4];
	sd_ang_  = data_[5];
	sd_dih_  = data_[6];
	energy_ = -log(prob_);
	norm_factor_ = 1.0;
}

} // namespace cenrot
} // namespace dunbrack
} // namespace scoring
} // namespace core

