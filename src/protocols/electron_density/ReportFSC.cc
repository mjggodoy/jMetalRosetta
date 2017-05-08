// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file

#include <protocols/electron_density/ReportFSC.hh>
#include <protocols/electron_density/ReportFSCCreator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/cryst/util.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.ReportFSC" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;


ReportFSC::ReportFSC() : moves::Mover() {
	init();
}

void ReportFSC::init() {
	res_low_=1000.0;
	res_high_=0.0;
	nresbins_=50;
	asymm_only_=false;
	mask_=true;
	testmap_ = nullptr;
	bin_squared_=false;
}

void ReportFSC::apply(core::pose::Pose & pose) {
	// pose->poseCoords
	poseCoords litePose;
	core::conformation::symmetry::SymmetryInfoOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( asymm_only_ && symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for ( core::Size j = 1; j <= natoms; ++j ) {
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = pose.pdb_info()->temperature( i, j );
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

			litePose.push_back( coord_j );
		}
	}

	utility::vector1< core::Real > modelmap1FSC(nresbins_,1.0), modelmap2FSC(nresbins_,1.0);
	utility::vector1< core::Real > modelmap1Error(nresbins_,1.0), modelmap2Error(nresbins_,1.0);
	core::Real fsc1=0.0, fsc2=0.0;
	core::Size bincountsum=0;

	// train map
	ObjexxFCL::FArray3D< double > rhoC, rhoMask;
	ObjexxFCL::FArray3D< std::complex<double> > FrhoC, FrhoO;
	core::scoring::electron_density::getDensityMap().calcRhoC( litePose, 0, rhoC, rhoMask );
	numeric::fourier::fft3(rhoC, FrhoC);
	numeric::fourier::fft3(core::scoring::electron_density::getDensityMap().get_data(), FrhoO);

	utility::vector1< core::Size > resobin_counts;
	utility::vector1< core::Real > resobins;
	core::scoring::electron_density::getDensityMap().getResolutionBins(nresbins_, 1.0/res_low_, 1.0/res_high_, resobins, resobin_counts, bin_squared_);

	core::scoring::electron_density::getDensityMap().getFSC( FrhoC, FrhoO, nresbins_, 1.0/res_low_, 1.0/res_high_, modelmap1FSC, bin_squared_ );
	for ( Size i=1; i<=modelmap1FSC.size(); ++i ) {
		fsc1 += resobin_counts[i] * modelmap1FSC[i];
		bincountsum += resobin_counts[i];
	}
	fsc1 /= bincountsum;

	numeric::xyzVector<core::Real> apix = core::scoring::electron_density::getDensityMap().get_voxel_spacing(  );
	numeric::xyzVector<core::Real> origin = core::scoring::electron_density::getDensityMap().getOrigin(  );

	if ( testmap_ && testmap_->isMapLoaded() ) {
		// set voxel spacing and origin of testmap
		testmap_->set_voxel_spacing( apix );
		testmap_->setOrigin( origin );

		numeric::fourier::fft3(testmap_->get_data(), FrhoO);
		testmap_->getFSC( FrhoC, FrhoO, nresbins_, 1.0/res_low_, 1.0/res_high_, modelmap2FSC, bin_squared_ );
		for ( Size i=1; i<=modelmap2FSC.size(); ++i ) {
			fsc2 += resobin_counts[i] * modelmap2FSC[i];
		}
		fsc2 /= bincountsum;
	}

	// tag
	core::io::RemarkInfo remark;
	std::ostringstream oss;
	core::Real maskwidth =  core::scoring::electron_density::getDensityMap().getAtomMask();

	if ( mask_ ) {
		oss << "FSC[mask=" << maskwidth << "](" << res_low_ << ":" << res_high_ << ") = " << fsc1;
	} else {
		oss << "FSC(" << res_low_ << ":" << res_high_ << ") = " << fsc1;
	}

	if ( testmap_ && testmap_->isMapLoaded() ) {
		oss << " / " << fsc2;
	}

	TR << oss.str() << std::endl;
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	// now tag map parameters (which may have been refined)
	oss.str(""); oss.clear();
	oss << "voxel spacing = " << apix[0] << ", " << apix[1] << ", " << apix[2];
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	oss.str(""); oss.clear();
	oss << "origin = " << origin[0] << ", " << origin[1] << ", " << origin[2];
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
ReportFSC::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( tag->hasOption("res_low") ) {
		res_low_ = tag->getOption<core::Real>("res_low");
	}
	if ( tag->hasOption("res_high") ) {
		res_high_ = tag->getOption<core::Real>("res_high");
	}
	if ( tag->hasOption("nresbins") ) {
		nresbins_ = tag->getOption<core::Size>("nresbins");
	}
	if ( tag->hasOption("asymm_only") ) {
		asymm_only_ = tag->getOption<bool>("asymm_only");
	}
	if ( tag->hasOption("mask") ) {
		mask_ = tag->getOption<bool>("mask");
	}
	if ( tag->hasOption("bin_squared") ) {
		bin_squared_ = tag->getOption<bool>("bin_squared");
	}
	if ( tag->hasOption("testmap") ) {
		std::string mapfile = tag->getOption<std::string>("testmap");
		core::Real mapreso = option[ edensity::mapreso ]();
		core::Real mapsampling = option[ edensity::grid_spacing ]();
		std::cerr << "Loading alternate density map " << mapfile << std::endl;
		testmap_ = core::scoring::electron_density::ElectronDensityOP( new core::scoring::electron_density::ElectronDensity() );
		testmap_->readMRCandResize( mapfile , mapreso , mapsampling );
	}
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ReportFSCCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ReportFSC );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ReportFSCCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ReportFSC::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ReportFSC::mover_name()
// XRW TEMP {
// XRW TEMP  return "ReportFSC";
// XRW TEMP }

std::string ReportFSC::get_name() const {
	return mover_name();
}

std::string ReportFSC::mover_name() {
	return "ReportFSC";
}

void ReportFSC::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute("res_low", xsct_real, "Low resolution cut off (Angstroms)")
		+ XMLSchemaAttribute("res_high", xsct_real, "High resultion cut off (Angstroms)")
		+ XMLSchemaAttribute("nresbins", xsct_non_negative_integer, "XRW TO DO, default=50")
		+ XMLSchemaAttribute("asymm_only", xsct_rosetta_bool, "XRW TO DO, default=false")
		+ XMLSchemaAttribute("mask", xsct_rosetta_bool, "XRW TO DO, default=true")
		+ XMLSchemaAttribute("bin_squared", xsct_rosetta_bool, "XRW TO DO, default=false")
		+ XMLSchemaAttribute("testmap", xs_string, "XRW TO DO");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Calculate a model-map FSC (Fourier Shell Correlation), and integrate over the input resolution ranges.", attlist );
}

std::string ReportFSCCreator::keyname() const {
	return ReportFSC::mover_name();
}

protocols::moves::MoverOP
ReportFSCCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReportFSC );
}

void ReportFSCCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReportFSC::provide_xml_schema( xsd );
}


}
}
