// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit Headers
#include <protocols/topology_broker/util.hh>

// Package Headers
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/TopologyClaimerFactory.hh>

// Project Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/util.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/LoopFragmentClaimer.hh>
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/pose/PDBInfo.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <sstream>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/topology_broker/weights/LargeFragWeight.hh>
#include <protocols/topology_broker/weights/SmallFragWeight.hh>
#include <protocols/topology_broker/weights/SmoothFragWeight.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#undef small
#endif

//Auto Headers
static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace utility::excn;
using namespace core;
using namespace core::fragment;

struct CmdLineData {
	CmdLineData() : b_has_constraint_claimer( false ) {};
	bool b_has_constraint_claimer;
};

class FragmentContainer {

public:

	FragmentContainer( FragSetOP large, FragSetOP small, std::string label ):
		frags_large_ (std::move( large )),
		frags_small_ (std::move( small )),
		label_ (std::move( label ))
	{}

	FragSetOP frags_large(){
		return frags_large_;
	}

	FragSetOP frags_small(){
		return frags_small_;
	}

	std::string label(){
		return label_;
	}

	void set_label( std::string label ){
		label_ = label;
	}

private:
	FragSetOP frags_large_;
	FragSetOP frags_small_;
	std::string label_;


};

//some modifiers for the fragment reading can be read from the stream
core::fragment::FragSetOP read_frags( std::istream& is, core::fragment::FragmentIO& io ) {
	std::string file;
	while ( is >> file ) {
		if ( file[0] == '#' ) {
			getline( is, file );
			continue;
		}
		if ( file == "NTOP" ) {
			Size ntop;
			is >> ntop;
			io.set_top_frag_num( ntop );
		} else if ( file == "NCOPY" ) {
			Size ncopy;
			is >> ncopy;
			io.set_ncopies( ncopy );
		} else if ( file == "ANNOTATE" ) {
			bool yesno;
			is >> yesno;
			io.set_read_annotation( yesno );
		} else {
			return io.read_data( file );
		}
	}
	runtime_assert( 0 ); // shouldn't be able to get here
	return nullptr;
}

void add_claims_from_stream( TopologyBroker& broker, std::istream& is ,  CmdLineData& cmdline_data, utility::vector1< FragmentContainer> & fragment_list ) {
	std::string tag;
	using namespace core::fragment;
	using namespace basic::options;

	while ( is >> tag ) {
		if ( tag[ 0 ]=='#' ) {
			getline( is, tag );
			continue;
		} else if ( tag == "USE_INPUT_POSE" ) {
			broker.use_job_pose( true );
		} else if ( tag == "NO_USE_INPUT_POSE" ) {
			broker.use_job_pose( false );
		} else if ( tag == "CLAIMER" ) {
			std::string claim_type;
			is >> claim_type;
			tr.Debug << "instantiate CLAIMER: " << claim_type << std::endl;
			if ( claim_type == ConstraintClaimer::_static_type_name() ) cmdline_data.b_has_constraint_claimer = true;
			TopologyClaimerOP claim = TopologyClaimerFactory::get_instance()->newTopologyClaimer( claim_type );
			broker.add( claim ); //register first, so that messages work immediatly!
			is >> *claim;
		} else if ( tag == "ABINITIO_FRAGS" ) {
			tr.Debug << "install Abinitio fragments " << std::endl;

			FragSetOP frags_large;
			FragSetOP frags_small;
			std::string frags_label;

			//short cut for the FragmentClaimers need for abinitio runs:
			//std::string file;
			while ( is >> tag && tag.substr(0,4) != "END_" ) {
				if ( tag == "LARGE" ) {
					FragmentIO frag_io(option[ OptionKeys::abinitio::number_9mer_frags ](),
						option[ OptionKeys::frags::nr_large_copies ](),
						option[ OptionKeys::frags::annotate ]() );
					frags_large = read_frags( is, frag_io );

				} else if ( tag == "SMALL" ) {
					FragmentIO frag_io(option[ OptionKeys::abinitio::number_3mer_frags ](),
						1,
						option[ OptionKeys::frags::annotate ]() );
					frags_small = read_frags( is, frag_io );
				} else if ( tag == "LABEL" ) {
					is >> frags_label;
				} else if ( tag[0] == '#' ) {
					getline( is, tag );
				}
			}
			if ( tr.Trace.visible() ) {
				tr.Trace << "LARGE Fragments read for label: " << frags_label << std::endl;// << *fragment.frags_large() << std::endl;
				tr.Trace << "SMALL Fragments read for label: " << frags_label << std::endl;// << *fragment.frags_small_ << std::endl;
			}

			if ( !frags_large || !frags_small ) {
				throw utility::excn::EXCN_BadInput("ABINITIO_FRAGS macro used in topology broker setup without both 'LARGE' and 'SMALL' tags specifying 9-mer and 3-mer frags.");
			}

			fragment_list.add_back( FragmentContainer(frags_large, frags_small, frags_label) );
		} else {
			throw utility::excn::EXCN_BadInput(" unrecognized tag " + tag );
		}
	} //while is tag


}

void add_claims_from_file( TopologyBroker& broker, std::string const& file , CmdLineData& cmdline_data, utility::vector1< FragmentContainer> & fragment_list ) {
	utility::io::izstream is;
	if ( file != "NO_SETUP_FILE" ) {
		is.open( file );
		if ( !is.good() ) throw EXCN_FileNotFound( file );
	}

	try {
		add_claims_from_stream( broker, is , cmdline_data, fragment_list );
	} catch ( EXCN_BadInput &excn ) {
		throw EXCN_BadInput( excn.msg() + " occurred when reading file "+file ); //of course I loose the speciality of the EXCEPTION
	}
//that might be just eof check for is.fail() ??? don't check...not my problem ?
}

void add_cmdline_claims( TopologyBroker& broker, bool const do_I_need_frags ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using simple_moves::ClassicFragmentMover;
	using simple_moves::SmoothFragmentMover;
	using simple_moves::GunnCost;
	using weights::LargeFragWeight;
	using weights::SmallFragWeight;
	using weights::SmoothFragWeight;

	CmdLineData cmdline_data;

	utility::vector1< FragmentContainer > input_fragments;


	if ( option[ OptionKeys::broker::setup ].user() ) {
		FileVectorOption& files( option[ OptionKeys::broker::setup ] );
		for ( Size i = 1; i<= files.size(); ++i ) {
			add_claims_from_file( broker, files[ i ], cmdline_data, input_fragments );
		}
	}

	// If no fragments have been specified in the input file, read from command line.
	if ( input_fragments.size() == 0 && do_I_need_frags ) {
		tr.Debug << "Reading FragmentFiles from cmd." << std::endl;

		FragSetOP large;
		FragSetOP small;

		core::fragment::read_std_frags_from_cmd( large, small );
		input_fragments.add_back( FragmentContainer( large, small, "DEFAULT" ) );
	}
	if ( input_fragments.size() == 0 ) {
		if ( do_I_need_frags ) {
			throw utility::excn::EXCN_BadInput( "expected LARGE and SMALL fragment sets in ABINITIO_FRAGS section or via command line ");
		} else {
			return;
		}
	}


	//Add 5 claimers for each FragmentContainer in <input_fragments>
	for ( auto & input_fragment : input_fragments ) {

		if ( input_fragment.label() == "" ) {
			input_fragment.set_label( "DEFAULT" );
		}

		tr.Info << " Adding claimers for fragments with label '" << input_fragment.label() << "'." <<std::endl;

		simple_moves::ClassicFragmentMoverOP
			bml( new ClassicFragmentMover(input_fragment.frags_small()) ),
			bms( new ClassicFragmentMover(input_fragment.frags_small()) ),
			sms( new SmoothFragmentMover (input_fragment.frags_small(), simple_moves::FragmentCostOP( new GunnCost ) ) );

		broker.add( TopologyClaimerOP( new FragmentClaimer( bml, "LargeFrags", weights::AbinitioMoverWeightOP( new LargeFragWeight ), input_fragment.label(), input_fragment.frags_large() ) ) );
		broker.add( TopologyClaimerOP( new FragmentClaimer( bms, "SmallFrags", weights::AbinitioMoverWeightOP( new SmallFragWeight ), input_fragment.label(), input_fragment.frags_small() ) ) );
		broker.add( TopologyClaimerOP( new FragmentClaimer( sms, "SmoothFrags", weights::AbinitioMoverWeightOP( new SmoothFragWeight ), input_fragment.label(), input_fragment.frags_small() ) ) );
		broker.add( TopologyClaimerOP( new LoopFragmentClaimer( input_fragment.frags_small(), input_fragment.label() ) ));
		core::fragment::SecondaryStructureOP ss_def( new core::fragment::SecondaryStructure(*(input_fragment.frags_small()), false /*no JustUseCentralResidue */ ) );
		broker.add( TopologyClaimerOP( new CutBiasClaimer( *ss_def, input_fragment.label() ) ) );

		bms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] ); //default is 30.0
		bml->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );
		sms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );

	}

	//check if there are any SequenceClaimers: if not make at least one from FASTA file
	if ( !broker.has_sequence_claimer() ) {
		using namespace basic::options::OptionKeys;

		//std::string label;
		std::string sequence;
		if ( option[ in::file::fasta ].user() ) {
			std::string filename = option[ in::file::fasta ]()[1];
			sequence = core::sequence::read_fasta_file( filename  )[1]->sequence();
			// label = core::sequence::read_fasta_file( filename )[1]->id();
			broker.add( TopologyClaimerOP( new SequenceClaimer ( sequence, "DEFAULT", core::chemical::CENTROID ) ) );

			tr.Debug << "Read command-line to instantiate a default SequenceClaimer. Read sequence is '"
				<< sequence << "' from command-line (-in:file:fasta or similar) specified file "
				<< filename << std::endl;
		} else if ( option[ in::file::native ].user() ) {
			pose::PoseOP native_pose( new pose::Pose );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
			sequence = native_pose->annotated_sequence();
			utility::vector1< core::pose::PoseOP > chain_poses = native_pose->split_by_chain();

			broker.add( TopologyClaimerOP( new SequenceClaimer ( sequence, "DEFAULT", core::chemical::CENTROID ) ) );

		} else {
			throw utility::excn::EXCN_BadInput("Error: can't read sequence! Use -in::file::fasta sequence.fasta or -in::file::native native.pdb!");
		}

		//broker.add( new SequenceClaimer ( sequence, core::chemical::CENTROID ) );
	}

	if ( !cmdline_data.b_has_constraint_claimer && basic::options::option[ constraints::cst_file ].user() ) {
		tr.Info << "add ConstraintClaimer to account for cmd-line constraints" << std::endl;
		broker.add( TopologyClaimerOP( new ConstraintClaimer( true /* read from cmd-line */ ) ) );
	}
	if ( !cmdline_data.b_has_constraint_claimer && basic::options::option[ constraints::cst_fa_file ].user() ) {
		tr.Info << "add fa ConstraintClaimer to account for cmd-line constraints" << std::endl;
		broker.add( TopologyClaimerOP( new ConstraintClaimer( true /* read from cmd-line */, false /*centroid*/, true /*true*/ ) ) );
	}
}

}
}
