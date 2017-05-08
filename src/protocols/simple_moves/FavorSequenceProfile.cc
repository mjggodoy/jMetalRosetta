// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FavorSequenceProfile.cc
/// @brief Add a SequenceProfileConstraint to a pose.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Unit Headers
#include <protocols/simple_moves/FavorSequenceProfile.hh>
#include <protocols/simple_moves/FavorSequenceProfileCreator.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/selection.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/id/SequenceMapping.hh>


#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/exit.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace std;
using namespace core::scoring;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.FavorSequenceProfile" );

// XRW TEMP std::string FavorSequenceProfileCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FavorSequenceProfile::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FavorSequenceProfileCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FavorSequenceProfile );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FavorSequenceProfile::mover_name() {
// XRW TEMP  return "FavorSequenceProfile";
// XRW TEMP }


FavorSequenceProfile::FavorSequenceProfile( ) :
	protocols::moves::Mover( "FavorSequenceProfile" ),
	weight_( 1.0 ),
	use_current_(false),
	matrix_("BLOSUM62"),
	scaling_("prob"),
	chain_(0),
	string_exclude_resnums_("")
{}

void
FavorSequenceProfile::set_weight( core::Real weight ) {
	weight_ = weight;
}

void
FavorSequenceProfile::set_sequence( core::sequence::Sequence & seq, std::string matrix) {
	if ( ref_profile_ ) {
		TR.Warning << "Overwriting existing profile in FavorSequenceProfile." << std::endl;
	}
	ref_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile );
	ref_profile_->generate_from_sequence(seq, matrix);
}

void
FavorSequenceProfile::set_profile( core::sequence::SequenceProfile & profile) {
	if ( ref_profile_ ) {
		TR.Warning << "Overwriting existing profile in FavorSequenceProfile." << std::endl;
	}
	ref_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile( profile ) );
}

void
FavorSequenceProfile::set_scaling( std::string const & scaling ) {
	if ( scaling != "prob" && scaling != "none" && scaling != "global" ) {
		utility_exit_with_message("Scaling in FavorSequenceProfile must be one of 'prob', 'none', or 'global'.");
	}
	scaling_ = scaling;
}

void
FavorSequenceProfile::apply( core::pose::Pose & pose )
{
	core::sequence::SequenceProfileOP profile;
	if ( use_current_ ) {
		core::sequence::Sequence seq(pose);
		profile = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile );
		profile->generate_from_sequence(seq, matrix_);
	} else {
		runtime_assert( ref_profile_ != nullptr );
		profile = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile( *ref_profile_) );
	}

	if ( scaling_ == "prob" ) {
		profile->convert_profile_to_probs( 1.0 );
	} else if ( scaling_ == "none" ) {
		// pass
	} else if ( scaling_ == "global" ) {
		profile->global_auto_rescale();
	} else {
		utility_exit_with_message("Unrecognized scaling type '" + scaling_ + "' in FavorSequenceProfile.");
	}

	if ( weight_ != 1.0 ) {
		profile->rescale(weight_);
	}
	//using varibles for start/stop in case a sequence for only one chain was specified
	core::Size start_seq = 1;
	core::Size stop_seq = pose.size();
	core::id::SequenceMappingOP smap( new core::id::SequenceMapping() );
	if ( chain_ > 1 ) {
		*smap = core::id::SequenceMapping::identity( pose.size() );
		smap->set_offset(pose.conformation().chain_begin( chain_ )-1);
		smap->show();
	} else { //set smap to identity
		*smap = core::id::SequenceMapping::identity( pose.size() );
	}

	utility::vector1<bool> use_all_residues;
	use_all_residues.resize(pose.size(),true);

	if ( string_exclude_resnums_.length()!=0 ) {
		set< core::Size > const res_vec( core::pose::get_resnum_list( string_exclude_resnums_, pose ) );
		for ( core::Size const res : res_vec ) {
			TR<<"Turning off res_type_constraint weight for "<<res<<std::endl;
			use_all_residues[res]=false;
		}
	}

	for ( core::Size seqpos( start_seq ), end( stop_seq ); seqpos <= end; ++seqpos ) {
		if ( use_all_residues[seqpos] ) {
			pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( pose, seqpos, profile, smap ) ) ) );
		}
	}
}

void
FavorSequenceProfile::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose)
{
	weight_ = tag->getOption<core::Real>( "weight", 1 );

	if ( tag->hasOption("scorefxns") ) {
		std::string const sf_val( tag->getOption<std::string>("scorefxns") );
		typedef utility::vector1< std::string > StringVec;
		StringVec const sf_keys( utility::string_split( sf_val, ',' ) );
		for ( auto const & sf_key : sf_keys ) {
			ScoreFunctionOP scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", sf_key ) );
			if ( scorefxn->get_weight( res_type_constraint ) == 0.0 ) {
				scorefxn->set_weight( res_type_constraint, 1 );
				TR<<"Turning on res_type_constraint weight in scorefxn "<<sf_key<<std::endl;
			}
		}
	}

	core::Size num_struct(0);
	if ( tag->getOption< bool >( "use_native", false ) ) ++num_struct;
	if ( tag->getOption< bool >( "use_fasta", false ) ) ++num_struct;
	if ( tag->getOption< bool >( "use_starting", false ) ) ++num_struct;
	if ( tag->getOption< bool >( "use_current", false ) ) ++num_struct;
	if ( tag->hasOption("pdbname") ) ++num_struct;

	if ( ! num_struct &&  ! tag->hasOption("pssm") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Must set one of 'pssm', 'use_native', 'use_fasta', 'use_starting', 'use_current', or 'pdbname' in FavorSequenceProfile");
	}
	if ( num_struct && tag->hasOption("pssm") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Cannot set both 'pssm' and one of 'use_native', 'use_fasta', 'use_starting', 'use_current', or 'pdbname' in FavorSequenceProfile");
	}
	if ( num_struct > 1 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can only set one of 'use_native', 'use_fasta', 'use_starting', 'use_current', or 'pdbname' in FavorSequenceProfile");
	}
	if ( tag->hasOption("matrix") && tag->hasOption("pssm")  ) {
		TR.Warning << "WARNING In option matrix not used with pssm specification." << std::endl;
	}
	if ( tag->hasOption("chain") ) {
		chain_ = tag->getOption<core::Size>("chain", 0 );
	}

	set_scaling( tag->getOption< std::string >( "scaling", "prob" ) );

	matrix_ = tag->getOption< std::string >( "matrix", "BLOSUM62" );
	if ( tag->getOption< bool >( "use_native", false ) ) {
		core::pose::Pose nat_pose;
		core::import_pose::pose_from_file( nat_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		core::sequence::Sequence seq(nat_pose.sequence(), basic::options::option[ basic::options::OptionKeys::in::file::native ]);
		set_sequence( seq, matrix_ );
	}
	if ( tag->getOption< bool >( "use_fasta", false ) ) {
		std::string fasta_file( core::sequence::read_fasta_file_str( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1] );
		std::string name("unknown");
		core::sequence::Sequence seq( fasta_file, name );
		std::cout << seq << std::endl;
		set_sequence( seq, matrix_ );
	}
	if ( tag->getOption< bool >( "use_starting", false ) ) {
		core::sequence::Sequence seq(pose);
		set_sequence( seq, matrix_ );
	}
	if ( tag->getOption< bool >( "use_current", false ) ) {
		use_current_ = true;
	}
	if ( tag->hasOption("pdbname") ) {
		core::pose::Pose ref_pose;
		core::import_pose::pose_from_file( ref_pose, tag->getOption<std::string>( "pdbname" ) , core::import_pose::PDB_file);
		core::sequence::Sequence seq(ref_pose.sequence(), tag->getOption<std::string>( "pdbname" ) );
		set_sequence( seq, matrix_ );
	}

	if ( tag->hasOption("pssm") ) {
		ref_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile );
		ref_profile_->read_from_file( tag->getOption< std::string >( "pssm" ) );
	}

	if ( tag->hasOption("exclude_resnums") ) {
		string_exclude_resnums_ = tag->getOption<std::string>( "exclude_resnums" );
	}
}

std::string FavorSequenceProfile::get_name() const {
	return mover_name();
}

std::string FavorSequenceProfile::mover_name() {
	return "FavorSequenceProfile";
}

void FavorSequenceProfile::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;

	XMLSchemaRestriction scaling_type;
	scaling_type.name( "favor_seqprof_scaling_type" );
	scaling_type.add_restriction( xsr_enumeration, "prob" );
	scaling_type.add_restriction( xsr_enumeration, "none" );
	scaling_type.add_restriction( xsr_enumeration, "global" );
	xsd.add_top_level_element( scaling_type );

	AttributeList attlist;
	attlist
		+ Attr::attribute_w_default( "weight", xsct_real, "Adjust the post-scaling strength of the constraints.", "1" )
		+ Attr( "scorefxns", xs_string,
		"Convenience feature to automatically set res_type_constraint to 1"
		"in the listed functions where it is currently turned off." )
		+ Attr( "use_native", xsct_rosetta_bool, "use the structure specified by -in:file:native as reference" )
		+ Attr( "use_fasta", xsct_rosetta_bool, "use a native FASTA sequence specified by the -in:file:fasta as reference" )
		+ Attr( "use_starting", xsct_rosetta_bool, "use the starting input structure (e.g. one passed to -s) as reference" )
		+ Attr( "use_current", xsct_rosetta_bool, "use the current structure (the one passed to apply) as the reference" )
		+ Attr( "pdbname", xsct_rosetta_bool, "use the structure specified by the filename as the reference" )
		+ Attr( "pssm", xs_string, "a filename of a blast formatted pssm file containing the sequence profile to use" )
		+ Attr( "chain", xsct_non_negative_integer, "0 is all chains, otherwise if a sequence is added, align it to the specified chain" )
		+ Attr::attribute_w_default( "scaling", "favor_seqprof_scaling_type",
		"Set how to scale the given values."
		"\"prob\"=Boltzmann-weighted probability based on the profile score"
		"\"global\"= global linear fixed-zero rescaling"
		"such that all (pre-weighted) values fall in the range of -1.0 to 1.0"
		"\"none\" does no adjustment of values.", "prob" )
		+ Attr::attribute_w_default( "matrix", xs_string, "Set substitution matrix; valid:  BLOSUM62,MATCH,IDENTITY", "BLOSUM62" )
		+ Attr( "exclude_resnums", xs_string, "Exclude residues from being contrained." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Sets residue type constraints (SequenceProfileConstraint) on the pose according to the given profile and weight.", attlist );
}

std::string FavorSequenceProfileCreator::keyname() const {
	return FavorSequenceProfile::mover_name();
}

protocols::moves::MoverOP
FavorSequenceProfileCreator::create_mover() const {
	return protocols::moves::MoverOP( new FavorSequenceProfile );
}

void FavorSequenceProfileCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FavorSequenceProfile::provide_xml_schema( xsd );
}



} //moves
} //protocols
