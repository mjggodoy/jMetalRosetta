// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/EnvMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/EnvMover.hh>
#include <protocols/environment/EnvMoverCreator.hh>

// Package headers
#include <protocols/environment/Environment.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/tag/Tag.hh>

// tracer
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.EnvMover", basic::t_info );

namespace protocols {
namespace environment {

// creator
// XRW TEMP std::string
// XRW TEMP EnvMoverCreator::keyname() const {
// XRW TEMP  return EnvMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP EnvMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new EnvMover() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EnvMover::mover_name() {
// XRW TEMP  return "Environment";
// XRW TEMP }

EnvMover::EnvMover():
	Mover(),
	movers_( new moves::SequenceMover )
{
	env_ = EnvironmentOP( new Environment( "env" ) );
	movers_->use_mover_status( true );
}

EnvMover::~EnvMover() = default;

void EnvMover::apply( Pose& pose ) {
	if ( movers_->size() == 0 ) {
		tr.Warning << "[WARNING] The environment " << env_->name()
			<< " is being run without any registrant movers." << std::endl;
	}

	EnvironmentOP env( new Environment( env_->name() ) );
	env->auto_cut( env_->auto_cut() );
	env->inherit_cuts( env_->inherit_cuts() );
	env->allow_pure_movers( env_->allow_pure_movers() );

	env->register_mover( movers_ );

	for ( auto const & reg_only_mover : reg_only_movers_ ) {
		env->register_mover( reg_only_mover );
	}

	core::pose::Pose ppose;
	try {
		ppose = env->start( pose );
	} catch( ... ) {
		tr.Error << "[ERROR] Error during broking in environment '" << get_name() << "'. " << std::endl;
		throw;
	}

	try {
		movers_->apply( ppose );
	} catch( EXCN_Env_Security_Exception& e ) {
		std::ostringstream ss;
		ss << "Error ocurred within " << this->get_name() << ". ";

		core::pose::PDBInfoCOP info( ppose.pdb_info() );

		if ( info &&
				e.id() != core::id::BOGUS_DOF_ID &&
				e.id().atom_id().rsd() <= ppose.size() ) {
			ss << "According to the PDBInfo object, the violating residue was "
				<< info->number( e.id().atom_id().rsd() )
				<< info->chain( e.id().atom_id().rsd() ) << ". ";
		}

		e.add_msg( ss.str() );

		throw;// e;
	}

	set_last_move_status( movers_->get_last_move_status() );
	tr.Debug << this->get_name() << " exited with status " << get_last_move_status() << std::endl;

	pose = env->end( ppose );
}

void EnvMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	core::pose::Pose const& pose ) {
	typedef utility::vector0< TagCOP > TagCOPs;

	env_ = EnvironmentOP( new Environment( tag->getOption<std::string>( "name" ) ) );

	env_->auto_cut( tag->getOption< bool >( "auto_cut", env_->auto_cut() ) );
	env_->inherit_cuts( tag->getOption< bool >( "inherit_cuts", env_->inherit_cuts() ) );
	env_->allow_pure_movers( tag->getOption< bool >( "allow_pure_movers", env_->allow_pure_movers() ) );

	TagCOPs const& subtags = tag->getTags();
	for ( auto const & subtag : subtags ) {
		parse_subtag( subtag, data, filters, movers, pose );
	}
}

moves::MoverOP find_mover( utility::tag::TagCOP tag,
	moves::Movers_map movers ){

	std::string mover_name("[NOT_SET]");
	if ( tag->hasOption( "name" ) ) {
		mover_name = tag->getOption< std::string >( "name" );
	} else if ( tag->hasOption( "mover" ) ) {
		mover_name = tag->getOption< std::string >( "mover" );
	} else {
		std::string err = "An environment was provided a subtag without the appropriate options. Either 'name' or 'mover' must be specified.";
		throw utility::excn::EXCN_RosettaScriptsOption( err );
	}
	moves::Movers_map::const_iterator mv_it = movers.find( mover_name );

	if ( mv_it != movers.end() ) {
		assert( mv_it->second );
		return mv_it->second;
	} else {
		std::string err = "The mover " + mover_name + " could not be found. Check your spelling in the xml script.";
		throw utility::excn::EXCN_RosettaScriptsOption( err );
	}

}

void EnvMover::parse_subtag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const & movers,
	core::pose::Pose const& ) {
	try {
		if ( tag->getName() == "Apply" ) {
			moves::MoverOP mover = find_mover( tag, movers );
			this->add_apply_mover( mover );
			if ( reg_only_movers_.find( utility::pointer::static_pointer_cast< Mover >( mover ) ) != reg_only_movers_.end() ) {
				tr.Warning << "[TIP] You don't need to register the mover "
					<< tag->getOption< std::string >( "name " )
					<< " with the 'Apply' tag ahead of time. The 'Apply' Environment tag does that automatically." << std::endl;
			}
		} else if ( tag->getName() == "Register" ) {
			moves::MoverOP mover = find_mover( tag, movers );
			this->add_registered_mover( mover );
		} else {
			std::ostringstream err;
			err << "The Environment cannot be used with the tag '" << *tag << "'.";
			throw utility::excn::EXCN_RosettaScriptsOption( err.str() );
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "In Environment '" + get_name() + "': " + e.msg() );
	}
}

void EnvMover::add_apply_mover( protocols::moves::MoverOP mover_in ) {
	if ( !mover_in ) {
		std::ostringstream ss;
		ss << "Mover '" << this->get_name() << "' recieved a null pointer in " << __FUNCTION__ << ".";
		throw utility::excn::EXCN_NullPointer( ss.str() );
	}
	movers_->add_mover( mover_in );
}

void EnvMover::add_registered_mover( protocols::moves::MoverOP mover_in ) {
	if ( !mover_in ) {
		std::ostringstream ss;
		ss << "Mover '" << this->get_name() << "' recieved a null pointer in " << __FUNCTION__ << ".";
		throw utility::excn::EXCN_NullPointer( ss.str() );
	}
	reg_only_movers_.insert( mover_in );
}


// XRW TEMP std::string EnvMover::get_name() const {
// XRW TEMP  return "EnvMover("+env_->name()+")";
// XRW TEMP }

moves::MoverOP EnvMover::clone() const {
	return moves::MoverOP( new EnvMover( *this ) );
}

std::string EnvMover::get_name() const {
	return mover_name();
}

std::string EnvMover::mover_name() {
	return "Environment";
}

void EnvMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"auto_cut", xsct_rosetta_bool,
		"XSD_XRW: TO DO");
	attlist + XMLSchemaAttribute(
		"inherit_cuts", xsct_rosetta_bool,
		"XSD_XRW: TO DO");
	attlist + XMLSchemaAttribute(
		"allow_pure_movers", xsct_rosetta_bool,
		"XSD_XRW: TO DO");

	XMLSchemaSimpleSubelementList ssl;
	AttributeList ssl_attlist;
	// Either name or mover must be provided
	ssl_attlist + XMLSchemaAttribute(
		"name", xs_string,
		"XSD_XRW: TO DO");
	ssl_attlist + XMLSchemaAttribute(
		"mover", xs_string,
		"XSD_XRW: TO DO");
	ssl.add_simple_subelement("Register", ssl_attlist, "XSD_XRW: TO DO");

	AttributeList ssl_attlist2;
	ssl_attlist2 + XMLSchemaAttribute(
		"name", xs_string, "XSD XRW TO DO");
	ssl_attlist2 + XMLSchemaAttribute(
		"mover", xs_string,
		"XSD_XRW: TO DO");
	ssl.add_simple_subelement("Apply", ssl_attlist2, "XSD_XRW: TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"XSD_RW: TO DO",
		attlist, ssl);
}

std::string EnvMoverCreator::keyname() const {
	return EnvMover::mover_name();
}

protocols::moves::MoverOP
EnvMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new EnvMover );
}

void EnvMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnvMover::provide_xml_schema( xsd );
}


} // environment
} // protocols
