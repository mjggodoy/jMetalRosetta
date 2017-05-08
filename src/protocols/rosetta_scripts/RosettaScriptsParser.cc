// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/RosettaScriptsParser.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Interface base class Parser
/// @author Sarel Fleishman sarelf@u.washington.edu
/// @author Luki Goldschmidt lugo@uw.edu

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

// Package headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/parser/DataLoader.hh>
#include <protocols/jd2/parser/DataLoaderFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Project headers
#include <basic/options/option.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// movers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/moves/mover_schemas.hh>

// Pose Metric Calculators for filters
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

// option key includes
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>

namespace protocols {
namespace rosetta_scripts {

using namespace core;
using namespace basic::options;
using namespace scoring;
using namespace moves;
using namespace jd2;

static THREAD_LOCAL basic::Tracer TR( "protocols.rosetta_scripts.RosettaScriptsParser" );

RosettaScriptsParser::RosettaScriptsParser() :
	Parser(),
	validator_( new utility::tag::XMLValidator ),
	filenames_already_validated_()  //Empty list initially
{
	register_factory_prototypes();
}

RosettaScriptsParser::~RosettaScriptsParser()= default;

using utility::tag::TagOP;
using utility::tag::TagCOP;
using utility::tag::TagCAP;

typedef utility::vector0< TagCOP > TagCOPs;

/// @brief Actually read in the XML file.  Called recursively to read in XML files that
/// this XML file includes.  At the end of this operation, fin contains the contents
/// of the XML file, with all xi:includes replaced with the contents of included XML
/// files.  Files are opened and closed here.
///
/// @details Note that filenames_encountered is passed by copy rather than by reference
/// DELIBERATELY.  This is so that it remains a list of parent files, so that only
/// circular dependencies (attempts to include one's own parent, grandparent, etc.) are
/// detected.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
RosettaScriptsParser::read_in_and_recursively_replace_includes(
	std::string const &filename,
	std::stringstream &fin,
	utility::vector1 < std::string > filenames_encountered
) const {
	using namespace utility::io;

	//Check whether we've already encountered this filename
	for ( core::Size i=1, imax=filenames_encountered.size(); i<=imax; ++i ) {
		if ( filenames_encountered[i] == filename ) {
			throw utility::excn::EXCN_BadInput( "Error in protocols::rosetta_scipts::RosettaScriptsParser::read_in_and_recursively_replace_includes(): Circular inclusion pattern detected when reading \"" + filename + "\"." );
		}
	}
	if ( filenames_encountered.size() != 0 && TR.visible() ) TR << "Including RosettaScripts XML file " << filename << "." << std::endl;
	filenames_encountered.push_back(filename);

	//Actually read in the file:
	izstream inputstream;
	inputstream.open(filename.c_str() );
	if ( !inputstream.good() ) {
		utility_exit_with_message("Unable to open RosettaScripts XML file: \"" + filename + "\".");
	}
	std::string line;
	while ( getline( inputstream, line ) ) {
		fin << line << '\n';
	}
	inputstream.close();

	//Find xi:include tags:
	//Will be an innermost set of <> brackets.
	std::string finstr;
	core::Size bracketcount(0);
	core::Size innerbracketstart(0);
	core::Size last_chunk_end(0);
	for ( core::Size i=0, imax=fin.str().length(); i<imax; ++i ) { //Loop through the XML
		if ( fin.str()[i] == '<' ) {
			++bracketcount;
			innerbracketstart=i;
		} else if ( fin.str()[i] == '>' &&  bracketcount > 0 ) {
			--bracketcount;
			if ( has_slash_at_end (fin.str(), i) && !has_double_percent_signs(fin.str(), i) ) {
				TagCOP tag( utility::tag::Tag::create( fin.str().substr( innerbracketstart, i-innerbracketstart+1 ) ) ); //Parse the current inner tag.
				if ( tag->getName() == "xi:include" ) {
					finstr+=fin.str().substr(last_chunk_end, innerbracketstart-1-last_chunk_end)+'\n'; //Copy all of the text up to the start of the bracket.
					last_chunk_end = i+1;
					runtime_assert_string_msg( tag->hasOption("href"), "Error in protocols::rosetta_scipts::RosettaScriptsParser::read_in_and_recursively_replace_includes(): An \"xi:include\" tag must include an \"href=...\" statement." );
					runtime_assert_string_msg( !tag->hasOption("parse") || tag->getOption<std::string>("parse") == "XML", "Error in protocols::rosetta_scipts::RosettaScriptsParser::read_in_and_recursively_replace_includes(): An \"xi:include\" tag can ONLY be used to include XML for parsing.  Other parse types are unsupporrted in RosettaScripts." );
					std::string tempfilename( tag->getOption<std::string>("href") );
					std::stringstream tempstream;
					read_in_and_recursively_replace_includes( tempfilename, tempstream, filenames_encountered );  //Recursively call this function to read in included XML.
					finstr+=tempstream.str();
				}
			}
		}
	}
	if ( last_chunk_end < fin.str().length() ) {
		finstr+=fin.str().substr(last_chunk_end, fin.str().length()-last_chunk_end); //Copy all of the rest of the text.
	}
	fin.str( finstr );

	TR.flush();

}

/// @details Uses the Tag interface to the xml reader library in boost to parse
/// an xml file that contains design protocol information. A sample protocol
/// file can be found in src/pilot/apps/sarel/dock_design.protocol.
///
/// SCOREFXNS provides a way to define scorefunctions as they are defined in the
/// rosetta database, using the weights/patch convenctions. Several default
/// scorefunctions are preset and can be used without defining them explicitly.
///
/// FILTERS defines a set of filters that can be used together with the
/// dockdesign movers to prune out poses that don't meet certain criteria
///
/// MOVERS defines the movers that will be used
///
/// PROTOCOLS is the only order-sensitive section where subsequent movers and
/// filters are expected to be defined. These movers and filters were defined
/// in the previous two sections. The order in which the protocols are specified
/// by the user will be maintained by the DockDesign mover.
///
/// APPLY_TO_POSE This section allows for certain movers to be applied to the
/// pose prior to starting the DockDesign protocol. For instance, these movers
/// could set constraints, such as favor_native_residue. In this case, for
/// example, the weights of res_type_constraint in all of the scorefunctions
/// that are defined in SCOREFXNS or by default are set to 1.0, but can be
/// changed by the user globally (in the definition of the weight on the
/// constraint), or in particular for each of the scorefunctions by changing
/// the relevant term (which is set by default to the global value).
///
/// OUTPUT is a section which allows the XML control of how things are output.
/// Notice that the order of the sections by which the protocol is written
/// doesn't matter, BUT the order of the mover-filter pairs in PROTOCOLS section
/// does matter.
///
/// xi:include statements can be placed anywhere to trigger read-in of another
/// RosettaScripts XML file.  Contents of the file replace the xi:include
/// statement.
bool
RosettaScriptsParser::generate_mover_from_pose(
	JobCOP,
	Pose & pose,
	MoverOP & in_mover,
	bool new_input,
	std::string const & xml_fname,
	bool guarantee_new_mover
)
{

	bool modified_pose( false );

	if ( !new_input && !guarantee_new_mover ) return modified_pose;

	std::string const dock_design_filename( xml_fname == "" ? option[ OptionKeys::parser::protocol ] : xml_fname );
	TR << "dock_design_filename=" << dock_design_filename << std::endl;

	// The input string.  Input files are opened and closed by read_in_and_recursively_replace_includes().
	// After calling this function, fin will contain the XML with includes replaced by whatever XML they include.
	std::stringstream fin;
	// A list of all of the files that have been included, which is checked whenever a
	// new file is included to prevent circular dependencies.
	utility::vector1 < std::string > filenames_encountered;
	read_in_and_recursively_replace_includes( dock_design_filename, fin, filenames_encountered ); //Actually read the input XML, reading in all included XML files as well.  No interpretation is done yet (except for recognizing xi:include).

	//fin.open(dock_design_filename.c_str() );
	//if ( !fin.good() ) {
	// utility_exit_with_message("Unable to open Rosetta Scripts XML file: '" + dock_design_filename + "'.");
	//}

	if ( option[ OptionKeys::parser::script_vars ].user() ) {
		std::stringstream fin_sub;
		substitute_variables_in_stream(fin, option[ OptionKeys::parser::script_vars ], fin_sub);
		fin.clear();
		fin.str( fin_sub.str() );
	}

	// Validate the input script against the XSD for RosettaScripts, if it hasn't yet been validated.
	if ( !was_already_validated( dock_design_filename ) ) {
		validate_input_script_against_xsd( dock_design_filename, fin );
		set_validated( dock_design_filename );
	}

	TagCOP tag = utility::tag::Tag::create(fin);

	//fin.close();
	TR << "Parsed script:" << "\n";
	TR << tag;
	TR.flush();
	runtime_assert( tag->getName() == "dock_design" || tag->getName() == "ROSETTASCRIPTS" );

	if ( tag->getName() == "dock_design" ) {
		TR << "<dock_design> tag is deprecated; please use <ROSETTASCRIPTS> instead." << std::endl;
	}

	in_mover = generate_mover_for_protocol( pose, modified_pose, tag );

	return modified_pose;
}

MoverOP RosettaScriptsParser::parse_protocol_tag(Pose & pose, utility::tag::TagCOP protocol_tag)
{
	bool modified_pose = false;

	MoverOP mover =  generate_mover_for_protocol(pose, modified_pose, protocol_tag);

	if ( modified_pose ) {
		throw utility::excn::EXCN_RosettaScriptsOption("RosettaScriptsParser::parse_protocol_tag resulted in modified_pose");
	}

	return mover;
}

MoverOP RosettaScriptsParser::parse_protocol_tag(TagCOP protocol_tag)
{
	Pose temp_pose;
	bool modified_pose = false;

	MoverOP mover =  generate_mover_for_protocol(temp_pose, modified_pose, protocol_tag);

	if ( modified_pose ) {
		throw utility::excn::EXCN_RosettaScriptsOption("RosettaScriptsParser::parse_protocol_tag resulted in modified_pose");
	}

	return mover;
}

MoverOP RosettaScriptsParser::generate_mover_for_protocol(
	Pose & pose,
	bool & modified_pose,
	TagCOP tag
)
{
	protocols::rosetta_scripts::ParsedProtocolOP protocol( new protocols::rosetta_scripts::ParsedProtocol );

	Movers_map movers;
	protocols::filters::Filters_map filters;
	basic::datacache::DataMap data; // abstract objects, such as scorefunctions, to be used by filter and movers
	std::set< ImportTagName > import_tag_names;

	MoverOP mover;

	typedef std::pair< std::string const, MoverOP > StringMover_pair;
	typedef std::pair< std::string const, protocols::filters::FilterOP > StringFilter_pair;

	//setting up some defaults
	protocols::filters::FilterOP true_filter( new protocols::filters::TrueFilter );
	protocols::filters::FilterOP false_filter( new protocols::filters::FalseFilter );
	filters.insert( StringFilter_pair( "true_filter", true_filter ) );
	filters.insert( StringFilter_pair( "false_filter", false_filter ) );

	MoverOP null_mover( new protocols::moves::NullMover() );
	movers.insert( StringMover_pair( "null", null_mover) );

	// default scorefxns
	ScoreFunctionOP commandline_sfxn  = core::scoring::get_score_function();
	ScoreFunctionOP talaris2014       = ScoreFunctionFactory::create_score_function(TALARIS_2014);
	ScoreFunctionOP talaris2013       = ScoreFunctionFactory::create_score_function(TALARIS_2013);
	ScoreFunctionOP score12           = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH );
	ScoreFunctionOP docking_score     = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, DOCK_PATCH );
	ScoreFunctionOP soft_rep          = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS );
	ScoreFunctionOP docking_score_low = ScoreFunctionFactory::create_score_function( "interchain_cen" );
	ScoreFunctionOP score4L           = ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );

	data.add( "scorefxns", "commandline", commandline_sfxn );
	data.add( "scorefxns", "talaris2014", talaris2014 );
	data.add( "scorefxns", "talaris2013", talaris2013 );
	data.add( "scorefxns", "score12", score12 );
	data.add( "scorefxns", "score_docking", docking_score );
	data.add( "scorefxns", "soft_rep", soft_rep );
	data.add( "scorefxns", "score_docking_low", docking_score_low );
	data.add( "scorefxns", "score4L", score4L );
	//default scorefxns end

	// Data Loaders -- each data loader handles one of the top-level blocks of a
	// rosetta script.  These blocks are handled by the RosettaScriptsParser itself;
	// other data loaders may be registered with the DataLoaderFactory at load time
	// so that arbitrary data, living in any library, can be loaded into the DataMap
	// through the XML interface.
	std::set< std::string > non_data_loader_tags;
	non_data_loader_tags.insert( "MOVERS" );
	non_data_loader_tags.insert( "APPLY_TO_POSE" );
	non_data_loader_tags.insert( "FILTERS" );
	non_data_loader_tags.insert( "PROTOCOLS" );
	non_data_loader_tags.insert( "OUTPUT" );
	non_data_loader_tags.insert( "IMPORT" );

	// Load in data into the basic::datacache::DataMap object.  All tags beside those listed
	// in the non_data_loader_tags set are considered DataLoader tags.
	TagCOPs const all_tags = tag->getTags();
	for ( Size ii = 0; ii < all_tags.size(); ++ii ) {
		using namespace parser;
		TagCOP iitag = all_tags[ ii ];
		if ( non_data_loader_tags.find( iitag->getName() ) != non_data_loader_tags.end() ) continue;
		DataLoaderOP loader = DataLoaderFactory::get_instance()->newDataLoader( iitag->getName() );
		loader->load_data( pose, iitag, data );
	}

	if ( !tag->hasTag("PROTOCOLS") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("parser::protocol file must specify PROTOCOLS section");
	}

	// Round 1: Import previously defined filters and movers
	for ( TagCOP const curr_tag : all_tags ) {

		////// IMPORT
		if ( curr_tag->getName() != "IMPORT" ) {
			continue;
		}

		if ( curr_tag->hasOption("taskoperations") ) {
			// Import task operations
			std::string taskoperations_str( curr_tag->getOption<std::string>("taskoperations") );
			std::istringstream taskoperation_ss(taskoperations_str);
			std::string taskoperation_name;
			while ( std::getline(taskoperation_ss, taskoperation_name, ',') ) {
				import_tag_names.insert(std::make_pair("TASKOPERATIONS", taskoperation_name));
			}
		}//fi taskoperations

		if ( curr_tag->hasOption("filters") ) {
			// Import filters
			std::string filters_str( curr_tag->getOption<std::string>("filters") );
			std::istringstream filters_ss(filters_str);
			std::string filter_name;
			while ( std::getline(filters_ss, filter_name, ',') ) {
				import_tag_names.insert(std::make_pair("FILTERS", filter_name));
			}
		}//fi filters

		if ( curr_tag->hasOption("movers") ) {
			// Import movers
			std::string movers_str( curr_tag->getOption<std::string>("movers") );
			std::istringstream movers_ss(movers_str);
			std::string mover_name;
			while ( std::getline(movers_ss, mover_name, ',') ) {
				import_tag_names.insert(std::make_pair("MOVERS", mover_name));
			}
		}//fi movers

	}// for curr_tag

	// Attempt to find and import requested objects; throws exception on failure
	if ( import_tag_names.size() > 0 ) {
		import_tags(import_tag_names, tag, data, filters, movers, pose);
	}

	// Round 2: Process definitions in this ROSETTASCRIPTS block
	for ( TagCOP const curr_tag : all_tags ) {

		///// APPLY TO POSE
		if ( curr_tag->getName() == "APPLY_TO_POSE" ) {
			// section is not mandatory
			// apply to pose may affect all of the scorefxn definitions below, so it is called first.
			TagCOPs const apply_tags( curr_tag->getTags() );

			for ( TagCOP apply_tag_ptr : apply_tags ) {
				std::string const mover_type( apply_tag_ptr->getName() );
				MoverOP new_mover( MoverFactory::get_instance()->newMover( apply_tag_ptr, data, filters, movers, pose ) );
				runtime_assert( new_mover != nullptr );
				new_mover->apply( pose );
				TR << "Defined and applied mover of type " << mover_type << std::endl;
				bool const name_exists( movers.find( mover_type ) != movers.end() );
				if ( name_exists ) {
					throw utility::excn::EXCN_RosettaScriptsOption("Can't apply_to_pose the same mover twice for " + mover_type);
				}

				modified_pose = true;
			} // done apply_tag_ptr

		}//  fi apply_to_pose
		TR.flush();

		////// FILTERS
		if ( curr_tag->getName() == "FILTERS" ) {
			for ( TagCOP const tag_ptr : curr_tag->getTags() ) {
				instantiate_filter(tag_ptr, data, filters, movers, pose);
			}
		}
		TR.flush();

		////// MOVERS
		if ( curr_tag->getName() == "MOVERS" ) {
			for ( TagCOP const tag_ptr : curr_tag->getTags() ) {
				instantiate_mover(tag_ptr, data, filters, movers, pose);
			}
		}
		TR.flush();

	}// for curr_tag

	////// ADD MOVER FILTER PAIRS
	TagCOP const protocols_tag( tag->getTag("PROTOCOLS") );
	protocol->parse_my_tag( protocols_tag, data, filters, movers, pose );
	TR.flush();

	////// Set Output options
	if ( tag->hasTag("OUTPUT") ) {
		TagCOP const output_tag( tag->getTag("OUTPUT") );
		protocol->final_scorefxn( rosetta_scripts::parse_score_function( output_tag, data, "commandline" ) );
	} else {
		protocol->final_scorefxn( commandline_sfxn );
	}

	tag->die_for_unaccessed_options_recursively();

	return protocol;
}

/// @brief Instantiate a new filter and add it to the list of defined filters for this ROSETTASCRIPTS block
void RosettaScriptsParser::instantiate_filter(
	TagCOP const & tag_ptr,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	Movers_map & movers,
	core::pose::Pose & pose
) {
	std::string const type( tag_ptr->getName() );
	if ( ! tag_ptr->hasOption("name") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't define unnamed Filter of type " + type);
	}

	std::string const user_defined_name( tag_ptr->getOption<std::string>("name") );
	bool const name_exists( filters.find( user_defined_name ) != filters.end() );
	if ( name_exists ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Filter of name \"" + user_defined_name + "\" (of type " + type + ") already exists.");
	}

	protocols::filters::FilterOP new_filter( protocols::filters::FilterFactory::get_instance()->newFilter( tag_ptr, data, filters, movers, pose ) );
	runtime_assert( new_filter != nullptr );
	filters.insert( std::make_pair( user_defined_name, new_filter ) );
	TR << "Defined filter named \"" << user_defined_name << "\" of type " << type << std::endl;
}

/// @brief Instantiate a new mover and add it to the list of defined movers for this ROSETTASCRIPTS block
void RosettaScriptsParser::instantiate_mover(
	TagCOP const & tag_ptr,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	Movers_map & movers,
	core::pose::Pose & pose
) {
	std::string const type( tag_ptr->getName() );
	if ( ! tag_ptr->hasOption("name") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't define unnamed Mover of type " + type);
	}

	std::string const user_defined_name( tag_ptr->getOption<std::string>("name") );
	bool const name_exists( movers.find( user_defined_name ) != movers.end() );
	if ( name_exists ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Mover of name \"" + user_defined_name + "\" (of type " + type + ") already exists.");
	}

	MoverOP new_mover( MoverFactory::get_instance()->newMover( tag_ptr, data, filters, movers, pose ) );
	runtime_assert( new_mover != nullptr );
	movers.insert( std::make_pair( user_defined_name, new_mover ) );
	TR << "Defined mover named \"" << user_defined_name << "\" of type " << type << std::endl;
}


/// @brief Instantiate a new task operation (used in IMPORT tag)
void RosettaScriptsParser::instantiate_taskoperation(
	TagCOP const & tag_ptr,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & /*filters*/,
	Movers_map & /*movers*/,
	core::pose::Pose & /*pose*/
) {
	using namespace core::pack::task::operation;

	std::string const type( tag_ptr->getName() );
	if ( ! tag_ptr->hasOption("name") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't define unnamed TaskOperation of type " + type);
	}

	std::string const user_defined_name( tag_ptr->getOption<std::string>("name") );
	if ( data.has( "task_operations", user_defined_name ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation with name \"" + user_defined_name + "\" (of type " + type + ") already exists.");
	}

	TaskOperationOP new_t_o( TaskOperationFactory::get_instance()->newTaskOperation( type, data, tag_ptr ) );
	runtime_assert( new_t_o != nullptr );
	data.add("task_operations", user_defined_name, new_t_o );
	TR << "Defined TaskOperation named \"" << user_defined_name << "\" of type " << type << std::endl;
}

/// @brief Recursively find a specific tag by option name and valie in any ROSETTASCRIPTS tag that's a child of rootTag
TagCOP RosettaScriptsParser::find_rosettascript_tag(
	TagCOP tag_in,
	const std::string & section_name,
	const std::string & option_name,
	const std::string & option_value
) {

	TagCAP tag_ap(tag_in);

	// Find the next higher ROSETTASCRIPTS block
	do {
		TagCOP tag = tag_ap.lock();
		tag_ap = tag ? tag->getParent() : TagCAP();
	} while(!tag_ap.expired() && tag_ap.lock()->getName() != "ROSETTASCRIPTS");

	TagCOP tag( tag_ap.lock() );
	if ( !tag ) {
		return nullptr;
	}

	// Look for section_name (MOVERS, FITLERS, ...) directly below ROSETTASCRIPTS
	if ( tag->hasTag(section_name) ) {
		TagCOP section_tag( tag->getTag(section_name) );
		if ( !option_name.length() ) {
			return section_tag;
		}
		for ( TagCOP const child_tag : section_tag->getTags() ) {
			if ( child_tag->getOption<std::string>(option_name) == option_value ) {
				return child_tag;
			}
		}
	}

	// No match, walk up another level
	TagCOP tagParent( tag->getParent().lock() );
	if ( tagParent ) {
		return find_rosettascript_tag(tagParent, section_name, option_name, option_value);
	}

	return nullptr;
}

/// @brief Import filters, movers, ... specified in the IMPORT tag
/// in the order they were originally defined elsewhere in the script
void RosettaScriptsParser::import_tags(
	std::set< ImportTagName > & import_tag_names,
	utility::tag::TagCOP & my_tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map & movers,
	core::pose::Pose & pose
) {
	// Process all parent ROSETTASCRIPTS tags, one at a time
	TagCAP curr_level_tag_ap(my_tag);

	while ( !curr_level_tag_ap.expired() ) {

		// Find direct parent ROSETTASCRIPTS tag
		do {
			TagCOP curr_level_tag = curr_level_tag_ap.lock();
			curr_level_tag_ap = curr_level_tag ? curr_level_tag->getParent() : TagCAP();
		} while(!curr_level_tag_ap.expired() && curr_level_tag_ap.lock()->getName() != "ROSETTASCRIPTS");

		TagCOP curr_level_tag( curr_level_tag_ap.lock() );
		if ( !curr_level_tag ) {
			break; // No ROSETTASCRIPTS tag to import from
		}

		TR.Debug << "Importing from tag: " << curr_level_tag << std::endl;

		// Check what we'd like to import from it
		for ( TagCOP tag : curr_level_tag->getTags() ) {

			if ( tag->getName() == "TASKOPERATIONS" ) {
				for ( TagCOP taskoperation_tag : tag->getTags() ) {
					std::string taskoperation_name( taskoperation_tag->getOption<std::string>("name") );
					ImportTagName key( std::make_pair( tag->getName(), taskoperation_name ) );
					bool const need_import( import_tag_names.find( key ) != import_tag_names.end() );
					if ( need_import ) {
						instantiate_taskoperation(taskoperation_tag, data, filters, movers, pose);
						import_tag_names.erase(key);
					}
				}
			}
			if ( tag->getName() == "FILTERS" ) {
				for ( TagCOP filter_tag : tag->getTags() ) {
					std::string filter_name( filter_tag->getOption<std::string>("name") );
					ImportTagName key( std::make_pair( tag->getName(), filter_name ) );
					bool const need_import( import_tag_names.find( key ) != import_tag_names.end() );
					if ( need_import ) {
						instantiate_filter(filter_tag, data, filters, movers, pose);
						import_tag_names.erase(key);
					}
				}
			}

			if ( tag->getName() == "MOVERS" ) {
				for ( TagCOP mover_tag : tag->getTags() ) {
					std::string mover_name( mover_tag->getOption<std::string>("name") );
					ImportTagName key( std::make_pair( tag->getName(), mover_name ) );
					bool const need_import( import_tag_names.find( key ) != import_tag_names.end() );
					if ( need_import ) {
						TagCAP parent = my_tag->getParent();
						if ( utility::pointer::equal(parent, mover_tag) ) {
							// Current mover_tag is our parent, i.e. same ROSETTASCRIPTS tag
							throw utility::excn::EXCN_RosettaScriptsOption("Cannot import mover " + mover_name + " into itself; recursion detected");
						}
						instantiate_mover(mover_tag, data, filters, movers, pose);
						import_tag_names.erase(key);
					}
				}
			}

		} //for curr_level_tag->getTags()

		if ( import_tag_names.size() <= 0 ) {
			break; // All done
		}

		// Go up the tree
		curr_level_tag_ap = curr_level_tag->getParent();

	} // while

	// Check if there are any remaining imports that could not be satisfied
	for ( ImportTagName import_tag : import_tag_names ) {
		std::string msg("Failed to import " + import_tag.second + " from " + import_tag.first);
		TR << msg << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption(msg);
	}
}


/// @brief Create a variable substituted version of the input stream, given a StringVectorOption formated list of variables
///to substitiute. Each item in script_vars should be in the form of "variable=value", where "value" is the string to substitiute
///into the input stream whereever the string "%%variable%%" is found in the input.
///
///Having the return value be passed through a parameter is to get around some copy-constructor limitations of C++ streams.
void
RosettaScriptsParser::substitute_variables_in_stream( std::istream & instream, utility::options::StringVectorOption const& script_vars, std::stringstream & out){
	using namespace std;
	//Parse variable substitutions
	map<string,string> var_map;
	for ( Size ii = 1; ii <= script_vars.size(); ii++ ) {
		Size eq_pos(script_vars[ii].find('='));
		if ( eq_pos != string::npos ) { //ignore ones without a '='
			var_map[ script_vars[ii].substr(0,eq_pos) ] = script_vars[ii].substr(eq_pos+1);
		}
	}
	//Print parsing for confirmation
	TR << "Variable substitution will occur with the following values: ";
	for ( auto const & elem : var_map ) {
		TR << "'%%" << elem.first << "%%'='" << elem.second << "';  ";
	}
	TR << std::endl;
	var_map[""] = "%%"; //for %%%%->%% substitutions
	//Perform actual variable substitution
	TR << "Substituted script:" << "\n";
	string line;
	while ( getline( instream, line ) ) {
		Size pos, end, last(0);
		string outline; // empty string
		while ( ( pos = line.find("%%", last) ) != string::npos ) {
			end = line.find("%%", pos+2);
			if ( end == string::npos ) break; //Unmatched %% on line - keep as is.
			outline += line.substr( last, pos-last );
			last = end+2;
			string var( line.substr( pos+2, end-(pos+2) ) );
			if ( var_map.find( var ) != var_map.end() ) {
				outline += var_map[var];
			} else {
				outline += "%%" + var + "%%"; // Silently ignore missing values - will probably cause error later.
			}
		}
		outline += line.substr( last ); // Add the rest of the line, starting at last
		TR << outline << "\n";
		out << outline << "\n";
	}
	TR.flush();
}

/// @brief Factories avoid requiring compiler dependencies for all possible constructible derived classes,
///by allowing runtime registration of derived class prototypes. However, this requires
///pre-registration of a derived type with the factory prior to asking the factory to return an
///instance of that type. This method registers those additional derived classes that are available for
///construction in the RosettaScriptsParser context.
/// TO-DO: replace this manual factory registration system with a load-time-initialized singleton scheme (see r32404 for example)
void
RosettaScriptsParser::register_factory_prototypes()
{
	// note: TaskOperations are now registered with a singleton factory at load time using apl's creator/registrator scheme

	// also register some constraint types with the ConstraintFactory (global singleton class)
	// this allows derived non-core constraints to be constructed from string definitions in constraints files
	//using namespace core::scoring::constraints;
	//ConstraintFactory & cstf( ConstraintIO::get_cst_factory() );
	//cstf.add_type( new core::scoring::constraints::SequenceProfileConstraint(
	// Size(), utility::vector1< id::AtomID >(), NULL ) );

	// register calculators
	core::Size const chain1( 1 ), chain2( 2 );
	using namespace core::pose::metrics;

	if ( !CalculatorFactory::Instance().check_calculator_exists( "sasa_interface" ) ) {
		PoseMetricCalculatorOP int_sasa_calculator( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator( chain1, chain2 ) );
		CalculatorFactory::Instance().register_calculator( "sasa_interface", int_sasa_calculator );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
		CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( "ligneigh" ) ) {
		PoseMetricCalculatorOP lig_neighbor_calc( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( chain1, chain2 ) );
		CalculatorFactory::Instance().register_calculator( "ligneigh", lig_neighbor_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( "liginterfE" ) ) {
		PoseMetricCalculatorOP lig_interf_E_calc( new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( "ligneigh" ) );
		CalculatorFactory::Instance().register_calculator( "liginterfE", lig_interf_E_calc );
	}
}

/// @brief Is the current tag one with a slash at the end?
/// @details Starting from endpos, crawl backward until we hit a
/// '/' or a non-whitespace character.  Return true if and only
/// if it's a slash, false otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
RosettaScriptsParser::has_slash_at_end(
	std::string const &str,
	core::Size const endpos
) const {
	debug_assert( str.length() >=endpos );
	debug_assert( endpos > 0 );
	core::Size i(endpos);

	do {
		--i;
		if ( str[i] == '/' ) return true;
		else if ( str[i] != ' ' && str[i] != '\t' && str[i] != '\n' ) return false;
	} while(i!=0);

	return false;
}

/// @brief Is the current tag one with something that will be replaced by variable replacement?
/// @details Starting from endpos, crawl backward until we hit a '<' or '%'.  Return 'false' if
/// the former and 'true' if the latter.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
RosettaScriptsParser::has_double_percent_signs(
	std::string const &str,
	core::Size const endpos
) const {
	debug_assert( str.length() >=endpos );
	debug_assert( endpos > 0 );
	core::Size i(endpos);

	do {
		--i;
		if ( str[i] == '<' ) return false;
		else if ( str[i] == '%' && i>0 && str[i-1] == '%' ) return true;
	} while(i!=0);

	return false;

}

/// @brief Is a filename in the list of files already validated?
/// @details Returns ture if it is, false if it is not.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
RosettaScriptsParser::was_already_validated(
	std::string const &filename
) const {
	for ( core::Size i=1, imax=filenames_already_validated_.size(); i<=imax; ++i ) {
		if ( !filenames_already_validated_[i].compare( filename ) ) return true;
	}
	return false;
}

/// @brief Add a filename to the list already validated.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
RosettaScriptsParser::set_validated(
	std::string const &filename
) {
	filenames_already_validated_.push_back(filename);
}

void
RosettaScriptsParser::validate_input_script_against_xsd( std::string const & fname, std::stringstream const & input_xml ) const
{
	using namespace utility::tag;

	// Step 1: Internal generation of the schema and initialization of the libxml2 data structure for schema validation
	if ( ! validator_->schema_has_been_set() ) {
		std::ostringstream oss;
		oss << "If you are seeing this message, the internally-generated XML Schema for rosetta_scripts could not be properly generated\nThis failure occurred before your XML script that provided was examined. The error has been compiled into Rosetta and will need to be fixed by a developer.";
		std::string schema;
		try {
			schema = xsd_for_rosetta_scripts();
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			oss << "An error was encountered while the string of the schema was being generated; this is before the schema is analyzed for whether it is correct or not.\n";
			oss << e.msg() << "\n";
			throw utility::excn::EXCN_Msg_Exception( oss.str() );
		}

		if ( basic::options::option[ basic::options::OptionKeys::parser::output_schema ].user() ) {
			std::ofstream ofs( basic::options::option[ basic::options::OptionKeys::parser::output_schema ] );
			ofs << schema << std::endl;
			exit( 0 );
		}

		XMLValidationOutput schema_valid_output;
		try {
			TR << "Initializing schema validator..." << std::endl;
			schema_valid_output = validator_->set_schema( schema );
			TR << "...done" << std::endl;
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			oss << e.msg() << "\n";
			throw utility::excn::EXCN_Msg_Exception( oss.str() );
		}
		if ( ! schema_valid_output.valid() ) {
			oss << "If there is an error message immediately above stating that the schema failed to validate, and errors below that look like 'real XML' with lots of <xs:something> tags and NOT like your XML input, then you have a global schema validation error not an XML input validation error.  Read the error message below and fix your schema in the C++ code. \n\n";
			oss << "Errors: " << schema_valid_output.error_messages() << "\n";
			oss << "Warnings: " << schema_valid_output.warning_messages() << "\n";
			throw utility::excn::EXCN_Msg_Exception( oss.str() );
		}
	}

	TR << "Validating input script..." << std::endl;

	std::string xml_file = input_xml.str();

	std::ostringstream oss;
	oss << "Input rosetta scripts XML file \"" << fname << "\" failed to validate against"
		" the rosetta scripts schema. Use the option -parser::output_schema <output filename> to output"
		" the schema to a file to see all valid options.\n";

	std::string const long_error_string("Your XML has failed validation.  The error message below will tell you where in your XML file the error occurred.  Here's how to fix it:\n\n1) If the validation fails on something obvious, like an illegal attribute due to a spelling error (perhaps you used scorefnction instead of scorefunction), then you need to fix your XML file.\n2) If you haven’t run the XML rewriter script and this might be pre-2017 Rosetta XML, run the rewriter script (tools/xsd_xrw/rewrite_rosetta_script.py) on your input XML first.  The attribute values not being in quotes (scorefunction=talaris2014 instead of scorefunction=\"talaris2014\") is a good indicator that this is your problem.\n3) If you are a developer and neither 1 nor 2 worked - email the developer’s mailing list or try Slack.\n4) If you are an academic or commercial user - try the Rosetta Forums https://www.rosettacommons.org/forum\n5)\n");

	oss << long_error_string << "\n\n";

	XMLValidationOutput validator_output;
	try {
		validator_output = validator_->validate_xml_against_schema( xml_file );
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		oss << e.msg() << "\n";
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
	}

	if ( ! validator_output.valid() ) {
		//I don't think this part actually works - these messages are always empty and there is a matching "does this work?" comment in XMLSchemaValidation.cc XRW TODO
		oss << "Error messages were:\n" << validator_output.error_messages();;
		oss << "------------------------------------------------------------\n";
		oss << "Warning messages were:\n" << validator_output.warning_messages();
		oss << "------------------------------------------------------------\n";
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
	}

	if ( basic::options::option[ basic::options::OptionKeys::parser::validate_and_exit ] ) {
		TR << "Validation succeeded. Exiting" << std::endl;
		exit( 0 );
	}

	TR << "...done" << std::endl;
}

std::string
RosettaScriptsParser::xsd_for_rosetta_scripts()
{
	TR << "Generating XML Schema for rosetta_scripts..." << std::endl;
	using namespace utility::tag;
	XMLSchemaDefinition xsd;

	write_ROSETTASCRIPTS_complex_type( xsd );

	// Finally, the root of the tree.
	// TO DO: separate the definition of the ROSETTA_SCRIPTS complex type and everything
	// that goes inside of it into its own function so that the ROSETTA_SCRIPTS element
	// becomes optional.
	XMLSchemaElement rosetta_scripts_element;
	rosetta_scripts_element.name( rosetta_scripts_element_name() )
		.type_name( rosetta_scripts_complex_type_naming_func( rosetta_scripts_element_name() ));
	xsd.add_top_level_element( rosetta_scripts_element );

	TR << "...done" << std::endl;

	return xsd.full_definition();

}

void
RosettaScriptsParser::write_ROSETTASCRIPTS_complex_type( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// early exit if we've been here already
	if ( xsd.has_top_level_element( rosetta_scripts_complex_type_naming_func( rosetta_scripts_element_name() ) ) ) return;

	// We have lots of complexTypes to define; we need them for the following elements:
	//
	// MOVERS
	// FILTERS
	// All the data loaders, (TASKOPERATIONS, SCOREFXNS, etc.)
	// IMPORT
	// PROTOCOLS -- PROTOCOLS subelement is just the complexType for the ParsedProtocol mover!
	// OUTPUT
	// APPLY_TO_POSE ??
	// ROSETTASCRIPTS
	//
	// Then we need complexTypes defined from all of the DataLoaders -- we'll ask the
	// DataLoaderFactory for these.

	// Do the ROSETTASCRIPTS complex type first, so that if we recurse here, we can
	// quit
	XMLSchemaSimpleSubelementList rosetta_scripts_initial_subelements;
	for ( auto data_loader_pair : parser::DataLoaderFactory::get_instance()->loader_map() ) {
		rosetta_scripts_initial_subelements.add_already_defined_subelement(
			data_loader_pair.first, data_loader_pair.second->schema_ct_naming_function() );
	}
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"MOVERS", & parser::DataLoaderFactory::data_loader_ct_namer );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"FILTERS", & parser::DataLoaderFactory::data_loader_ct_namer );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"APPLY_TO_POSE", & rosetta_scripts_complex_type_naming_func );
	rosetta_scripts_initial_subelements.add_already_defined_subelement(
		"IMPORT", & rosetta_scripts_complex_type_naming_func );

	XMLSchemaSimpleSubelementList rosetta_scripts_protocols_subelement;
	rosetta_scripts_protocols_subelement.add_already_defined_subelement_w_alt_element_name(
		"PROTOCOLS", ParsedProtocol::mover_name(), & moves::complex_type_name_for_mover );

	XMLSchemaSimpleSubelementList rosetta_scripts_output_subelement;
	rosetta_scripts_output_subelement.add_already_defined_subelement(
		"OUTPUT", & rosetta_scripts_complex_type_naming_func );

	XMLSchemaComplexTypeGenerator rosetta_scripts_ct;
	rosetta_scripts_ct.element_name( rosetta_scripts_element_name() )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The main ROSETTASCRIPTS block allows the definition of many different types of objects which can be generated"
		" from a text file. The combination of Movers and Filters can define an end-to-end protocol, and the use of the MultiplePoseMover"
		" can define pseudo-parallel diversification and pruning of a set of structures." )
		.add_ordered_subelement_set_as_repeatable( rosetta_scripts_initial_subelements )
		.add_ordered_subelement_set_as_required( rosetta_scripts_protocols_subelement )
		.add_ordered_subelement_set_as_optional( rosetta_scripts_output_subelement )
		.write_complex_type_to_schema( xsd );


	// Make sure the schemas for the DataLoaders, for all of the Movers, and all of the Filters are
	// written to the XSD.
	for ( auto data_loader_pair : parser::DataLoaderFactory::get_instance()->loader_map() ) {
		data_loader_pair.second->provide_xml_schema( xsd );
	}
	moves::MoverFactory::get_instance()->define_mover_xml_schema( xsd );
	filters::FilterFactory::get_instance()->define_filter_xml_schema( xsd );

	// MOVERS
	XMLSchemaSimpleSubelementList movers_subelements;
	movers_subelements.add_group_subelement( & moves::MoverFactory::mover_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator movers_ct_gen;
	movers_ct_gen.element_name( "MOVERS" )
		.complex_type_naming_func( & parser::DataLoaderFactory::data_loader_ct_namer )
		.description( "The set of all of the Movers that are to be used in the script" )
		.set_subelements_repeatable( movers_subelements )
		.write_complex_type_to_schema( xsd );

	// FILTERS
	XMLSchemaSimpleSubelementList filters_subelements;
	filters_subelements.add_group_subelement( & filters::FilterFactory::filter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator filters_ct_gen;
	filters_ct_gen.element_name( "FILTERS" )
		.complex_type_naming_func( & parser::DataLoaderFactory::data_loader_ct_namer )
		.description( "The set of all of the Filters that are to be used in the script" )
		.set_subelements_repeatable( filters_subelements )
		.write_complex_type_to_schema( xsd );

	// IMPORT
	AttributeList import_attributes;
	import_attributes
		+ XMLSchemaAttribute( "taskoperations", xs_string, "A comma separated list of TaskOperations that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" )
		+ XMLSchemaAttribute( "movers", xs_string, "A comma separated list of Movers that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" )
		+ XMLSchemaAttribute( "filters", xs_string, "A comma separated list of Filters that have been"
		" defined at a higher level than the MultiplePoseMover that this IMPORT statement is inside of" );
	XMLSchemaComplexTypeGenerator import_ct_gen;
	import_ct_gen.element_name( "IMPORT" )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The IMPORT statement is meant to be used in conjuction with the MultiplePoseMover. It allows"
		" users who have defined Movers, Filters, and TaskOperations at a higher level of the protocol to avoid"
		" having to redefine them in order to use them within the MultiplePoseMover. The classes that are imported"
		" will be re-parsed within the inner context. This parsing occurs before any of sub-elements of the"
		" ROSETTASCRIPTS block are loaded, reguardless of the order in which the IMPORT statement appears. Multiple"
		" IMPORT statements may be needed, e.g. if certain TaskOperations must be loaded before other Movers which"
		" require those TaskOperations are parsed." )
		.add_attributes( import_attributes )
		.write_complex_type_to_schema( xsd );

	// OUTPUT
	// Output element only contains attributes defined by parse_score_function
	XMLSchemaComplexTypeGenerator output_ct;
	AttributeList attributes;
	rosetta_scripts::attributes_for_parse_score_function( attributes );
	output_ct
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.element_name( "OUTPUT" )
		.add_attributes( attributes )
		.description( "The OUTPUT element controls the final score function that is used by RosettaScripts before outputting a Pose."
		" If for example you have been using a score function with constraints in your protocol, but at the end, want to output the"
		" unconstrained scores, then you can set the score function in the OUTPUT element to some ScoreFunction that has constraint"
		" weights of zero and then RosettaScripts will rescore the Pose with that ScoreFunction right before writing the scores to"
		" disk." )
		.write_complex_type_to_schema( xsd );

	// APPLY_TO_POSE
	// This section is just a bunch of movers
	XMLSchemaSimpleSubelementList apply_to_pose_subelements;
	apply_to_pose_subelements.add_group_subelement( & moves::MoverFactory::mover_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator apply_to_pose_ct;
	apply_to_pose_ct.element_name( "APPLY_TO_POSE" )
		.complex_type_naming_func( & rosetta_scripts_complex_type_naming_func )
		.description( "The APPLY_TO_POSE block is an odd duck, and is deprecated, and you should not use it. What it allows you to do is"
		" define a set of movers, which the MOVERS block also allows you to do, and then to have the apply() method for these movers"
		" called immediately. As the Pose is handed to Movers and Filters in their parse_my_tag() methods, the Movers in this block"
		" have the ability to change the Pose before it is handed to the the parse_my_tag() methods of other objects. The upshot,"
		" however, is that the movers in the APPLY_TO_POSE block are going to run in the master thread, and will represent a bottleneck"
		" in multithreaded protocols. It also is just weird that the APPLY_TO_POSE block is taking over some of the functionality of"
		" the PROTOCOLS block. The original justification for why the APPLY_TO_POSE block is needed was so that you could add constraints"
		" to a Pose or change data in the DataMap, but both of these things are possible by adding Movers to the PROTOCOLS block." )
		.set_subelements_repeatable( apply_to_pose_subelements )
		.write_complex_type_to_schema( xsd );

}

std::string
RosettaScriptsParser::rosetta_scripts_element_name()
{
	return "ROSETTASCRIPTS";
}

std::string
RosettaScriptsParser::rosetta_scripts_complex_type_naming_func( std::string const & element_name )
{
	return "rosetta_scripts_parser_" + element_name + "_type";
}

} // namespace rosetta_scripts
} // namespace protocols
