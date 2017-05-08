// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/CDRSetOptionsParser.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/database/CDRSetOptionsParser.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>

#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/string_util.hh>
#include <utility/py/PyAssert.hh>
#include <basic/Tracer.hh>
//#include <utility/io/izstream.hh>

#include <basic/database/open.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <iostream>
#include <fstream>
#include <cctype>
#include <utility/io/izstream.hh>

#include <boost/algorithm/string.hpp>


static THREAD_LOCAL basic::Tracer TR("protocols.antibody.design.CDRSetOptionsParser");

namespace protocols {
namespace antibody {

using namespace core;
using namespace boost;
using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using utility::io::izstream;
using std::string;
using utility::vector1;

CDRSetOptionsParser::CDRSetOptionsParser():
	utility::pointer::ReferenceCount(),
	default_and_user_(false)
{
	ab_manager_ = AntibodyEnumManagerCOP( AntibodyEnumManagerOP( new AntibodyEnumManager() ) );
	cluster_manager_ = clusters::CDRClusterEnumManagerCOP( clusters::CDRClusterEnumManagerOP( new CDRClusterEnumManager() ) );
}


CDRSetOptionsParser::~CDRSetOptionsParser() {}

CDRSetOptionsParser::CDRSetOptionsParser( CDRSetOptionsParser const & src ):
	instructions_path_(src.instructions_path_),
	default_and_user_( src.default_and_user_)
{
	using namespace clusters;
	
	if ( src.cdr_options_) cdr_options_ = CDRSetOptionsOP( new CDRSetOptions( *src.cdr_options_ ));
	if ( src.ab_manager_ ) ab_manager_ = AntibodyEnumManagerOP( new AntibodyEnumManager( * src.ab_manager_ ));
	if ( src.cluster_manager_) cluster_manager_ = CDRClusterEnumManagerOP( new CDRClusterEnumManager( *src.cluster_manager_ ));
}

CDRSetOptionsParserOP
CDRSetOptionsParser::clone() const {
	return CDRSetOptionsParserOP( new CDRSetOptionsParser( *this ));
}

utility::vector1<CDRSetOptionsOP>
CDRSetOptionsParser::parse_default_and_user_options(std::string const & filename) {
	utility::vector1<CDRSetOptionsOP> antibody_options;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		antibody_options.push_back(parse_default_and_user_options(cdr, filename));
	}
	return antibody_options;
}

CDRSetOptionsOP
CDRSetOptionsParser::parse_default_and_user_options(CDRNameEnum cdr, std::string const & filename) {

	cdr_options_ = CDRSetOptionsOP( new CDRSetOptions(cdr) );
	std::string path = basic::options::option [basic::options::OptionKeys::antibody::design::base_cdr_instructions]();
	default_and_user_ = true;
	parse_options(cdr, path);
	parse_options(cdr, filename);
	default_and_user_ = false;
	return cdr_options_->clone();

}

utility::vector1<CDRSetOptionsOP>
CDRSetOptionsParser::parse_options(std::string const & filename) {
	utility::vector1<CDRSetOptionsOP> antibody_options;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		antibody_options.push_back(parse_options(cdr, filename));
	}
	return antibody_options;
}

CDRSetOptionsOP
CDRSetOptionsParser::parse_options(CDRNameEnum cdr, std::string const & path) {

	using namespace utility;
	using namespace std;


	if ( default_and_user_ ) {
		cdr_options_->set_cdr(cdr);
	} else {
		cdr_options_ = CDRSetOptionsOP( new CDRSetOptions(cdr) );
	}



	instructions_path_ = path;

	check_path();
	//This is straight from C++ tutorials.
	string line;
	izstream instruction_file(instructions_path_);
	if ( instruction_file.bad() ) {
		utility_exit_with_message("Unable to open grafting instruction file.");
	}
	//TR <<"Reading "<<path << " for "<< ab_manager_->cdr_name_enum_to_string(cdr) << std::endl;
	while ( getline(instruction_file, line) ) {

		//Skip any comments + empty lines
		utility::trim(line, "\n"); //Remove trailing line break
		boost::algorithm::trim(line); //Remove any whitespace on either side of the string

		//Continue to next line on empty string, comment
		if ( startswith(line, "#") || startswith(line, "\n") || line.empty()  ||  (line.find_first_not_of(' ') == std::string::npos) ) {
			continue;
		}

		boost::to_upper(line); //Capitalize entire line.

		vector1< string > lineSP = string_split_multi_delim(line); //Split on space or tab

		check_line_len(lineSP, 2);
		//TR << utility::to_string(lineSP) << std::endl;

		//Everything besides comments needs to have a CDR or ALL associated with it.
		std::string cdr_type = lineSP[1];
		std::string mode = lineSP[2];


		if ( cdr_type == "ALL" && !(cdr == l4 || cdr == h4) )  {
			parse_cdr_option(mode, lineSP);
		} else if ( ( cdr_type == "DE" || cdr_type == "CDR4") && (cdr == l4 || cdr == h4) ) {
			parse_cdr_option(mode, lineSP);
		} else if ( ab_manager_->cdr_name_is_present(cdr_type) ) {
			if ( ab_manager_->cdr_name_string_to_enum(cdr_type) == cdr ) {
				parse_cdr_option(mode, lineSP);
			}
		} else {
			//If expansion to chains, frameworks, etc.  Do it here.
			//We may have separate parsers for framework or L2.5 or whatever.
			//If its not a CDR, just skip it for now so we can have
			//TR << "Unrecognized CDR: "<<cdr_type <<" skipping...."<<std::endl;
			continue;
		}
	}
	instruction_file.close();
	//TR << "Instructions read successfully" <<std::endl;

	if ( cdr == l4 || cdr == h4 ) {
		cdr_options_->load( false ); //Can't do any DE grafting yet.  Disable, but have it in our list!
	}


	return cdr_options_->clone();
}

void
CDRSetOptionsParser::check_path() {
	using namespace std;
	izstream check( instructions_path_, ifstream::in);
	if ( check.good() ) { return;}
	else {
		ifstream check2((basic::database::full_name(instructions_path_, false)).c_str(), ifstream::in);
		if ( check2.good() ) {
			instructions_path_ = basic::database::full_name(instructions_path_);
			return;
		} else {
			utility_exit_with_message("Instructions file path not good.  Please check path.");
		}
	}
}

void
CDRSetOptionsParser::parse_cdr_option(std::string const & mode, vector1<string> const & lineSP) {



	if ( mode == "CDR_SET" || mode == "CDRS" || mode == "GRAFT_SET" || mode == "CDRSET" ) {
		check_line_len(lineSP, 3);
		std::string adjective = lineSP[3];
		parse_cdr_set_option(adjective, lineSP);
	}
}

void
CDRSetOptionsParser::check_line_len(const vector1<string> & lineSP, const Size len_check) const {
	if ( lineSP.size() < len_check ) {
		utility_exit_with_message("Could not parse cdr_set instructions. Line not long enough: "+utility::to_string(len_check)+" "+utility::to_string(lineSP));
	}
}

void
CDRSetOptionsParser::parse_cdr_set_option(std::string const & setting, vector1<string> const & lineSP) {

	//Here we match.  This is rather ugly, as I don't have much C++ expereince in this.  Python however...

	if ( lineSP.size() == 3 ) {
		std::string option = lineSP[3];
		set_cdr_set_general_option(option);
		return;
	}

	check_line_len(lineSP, 4);
	std::string type = lineSP[4];

	if ( setting == "INCLUDEONLY" ) {

		set_cdr_set_include_options(type, lineSP);
	} else if ( setting == "EXCLUDE" ) {

		set_cdr_set_exclude_options(type, lineSP);
	} else if ( setting == "LENGTH" ) {

		if ( type == "MIN" ) {
			check_line_len(lineSP, 5);
			cdr_options_->min_length(utility::string2int(lineSP[5]));
		} else if ( type == "MAX" ) {
			check_line_len(lineSP, 5);
			cdr_options_->max_length(utility::string2int(lineSP[5]));
		} else {
			utility_exit_with_message("Cannot parse cdr_set instructions. " + type);
		}
	} else if ( setting == "CLUSTER_CUTOFF" ||
			setting == "CLUSTER_CUTOFFS"||
			setting == "USE_CLUSTER_CUTOFF" ||
			setting == "USE_CLUSTER_CUTOFFS" ||
			setting == "LIMIT_CLUSTERS" ||
			setting == "LARGE_CLUSTERS_ONLY" ) {
		check_line_len(lineSP, 4);
		std::string option = lineSP[4];

		if ( option == "FALSE" || option == "NO" || option == "NONE" ) {
			cdr_options_->cluster_sampling_cutoff(0);
		} else {
			cdr_options_->cluster_sampling_cutoff(utility::string2int(option));
		}
	}
}

void
CDRSetOptionsParser::set_cdr_set_general_option(std::string const & option) {

	if ( option == "LOAD" ) {
		cdr_options_->load(true);
	} else if ( option == "NONE" ) {
		cdr_options_->load(false);
	} else if ( option == "STAYNATIVECLUSTER" || option == "STAY_NATIVE_CLUSTER" || option == "ONLY_CURRENT_CLUSTER" ) {
		cdr_options_->include_only_current_cluster(true);
	} else if ( option == "CENTERSONLY" || option == "CENTERS_ONLY" || option == "CENTER_CLUSTERS_ONLY" ) {
		cdr_options_->include_only_center_clusters(true);

	} else { utility_exit_with_message("Unrecognized option: "+option); }
}

void
CDRSetOptionsParser::set_cdr_set_include_options(std::string const & type, vector1<string> const & lineSP) {

	this->clear_cdr_set_include_options(type);
	for ( Size i=5; i<=lineSP.size(); ++i ) {

		std::string item = lineSP[i];

		if ( type == "CLUSTERS" || type == "CLUSTER" ) {
			TR << item<<std::endl;
			cdr_options_->include_only_clusters_add(cluster_manager_->cdr_cluster_string_to_enum(item));
		} else if ( type == "PDBIDS" || type == "PDBID" || type == "PDB" ) {
			cdr_options_->include_only_pdbs_add(item);
		} else if ( type == "TYPES" || type == "LENGTH_TYPES" ) {
			cdr_options_->length_type(utility::string2int(item), true);
		} else if ( type == "SPECIES" ) {
			if ( item.length() != 2 ) {
				utility_exit_with_message("Unknown species - parser requires two letter abbreviation: " + item); //Better to crash then to continue and not realize your results are Wrong.
			}
			item[1] = tolower(item[1]);
			cdr_options_->include_only_species_add(item);
		} else if ( type == "GERMLINE" || type == "GERMLINES" ) {
			cdr_options_->include_only_germlines_add(item);
		} else { utility_exit_with_message("Unrecognized cluster option"); }

	}
}

void
CDRSetOptionsParser::clear_cdr_set_include_options( std::string const & type ) {

	if ( type == "CLUSTERS" || type == "CLUSTER" ) {
		cdr_options_->include_only_clusters_clear();
	} else if ( type == "PDBIDS" || type == "PDBID" || type == "PDB" ) {
		cdr_options_->include_only_pdbs_clear();
	} else if ( type == "TYPES" || type == "LENGTH_TYPES" ) {
		for ( Size i = 1; i <=3; ++i ) {
			cdr_options_->length_type(i, false);
		}
	} else if ( type == "SPECIES" ) {
		cdr_options_->include_only_species_clear();
	} else if ( type == "GERMLINE" || type == "GERMLINES" ) {
		cdr_options_->include_only_germlines_clear();
	}
}

void
CDRSetOptionsParser::set_cdr_set_exclude_options(std::string const & type, vector1<string> const & lineSP){

	this->clear_cdr_set_exclude_options(type);
	for ( Size i=5; i<=lineSP.size(); ++i ) {

		std::string item = lineSP[i];

		if ( type == "CLUSTERS" || type == "CLUSTER" ) {
			cdr_options_->exclude_clusters_add(cluster_manager_->cdr_cluster_string_to_enum(item));
		} else if ( type == "PDBIDS" || type == "PDBID" || type == "PDB" ) {
			cdr_options_->exclude_pdbs_add(item);
		} else if ( type == "TYPES" || type == "LENGTH_TYPES" ) {
			cdr_options_->length_type(utility::string2int(item),false);
		} else if ( type == "SPECIES" ) {
			if ( item.length() != 2 ) {
				utility_exit_with_message("Unknown species - parser requires two letter abbreviation: " + item); //Better to crash then to continue and not realize your results are Wrong.
			}
			item[1] = tolower(item[1]);
			cdr_options_->exclude_species_add(item);
		} else if ( type == "GERMLINE" || type == "GERMLINES" ) {
			cdr_options_->exclude_germlines_add(item);
		} else {
			utility_exit_with_message("Unrecognized cluster option");
		}
	}
}

void
CDRSetOptionsParser::clear_cdr_set_exclude_options( std::string const & type ) {

	if ( type == "CLUSTERS" || type == "CLUSTER" ) {
		cdr_options_->exclude_clusters_clear();
	} else if ( type == "PDBIDS" || type == "PDBID" || type == "PDB" ) {
		cdr_options_->exclude_pdbs_clear();
	} else if ( type == "SPECIES" ) {
		cdr_options_->exclude_species_clear();
	} else if ( type == "GERMLINE" || type == "GERMLINES" ) {
		cdr_options_->exclude_germlines_clear();
	}
}

}
}












