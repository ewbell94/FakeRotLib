// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobOutputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class JobOutputter
/// @author Steven Lewis smlewi@gmail.com

// Unit headers
#include <protocols/jd2/JobOutputter.hh>

// package headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.fwd.hh>

// Project headers
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>


#include <basic/options/option.hh>
#include <basic/MemTracer.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <cstdlib>

///C++ headers
//#include <string> //in the .hh anyway

namespace protocols {
namespace jd2 {

static basic::MemTracer mem_tr;

JobOutputter::JobOutputter() : evaluators_(utility::pointer::make_shared< protocols::evaluation::PoseEvaluators >()) {
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluators_);
	mem_tr << "JobOutputter CSTOR" << std::endl;
	set_defaults();
}

JobOutputter::~JobOutputter() = default;

///base implementation does nothing
void JobOutputter::flush() {}

/// @details this is copied from protocols/jobdist/Jobs.cc, r24761
/// Checks the current JobInputter input source from the singleton instance of
/// JobDistributor.  If a PDB_FILE, it will deliberately discard any path
/// information in the input tag as well as any file name extension (since
/// input tags are usually file names).  Otherwise, it will keep the entire input
/// tag as-is.
std::string JobOutputter::affixed_numbered_name( JobCOP job ){
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	if ( option[out::mover_controlled_output_name]() ) {
		std::string pose_name = job->inner_job()->optional_output_name();
		if ( pose_name.size()>0 ) {
			return(pose_name);
		}
	}
	// Use at least 4 digits in number to match Rosetta++
	core::Size nstruct_width = 0;
	core::Size const nstruct = job->nstruct_max();
	for ( core::Size i = 1; i <= nstruct || nstruct_width < 4; i *= 10 ) {
		nstruct_width += 1;
	}

	// Construct the base name.  For pdb files we need to hack away the paths
	// and extensions.  For all other types of input sources, assume that the
	// string that comes in is exactly the tag we want.  This is important
	// for silent file <-> BOINC usage to avoid problems with duplicate tags
	// caused by path/ext removal.
	std::string base_name = job->input_tag();
	switch ( JobDistributor::get_instance()->job_inputter_input_source() ) {
	case JobInputterInputSource::SERIALIZED_POSE_FILE:
	case JobInputterInputSource::PDB_FILE : { // remove paths and ext.

		// Treat tags as file names so that we put the number before the extension.
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = out_name.base();

		break;
	}
	case JobInputterInputSource::SCREENING_FILE : {

		//combine file names and group name.
		Job::StringStringPairs string_pairs =job->get_string_string_pairs();

		auto group_name_it(string_pairs.find("input_group_name"));
		debug_assert(group_name_it != string_pairs.end()); //You shouldn't be able to get this far without a properly set group name
		std::string group_name(group_name_it->second);
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = group_name + "_" + out_name.base();

		break;
	}
	default : { // do nothing
		break;
	}
	}

	// now construct the full name
	std::string prefix( prefix_ );
	if ( job->status_prefix().size() ) {
		prefix = job->status_prefix()+"_";
	}

	//If the user requests it, mangle the job with the $HOSTNAME env var.
	//This is useful when job distribution is being done outside of Rosetta (for
	//example in kubernetes pods) and the hostnames are distinct hashes, as a
	//control for the fact that all the nstructs are 0001.
	//Note that currently this looks for a particular env var and
	//mangles thereby, but it would be easy to take a list of env vars and
	//mangle with all of them.  $HOSTNAME is common on containerized systems.
	if ( option[basic::options::OptionKeys::jd2::HOSTNAME_in_jobname].value() ) {
		//getenv returns a char* if found, or a null pointer if the env var is not found
		//we need to catch the nullptr case before trying to make a string out of it
		const char* hostname_cstar(std::getenv("HOSTNAME"));
		if ( hostname_cstar != nullptr ) {
			prefix = std::string(hostname_cstar) + "_" + prefix;
		} // do nothing if else; we could append the _ anyway I guess.
	}

	std::ostringstream oss;
	oss << prefix << base_name << suffix_;
	if ( ! no_nstruct_label_ || job->nstruct_index() != 1 ) {
		oss << '_' << std::setfill('0') << std::setw(nstruct_width) << job->nstruct_index();
	}
	return oss.str();
}

/// @details run the PoseEvaluators on the pose
/// evaluation creates string value pairs which end up in the SilentStruct energy object
/// instead of filling things into a SilentStruct we could provide a different interface...
///  wait until Steven has finished his string/value pair output
void JobOutputter::evaluate(
	core::pose::Pose & pose,
	std::string tag,
	core::io::silent::SilentStruct & pss
) const {
	if ( evaluators()->size() ) {
		evaluators()->apply( pose, tag, pss );
	}
}

/// @brief optionally pass a starting (reference) pose to a JobOutputter for comparison purposes and/or as interface for initializing evaluators. (Currently does nothing in this base class.)
void JobOutputter::starting_pose( core::pose::Pose const & ){}

/// @details add another PoseEvaluator to the list of evaluations
void JobOutputter::add_evaluation( evaluation::PoseEvaluatorOP ev_in ) {
	evaluators_->add_evaluation( ev_in );
}

/// @details set a list of Evaluations
/// ( the list will be copied, the evaluations are OPs )
///
void JobOutputter::set_evaluators( evaluation::PoseEvaluators const& ev_in ) {
	evaluators_ = utility::pointer::make_shared< protocols::evaluation::PoseEvaluators >(ev_in);
}

/// @brief clear the list of evaluators
void JobOutputter::clear_evaluators() {
	evaluators_->clear();
}

/// @brief clear the list of evaluators
void JobOutputter::set_defaults() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	prefix_ = option[ out::prefix ]();
	suffix_ = option[ out::suffix ]();
	no_nstruct_label_ = option[ out::no_nstruct_label ]();
}

/// @details return the list of PoseEvaluators
evaluation::PoseEvaluatorsCOP JobOutputter::evaluators() const {
	return evaluators_;
}

void JobOutputter::call_output_observers( core::pose::Pose const& pose, JobOP job  ) const {
	if ( !job ) return;
	job->call_output_observers( pose );
}

} // jd2
} // protocols
