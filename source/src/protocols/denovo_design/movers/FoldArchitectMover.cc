// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/FoldArchitectMover.cc
/// @brief Mover that builds and folds a structure via fragment insertion
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/FoldArchitectMover.hh>
#include <protocols/denovo_design/movers/FoldArchitectMoverCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/DivideAndConqueror.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/NullPoseFolder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/StructureDataPerturber.hh>
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh>
#include <protocols/denovo_design/movers/SealFoldTreeMover.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/filters/CalculatorFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/rosetta_scripts/util.hh>
// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost headers

// C++ headers
#include <fstream>
#include <ctime>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "protocols.denovo_design.movers.FoldArchitectMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

FoldArchitectMover::FoldArchitectMover():
	protocols::moves::Mover( FoldArchitectMover::mover_name() ),
	architect_(),
	builder_( new components::ExtendedPoseBuilder ),
	folder_( new components::RemodelLoopMoverPoseFolder ),
	perturber_( new components::NullPerturber ),
	prefold_movers_(),
	postfold_movers_(),
	filters_(),
	score_filter_(),
	id_(),
	dry_run_( false ),
	dump_pdbs_( false ),
	debug_( false ),
	build_overlap_( 1 ),
	iterations_per_phase_( 0 ),
	start_segments_(),
	stop_segments_()
{
	set_debug( debug_ );
}

FoldArchitectMover::~FoldArchitectMover() = default;

void
FoldArchitectMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	set_id( tag->getOption< std::string >( "name" ) );

	set_build_overlap( tag->getOption< core::Size >( "build_overlap", build_overlap_ ) );

	std::string const segments_csv = tag->getOption< std::string >( "start_segments", "" );
	if ( tag->hasOption( "start_segments" ) ) set_start_segments( segments_csv );

	std::string const stop_segments_csv = tag->getOption< std::string >( "stop_segments", "" );
	if ( tag->hasOption( "stop_segments" ) ) set_stop_segments( stop_segments_csv );

	set_iterations_per_phase( tag->getOption< core::Size >( "iterations_per_phase", iterations_per_phase_ ) );
	set_dump_pdbs( tag->getOption< bool >( "dump_pdbs", dump_pdbs_ ) );
	set_debug( tag->getOption< bool >( "debug", debug_ ) );

	if ( tag->hasOption( "xml" ) ) {
		setup_from_xml_file( tag->getOption< std::string >( "xml" ), data );
	} else {
		setup_from_xml_tag( tag, data );
	}
}

void
FoldArchitectMover::setup_from_xml_file( std::string const & xml_file, basic::datacache::DataMap & data )
{
	std::ifstream is( xml_file );
	utility::tag::TagCOP tag = utility::tag::Tag::create( is );
	setup_from_xml_tag( tag, data );
}

void
FoldArchitectMover::setup_from_xml_string( std::string const & xml_string, basic::datacache::DataMap & data )
{
	utility::tag::TagCOP tag = utility::tag::Tag::create( xml_string );
	setup_from_xml_tag( tag, data );
}

void
FoldArchitectMover::setup_from_xml_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	clear_prefold_movers();
	clear_postfold_movers();
	clear_filters();

	core::Size idx = 1;
	for ( utility::tag::TagCOP subtag : tag->getTags() ) {
		if ( subtag->getName() == "PreFoldMovers" ) parse_prefold_movers( subtag, data );
		else if ( subtag->getName() == "PostFoldMovers" ) parse_postfold_movers( subtag, data );
		else if ( subtag->getName() == "Filters" ) parse_filters( subtag, data );
		else if ( idx == 1 ) {
			parse_architect( subtag, data );
			++idx;
		} else if ( idx == 2 ) {
			parse_folder( subtag, data );
			++idx;
		} else if ( idx == 3 ) {
			parse_perturber( subtag, data );
			++idx;
		} else {
			std::stringstream msg;
			msg << get_name() << "::parse_my_tag(): Invalid subtag found:" << *subtag << std::endl;
			msg << "Valid subtags are \"PreFoldMovers\", \"PostFoldMovers\", a DeNovo architect, and a Pose Folder."
				<< std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
	}
	TR << "Finished parsing tag: Found " << prefold_movers_.size() << " prefold movers, "
		<< postfold_movers_.size() << " postfold movers, and " << filters_.size()
		<< " filters." << std::endl;
}

protocols::moves::MoverOP
FoldArchitectMover::clone() const
{
	return utility::pointer::make_shared< FoldArchitectMover >( *this );
}


protocols::moves::MoverOP
FoldArchitectMover::fresh_instance() const
{
	return utility::pointer::make_shared< FoldArchitectMover >();
}



std::string const &
FoldArchitectMover::id() const
{
	return id_;
}

void
FoldArchitectMover::set_id( std::string const & id_value )
{
	id_ = id_value;
}

void
FoldArchitectMover::set_iterations_per_phase( core::Size const niter )
{
	iterations_per_phase_ = niter;
}

void
FoldArchitectMover::set_debug( bool const debug )
{
	debug_ = debug;
	builder_->set_debug( debug );
}

protocols::filters::FilterOP
FoldArchitectMover::default_score_filter()
{
	protocols::filters::CalculatorFilterOP calc_filt( new protocols::filters::CalculatorFilter( "-SS" ) );
	protocols::filters::FilterOP ss_filt( new protocols::fldsgn::filters::SecondaryStructureFilter );
	calc_filt->add_filter( "SS", ss_filt );
	return calc_filt;
}

void
FoldArchitectMover::apply( core::pose::Pose & pose )
{
	if ( !architect_ ) {
		std::stringstream msg;
		msg << get_name() << "::apply(): No architect is set!"
			<< get_name() << " requires an architect to build and modify poses." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// create StructureData "blueprint"
	components::StructureDataOP sd = architect_->apply( pose );
	if ( !sd ) {
		std::stringstream msg;
		msg << "FoldArchitectMover::architect " << architect_->id()
			<< " failed to generate anything." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// Break up into managable pieces
	components::DivideAndConqueror divider;
	divider.set_start_segments( start_segments_ );
	divider.set_stop_segments( stop_segments_ );

	components::BuildPhases const build_phases = divider.divide_and_conquer( *sd );

	TR << "Trying to build the following structure:" << std::endl;
	TR << *sd << std::endl;

	core::pose::PoseOP pose_ptr = build_in_phases( *sd, build_phases, SegmentNameSet(), 1 );
	if ( pose_ptr ) pose = *pose_ptr;
	else set_last_move_status( protocols::moves::FAIL_RETRY );
}

/// @brief given an build phase number, returns a pdb filename for checkpointing
///        Example filename: iter_01_20160718133800.pdb
/// @param[in] phase_num  Build phase number
std::string
pdb_filename( core::Size const phase_num )
{
	std::stringstream filename;
	filename << time( nullptr ) << "_iter_" << std::setw( 2 )
		<< std::setfill('0') << phase_num << ".pdb";
	return filename.str();
}

/// @brief builds/folds pose in phases using recursive algorithm
/// @throws EXCN_Fold if we couldn't fold the pose
core::pose::PoseOP
FoldArchitectMover::build_in_phases(
	components::StructureData const & full_sd,
	components::BuildPhases const & phases,
	SegmentNameSet const & finished,
	core::Size const phase_num ) const
{
	core::Size const max_iter = iterations_per_phase_ == 0 ? 50 : iterations_per_phase_;

	SegmentNames const & phase = *phases.begin();
	TR << "Building phase " << phase_num << " " << phase << std::endl;

	// compute segments to slice out of StructureData
	SegmentNameSet all_phase_segments( finished.begin(), finished.end() );
	all_phase_segments.insert( phase.begin(), phase.end() );

	components::StructureData working_sd = full_sd;
	components::StructureDataPerturberOP phase_perturber = perturber_->clone();

	// compute segments not in this phase
	SegmentNameSet const in_this_phase( phase.begin(), phase.end() );
	SegmentNameSet const in_sd( full_sd.segments_begin(), full_sd.segments_end() );
	SegmentNames const ignore_vec = set_difference< SegmentNames, SegmentName >( in_sd.begin(), in_sd.end(), in_this_phase.begin(), in_this_phase.end() );
	SegmentNameSet const ignore_set( ignore_vec.begin(), ignore_vec.end() );

	phase_perturber->set_ignore_segments( ignore_set );
	TR << "Set ignored segments to " << ignore_set << std::endl;

	for ( core::Size iter=1; iter<=max_iter; ++iter ) {
		TR << "Attempting to fold phase " << phase_num << " " << phase << " "
			<< " iteration " << iter << " / " << max_iter << std::endl;

		// perturb
		phase_perturber->apply( working_sd );

		// slice
		TR << "Slicing " << all_phase_segments << " out of StructureData" << std::endl;
		components::StructureData const phase_sd = working_sd.slice( all_phase_segments, true );
		TR << "Done slicing. New SD=" << phase_sd << std::endl;

		// build pose
		core::pose::PoseOP built = builder_->apply( phase_sd );
		try {
			fold_attempt( *built );

			/// dump pdb
			if ( dump_pdbs_ ) built->dump_pdb( pdb_filename( phase_num ) );

			/// add template information based on the folded segments, so that the next phase can use it
			/// this will never be reached if check_pose() throws an exception
			components::StructureData new_full_sd = working_sd;
			components::StructureData const & built_sd = components::StructureDataFactory::get_instance()->get_from_const_pose( *built );
			for ( auto s=built_sd.segments_begin(); s!=built_sd.segments_end(); ++s ) {
				new_full_sd.set_template_pose( *s, *built, built_sd.segment( *s ).start(), built_sd.segment( *s ).stop() );
			}

			if ( ++phases.begin() == phases.end() ) {
				return built;
			} else {
				return build_in_phases(
					new_full_sd,
					components::BuildPhases( ++phases.begin(), phases.end() ),
					SegmentNameSet( all_phase_segments.begin(), all_phase_segments.end() ),
					phase_num + 1 );
			}
		} catch ( components::EXCN_Fold const & e ) {
			e.show( TR );
			TR << std::endl << "Failing phase " << phase_num << " " << phase << std::endl;
		} catch ( EXCN_NothingToFold const & e ) {
			e.show( TR );
			TR.flush();
			TR.Warning << "not folding anything!" << std::endl;
		} catch (EXCN_FilterFailed const & e ) {
			e.show( TR );
			TR.flush();
		}
	}

	// if we reach here without returning something, folding failed.
	std::stringstream msg;
	msg << mover_name() << "build_in_phases():  Failed to fold anything passing filters after "
		<< max_iter << " attempts." << std::endl;
	throw CREATE_EXCEPTION(components::EXCN_Fold, msg.str() );
	return core::pose::PoseOP();
}

/// @brief checks the pose for whether it folded correctly, according to user's filters
void
FoldArchitectMover::check_pose( core::pose::Pose const & pose ) const
{
	core::Size filter_num = 1;
	for ( auto f=filters_.begin(); f!=filters_.end(); ++f, ++filter_num ) {
		if ( !(*f)->apply( pose ) ) throw CREATE_EXCEPTION(EXCN_FilterFailed, (*f)->get_type(), filter_num );
	}
}

SegmentNames
FoldArchitectMover::segments_to_fold( components::StructureData const & sd ) const
{
	SegmentNames segs;
	for ( auto s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		if ( !sd.segment( *s ).template_pose() ) segs.push_back( *s );
	}
	return segs;
}

SegmentNames
FoldArchitectMover::find_roots(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops ) const
{
	SegmentNames roots;

	components::StructureData const & sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	for ( auto const & loop : loops ) {
		core::Size const startres =
			loop_start_without_overlap( pose, loop.start(), build_overlap_ );
		std::string const & start_comp = sd.segment_name( startres );
		std::string const & prevcomp = sd.segment( start_comp ).lower_segment();
		if ( prevcomp != "" ) {
			roots.push_back( prevcomp );
		}
		core::Size const stopres =
			loop_stop_without_overlap( pose, loop.stop(), build_overlap_ );
		std::string const & stop_comp = sd.segment_name( stopres );
		std::string const & nextcomp = sd.segment(stop_comp).upper_segment();
		if ( nextcomp != "" ) {
			roots.push_back( nextcomp );
		}
		TR.Debug << id() << ": Loop: " << startres << " " << prevcomp << "__" << start_comp << " --> " << stopres << " " << stop_comp << "__" << nextcomp << std::endl;
	}
	if ( roots.empty() ) {
		roots.push_back( *sd.segments_begin() );
	}
	return roots;
}

void
FoldArchitectMover::set_perturber( components::StructureDataPerturber const & perturber )
{
	set_perturber( perturber.clone() );
}

void
FoldArchitectMover::set_perturber( components::StructureDataPerturberOP perturber )
{
	perturber_ = perturber;
}

void
FoldArchitectMover::parse_perturber( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	set_perturber( components::StructureDataPerturber::create( *tag, data ) );
}

void
FoldArchitectMover::parse_prefold_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	MoverOPs const & movers = parse_movers( tag, data );
	for ( auto const & mover : movers ) {
		add_prefold_mover( *mover );
	}
}

void
FoldArchitectMover::parse_postfold_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	MoverOPs const & movers = parse_movers( tag, data );
	for ( auto const & mover : movers ) {
		add_postfold_mover( *mover );
	}
}

FoldArchitectMover::MoverOPs
FoldArchitectMover::parse_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) const
{
	MoverOPs movers;
	for ( auto subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() != "Add" ) {
			std::stringstream msg;
			msg << type() << ": Invalid xml tag name (" << (*subtag)->getName() << ") found in tag " << *tag << std::endl;
			msg << "Valid tags are: \"Add\"." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
		std::string const mover_name = (*subtag)->getOption< std::string >( "mover" );
		protocols::moves::MoverOP mover = protocols::rosetta_scripts::parse_mover_or_null( mover_name, data );
		if ( ! mover ) {
			std::stringstream msg;
			msg << type() << "::parse_movers(): ERROR !! mover not found in map: \n" << **subtag << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
		movers.push_back( mover );
		TR.Debug << "found mover " << mover_name << std::endl;
	}
	return movers;
}

void
FoldArchitectMover::parse_filters( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	clear_filters();
	for ( auto subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() != "Add" ) {
			std::stringstream msg;
			msg << type() << ": Invalid xml tag name (" << (*subtag)->getName() << ") found in tag " << *tag << std::endl;
			msg << "Valid tags are: \"Add\"." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
		std::string const filter_name = (*subtag)->getOption< std::string >( "filter" );
		protocols::filters::FilterOP filter_ptr = protocols::rosetta_scripts::parse_filter_or_null( filter_name, data );
		if ( ! filter_ptr ) {
			std::stringstream msg;
			msg << type() << "::parse_filters(): ERROR !! filter not found in map: \n" << **subtag << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
		add_filter( *filter_ptr );
		TR.Debug << "found filter " << filter_name << std::endl;
	}
}

void
FoldArchitectMover::add_filter( protocols::filters::Filter const & filter )
{
	filters_.push_back( filter.clone() );
}

void
FoldArchitectMover::clear_filters()
{
	filters_.clear();
}

FoldArchitectMover::MoverOPs
FoldArchitectMover::prefold_movers(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	ConstraintGeneratorCOPs const & generators ) const
{
	MoverOPs movers;

	// 1. switch to centroid
	protocols::moves::MoverOP switch_cen( new protocols::simple_moves::SwitchResidueTypeSetMover( "centroid" ) );
	movers.push_back( switch_cen );

	// 2. Set secondary structure in pose
	SetPoseSecstructFromStructureDataMover set_ss;
	movers.push_back( set_ss.clone() );

	// 3. Run user-provided movers
	for ( auto const & prefold_mover : prefold_movers_ ) {
		movers.push_back( prefold_mover );
	}

	// 4. Add coordinate constraints for overlapping regions
	if ( !generators.empty() ) {
		protocols::moves::MoverOP add_csts( new protocols::constraint_generator::AddConstraints( generators ) );
		if ( add_csts ) movers.push_back( add_csts );
	}

	// 4. Build proper fold tree for remodel
	SegmentNames const roots = find_roots( pose, loops );
	protocols::moves::MoverOP set_ft( new FoldTreeFromFoldGraphMover( roots, loops ) );
	movers.push_back( set_ft );

	return movers;
}

FoldArchitectMover::MoverOPs
FoldArchitectMover::postfold_movers(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	ConstraintGeneratorCOPs const & generators ) const
{
	MoverOPs movers;

	components::StructureData const sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	// 1. Remove coodinate constraints for overlapping regions
	//    if they were added
	if ( !generators.empty() ) {
		protocols::moves::MoverOP rm_csts( new protocols::constraint_generator::RemoveConstraints( generators, false ) );
		if ( rm_csts ) movers.push_back( rm_csts );
	}

	// 2. Seal fold tree
	protocols::moves::MoverOP seal_ft( new SealFoldTreeMover( sd, loops ) );
	movers.push_back( seal_ft );

	// 2. Run user-provided movers
	for ( auto const & postfold_mover : postfold_movers_ ) {
		movers.push_back( postfold_mover );
	}

	// 3. Switch to full atom
	protocols::moves::MoverOP switch_fa( new protocols::simple_moves::SwitchResidueTypeSetMover( "fa_standard" ) );
	movers.push_back( switch_fa );

	// 4. Replace sidechains
	protocols::moves::MoverOP restore_sc_mover(
		new simple_moves::ReturnSidechainMover( pose, 1, pose.size() ) );
	//protocols::moves::MoverOP sars(
	// new simple_moves::SaveAndRetrieveSidechains(
	// pose,
	// true, // allsc -- replace all sidechains, as opposed to only ALAs
	// false, // ensure_variant_matching
	// 0 ) ); // jumpid -- no idea why this needs to be set
	movers.push_back( restore_sc_mover );

	// 5. Re-seal fold tree for full-atom mode
	//    This is necessary if there are sidechain-sidechain covalent bonds
	movers.push_back( seal_ft );

	return movers;
}

void
FoldArchitectMover::clear_prefold_movers()
{
	prefold_movers_.clear();
}

void
FoldArchitectMover::add_prefold_mover( protocols::moves::Mover const & prefold_mover )
{
	prefold_movers_.push_back( prefold_mover.clone() );
}

void
FoldArchitectMover::clear_postfold_movers()
{
	postfold_movers_.clear();
}

void
FoldArchitectMover::add_postfold_mover( protocols::moves::Mover const & postfold_mover )
{
	postfold_movers_.push_back( postfold_mover.clone() );
}

protocols::loops::Loops
FoldArchitectMover::create_loops(
	core::pose::Pose const & pose,
	ResidueVector & overlap_residues,
	SegmentNames & movable_segments ) const
{
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	components::StructureData const & sd = factory.get_from_const_pose( pose );
	components::FoldGraph fg( sd );
	movable_segments = segments_to_fold( sd );
	protocols::loops::LoopsOP loops = fg.create_loops( movable_segments );
	if ( !loops ) {
		std::stringstream msg;
		msg << type() << "::create_loops(): Loops object could not be created! Movable segments="
			<< movable_segments << " SD=" << sd << std::endl;
		throw CREATE_EXCEPTION(EXCN_NothingToFold, msg.str() );
	}

	// add requested connection overlap
	overlap_residues = add_overlap_to_loops( *loops, build_overlap_, pose );
	return *loops;
}

core::scoring::ScoreFunctionOP
get_score_function( core::pose::Pose const & )
{
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	return scorefxn;
}

void
FoldArchitectMover::fold_attempt( core::pose::Pose & pose ) const
{
	// Get StructureData object from pose
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	if ( !factory.has_cached_data( pose ) ) {
		utility_exit_with_message( "No StructureData at start of fold attempt" );
	}

	// make Loops object to pass to PoseFolder
	ResidueVector overlap_residues;
	SegmentNames movable_segments;
	protocols::loops::Loops const loops = create_loops( pose, overlap_residues, movable_segments );
	if ( loops.empty() ) {
		TR.Warning << "Loops object is empty.  Not doing anything." << std::endl;
		return;
	}

	// Make movable residue subset for PoseFolder
	core::select::residue_selector::ResidueSubset const movable =
		select_movable_residues( pose, movable_segments, overlap_residues );

	// Create pre-fold and post-fold movers
	ConstraintGeneratorCOPs const generators = create_constraint_generators( overlap_residues );

	components::StructureData const sd = factory.get_from_pose( pose );

	// apply user-provided movers
	MoverOPs const prefold = prefold_movers( pose, loops, generators );
	MoverOPs const postfold = postfold_movers( pose, loops, generators );

	// apply user-provided movers
	if ( !factory.has_cached_data( pose ) ) {
		utility_exit_with_message( "No data before applying premovers" );
	}
	apply_movers( prefold, pose );

	// fold structure
	// folders need to keep sd and pose in sync
	folder_->apply( pose, movable, loops );

	try {
		factory.get_from_pose( pose ).check_pose_consistency( pose );
	} catch ( components::EXCN_PoseInconsistent const & e ) {
		std::stringstream msg;
		msg << "FoldArchitectMover::fold_attempt(): After pose folder, StructureData no longer agrees with pose."
			<< " PoseFolders need to ensure that StructureData and pose are synced.  Error =" << std::endl;
		e.show( msg );
		utility_exit_with_message( msg.str() );
	}

	// apply user-provided movers
	apply_movers( postfold, pose );

	// set Energies object -- sets hbonds info in the Energies object in the pose
	pose.energies().clear();
	(*get_score_function(pose))( pose );

	// check pose against user-provided filters
	check_pose( pose );
}

void
FoldArchitectMover::remove_cutpoints( components::StructureData & sd, protocols::loops::Loops const & loops ) const
{
	for ( auto const & loop : loops ) {
		if ( ! loop.cut() ) continue;
		TR.Debug << "Removing cutpoint from StructureData: " << loop.cut() << std::endl;
		sd.set_cutpoint( sd.segment_name( loop.cut() ), 0 );
	}
}

void
FoldArchitectMover::apply_movers( MoverOPs const & movers, core::pose::Pose & pose ) const
{
	//pose.dump_pdb( "prefold_0.pdb" );
	core::Size count = 1;
	for ( auto m=movers.begin(); m!=movers.end(); ++m, ++count ) {
		TR.Debug << "Running mover " << (*m)->get_name() << std::endl;
		(*m)->apply( pose );
		//pose.dump_pdb( "prefold_" + (*m)->get_name() + boost::lexical_cast< std::string >( count ) + ".pdb" );
	}
}

core::select::residue_selector::ResidueSubset
FoldArchitectMover::select_movable_residues(
	core::pose::Pose const & pose,
	SegmentNames const & movable_segments,
	ResidueVector const & overlap_residues ) const
{
	components::StructureData const & sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	core::select::residue_selector::ResidueSubset subset( sd.pose_length(), false );
	for ( auto const & movable_segment : movable_segments ) {
		for ( core::Size resid=sd.segment( movable_segment ).lower(); resid<=sd.segment( movable_segment ).upper(); ++resid ) {
			subset[ resid ] = true;
		}
	}
	for ( unsigned long overlap_residue : overlap_residues ) {
		debug_assert( overlap_residue > 0 );
		debug_assert( overlap_residue <= pose.size() );
		subset[ overlap_residue ] = true;
	}
	return subset;
}

FoldArchitectMover::ConstraintGeneratorCOPs
FoldArchitectMover::create_constraint_generators( ResidueVector const & residues ) const
{
	if ( residues.empty() ) {
		return FoldArchitectMover::ConstraintGeneratorCOPs();
	}
	return { create_overlap_constraint_generator( residues ) };
}

/// @brief builds a constraint generator for residues that overlap between phases
/// @param[in] residues Vector of resids corresponding to the overlap residues
FoldArchitectMover::ConstraintGeneratorOP
FoldArchitectMover::create_overlap_constraint_generator( ResidueVector const & residues ) const
{
	using core::select::residue_selector::ResidueIndexSelector;

	protocols::constraint_generator::AtomPairConstraintGeneratorOP dist_csts(
		new protocols::constraint_generator::AtomPairConstraintGenerator );

	dist_csts->set_id( "FoldArchitectMover_constrain_overlap_residues" );
	dist_csts->set_ca_only( false );
	dist_csts->set_min_seq_sep( 0 );
	dist_csts->set_sd( 1.0 );

	std::stringstream index_ss;
	for ( auto r=residues.begin(); r!=residues.end(); ++r ) {
		if ( r != residues.begin() ) index_ss << ',';
		index_ss << *r;
	}

	TR << "Creating index selector from " << index_ss.str() << std::endl;
	ResidueIndexSelector const selector( index_ss.str() );
	dist_csts->set_residue_selector( selector );
	return dist_csts;
}

FoldArchitectMover::FoldScore
FoldArchitectMover::folding_score( core::pose::Pose const & pose ) const
{
	if ( !score_filter_ ) return 0.0;

	TR << "Computing score" << std::endl;
	core::Real const score = score_filter_->report_sm( pose );
	TR << "Score is " << score << std::endl;
	return score;
}

void
FoldArchitectMover::set_score_filter( protocols::filters::Filter const & filter )
{
	score_filter_ = filter.clone();
}

void
FoldArchitectMover::parse_architect( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	set_architect( *architects::DeNovoArchitectFactory::get_instance()->create_from_tag( tag, data ) );
}

void
FoldArchitectMover::parse_folder( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	components::PoseFolderOP folder;
	if ( tag->getName() == components::RemodelLoopMoverPoseFolder::class_name() ) {
		folder = utility::pointer::make_shared< components::RemodelLoopMoverPoseFolder >();
	} else if ( tag->getName() == components::RandomTorsionPoseFolder::class_name() ) {
		folder = utility::pointer::make_shared< components::RandomTorsionPoseFolder >();
	} else if ( tag->getName() == components::NullPoseFolder::class_name() ) {
		folder = utility::pointer::make_shared< components::NullPoseFolder >();
	} else {
		std::stringstream msg;
		msg << "FoldArchitectMover::parse_folder(): Unknown folder type: "
			<< tag->getName() << std::endl;
		utility_exit_with_message( msg.str() );
	}

	folder->parse_my_tag( tag, data );
	set_folder( folder );
}

void
FoldArchitectMover::set_architect( architects::DeNovoArchitect const & architect )
{
	architect_ = architect.clone();
}

void
FoldArchitectMover::set_folder( components::PoseFolder const & folder )
{
	set_folder( folder.clone() );
}

void
FoldArchitectMover::set_folder( components::PoseFolderOP folder )
{
	folder_ = folder;
}

void
FoldArchitectMover::set_builder( components::PoseBuilder const & builder )
{
	set_builder( builder.clone() );
}

void
FoldArchitectMover::set_builder( components::PoseBuilderOP builder )
{
	builder_ = builder;
}

void
FoldArchitectMover::set_build_overlap( core::Size const overlap_val )
{
	build_overlap_ = overlap_val;
}

/// @brief sets names of segments to be included in the starting build phase
/// @param[in] segments_csv Comma-separated string containing segment names
void
FoldArchitectMover::set_start_segments( std::string const & segments_csv )
{
	set_start_segments( csv_to_container< SegmentNameSet >( segments_csv ) );
}

/// @brief sets names of segments to be included in the starting build phase
/// @param[in] segments Set of segment names
void
FoldArchitectMover::set_start_segments( SegmentNameSet const & segments )
{
	start_segments_ = segments;
}

/// @brief sets names of segments to be included in the final build phase
/// @param[in] segments_csv Comma-separated string containing segment names
void
FoldArchitectMover::set_stop_segments( std::string const & segments_csv )
{
	set_stop_segments( csv_to_container< SegmentNameSet >( segments_csv ) );
}

/// @brief sets names of segments to be included in the final build phase
/// @param[in] segments Set of segment names
void
FoldArchitectMover::set_stop_segments( SegmentNameSet const & segments )
{
	stop_segments_ = segments;
}


/////////////// Creator ///////////////

std::string FoldArchitectMover::get_name() const {
	return mover_name();
}

std::string FoldArchitectMover::mover_name() {
	return "FoldArchitect";
}

std::string FoldArchitectMover::folder_ct_namer( std::string const & folder_name ){
	return "denovo_folder_" + folder_name + "_complex_type";
}

std::string FoldArchitectMover::perturber_ct_namer( std::string const & perturber_name ){
	return "denovo_perturber_" + perturber_name + "_complex_type";
}

std::string FoldArchitectMover::folder_group_name(){
	return "denovo_folder";
}

std::string FoldArchitectMover::perturber_group_name(){
	return "denovo_perturber";
}


std::string FoldArchitectMover::build_denovo_backbone_ct_naming_func( std::string const & subelement_name ){
	return "build_denovo_backbone_" + subelement_name + "_complex_type";
}

std::string FoldArchitectMover::prefold_ct_naming_func( std::string const & subelement_name ){
	return "denovo_prefold_" + subelement_name + "_complex_type";
}

std::string FoldArchitectMover::postfold_ct_naming_func( std::string const & subelement_name ){
	return "denovo_postfold_" + subelement_name + "_complex_type";
}

std::string FoldArchitectMover::filters_ct_naming_func( std::string const & subelement_name ){
	return "denovo_filters_" + subelement_name + "_complex_type";
}

void FoldArchitectMover::define_perturber_group( utility::tag::XMLSchemaDefinition & xsd ){

	using namespace utility::tag;
	AttributeList helix_or_connection_perturber_attlist;
	helix_or_connection_perturber_attlist + XMLSchemaAttribute::required_attribute( "architect", xs_string, "Architect to use for this mover" );
	AttributeList null_perturber_attlist; //empty

	//Define the complex types for all of these
	//Helix and connection have the same attributes
	XMLSchemaComplexTypeGenerator helix_perturber_ct;
	helix_perturber_ct.element_name( "HelixPerturber" )
		.complex_type_naming_func ( & perturber_ct_namer )
		.description( "Used to perturb helices" )
		.add_attributes( helix_or_connection_perturber_attlist )
		.write_complex_type_to_schema( xsd );

	XMLSchemaComplexTypeGenerator connection_perturber_ct;
	connection_perturber_ct.element_name( "ConnectionPerturber" )
		.complex_type_naming_func ( & perturber_ct_namer )
		.description( "Used to perturb connections" )
		.add_attributes( helix_or_connection_perturber_attlist )
		.write_complex_type_to_schema( xsd );

	//Null has no attributes
	XMLSchemaComplexTypeGenerator null_perturber_ct;
	null_perturber_ct.element_name( "NullPerturber" )
		.complex_type_naming_func ( & perturber_ct_namer )
		.description( "Doesn't do anything" )
		.write_complex_type_to_schema( xsd );

	//Compound takes other perturbers as subtags

	XMLSchemaRestriction mode;
	mode.name( "compound_perturber_mode" );
	mode.base_type( xs_string );
	mode.add_restriction( xsr_enumeration, "AND" );
	mode.add_restriction( xsr_enumeration, "OR" );
	xsd.add_top_level_element( mode );

	AttributeList compound_perturber_attlist;
	compound_perturber_attlist
		+ XMLSchemaAttribute::attribute_w_default( "mode", "compound_perturber_mode", "Which logical operator to apply to these perturbers?", "AND" );
	XMLSchemaSimpleSubelementList compound_perturber_subelements;
	compound_perturber_subelements
		.add_group_subelement( & perturber_group_name );
	XMLSchemaComplexTypeGenerator compound_perturber_ct;
	compound_perturber_ct.element_name( "CompoundPerturber" )
		.complex_type_naming_func ( & perturber_ct_namer )
		.description( "Used to combine other perturbers" )
		.set_subelements_repeatable( compound_perturber_subelements)
		.add_attributes( compound_perturber_attlist )
		.write_complex_type_to_schema( xsd );

	XMLSchemaElementOP helix_perturber_subelement( new XMLSchemaElement );
	helix_perturber_subelement->name( "HelixPerturber" );
	helix_perturber_subelement->type_name( perturber_ct_namer( "HelixPerturber" ) );

	XMLSchemaElementOP connection_perturber_subelement( new XMLSchemaElement );
	connection_perturber_subelement->name( "ConnectionPerturber" );
	connection_perturber_subelement->type_name( perturber_ct_namer( "ConnectionPerturber" ) );

	XMLSchemaElementOP null_perturber_subelement( new XMLSchemaElement );
	null_perturber_subelement->name( "NullPerturber" );
	null_perturber_subelement->type_name( perturber_ct_namer( "NullPerturber" ) );

	XMLSchemaElementOP compound_perturber_subelement( new XMLSchemaElement );
	compound_perturber_subelement->name( "CompoundPerturber" );
	compound_perturber_subelement->type_name( perturber_ct_namer( "CompoundPerturber" ) );

	XMLSchemaModelGroupOP perturber_choice( new XMLSchemaModelGroup );
	perturber_choice->type( xsmgt_choice );
	perturber_choice->append_particle( helix_perturber_subelement );
	perturber_choice->append_particle( connection_perturber_subelement );
	perturber_choice->append_particle( null_perturber_subelement );
	perturber_choice->append_particle( compound_perturber_subelement );

	XMLSchemaModelGroup perturber_group;
	perturber_group.group_name( perturber_group_name() );
	perturber_group.append_particle( perturber_choice );
	xsd.add_top_level_element( perturber_group );

}

void
FoldArchitectMover::define_filters_ct( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList filters_subelements_attlist;
	filters_subelements_attlist
		+ XMLSchemaAttribute::required_attribute( "filter", xs_string, "Filter to add" );

	XMLSchemaSimpleSubelementList filters_subelements;
	filters_subelements.add_simple_subelement( "Add", filters_subelements_attlist, "Specifies a filter to apply to the backbones" );

	XMLSchemaComplexTypeGenerator filters_ct_gen;
	filters_ct_gen
		.set_subelements_repeatable( filters_subelements )
		.complex_type_naming_func( & filters_ct_naming_func )
		.element_name( "Filters" )
		.description( "Filters to apply to generated backbones" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

void FoldArchitectMover::define_folder_group( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList remodel_folder_attlist;
	rosetta_scripts::attributes_for_parse_score_function( remodel_folder_attlist );

	//Define the complex types for all of these
	XMLSchemaComplexTypeGenerator remodel_folder_ct;
	remodel_folder_ct.element_name( "RemodelLoopMoverPoseFolder" )
		.complex_type_naming_func ( & folder_ct_namer )
		.description( "Folds residues in a pose using RemodelLoopMover" )
		.add_attributes( remodel_folder_attlist )
		.write_complex_type_to_schema( xsd );


	//random and null have no attributes
	XMLSchemaComplexTypeGenerator random_folder_ct;
	random_folder_ct.element_name( "RandomTorsionPoseFolder" )
		.complex_type_naming_func ( & folder_ct_namer )
		.description( "Folds pose using random phi/psi torsions" )
		.write_complex_type_to_schema( xsd );


	XMLSchemaComplexTypeGenerator null_folder_ct;
	null_folder_ct.element_name( "NullPoseFolder" )
		.complex_type_naming_func ( & folder_ct_namer )
		.description( "Does nothing" )
		.write_complex_type_to_schema( xsd );


	XMLSchemaElementOP remodel_folder_subelement( new XMLSchemaElement );
	remodel_folder_subelement->name( "RemodelLoopMoverPoseFolder" );
	remodel_folder_subelement->type_name( folder_ct_namer( "RemodelLoopMoverPoseFolder" ) );

	XMLSchemaElementOP random_folder_subelement( new XMLSchemaElement );
	random_folder_subelement->name( "RandomTorsionPoseFolder" );
	random_folder_subelement->type_name( folder_ct_namer( "RandomTorsionPoseFolder" ) );

	XMLSchemaElementOP null_folder_subelement( new XMLSchemaElement );
	null_folder_subelement->name( "NullPoseFolder" );
	null_folder_subelement->type_name( folder_ct_namer( "NullPoseFolder" ) );

	XMLSchemaModelGroupOP folder_choice( new XMLSchemaModelGroup );
	folder_choice->type( xsmgt_choice );
	folder_choice->append_particle( remodel_folder_subelement );
	folder_choice->append_particle( random_folder_subelement );
	folder_choice->append_particle( null_folder_subelement );

	XMLSchemaModelGroup folder_group;
	folder_group.group_name( folder_group_name() );
	folder_group.append_particle( folder_choice );
	xsd.add_top_level_element( folder_group );
}

utility::tag::XMLSchemaSimpleSubelementList
pre_and_post_fold_mover_subelements()
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "mover", xs_string, "Mover to add to this step step" );

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_simple_subelement( "Add", attlist, "Specify a mover to add to this step" );

	return subelements;
}

void
FoldArchitectMover::define_prefold_movers_ct( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// prefold and postfold take the same subelements
	XMLSchemaSimpleSubelementList const prefold_subelements =
		pre_and_post_fold_mover_subelements();

	XMLSchemaComplexTypeGenerator prefold_ct_gen;
	prefold_ct_gen
		.set_subelements_repeatable( prefold_subelements )
		.complex_type_naming_func( & prefold_ct_naming_func )
		.element_name( "PreFoldMovers" )
		.description( "Prefold protocol" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

void
FoldArchitectMover::define_postfold_movers_ct( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// prefold and postfold take the same subelements
	XMLSchemaSimpleSubelementList const postfold_subelements =
		pre_and_post_fold_mover_subelements();

	XMLSchemaComplexTypeGenerator postfold_ct_gen;
	postfold_ct_gen
		.set_subelements_repeatable( postfold_subelements )
		.complex_type_naming_func( & postfold_ct_naming_func )
		.element_name( "PostFoldMovers" )
		.description( "Postfold protocol" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}


void
FoldArchitectMover::provide_xml_schema_with_name(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & mover_name )
{
	using namespace utility::tag;

	///
	/// PreFoldMovers
	///
	define_prefold_movers_ct( xsd );
	XMLSchemaSimpleSubelementList prefold_subelements;
	prefold_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_already_defined_subelement( "PreFoldMovers", & prefold_ct_naming_func );

	///
	/// PostFoldMovers
	///
	define_postfold_movers_ct( xsd );
	XMLSchemaSimpleSubelementList postfold_subelements;
	postfold_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_already_defined_subelement( "PostFoldMovers", & postfold_ct_naming_func );

	///
	/// Filters
	///
	define_filters_ct( xsd );
	XMLSchemaSimpleSubelementList filters_subelements;
	filters_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_already_defined_subelement( "Filters", & filters_ct_naming_func );

	///
	/// Architects
	///
	architects::DeNovoArchitectFactory::get_instance()->define_architect_group( xsd );
	XMLSchemaSimpleSubelementList architect_subelements;
	architect_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_group_subelement( & architects::DeNovoArchitectFactory::architect_group_name );

	///
	/// Folders
	///
	define_folder_group( xsd );
	XMLSchemaSimpleSubelementList folder_subelements;
	folder_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_group_subelement( & folder_group_name );

	///
	/// Perturbers
	///
	define_perturber_group( xsd );
	XMLSchemaSimpleSubelementList perturber_subelements;
	perturber_subelements.complex_type_naming_func( & build_denovo_backbone_ct_naming_func )
		.add_group_subelement( & perturber_group_name );


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "Unique ID for this mover")
		+ XMLSchemaAttribute( "build_overlap", xsct_non_negative_integer, "Overlap to use when building loops")
		+ XMLSchemaAttribute( "start_segments", xs_string, "Set names of segments to include in the first build phase")
		+ XMLSchemaAttribute( "stop_segments", xs_string, "Set names of segments to be included in the final build phase")
		+ XMLSchemaAttribute( "iterations_per_phase", xsct_non_negative_integer, "Number of iterations per build phase")
		+ XMLSchemaAttribute( "dump_pdbs", xsct_rosetta_bool, "Dump output to PDB files?")
		+ XMLSchemaAttribute( "debug", xsct_rosetta_bool, "Dump debugging PDBs?");

	XMLSchemaComplexTypeGenerator complex_type_generator;
	complex_type_generator
		.element_name( mover_name )
		.description( "Mover to generate backbones for de novo design" )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.add_ordered_subelement_set_as_required( architect_subelements )
		.add_ordered_subelement_set_as_required( folder_subelements )
		.add_ordered_subelement_set_as_optional( perturber_subelements )
		.add_ordered_subelement_set_as_optional( prefold_subelements )
		.add_ordered_subelement_set_as_optional( postfold_subelements )
		.add_ordered_subelement_set_as_optional( filters_subelements )
		.write_complex_type_to_schema( xsd );
}

void FoldArchitectMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	provide_xml_schema_with_name( xsd, mover_name() );
}

std::string FoldArchitectMoverCreator::keyname() const {
	return FoldArchitectMover::mover_name();
}

protocols::moves::MoverOP
FoldArchitectMoverCreator::create_mover() const {
	return utility::pointer::make_shared< FoldArchitectMover >();
}

void FoldArchitectMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FoldArchitectMover::provide_xml_schema( xsd );
}

std::string BuildDeNovoBackboneMover::get_name() const {
	return mover_name();
}

std::string BuildDeNovoBackboneMover::mover_name() {
	return "BuildDeNovoBackboneMover";
}

void BuildDeNovoBackboneMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	provide_xml_schema_with_name( xsd, mover_name() );
}

std::string BuildDeNovoBackboneMoverCreator::keyname() const {
	return BuildDeNovoBackboneMover::mover_name();
}

protocols::moves::MoverOP
BuildDeNovoBackboneMoverCreator::create_mover() const {
	return utility::pointer::make_shared< BuildDeNovoBackboneMover >();
}

void BuildDeNovoBackboneMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BuildDeNovoBackboneMover::provide_xml_schema( xsd );
}


////////////// Helper movers /////////////////
std::string
SetPoseSecstructFromStructureDataMover::get_name() const
{
	return class_name();
}

protocols::moves::MoverOP
SetPoseSecstructFromStructureDataMover::clone() const
{
	return utility::pointer::make_shared< SetPoseSecstructFromStructureDataMover >(*this);
}

void
SetPoseSecstructFromStructureDataMover::apply( core::pose::Pose & pose )
{
	components::StructureData const & sd = components::StructureDataFactory::get_instance()->get_from_const_pose( pose );
	std::string target_ss = "";

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		target_ss = symmetric_secstruct( pose, sd.ss() );
	} else {
		target_ss = sd.ss();
	}

	if ( target_ss.size() != pose.size() ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): StructureData ss size ("
			<< target_ss.size() << ") does not match pose size ("
			<< pose.size() << ") " << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size resid = 1;
	for ( std::string::const_iterator sschar=target_ss.begin(); sschar!=target_ss.end(); ++sschar, ++resid ) {
		pose.set_secstruct( resid, *sschar );
	}
	TR << "Set Pose secstruct to " << target_ss << std::endl;

}

////////////// Helper functions /////////////////

/// @brief goes through loops and adds overlapping residues
FoldArchitectMover::ResidueVector
add_overlap_to_loops(
	protocols::loops::Loops & loops,
	core::Size const overlap,
	core::pose::Pose const & pose )
{
	FoldArchitectMover::ResidueVector residues;
	if ( !overlap ) return residues;

	// go through loops and add proper overlap
	for ( auto l=loops.v_begin(); l!=loops.v_end(); ++l ) {
		for ( core::Size i=1; i<=overlap; ++i ) {
			// see if this loop is lower-terminal
			if ( core::pose::is_lower_terminus( pose, l->start() ) ) {
				break;
			}
			l->set_start( l->start()-1 );
			residues.push_back( l->start() );
		}
		for ( core::Size i=1; i<=overlap; ++i ) {
			// check if this loop is upper-terminal
			if ( core::pose::is_upper_terminus( pose, l->stop() ) ) {
				break;
			}
			l->set_stop( l->stop()+1 );
			residues.push_back( l->stop() );
		}
		TR << "LOOP is now " << l->start() << " " << l->stop() << " " << l->cut() << std::endl;
	}
	return residues;
}

} //protocols
} //denovo_design
} //movers
