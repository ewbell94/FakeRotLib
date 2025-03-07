// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/LRmsdFilter.hh>
#include <protocols/protein_interface_design/filters/LRmsdFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

LRmsdFilter::LRmsdFilter() :
	protocols::filters::Filter( "LRmsd" ),
	threshold_( 5.0 ),
	reference_pose_( /* NULL */ )
{}

LRmsdFilter::LRmsdFilter(protocols::docking::DockJumps const movable_jumps,
	core::Real const threshold,
	core::pose::PoseOP reference_pose)
: protocols::filters::Filter( "LRmsd" ),
	threshold_(threshold),
	reference_pose_(std::move(reference_pose)),
	movable_jumps_(movable_jumps)
{}

LRmsdFilter::~LRmsdFilter() = default;

protocols::filters::FilterOP
LRmsdFilter::clone() const {
	return utility::pointer::make_shared< LRmsdFilter >( *this );
}

static basic::Tracer TR( "protocols.protein_interface_design.filters.LRmsdFilter" );
core::Real
LRmsdFilter::compute( core::pose::Pose const & pose ) const
{
	return protocols::docking::calc_Lrmsd(pose, *reference_pose_, movable_jumps_);
}

bool
LRmsdFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const rmsd( compute( pose ));
	TR << "L_rmsd: " << rmsd ;
	if ( rmsd <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else TR<<" failing." << std::endl;
	return( false );
}

void
LRmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd <<'\n';
}

core::Real
LRmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
LRmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map )
{
	reference_pose_ = protocols::rosetta_scripts::legacy_saved_pose_or_input( tag, data_map, class_name() );

	threshold_ = tag->getOption<core::Real>( "threshold", 5.0 );

	//TODO: support multiple jumps
	auto jump_num = tag->getOption<core::Size>( "jump", 1);

	//TODO: convert jump_num to movable_jumps_ (vector0?)
	movable_jumps_.push_back(jump_num);
}



std::string LRmsdFilter::name() const {
	return class_name();
}

std::string LRmsdFilter::class_name() {
	return "LRmsd";
}

void LRmsdFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_saved_reference_pose( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_non_negative_integer, "Holes score threshold above which we fail the filter", "200" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_non_negative_integer, "Jump across which to evaluate the holes score, numbered sequentially from 1", "1" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string LRmsdFilterCreator::keyname() const {
	return LRmsdFilter::class_name();
}

protocols::filters::FilterOP
LRmsdFilterCreator::create_filter() const {
	return utility::pointer::make_shared< LRmsdFilter >();
}

void LRmsdFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LRmsdFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // devel


