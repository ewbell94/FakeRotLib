// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/ExtractSubposeMover.hh
/// @brief  Extract primary component associated with symdofs and all neighboring components.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_ExtractSubposeMover_hh
#define INCLUDED_protocols_matdes_ExtractSubposeMover_hh

// Unit Headers
#include <protocols/matdes/ExtractSubposeMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/types.hh>

// Utility Headers

// C++ Headers

namespace protocols {
namespace matdes {

class ExtractSubposeMover : public protocols::moves::Mover {
public:

	// default constructor
	ExtractSubposeMover();

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	std::string sym_dof_names_;
	std::string prefix_;
	std::string suffix_;
	core::Real contact_dist_;
	bool extras_;
};

} //namespace matdes
} //namespace protocols

#endif
