// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/ISecretionPredictionFilterCreator.hh
/// @brief  Filter for greasy helices  https://www.nature.com/articles/nature06387
/// @author Designed by John Wang(jyjwang@uw.edu) converted to a filter by TJ Brunette(tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_SecretionPredictionFilterCreator_hh
#define INCLUDED_protocols_simple_filters_SecretionPredictionFilterCreator_hh


// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers

// c++ headers
#include <string>

namespace protocols {
namespace simple_filters {

class SecretionPredictionFilterCreator : public protocols::filters::FilterCreator
{
public:
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


} //namespace filters
} //namespace protocols

#endif
