// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.fwd.hh
/// @brief  forward declaration for InterfaceHydrophobicResidueContactsFilter
/// @author  Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_protocols_simple_filters_InterfaceHydrophobicResidueContactsFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_InterfaceHydrophobicResidueContactsFilter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_filters {
// Forward
class InterfaceHydrophobicResidueContactsFilter;


typedef utility::pointer::shared_ptr< InterfaceHydrophobicResidueContactsFilter > InterfaceHydrophobicResidueContactsFilterOP;
typedef utility::pointer::shared_ptr< InterfaceHydrophobicResidueContactsFilter const > InterfaceHydrophobicResidueContactsFilterCOP;

} // namespace simple_filters
} // namespace protocols

#endif
