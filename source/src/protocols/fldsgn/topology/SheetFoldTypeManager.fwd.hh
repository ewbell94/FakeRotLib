// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/SheetFoldTypeManager.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )


#ifndef INCLUDED_protocols_fldsgn_topology_SheetFoldTypeManager_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_SheetFoldTypeManager_fwd_hh
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace fldsgn {
namespace topology {

enum SheetFoldType : char;

typedef utility::vector1< SheetFoldType > SheetFoldTypes;
class SheetFoldTypeManager;

} // namespace topology
} // namespace fldsgn
} // namespace protocols


#endif
