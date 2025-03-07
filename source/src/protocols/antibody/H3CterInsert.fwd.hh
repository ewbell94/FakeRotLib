// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/H3CterInsert.fwd.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_H3CterInsert_fwd_hh
#define INCLUDED_protocols_antibody_H3CterInsert_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace antibody {

class H3CterInsert;
typedef utility::pointer::shared_ptr< H3CterInsert > H3CterInsertOP;
typedef utility::pointer::shared_ptr< H3CterInsert const > H3CterInsertCOP;

}
}


#endif


