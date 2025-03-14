// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/GapCloser.fwd.hh
/// @brief GapCloser closes the gaps after moving the free peptide.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_fwd_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

class GapCloser;

typedef utility::pointer::shared_ptr< GapCloser > GapCloserOP;
typedef utility::pointer::shared_ptr< GapCloser const > GapCloserCOP;



} //protocols
} //backbone_moves
} //local_backbone_mover


#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_fwd_hh





