// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PeakCalibratorList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_StructureIndependentPeakCalibrator_hh
#define INCLUDED_protocols_noesy_assign_StructureIndependentPeakCalibrator_hh


// Unit Header
//#include <devel/NoesyAssign/StructureIndependentPeakCalibrator.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.hh>

// Package Headers
// #include <devel/NoesyAssign/StructureIndependentPeakCalibratorInfo.hh>
// #include <devel/NoesyAssign/PeakAssignment.hh>
// #include <devel/NoesyAssign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>
//#include <core/chemical/AA.hh>

// Utility headers
//#include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
//#include <utility/vector1.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>

//// C++ headers
//#include <cstdlib>
// #include <string>
// #include <list>



namespace protocols {
namespace noesy_assign {

class StructureIndependentPeakCalibrator : public PeakCalibrator {
public:

	StructureIndependentPeakCalibrator() : PeakCalibrator( 1 /*normal sign*/ ) {};
	PeakCalibratorOP fresh_instance() override {
		return utility::pointer::make_shared< StructureIndependentPeakCalibrator >();
	}

	void collect_upperbound_statistics( core::Size  /*peak*/, TypeCumulator const& /*types*/ ) override;
private:
};


}
}

#endif
