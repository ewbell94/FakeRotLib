// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/simple_metrics/composite_metrics/ElectrostaticSimilarityMetric.fwd.hh
/// @brief   Evalute electrostatic similarity of structures
/// @author  Andreas Scheck (andreas.scheck@epfl.ch)


#ifndef INCLUDED_core_simple_metrics_ElectrostaticSimilarityMetric_fwd_hh
#define INCLUDED_core_simple_metrics_ElectrostaticSimilarityMetric_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace simple_metrics {
namespace composite_metrics {

// Forward
class ElectrostaticSimilarityMetric;

// Types
typedef utility::pointer::shared_ptr< ElectrostaticSimilarityMetric >  ElectrostaticSimilarityMetricOP;
typedef utility::pointer::shared_ptr< ElectrostaticSimilarityMetric const >  ElectrostaticSimilarityMetricCOP;

} // namespace core
} // namespace simple_metrics
} // namespace composite_metrics

#endif
