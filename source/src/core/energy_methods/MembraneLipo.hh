// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/MembraneLipo.hh
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


#ifndef INCLUDED_core_energy_methods_MembraneLipo_hh
#define INCLUDED_core_energy_methods_MembraneLipo_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// Project headers

// Utility headers


namespace core {
namespace energy_methods {



class MembraneLipo : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:


	MembraneLipo();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;
	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;
	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:
	// const-ref to scoring database
	core::scoring::MembranePotential const & potential_;
	core::Size version() const override;
};


} // scoring
} // core

#endif
