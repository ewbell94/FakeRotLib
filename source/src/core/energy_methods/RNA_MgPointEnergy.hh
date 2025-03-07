// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RNA_MgPointEnergy.hh
/// @brief  Statistically derived Mg(2+) binding potential for RNA.
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_MgPointEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_MgPointEnergy_hh

// Unit Headers
#include <core/scoring/magnesium/MgKnowledgeBasedPotential.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace energy_methods {


class RNA_MgPointEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;

public:


	RNA_MgPointEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override {}


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

private:

	Size
	get_vdw_atom_number(
		utility::vector1< utility::vector1< Size > > const & atom_numbers_for_vdw_calculation,
		Size const & pos1,
		Size const & i ) const;

	Size
	get_vdw_atom_number(
		char const which_nucleotide,
		Size const & i ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	//Helper functions just to get things set up.
	void
	setup_info_for_mg_calculation( pose::Pose & pose ) const;

	void
	residue_pair_energy_one_way(
		conformation::Residue const & rsd1, // The RNA residue
		conformation::Residue const & rsd2, // The Mg(2+)
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	core::Size version() const override;

	core::scoring::magnesium::MgKnowledgeBasedPotentialOP rna_mg_knowledge_based_potential_;

	bool const verbose_;

};


} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
