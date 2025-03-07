// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/StackElecEnergy.hh
/// @brief  FA_Elec, but just 'perpendicular' to bases... trying to separate from H-bonds & geom_sol.
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_StackElecEnergy_HH
#define INCLUDED_core_scoring_rna_StackElecEnergy_HH

// Unit Headers
#include <core/energy_methods/StackElecEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh> // AUTO IWYU For CountPairFunctionCOP


namespace core {
namespace energy_methods {

typedef  numeric::xyzMatrix< Real > Matrix;


class StackElecEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;

public:


	StackElecEnergy( core::scoring::methods::EnergyMethodOptions const & options );


	StackElecEnergy( StackElecEnergy const & src );

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

	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const override;

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
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// totals
	) const override;

	void
	setup_for_minimizing(
		pose::Pose & pose,
		core::scoring::ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const override;

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool res_moving_wrt_eachother
	) const override;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const res1,
		Size const res2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &
	) const;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &
	) const;

	bool
	use_extended_residue_pair_energy_interface() const override;

	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData & pair_data
	) const override;

	void
	setup_for_packing(
		pose::Pose  & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const & designing_residues
	) const override;


	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

protected:

	inline
	core::scoring::etable::coulomb::Coulomb const &
	coulomb() const { return coulomb_; }

private:

	// Real
	// get_fa_stack_score( Distance const dist, Real const cos_kappa ) const;

	Real
	get_stack_elec_score( Vector const & r_i,
		Vector const & r_j,
		Real const & i_charge,
		Real const & j_charge,
		Matrix const & M_i,
		Real & cos_kappa2 ) const;


	Vector
	get_stack_elec_deriv( Vector const & r_i,
		Vector const & r_j,
		Real const & i_charge,
		Real const & j_charge,
		Matrix const & M_i ) const;

	Real
	residue_pair_energy_one_way(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		Real & score_base_base,
		Real & score_base_bb
	) const;

	bool
	check_base_base_OK(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Size const & m, Size const & n ) const;

	bool
	is_rna_base(
		conformation::Residue const & rsd1,
		Size const & m ) const;

	core::Size version() const override;

	core::scoring::etable::coulomb::Coulomb coulomb_;

	bool const base_base_only_;

	bool const verbose_;

	mutable bool might_be_designing_;

	mutable bool using_extended_method_;

};


} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
