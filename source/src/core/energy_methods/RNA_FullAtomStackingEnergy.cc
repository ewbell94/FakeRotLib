// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_FullAtomStacking.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das


// Unit headers
#include <core/energy_methods/RNA_FullAtomStackingEnergy.hh>
#include <core/energy_methods/RNA_FullAtomStackingEnergyCreator.hh>

// Package headers
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MinimizerMapBase.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a simple extension of the van der Waals attraction to give extra attraction between
//  delocalized electron clouds on aromatic residues.
//
// The attraction is full-strength for r < 4.0 Angstroms, and then cutoff smoothly at 6.0 A.
//  (full strength at neighboring base stacks; but no next-nearest base stacks in nucleic acids).
//
// There is also an orientation dependence of cos( kappa )^2, where kappa is the angle of the
//  the r vector to the base axis. This limits the calculation to atoms that are right 'above'
//  each other in the stacks.
//
// This is basically arbitrary and has not yet been compared to, e.g., quantum computations that
//  take into account electron correlations.
//
// Note that base-base contributions count twice, and by default, only base-base interactions
//  are currently computed.
//
// Some of the code below is currently RNA specific, but could easily be generalized to DNA,
//  and with a tiny bit of extra effort to general aromatics. This physics may be better
//  modeled in the orbitals code, however.
//
// Currently, intrares not supported.
// And... probably could accelerate this with the nblist & residue_pair_energy_ext() formalism.
//   (see, e.g., StackElecEnergy.cc). Seems like someone (arvind?) checked in part of this.
//
//    -- rhiju, 2014
//
// Extension to penalize desolvation of waters stacked on nucleobase ('sol_stack') and
//  longer-range water-mediated stacking ('lr_stack') added in may 2015.
//
// Parameters set to match water potential-of-mean-force around adenosine, estimated from MD (TIP3P) simulations and
//  compared to rosetta calculations with 'nucleobase_sample_around', available at:
//
//  https://drive.google.com/file/d/0By0BpYZBGuK-R2dhcEpKN0E4d28/view
//
//    -- rhiju, 2015
//
///////////////////////////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "core.energy_methods.RNA_FullAtomStackingEnergy" );
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace energy_methods {

/// @details This must return a fresh instance of the RNA_FullAtomStackingEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RNA_FullAtomStackingEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< RNA_FullAtomStackingEnergy >( options );
}

core::scoring::ScoreTypes
RNA_FullAtomStackingEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( fa_stack );
	sts.push_back( fa_stack_lower ); // just lower residue, used by free_base entropy estimator
	sts.push_back( fa_stack_upper ); // just upper residue, used by free_base entropy estimator
	sts.push_back( fa_stack_aro  );  // just component where occluding atom is aromatic.

	sts.push_back( fa_stack_rna_protein );
	sts.push_back( fa_stack_rna_protein_aro );

	sts.push_back( fa_stack_ext );
	sts.push_back( fa_stack_sol );
	sts.push_back( fa_stack_lr  );
	return sts;
}


/// c-tor
RNA_FullAtomStackingEnergy::RNA_FullAtomStackingEnergy( core::scoring::methods::EnergyMethodOptions const & options ) :
	parent( utility::pointer::make_shared< RNA_FullAtomStackingEnergyCreator >() ),
	//Parameters are totally arbitrary and made up!
	prefactor_    ( -0.2 ),
	stack_cutoff_ ( 4.0 ), //Encompass next base stack
	dist_cutoff_  ( 6.0 ), //  (Do not go to next-nearest base stack!)
	/////////////////////////////////////////////////
	sol_prefactor_      ( option[ fa_stack_sol_prefactor ]()    /*+0.1*/ ), // Penalize displacement of water that was stacked
	sol_stack_cutoff_   ( option[ fa_stack_sol_stack_cutoff ]() /* 5.5*/ ), //
	sol_dist_cutoff_    ( option[ fa_stack_sol_dist_cutoff ]()  /* 6.5*/ ), // how far to penalize displacement
	/////////////////////////////////////////////////
	lr_prefactor_       ( option[ fa_stack_lr_prefactor ]()    /*-0.05*/ ), // water-mediated stacking
	lr_stack_cutoff_    ( option[ fa_stack_lr_stack_cutoff ]() /*6.5*/   ), //
	lr_dist_cutoff_     ( option[ fa_stack_lr_dist_cutoff ]()  /*7.5*/  ),//
	base_base_only_( !options.fa_stack_base_all() /* !false --> true */ ),
	include_rna_prot_( option[ fa_stack_include_rna_protein ]() /*true */ ) // for fa_stack_rna_protein(_aro)
{}

/// copy c-tor
RNA_FullAtomStackingEnergy::RNA_FullAtomStackingEnergy( RNA_FullAtomStackingEnergy const & /*src*/ ) = default;

//clone
core::scoring::methods::EnergyMethodOP
RNA_FullAtomStackingEnergy::clone() const
{
	return utility::pointer::make_shared< RNA_FullAtomStackingEnergy >( *this );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::setup_for_minimizing(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const {

	if ( !pose.energies().use_nblist() ) return;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// I think jyesselm put this in, but use of nblist is not activated -- need to define residue_pair_energy_ext(),
	//   eval_residue_pair_derivatives(), use_extended_residue_pair_energy_interface(), etc.
	//
	// See, e.g., src/core/scoring/geom_sol/ContextDependentGeometricSolEnergy or src/core/scoring/magnesium/MgEnergy.cc
	//
	// Putting this in correctly could significantly accelerate RNA and RNA/protein code if we continue
	//  use of fa_stack, fa_stack_sol, or fa_stack_lr.
	//
	// -- rhiju.
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	// stash our nblist inside the pose's energies object
	core::scoring::Energies & energies( pose.energies() );

	// untested -- rhiju
	Distance const dist_cutoff = sfxn.has_nonzero_weight( core::scoring::fa_stack_lr ) ? lr_dist_cutoff_ : ( sfxn.has_nonzero_weight( core::scoring::fa_stack_sol ) ? sol_dist_cutoff_ : dist_cutoff_ );

	// setup the atom-atom nblist
	core::scoring::NeighborListOP nblist;
	Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? basic::options::option[ basic::options::OptionKeys::run::nblist_autoupdate_narrow ] : 1.5;
	Real const XX = dist_cutoff + 2 * tolerated_motion;
	nblist = utility::pointer::make_shared< core::scoring::NeighborList >( min_map.domain_map(), XX*XX, XX*XX, XX*XX );
	if ( pose.energies().use_nblist_auto_update() ) {
		nblist->set_auto_update( tolerated_motion );
	}
	nblist->setup( pose, sfxn, *this );
	energies.set_nblist( core::scoring::EnergiesCacheableDataType::FA_STACK_NBLIST, nblist );
}


void
RNA_FullAtomStackingEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & scfxn ) const
{
	pose.update_residue_neighbors();

	core::scoring::rna::RNA_ScoringInfo  & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );

	pose.update_residue_neighbors();
	if ( pose.energies().use_nblist() ) {
		core::scoring::NeighborList const & nblist( pose.energies().nblist( core::scoring::EnergiesCacheableDataType::FA_STACK_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

void
RNA_FullAtomStackingEnergy::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	pose.update_residue_neighbors();

	core::scoring::rna::RNA_ScoringInfo  & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::EnergyMap & emap
) const {
	// This is ok because this interaction is defined as only attractive.'
	// (before I got here, we just 'continue'd for repulsive atoms!)
	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd2.has_variant_type( REPLONLY ) ) return;

	Real score_aro1( 0.0 ), score_aro2( 0.0 );
	Real const score1 = residue_pair_energy_one_way( rsd1, rsd2, pose, score_aro1, prefactor_, stack_cutoff_, dist_cutoff_ );
	Real const score2 = residue_pair_energy_one_way( rsd2, rsd1, pose, score_aro2, prefactor_, stack_cutoff_, dist_cutoff_ ) ;

	if ( (base_base_only_ && rsd1.is_RNA() && rsd2.is_RNA()) || !base_base_only_ ) {
		emap[ core::scoring::fa_stack ]       += score1 + score2;
		emap[ core::scoring::fa_stack_ext ]   += score1 + score2;
		emap[ core::scoring::fa_stack_aro ]   += score_aro1 + score_aro2;
	}

	if ( include_rna_prot_ ) {
		if ( (rsd1.is_RNA() && rsd2.is_protein()) || (rsd1.is_protein() && rsd2.is_RNA()) ) {
			emap[ core::scoring::fa_stack_rna_protein ] += score1 + score2; // one of these will be 0
			emap[ core::scoring::fa_stack_rna_protein_aro ] += score_aro1 + score_aro2; // one of these will be 0
		}
	}

	if ( rsd1.seqpos() <= rsd2.seqpos() ) {
		if ( (base_base_only_ && rsd1.is_RNA() && rsd2.is_RNA()) || !base_base_only_ ) {
			emap[ core::scoring::fa_stack_lower ]  += score1;
			emap[ core::scoring::fa_stack_upper ]  += score2;
		}
	} else {
		if ( (base_base_only_ && rsd1.is_RNA() && rsd2.is_RNA()) || !base_base_only_ ) {
			emap[ core::scoring::fa_stack_lower ]  += score2;
			emap[ core::scoring::fa_stack_upper ]  += score1;
		}
	}

	if ( sfxn.get_weight( core::scoring::fa_stack_sol ) != 0 || sfxn.get_weight( core::scoring::fa_stack_ext ) != 0 ) {
		Real const sol_score1 = residue_pair_energy_one_way( rsd1, rsd2, pose, score_aro1, sol_prefactor_, sol_stack_cutoff_, sol_dist_cutoff_);
		Real const sol_score2 = residue_pair_energy_one_way( rsd2, rsd1, pose, score_aro2, sol_prefactor_, sol_stack_cutoff_, sol_dist_cutoff_ ) ;
		emap[ core::scoring::fa_stack_sol ] += sol_score1 + sol_score2;
		emap[ core::scoring::fa_stack_ext ] += sol_score1 + sol_score2;
	}

	if ( sfxn.get_weight( core::scoring::fa_stack_lr ) != 0 || sfxn.get_weight( core::scoring::fa_stack_ext ) != 0 ) {
		Real const lr_score1 = residue_pair_energy_one_way( rsd1, rsd2, pose, score_aro1, lr_prefactor_, lr_stack_cutoff_, lr_dist_cutoff_);
		Real const lr_score2 = residue_pair_energy_one_way( rsd2, rsd1, pose, score_aro2, lr_prefactor_, lr_stack_cutoff_, lr_dist_cutoff_ ) ;
		emap[ core::scoring::fa_stack_lr  ] += lr_score1 + lr_score2;
		emap[ core::scoring::fa_stack_ext ] += lr_score1 + lr_score2;
	}

	// TR << rsd1.name3()  << rsd1.seqpos() << "---" << rsd2.name3() << rsd2.seqpos() << ": " << (score1+score2) << std::endl;
	debug_assert( emap.is_finite() );
}


//////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_FullAtomStackingEnergy::residue_pair_energy_one_way(
	conformation::Residue const & rsd1, // residue with nucleobase
	conformation::Residue const & rsd2, // occluding residue
	pose::Pose const & pose,
	Real & score_aro,
	Real     const & prefactor,
	Distance const & stack_cutoff,
	Distance const & dist_cutoff
) const {
	if ( rsd1.has_variant_type( REPLONLY ) ) return 0.0;
	if ( rsd2.has_variant_type( REPLONLY ) ) return 0.0;

	score_aro = 0.0;
	Real score( 0.0 );
	if ( !rsd1.is_RNA() && !rsd1.is_TNA() ) return score;
	if ( base_base_only_ ) {
		if ( !rsd2.is_protein() && !rsd2.is_RNA() && !rsd2.is_TNA() ) return score;
		if ( rsd2.is_protein() && !include_rna_prot_ ) return score;
	}

	core::scoring::rna::RNA_ScoringInfo  const & rna_scoring_info( core::scoring::rna::rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	// utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	Size const i( rsd1.seqpos() );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const M_i ( stub_i.M );

	// Loop over base heavy atoms.
	// If I want to generalize this to proteins, maybe could loop over "aromatic" atoms.
	Real const dist_cutoff2( dist_cutoff * dist_cutoff );
	for ( Size m = rsd1.first_sidechain_atom(); m <= rsd1.nheavyatoms(); ++m ) {
		//Need to be careful nheavyatoms count includes hydrogen! when the hydrogen is made virtual..

		if ( rsd1.is_virtual( m ) ) continue;
		if ( rsd1.is_repulsive( m ) ) continue;
		if ( m == rsd1.first_sidechain_atom() ) continue; // 2'-OH for RNA

		Vector const heavy_atom_i( rsd1.xyz( m ) );

		//Look for occlusion by other base atoms? Or all other heavy atoms?
		Size const atom_num_start = base_base_only_ ? rsd2.first_sidechain_atom() : 1 ;

		for ( Size n = atom_num_start; n <= rsd2.nheavyatoms(); ++n ) {

			if ( rsd2.is_virtual( n ) ) continue;
			if ( rsd2.is_repulsive( n ) ) continue;
			if ( base_base_only_ && n == rsd2.first_sidechain_atom() ) continue; // 2'-OH for RNA

			runtime_assert( check_base_base_OK( rsd1, rsd2, m, n ) );

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector r = heavy_atom_j - heavy_atom_i;
			Real const dist2 = r.length_squared();

			if ( dist2 >= dist_cutoff2 ) continue;
			Real const fa_stack_score = get_fa_stack_score( r, M_i, prefactor, stack_cutoff, dist_cutoff );
			score += fa_stack_score;

			if ( rsd1.atom_type( m ).is_aromatic() && rsd2.atom_type( n ).is_aromatic() ) score_aro += fa_stack_score;
		}
	}

	return score;
}

//////////////////////////////////////////////////////////////////////////////
// Depending on the option "base_base_only_" allow a stacking interaction
//  between an atom in a nucleobase and the nucleobase atoms or *any*
//  other atom.
//
bool
RNA_FullAtomStackingEnergy::check_base_base_OK(
	conformation::Residue const & rsd1, // harbors nucleobase
	conformation::Residue const & rsd2,
	Size const m, Size const n ) const {

	// must be in RNA base.
	if ( m < rsd1.first_sidechain_atom()  || m > rsd1.nheavyatoms() ) return false;
	if ( m == rsd1.first_sidechain_atom() ) return false;

	if ( n > rsd2.nheavyatoms() ) return false;
	if ( base_base_only_ ) {
		if ( !rsd2.is_protein() && !rsd2.is_RNA() && !rsd2.is_TNA() ) return false;
		if ( rsd2.is_protein() && !include_rna_prot_ ) return false;
		if ( rsd2.is_RNA() || rsd2.is_TNA() ) {
			if ( n < rsd2.first_sidechain_atom() ) return false;
			if ( n == rsd2.first_sidechain_atom() ) return false;  //2'-OH
		}
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	if ( pose.residue_type( atom_id.rsd() ).has_variant_type( REPLONLY ) ) return;

	// following repeats some calculations but is conceptually simplest.
	eval_atom_derivative(
		atom_id, pose, domain_map, F1, F2,
		weights[ core::scoring::fa_stack ] + weights[ core::scoring::fa_stack_ext ] + weights[ core::scoring::fa_stack_rna_protein ],
		//weights[ core::scoring::fa_stack ] + weights[ core::scoring::fa_stack_ext ],
		weights[ core::scoring::fa_stack_lower ],
		weights[ core::scoring::fa_stack_upper ],
		weights[ core::scoring::fa_stack_aro ] + weights[ core::scoring::fa_stack_rna_protein_aro ],
		//weights[ core::scoring::fa_stack_aro ] /* aro -- legacy */,
		prefactor_, stack_cutoff_, dist_cutoff_ );

	if ( sfxn.get_weight( core::scoring::fa_stack_sol ) != 0 || sfxn.get_weight( core::scoring::fa_stack_ext ) != 0 ) {
		eval_atom_derivative(
			atom_id, pose, domain_map, F1, F2,
			weights[ core::scoring::fa_stack_sol ] + weights[ core::scoring::fa_stack_ext ],
			0.0, 0.0, 0.0,
			sol_prefactor_, sol_stack_cutoff_, sol_dist_cutoff_ );
	}

	if ( sfxn.get_weight( core::scoring::fa_stack_lr ) != 0 || sfxn.get_weight( core::scoring::fa_stack_ext ) != 0 ) {
		eval_atom_derivative(
			atom_id, pose, domain_map, F1, F2,
			weights[ core::scoring::fa_stack_lr ] + weights[ core::scoring::fa_stack_ext ],
			0.0, 0.0, 0.0,
			lr_prefactor_, lr_stack_cutoff_, lr_dist_cutoff_ );
	}
}

//////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	Vector & F1,
	Vector & F2,
	Real const fa_stack_weight,
	Real const fa_stack_lower_weight,
	Real const fa_stack_upper_weight,
	Real const fa_stack_aro_weight,
	Real const prefactor,
	Distance const stack_cutoff,
	Distance const dist_cutoff
) const {
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( rsd1.has_variant_type( REPLONLY ) ) return;

	if ( base_base_only_ && !rsd1.is_RNA() && !rsd1.is_TNA() && !include_rna_prot_ ) return;
	if ( base_base_only_ && include_rna_prot_ && !(rsd1.is_protein() || rsd1.is_RNA() || rsd1.is_TNA() ) ) return;
	//if ( base_base_only_ && !rsd1.is_RNA() ) return;
	if ( m > rsd1.nheavyatoms() ) return;
	if ( rsd1.is_virtual( m ) ) return;
	if ( rsd1.is_repulsive( m ) ) return;

	core::scoring::rna::RNA_ScoringInfo  const & rna_scoring_info( core::scoring::rna::rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );
	kinematics::Stub const & stub_i( base_stubs[i] );

	Matrix const M_i ( stub_i.M );
	Vector const heavy_atom_i( rsd1.xyz( m ) );

	bool const pos1_fixed( domain_map( i ) != 0 );

	// cached energies object
	core::scoring::Energies const & energies( pose.energies() );

	// the neighbor/energy links
	core::scoring::EnergyGraph const & energy_graph( energies.energy_graph() );

	Real const dist_cutoff2( dist_cutoff * dist_cutoff );

	for ( utility::graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( ( *iter )->get_other_ind( i ) );

		if ( pos1_fixed && domain_map( i ) == domain_map( j ) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );
		//if ( base_base_only_ && !rsd2.is_RNA() ) continue;
		if ( base_base_only_ ) {
			if ( !rsd1.is_RNA() || !rsd2.is_RNA() || rsd2.is_TNA() ) {
				if ( !include_rna_prot_ ) continue;
				if ( !(( ( rsd1.is_RNA() || rsd1.is_TNA() ) && rsd2.is_protein()) || (rsd1.is_protein() && ( rsd2.is_RNA() || rsd2.is_TNA() ) )) ) continue;
			}
		}

		kinematics::Stub const & stub_j( base_stubs[j] );
		Matrix const M_j ( stub_j.M );

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			if ( rsd2.is_virtual( n ) ) continue;
			if ( rsd2.is_repulsive( n ) ) continue;

			//   if ( !check_base_base_OK( rsd1, rsd2, m, n ) && !check_base_base_OK( rsd2, rsd1, n, m ) ) continue;

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector r = heavy_atom_j - heavy_atom_i;
			Real const dist2 = r.length_squared();

			///////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////
			// Must update following to handle fa_stack_sol & fa_stack_ext derivatives
			///////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////

			if ( dist2 >= dist_cutoff2 ) continue;

			if ( check_base_base_OK( rsd1, rsd2, m, n ) ) {

				Vector const deriv = get_fa_stack_deriv( r, M_i, prefactor, stack_cutoff, dist_cutoff );

				Vector force_vector_i = fa_stack_weight * deriv;
				force_vector_i += fa_stack_lower_weight * deriv;
				if ( rsd1.atom_type( m ).is_aromatic() && rsd2.atom_type( n ).is_aromatic() ) force_vector_i += fa_stack_aro_weight * deriv;

				//Force/torque with which occluding atom j acts on "dipole" i.
				F1 += -1.0 * cross( force_vector_i, heavy_atom_j );
				F2 += -1.0 * force_vector_i;
			}

			if ( check_base_base_OK( rsd2, rsd1, n, m ) ) {
				//Force/torque with which occluding atom i acts on "dipole" j.
				// Note that this calculation is a repeat of some other calls to this function.
				// Might make more sense to do some (alternative) bookkeeping, e.g. with _ext interface.

				Vector const deriv = get_fa_stack_deriv( -r, M_j, prefactor, stack_cutoff, dist_cutoff );

				Vector force_vector_j = fa_stack_weight * deriv;
				force_vector_j += fa_stack_upper_weight * deriv;
				if ( rsd1.atom_type( m ).is_aromatic() && rsd2.atom_type( n ).is_aromatic() ) force_vector_j += fa_stack_aro_weight * deriv;

				F1 += cross( force_vector_j, heavy_atom_i );
				F2 += force_vector_j;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
//get_fa_stack_score evaluates the fa_stack_score between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by residue_pair_energy_one_way
Real
RNA_FullAtomStackingEnergy::get_fa_stack_score(
	Vector const & r_vec,
	Matrix const & M_i,
	Real const prefactor,
	Distance const stack_cutoff,
	Distance const dist_cutoff
) const {
	Vector const z_i = M_i.col_z();

	Real const r = r_vec.length();
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	Real score  = prefactor;

	//Orientation dependence
	score *= cos_kappa * cos_kappa ;

	Distance b = r - stack_cutoff;

	//Just use a simple cubic spline to fade from 1 to 0,
	// and to force derivatives to be continuous.
	if ( b > 0.0 ) {
		b /= ( dist_cutoff - stack_cutoff );
		Real const b2 = b*b;
		Real const b3 = b2*b;
		score *= ( 2 * b3 - 3 * b2 + 1 );
	}

	return score;
}


//////////////////////////////////////////////////////////////////////////////////////////
//get_fa_stack_deriv evaluates the fa_stack_deriv between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by eval_atom_derivative
Vector
RNA_FullAtomStackingEnergy::get_fa_stack_deriv(
	Vector const & r_vec,
	Matrix const & M_i,
	Real const prefactor,
	Distance const stack_cutoff,
	Distance const dist_cutoff
) const {
	//  Well, the energy is a function of the inter-atom distance and a special cos(angle)
	//            E = E ( r, cos(theta ) ).
	//  so
	//      dE/dx = dE/dr (x/r)   +   ( - x * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dy = dE/dr (y/r)   +   ( - y * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dz = dE/dr (z/r)   +   ( (r^2 - z^2) / r^3 ) ( dE/dcos(theta) )

	Vector const x_i = M_i.col_x();
	Vector const y_i = M_i.col_y();
	Vector const z_i = M_i.col_z();

	Real const r = r_vec.length();
	Real const x = dot( x_i, r_vec );
	Real const y = dot( y_i, r_vec );
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	/////////////////////////////////
	//dE_dcoskappa
	/////////////////////////////////
	Real dE_dcoskappa  = prefactor;

	//Orientation dependence
	dE_dcoskappa *= 2 * cos_kappa ;

	Distance b = r - stack_cutoff;
	Distance distance_scale = ( dist_cutoff - stack_cutoff );
	b /= distance_scale;

	Real b2( 0.0 ), b3( 0.0 );
	//Just use a simple cubic spline to fade from 1 to 0,
	// and to force derivatives to be continuous.
	if ( b > 0.0 && b < 1.0 ) {
		b2 = b*b;
		b3 = b2*b;
		dE_dcoskappa *= ( 2 * b3 - 3 * b2 + 1 );
	}

	/////////////////////////////////
	//dE_dr
	/////////////////////////////////
	Real dE_dr  = prefactor;

	//Orientation dependence
	dE_dr *= cos_kappa * cos_kappa ;

	//Just use a simple cubic spline to fade from 1 to 0,
	// and to force derivatives to be continuous.
	if ( b > 0.0 && b < 1.0 ) {
		dE_dr *= ( 6 * b2 - 6 * b );
		dE_dr /= distance_scale;
	} else {
		dE_dr = 0.0;
	}

	Real const dE_dx = ( dE_dr ) * ( x/r ) - ( dE_dcoskappa ) * ( x * z )/ ( r*r*r ) ;
	Real const dE_dy = ( dE_dr ) * ( y/r ) - ( dE_dcoskappa ) * ( y * z )/ ( r*r*r ) ;
	Real const dE_dz = ( dE_dr ) * ( z/r ) + ( dE_dcoskappa ) * ( x*x + y*y )/ ( r*r*r );

	return ( dE_dx * x_i + dE_dy * y_i + dE_dz * z_i );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const {

	core::scoring::rna::RNA_ScoringInfo  & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.calculated() = false;

}

/// @brief RNA_FullAtomStackingEnergy distance cutoff
Distance
RNA_FullAtomStackingEnergy::atomic_interaction_cutoff() const
{
	//  return dist_cutoff_ + 2 * chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH; //Similar to FA_ElecEnergy
	return lr_dist_cutoff_ + 2 * chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH; //Similar to FA_ElecEnergy
}

core::Size
RNA_FullAtomStackingEnergy::version() const
{
	return 1; // Initial versioning
}

core::scoring::etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_intrares_countpair(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &
) const
{
	utility_exit_with_message( "FA_ElecEnergy does not define intra - residue pair energies; do not call get_intrares_countpair()" );
	return nullptr;
}

core::scoring::etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &
) const
{
	using namespace core::scoring::etable::count_pair;
	if ( res1 == res2 ) {
		return utility::pointer::make_shared< CountPairNone >();
	}
	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}


core::scoring::etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace core::scoring::etable::count_pair;

	if ( !rsd1.is_RNA() && !rsd1.is_TNA() ) return utility::pointer::make_shared< CountPairNone >();
	if ( !rsd2.is_RNA() && !rsd2.is_TNA() ) return utility::pointer::make_shared< CountPairNone >();
	if ( rsd1.seqpos() == rsd2.seqpos() ) return utility::pointer::make_shared< CountPairNone >();
	return utility::pointer::make_shared< CountPairAll >();
}


} //scoring
} //core
