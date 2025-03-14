// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/ChemicalShiftAnisotropyEnergy.cc
/// @brief  CSA energy - Orientation dependent chemical shift
/// @author Lei Shi


//Unit headers
#include <core/energy_methods/ChemicalShiftAnisotropyEnergy.hh>
#include <core/energy_methods/ChemicalShiftAnisotropyEnergyCreator.hh>
#include <core/scoring/ChemicalShiftAnisotropy.hh>
#include <core/scoring/ChemicalShiftAnisotropy.fwd.hh>
#include <core/scoring/ScoreType.hh>
//Package headers

#include <core/conformation/Residue.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/datacache/CacheableDataType.hh>

//numeric headers
#include <numeric/xyzVector.hh>


//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

//Objexx headers


//#include <basic/options/keys/csa.OptionKeys.gen.hh>

//C++ headers

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>


static basic::Tracer tr( "core.energy_methods.ChemicalShiftAnisotropyEnergy" );

namespace core {
namespace energy_methods {


using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the ChemicalShiftAnisotropyEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
ChemicalShiftAnisotropyEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< ChemicalShiftAnisotropyEnergy >();
}

core::scoring::ScoreTypes
ChemicalShiftAnisotropyEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( csa );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ChemicalShiftAnisotropyEnergy::ChemicalShiftAnisotropyEnergy() :
	parent( utility::pointer::make_shared< ChemicalShiftAnisotropyEnergyCreator >() )
{}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
core::scoring::methods::EnergyMethodOP
ChemicalShiftAnisotropyEnergy::clone() const
{
	return utility::pointer::make_shared< ChemicalShiftAnisotropyEnergy >();
}

void ChemicalShiftAnisotropyEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &
) const
{
	csa_score_ = eval_csa( pose );
}

void ChemicalShiftAnisotropyEnergy::finalize_total_energy(
	pose::Pose &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const
{
	totals[ core::scoring::csa ] = csa_score_;
}

void ChemicalShiftAnisotropyEnergy::setup_for_minimizing(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	kinematics::MinimizerMapBase const &
) const
{
	using namespace core::scoring;
	ChemicalShiftAnisotropy const& csa_data( * retrieve_CSA_from_pose( pose ) );
	ChemicalShiftAnisotropy::CSA_lines const& All_CSA_lines( csa_data.get_CSA_data() );
	ChemicalShiftAnisotropy::CSA_lines::const_iterator it;
	Size ct = 0;
	for ( it = All_CSA_lines.begin(); it != All_CSA_lines.end(); ++it ) {
		id::AtomID atom1( pose.residue(it->res1()).atom_index(it->atom1()), it->res1());
		id::AtomID atom2( pose.residue(it->res2()).atom_index(it->atom2()), it->res2());
		id::AtomID atom3( pose.residue(it->res3()).atom_index(it->atom3()), it->res3());
		//tr.Trace << "method: it->res1(): " << it->res1() << " it->atom1() " << it->atom1() << std::endl;
		//tr.Trace << "method: it->res2(): " << it->res2() << " it->atom2() " << it->atom2() << std::endl;
		//tr.Trace << "method: it->res3(): " << it->res3() << " it->atom3() " << it->atom3() << std::endl;
		//tr.Trace << "insert in atom-map " << atom1 << " " << atom2 << " " << atom3 << std::endl;
		++ct;
		utility::vector1< core::Size > atm1_map = atom2csa_map_.get( atom1 );
		utility::vector1< core::Size > atm2_map = atom2csa_map_.get( atom2 );
		utility::vector1< core::Size > atm3_map = atom2csa_map_.get( atom3 );
		atm1_map.push_back( ct );
		atm2_map.push_back( ct );
		atm3_map.push_back( ct );
		atom2csa_map_.set( atom1, atm1_map );
		atom2csa_map_.set( atom2, atm2_map );
		atom2csa_map_.set( atom3, atm3_map );
	}
}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
core::scoring::ChemicalShiftAnisotropy &
ChemicalShiftAnisotropyEnergy::csa_from_pose(
	pose::Pose & pose
) const
{
	using namespace core::scoring;
	ChemicalShiftAnisotropyOP csa_info( retrieve_CSA_from_pose( pose ) );
	if ( !csa_info ) {
		csa_info = utility::pointer::make_shared< ChemicalShiftAnisotropy >();
		store_CSA_in_pose( csa_info, pose );
	}
	return *csa_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for CSA energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real ChemicalShiftAnisotropyEnergy::eval_csa(
	pose::Pose & pose
) const
{

	core::scoring::ChemicalShiftAnisotropy& csa_data( csa_from_pose( pose ) );
	Real score = csa_data.compute_csascore( pose );
	return score;
}

void
ChemicalShiftAnisotropyEnergy::eval_atom_derivative(
	id::AtomID const & aid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & score_weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace core::scoring;

	if ( !atom2csa_map_.has( aid ) ) return; //damn this "has" isn't correct at all
	utility::vector1< Size > const csa_nrs( atom2csa_map_[ aid ] );
	//tr.Trace << " aid " << aid << std::endl;

	if ( csa_nrs.size() == 0 ) {
		//  tr.Trace << "no CSA entry for " << aid << " skipping.. "<< std::endl;
		return;
	}
	//tr.Trace << "csa_nrs.size(): " << csa_nrs.size() << std::endl;

	Vector fij(0,0,0);

	for ( Size ii=1; ii<=csa_nrs.size(); ++ii ) {
		core::Size csa_nr = csa_nrs[ ii ];
		//tr.Trace << "csa_nr: " << csa_nr << std::endl;
		//tr.Trace << "aid.rsd(): " << aid.rsd() << std::endl;
		ChemicalShiftAnisotropy const& csa_cache( *retrieve_CSA_from_pose( pose ) );
		utility::vector1< core::scoring::CSA > All_CSA_lines( csa_cache.get_CSA_data() );
		runtime_assert( csa_nr <= All_CSA_lines.size() );
		CSA const& csa_data( All_CSA_lines[ csa_nr ] );
		conformation::Residue const& rsd1( pose.residue( csa_data.res1() ) );
		conformation::Residue const& rsd2( pose.residue( csa_data.res2() ) );
		conformation::Residue const& rsd3( pose.residue( csa_data.res3() ) );
		//tr.Trace << "rsd1.atom_name( aid.atomno() " << rsd1.atom_name( aid.atomno()) << std::endl;

		if ( aid.rsd() == csa_data.res1() && utility::trimmed_compare( rsd1.atom_name( aid.atomno() ), csa_data.atom1() ) ) {
			//tr.Trace << "aid.rsd(): " << aid.rsd() << " rsd1.atom_name( aid.atomno() ) "<< rsd1.atom_name( aid.atomno() ) << " csa_data.atom1() " << csa_data.atom1() <<  std::endl;
			fij += csa_data.f1ij();
		} else if ( aid.rsd() == csa_data.res2() && utility::trimmed_compare( rsd2.atom_name( aid.atomno() ), csa_data.atom2() ) ) {
			//tr.Trace << "aid.rsd(): " << aid.rsd() << " rsd2.atom_name( aid.atomno() ) "<< rsd2.atom_name( aid.atomno() ) << " csa_data.atom2() " << csa_data.atom2() << std::endl;
			fij += csa_data.f2ij();
		} else if ( aid.rsd() == csa_data.res3() && utility::trimmed_compare( rsd3.atom_name( aid.atomno() ), csa_data.atom3() ) ) {
			//tr.Trace << "aid.rsd(): " << aid.rsd() << " rsd3.atom_name( aid.atomno() ) "<< rsd3.atom_name( aid.atomno() ) << " csa_data.atom3() " << csa_data.atom3() << std::endl;
			fij += csa_data.f3ij();
		} else return;

	}

	//tr.Trace << "fij[0]: " << fij[0]<< " fij[1]: " << fij[1]<< " fij[2]: " << fij[2]<< std::endl;
	//tr.Trace << "torsion gradient: " << aid << std::endl;
	//tr.Trace << "score_weights[ core::scoring::csa ]: " << score_weights[ core::scoring::csa ] << std::endl;
	//thanks to Will Sheffler:
	numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
	numeric::xyzVector<core::Real> const f2( fij );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += score_weights[ core::scoring::csa ] * f1;
	F2 += score_weights[ core::scoring::csa ] * f2;

}

core::Size
ChemicalShiftAnisotropyEnergy::version() const
{
	return 1; // Initial versioning
}

} // scoring
} // core
