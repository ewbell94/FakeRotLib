// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/FaMPEnv.cxxtest.hh
///
/// @brief   Unit Test: FaMPEnv Object
/// @details This suite is for the FaMPEnv energy term.
///          It includes a test to make sure the analytical
///          and numeric derivatives are equal.
///
/// @author  Hope Woods (hope.woods@vanderbilt.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>


//Unit Headers

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>


#include <basic/Tracer.hh>


// Utility Headers

using namespace core;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;

static basic::Tracer TR("core.scoring.membrane.FaMPEnv.cxxtest");

class FaMPEnvTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {
		using namespace core::import_pose;
		using namespace core::pose;


		// Initialize core & options system
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior" );
		TR << "core_init" << std::endl;

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods ////////////////////////////////


	void test_fampenvenergy_deriv_pore() {
		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/scoring/membrane/2KSF_tr.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/scoring/membrane/2KSF_tr.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5);

	}

	void test_fampenvenergy_vesicle_deriv() {
		core_init_with_additional_options( "-mp:restore_lazaridis_imm_behavior -mp:geometry vesicle -mp:geo:vesicle_radius 100 -has_pore 0" );

		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1PY6_mp_coords.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1PY6.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( FaMPEnv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);

	}


private:



}; // FaMPEnv unit test
