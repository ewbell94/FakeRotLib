// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/benchmark/performance/Docking.bench.hh
///
/// @brief  Run a performance benchmark of docking
/// @author Monica Berrondo

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/docking/DockingProtocol.hh>

//#include <basic/Tracer.hh>
//#include <numeric/random/random.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>
#include <apps/benchmark/performance/init_util.hh>




enum DockType {High, Low};

template  <DockType dock, int TScale>
class DockingBenchmark : public PerformanceBenchmark
{
public:
	core::pose::PoseOP start_pose;
	protocols::docking::DockingProtocolOP docking;

	DockingBenchmark(std::string name) : PerformanceBenchmark(name) {};

	void setUp() override {
		if ( dock == Low ) core_init_with_additional_options( "-low_res_protocol_only -dock_pert 3 8 -run:constant_seed" );
		if ( dock == High ) core_init_with_additional_options( "-docking_local_refine -run:constant_seed" );
		start_pose = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file(*start_pose, "dock_in.pdb", core::import_pose::PDB_file);
		docking = utility::pointer::make_shared< protocols::docking::DockingProtocol >();
		docking->set_native_pose(start_pose);
		docking->set_input_pose(start_pose);
	};


	void run(core::Real scaleFactor) override {
		//for(int i=0; i<10; i++) {
		// std::cout << "i="<< i << " R=" << numeric::random::uniform() << std::endl;
		//}
		int reps( (int)(TScale*scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // do at least one repetition, regardless of scaling
		for ( int i=0; i<reps; i++ ) {
			core::pose::Pose pose;
			pose = *start_pose;
			docking->apply( pose );
		}
	};

	void tearDown() override {};
};

typedef DockingBenchmark<Low,  2> DockingBenchmark_low;
typedef DockingBenchmark<High, 1> DockingBenchmark_high;
