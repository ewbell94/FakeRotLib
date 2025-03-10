// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/EnableMultiCoolAnnealerCreator.cc
/// @brief  Creator for Task Operation for turning on the multi-cool annealer
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_task_operation_EnableMultiCoolAnnealerCreator_HH
#define INCLUDED_core_pack_task_operation_EnableMultiCoolAnnealerCreator_HH

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class EnableMultiCoolAnnealerCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};



} //core
} //pack
} //task
} //operation


#endif //INCLUDED_core/pack/task/operation_EnableMultiCoolAnnealerCreator_HH
