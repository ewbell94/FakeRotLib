// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperation.hh
/// @brief  Class that performs an operation on ResidueLevelTask,
///         usually by a PackerTaskFactory right after the task's construction.
/// @author ashworth

// only generalized base classes go here. ResLvlTaskOperations that actually do things do not belong here.

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperation_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperation_hh

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>

// Project Headers
#include <core/pack/task/PackerTask.fwd.hh>

// Basic headers

// Utility Headers
#include <utility/VirtualBase.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers


#include <basic/citation_manager/CitationCollectionBase.fwd.hh> // AUTO IWYU For CitationCollectionList

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResLvlTaskOperation : public utility::VirtualBase
{
public:
	typedef utility::tag::TagCOP TagCOP;
public:
	/// @brief Create another task operation of the type matching the most-derived
	/// version of the class.
	virtual ResLvlTaskOperationOP clone() const = 0;

	/// @brief Change a ResidueLevelTask in some way.  The input pose is the one to which the input
	/// task will be later applied.
	virtual void apply( ResidueLevelTask & ) const = 0;

	/// @brief parser xml-like tags to set class data/parameters
	virtual void parse_tag( TagCOP ) {}

	/// @brief Provide citations to the passed CitationCollectionList
	/// Subclasses should add the info for themselves and any other classes they use.
	/// @details The default implementation of this function does nothing.  It may be
	/// overriden by task operations wishing to provide citation information.
	virtual void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const {
		//GNDN.
	}
};

// only generalized base classes go here. ResLvlTaskOperations that actually do things do not belong here.

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
