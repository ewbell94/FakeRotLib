// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/STMStoredTask.hh
/// @brief
/// @author Neil King (neilking@uw.edu)


#ifndef INCLUDED_protocols_task_operations_STMStoredTask_hh
#define INCLUDED_protocols_task_operations_STMStoredTask_hh

#include <protocols/task_operations/STMStoredTask.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <map>
#include <string>

#include <core/pack/task/PackerTask.fwd.hh> // AUTO IWYU For PackerTaskOP

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace task_operations {

class STMStoredTask: public basic::datacache::CacheableData {
public:

	// constructors
	STMStoredTask();
	STMStoredTask( const STMStoredTask& rval );
	basic::datacache::CacheableDataOP clone() const override;
	virtual basic::datacache::CacheableDataOP fresh_instance() const;

	// setter
	void set_task( core::pack::task::PackerTaskOP const task, std::string const & task_name );

	// getters
	core::pack::task::PackerTaskOP get_task( std::string const & task_name ) const;
	bool has_task( std::string const & task_name ) const;

private:
	std::map< std::string, core::pack::task::PackerTaskOP > tasks_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // task_operations
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_task_operations_STMStoredTask )
#endif // SERIALIZATION


#endif
