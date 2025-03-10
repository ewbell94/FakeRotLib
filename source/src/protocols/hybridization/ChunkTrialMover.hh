// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Wrapper for InsertChunkMover. It can take a random template and steal coordinates of all chunks or a random one
/// @details
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_ChunkTrialMover_hh
#define INCLUDED_protocols_hybridization_ChunkTrialMover_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/ChunkTrialMover.fwd.hh>


#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>


#include <set>

namespace protocols {
namespace hybridization {

enum AlignOption { all_chunks, random_chunk };

class ChunkTrialMover: public protocols::moves::Mover
{

public:
	ChunkTrialMover(
		utility::vector1< core::pose::PoseCOP > const & template_poses,
		utility::vector1< protocols::loops::Loops > const & template_chunks,
		bool random_template = true,
		AlignOption align_option = all_chunks,
		utility::vector1<bool> residue_sample_template=utility::vector1<bool>() );

	void
	get_alignment_from_template(
		core::pose::PoseCOP template_pose,
		std::map <core::Size, core::Size> & seqpos_alignment );

	void apply(core::pose::Pose & pose) override;

	void set_template(core::Size const template_number);
	void set_max_registry_shift( core::Size max_registry_shift_in );
	core::Size template_number();
	void pick_random_template();
	void pick_random_chunk(core::pose::Pose & pose);
	core::Size trial_counter(core::Size ires);
	std::string get_name() const override;
	bool has_valid_moves( ) { return has_valid_moves_; }

	void set_templates_to_ignore( std::set< core::Size> template_indices_in ) {
		ignore_template_indices_ = template_indices_in;
	}

private:
	utility::vector1 < core::pose::PoseCOP > template_poses_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;

	bool random_template_, has_valid_moves_;

	AlignOption align_option_;
	InsertChunkMover align_chunk_;
	utility::vector1 < std::map <core::Size, core::Size> > sequence_alignments_;
	utility::vector1 < core::Size > residue_covered_by_template_;

	core::Size template_number_; // the jump to be realigned
	core::Size jump_number_; // the jump to be realigned
	core::Size highest_tmpl_resnum_; // the highest residue number from all templates

	std::set< core::Size > ignore_template_indices_;
	utility::vector1<bool> sampling_chunk_;
	core::Size max_registry_shift_global_;
}; //class ChunkTrialMover

} // hybridization
} // protocols

#endif
