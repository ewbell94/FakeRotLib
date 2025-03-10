#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_rmsd_dict = {}
cutoffs_score_dict = {}
failures = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "rmsBB_if"
y_label = "reweighted_sc"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_score = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_score = map( float, cutoffs_score )
cutoffs_rmsd_dict.update( dict( zip ( protein, cutoffs_rmsd )))
cutoffs_score_dict.update( dict( zip ( protein, cutoffs_score )))

# open results output file
f = open( outfile, "w" )

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):

	target_results = {}

	# read in score file, scores are sorted, first one is lowest
	x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

	# map values to floats (were strings)
	x = list( map( float, x ))
	y = list( map( float, y ))
	y, x = (list(t) for t in zip(*sorted(zip(y, x))))

	# only look at top10 scores
	y = y[0:10]
	x = x[0:10]

	# check for RMSDs below cutoff - 80% of the top10-scoring points should be below cutoff
	f.write( targets[i] + "\t" )
	val_cutoff = qm.check_xpercent_values_below_cutoff( x, cutoffs_rmsd_dict[targets[i]], "rmsd", f, 50)
	target_results.update( val_cutoff )

    # add to failues
	if val_cutoff['All rmsds < cutoff'] == False:
		failures.append( targets[i] )

	# check for scores below cutoff
	# f.write( targets[i] + "\t" )
	# val_cutoff = qm.check_xpercent_values_below_cutoff( y, cutoffs_score_dict[targets[i]], "score", f ,10/len(x))
	# target_results.update( val_cutoff )

    # add to failures
	# if val_cutoff['All scores < cutoff'] == False:
		# failures.append( targets[i] )

	# check for RMSD range
	f.write( targets[i] + "\t" )
	val_rms = qm.check_range( x, "rmsBB_if", f )
	target_results.update( val_rms )

	# check for score range
	f.write( targets[i] + "\t" )
	val_score = qm.check_range( y, "reweighted_sc", f )
	target_results.update( val_score )

	results.update( {targets[i] : target_results} )
	f.write( "\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
