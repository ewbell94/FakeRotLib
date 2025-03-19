#!/bin/bash

#SBATCH --exclude bf[49-51]
#SBATCH -p production
#SBATCH -t 160:00:00
#SBATCH --mem=20G
#SBATCH -o /dors/meilerlab/data/belle6/ncaaBenchmark/logs/%A_%a.log
#SBATCH -a 1-1000

ROTSET="$1"
echo $ROTSET

#file=$(ls /dors/meilerlab/data/belle6/ncaaBenchmark/nosc/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
#/dors/meilerlab/data/belle6/ncaaBenchmark/blueRotRecovery.sh /dors/meilerlab/data/belle6/ncaaBenchmark/$ROTSET $file
#basename=`basename $file`
#/dors/meilerlab/data/belle6/ncaaBenchmark/unMutPDB.sh /dors/meilerlab/data/belle6/ncaaBenchmark/$ROTSET/protres/$basename/*.pdb

file=$(ls /dors/meilerlab/data/belle6/ncaaBenchmark/nosc_sub/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
/dors/meilerlab/data/belle6/ncaaBenchmark/blueSeqRecovery.sh /dors/meilerlab/data/belle6/ncaaBenchmark/$ROTSET $file
