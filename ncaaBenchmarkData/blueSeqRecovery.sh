#!/bin/bash

ROTDIR=$1
TARGET=$2
TNAME=`echo $TARGET | sed "s/.*\///"`

mkdir -p $ROTDIR/protres/$TNAME
/dors/meilerlab/data/belle6/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease -s $TARGET -nstruct 1 -out:prefix $ROTDIR/protres/$TNAME/prd_ -parser:protocol seqRecovery.xml \
-pack_missing_sidechains false -out:level 300 -ignore_zero_occupancy false -extra_res_fa \
$ROTDIR/XXA/XXA.params \
$ROTDIR/XXC/XXC.params \
$ROTDIR/XXD/XXD.params \
$ROTDIR/XXE/XXE.params \
$ROTDIR/XXF/XXF.params \
$ROTDIR/XXG/XXG.params \
$ROTDIR/XXH/XXH.params \
$ROTDIR/XXI/XXI.params \
$ROTDIR/XXK/XXK.params \
$ROTDIR/XXL/XXL.params \
$ROTDIR/XXM/XXM.params \
$ROTDIR/XXN/XXN.params \
$ROTDIR/XXP/XXP.params \
$ROTDIR/XXQ/XXQ.params \
$ROTDIR/XXR/XXR.params \
$ROTDIR/XXS/XXS.params \
$ROTDIR/XXT/XXT.params \
$ROTDIR/XXV/XXV.params \
$ROTDIR/XXW/XXW.params \
$ROTDIR/XXY/XXY.params

