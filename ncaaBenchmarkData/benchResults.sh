#!/bin/bash
#SBATCH -A p_csb_meiler
#SBATCH -p production
#SBATCH -t 24:00:00
#SBATCH -o /dors/meilerlab/data/belle6/rtm.txt
#SBATCH --mem=10G
#SBATCH -n 32

TASK="rot"
export FORM="prd"

if [ $FORM == "rtm" ] || [ $FORM == "prm" ]; then
    cd /dors/meilerlab/data/belle6/ncaaBenchmark/orgpdb
    if [ $TASK == "atom" ]; then
	parallel -j 32 'i="{}"; echo $i Native:$(../rotRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatRandCharge:$(../rotRecovery.py ../NatRandCharge/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatBCLAtom:$(../rotRecovery.py ../NatBCLAtom/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatM2PPAtom:$(../rotRecovery.py ../NatM2PPAtom/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    elif [ $TASK == "geo" ]; then 
	parallel -j 32 'i="{}"; echo $i Native:$(../rotRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatBCLGeo:$(../rotRecovery.py ../NatBCLGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatGauGeo:$(../rotRecovery.py ../NatGauGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatUFFGeo:$(../rotRecovery.py ../NatUFFGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatMMFFGeo:$(../rotRecovery.py ../NatMMFFGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    elif [ $TASK == "rot" ]; then
	parallel -j 32 'i="{}"; echo $i Native:$(../rotRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) Parent:$(../rotRecovery.py ../Parent/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) RDKit:$(../rotRecovery.py ../RDKit/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) MakeRotLib:$(../rotRecovery.py ../MakeRotLib/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) BCL:$(../rotRecovery.py ../BCL/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) AutoRotLib:$(../rotRecovery.py ../AutoRotLib/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NoRot:$(../rotRecovery.py ../NoRot/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) FakeRotLib:$(../rotRecovery.py ../FakeCart/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) FRLNoExp:$(../rotRecovery.py ../FakeCartNon/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    fi
elif [ $FORM == "prd" ] || [ $FORM == "fsd" ]; then
    cd /dors/meilerlab/data/belle6/ncaaBenchmark/nosc_sub
    if [ $TASK == "atom" ]; then
	parallel -j 32 'i="{}"; echo $i Native:$(../seqRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb ../orgpdb/$i 2> /dev/null) NatRandCharge:$(../seqRecovery.py ../NatRandCharge/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatBCLAtom:$(../seqRecovery.py ../NatBCLAtom/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatM2PPAtom:$(../seqRecovery.py ../NatM2PPAtom/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    elif [ $TASK == "geo" ]; then
	parallel -j 32 'i="{}"; echo $i Native:$(../seqRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatBCLGeo:$(../seqRecovery.py ../NatBCLGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatGauGeo:$(../seqRecovery.py ../NatGauGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatUFFGeo:$(../seqRecovery.py ../NatUFFGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NatMMFFGeo:$(../seqRecovery.py ../NatMMFFGeo/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    elif [ $TASK == "rot" ]; then
	parallel -j 32 'i="{}"; echo $i Native:$(../seqRecovery.py ../Native/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) Parent:$(../seqRecovery.py ../Parent/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) RDKit:$(../seqRecovery.py ../RDKit/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) MakeRotLib:$(../seqRecovery.py ../MakeRotLib/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) BCL:$(../seqRecovery.py ../BCL/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) AutoRotLib:$(../seqRecovery.py ../AutoRotLib/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) NoRot:$(../seqRecovery.py ../NoRot/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) FakeRotLib:$(../seqRecovery.py ../FakeCart/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null) FRLNoExp:$(../seqRecovery.py ../FakeCartNon/protres/$i/$FORM\_$i\_0001.pdb $i 2> /dev/null)' ::: *
    fi
fi
