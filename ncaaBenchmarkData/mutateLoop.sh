#!/bin/bash

mutscript="Mutate_simple.xml"
rotset="$1"
prot="l_straight"

for i in C D E F H I L M N Q S T V W Y; do
    #for i in D E K R; do
    ./Mutate.sh $mutscript $prot.pdb $rotset/XX$i/XX$i.params 8 XX$i $rotset/XX$i/mut_ &
done
wait

for i in C D E F H I K L M N Q R S T V W Y; do
#for i in C D E F H I L M N Q S T V W Y; do
    ./compute_dihedrals.sh $rotset/XX$i/mut_$prot\_0001.pdb.gz $rotset/XX$i/mut_$prot\_0001_MH_Traj.pdb.gz 8 $rotset/XX$i/$prot\_dihedrals.txt
    ./plot.sh $rotset/XX$i/$prot\_dihedrals.txt
done
