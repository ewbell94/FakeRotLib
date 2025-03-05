#!/dors/meilerlab/data/belle6/miniforge3/envs/glypred/bin/python

import numpy as np
from sys import argv

def pdb2fasta(pdb):
    f = open(pdb)
    threetoone = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I",
                  "LYS":"K", "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S",
                  "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y", "XXC":"C", "XXD":"D", "XXE":"E", "XXF":"F",
                  "XXH":"H", "XXI":"I", "XXK":"K", "XXL":"L", "XXM":"M", "XXN":"N", "XXQ":"Q", "XXR":"R",
                  "XXS":"S", "XXT":"T", "XXV":"V", "XXW":"W", "XXY":"Y", "XXA":"A", "XXG":"G", "XXP":"P"}

    seq = ""
    resnum = ""
    letter = ""
    bba = [False, False, False]
    for line in f:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            if resnum != line[22:26]:
                if all(bba):
                    seq += letter
                bba = [False, False, False]
            resnum = line[22:26]
            if line[12:16].strip() == "N":
                bba[0] = True
            elif line[12:16].strip() == "CA":
                bba[1] = True
                try:
                    letter = threetoone[line[17:20]]
                except:
                    print("ERROR: AA not found %s"%line[17:20])
                    exit(2)
            elif line[12:16].strip() == "C":
                bba[2] = True
    f.close()
    if all(bba):
        seq += letter
    return seq

seqA = pdb2fasta(argv[1])
seqB = pdb2fasta(argv[2])

if len(seqA) != len(seqB):
    print("ERROR: sequence lengths don't match %d %d"%(len(seqA),len(seqB)))
    exit(1)

alphabet = "ACDEFGHIKLMNPQRSTVWY"
counts = np.zeros((20,20))
for i in range(len(seqA)):
    indexA = alphabet.index(seqA[i])
    indexB = alphabet.index(seqB[i])
    counts[indexA,indexB] += 1

np.save(argv[1]+".npy", counts)
print(np.trace(counts)/len(seqA))

