#!/dors/meilerlab/data/belle6/miniforge3/envs/glypred/bin/python

from Bio.PDB.PDBParser import PDBParser
from sys import argv

parser = PDBParser()
querypose=parser.get_structure("query",argv[1])[0][str(argv[2][4])]
refpose=parser.get_structure("ref",argv[2])[0][str(argv[2][4])]
querypose.atom_to_internal_coordinates()
refpose.atom_to_internal_coordinates()
qreslist = [r for r in querypose.get_residues() if r.internal_coord.pick_angle("tau") != None]
rreslist = [r for r in refpose.get_residues() if r.internal_coord.pick_angle("tau") != None]
if len(qreslist) != len(rreslist):
    print("ERROR: query and reference have different lengths %d %d"%(len(qreslist),len(rreslist)))
    exit(1)

recover = 0.
valid = 0.
for i in range(len(qreslist)):
    rres = rreslist[i].get_resname()
    qres = qreslist[i].get_resname()
    if rres != qres:
        print("ERROR: query and reference have residue mismatch at pos %d"%(i+1))
        exit(1)
    
    if rres == "ALA" or rres == "GLY":
        continue

    #print(rres)
    rchi = rreslist[i].internal_coord.pick_angle("chi1")
    if rchi == None:
        continue

    valid += 1.
    recovered = True
    j = 1
    while rchi != None:
        qchi = qreslist[i].internal_coord.pick_angle("chi"+str(j))
        #print(rchi,qchi)
        if qchi == None:
            print("ERROR: query does not have chi%d when reference does at pos %d:%s"%(j,i+1,rres))
            exit(1)
        else:
            qang = qchi.angle
            rang = rchi.angle
            anglediff = min(abs(qang-rang),360-abs(qang-rang))
            #print(i,j,qang,rang,anglediff)

            #Symmetry correction for terminal torsion angles
            if (rres == "PHE" and j == 2) or (rres == "TYR" and j == 2) or (rres == "ASP" and j == 2) or (rres == "GLU" and j == 3):
                anglediff = min(anglediff, 180-anglediff)

            if anglediff > 20.:
                #print("bigchi")
                recovered=False
                break
        j+=1
        rchi = rreslist[i].internal_coord.pick_angle("chi"+str(j))

    if recovered:
        recover += 1.
if valid > 30:
    print(recover/valid)
else:
    print("BAD")
