#!/usr/bin/env python

import math, string


def parsePDBMultiChains(infile, charge = 1, chargeFromInfile = False, bfactor = True, CG = False) :
    """ purpose: to parse a pdb file (infile)
        input: PDB file
        output: a dico dPDB which contains for each atom of each residue of
        each chain, its corresponding 3D coordinates. Please take a look to
        the code to understand the structure of the dico.
    """
    # lecture du fichier PDB
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var init
    dPDB = {}
    dPDB["reslist"] = []
    current_res = 1
    prev_res = None
    # parcoure le PDB
    for line in lines :
        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and ( (line[17:20].strip() == "MET") or  (line[17:20].strip() == "MSE") )) :
            code_res = str(current_res) + line[26:27].strip()
            cures = line[22:27].strip()
            if(cures != prev_res and prev_res != None) :  current_res += 1
            if not code_res in dPDB["reslist"] : # first time we encounter it
                dPDB["reslist"].append(code_res)
                dPDB[code_res] = {}
                dPDB[code_res]["resname"] = line[17:20].strip()
                dPDB[code_res]["atomlist"] = []
                #dPDB[curres]["atomlistTowrite"] = []
                alternateoccupancy = None #"%s"%(line[16:17])
                occupancy = "%s"%(line[16:17])
                if occupancy != " " :
                    alternateoccupancy = occupancy

            else: # this is not a new residue
                occupancy = "%s"%(line[16:17])

                if occupancy != " " and alternateoccupancy == None : # means we are in the first alternate location of that residue
                    alternateoccupancy = occupancy

            if CG : # means we are parsing a CG model so we have to treat the CSE atomtypes which can be redondant in terms of name the same res
                atomtype = "%s_%s"%(line[6:11].strip(), line[12:16].strip())
            else:
                atomtype = line[12:16].strip()

            #if not atomtype in dPDB[curres]["atomlist"] :
            if occupancy == alternateoccupancy  or occupancy == " " : # means this atom corresponds to the first rotamer found in the PDB for this residue

                dPDB[code_res]["atomlist"].append(atomtype)
                dPDB[code_res][atomtype] = {}
                dPDB[code_res][atomtype]["x"] = float(line[30:38])
                dPDB[code_res][atomtype]["y"] = float(line[38:46])
                dPDB[code_res][atomtype]["z"] = float(line[46:54])
                dPDB[code_res][atomtype]["id"] = line[6:11].strip()
                if bfactor == True :
                    dPDB[code_res][atomtype]["bfactor"] = float(line[60:67].strip())

            dPDB[code_res]["resnum"] = current_res
            prev_res = cures


    return dPDB


def writePDB(dPDB, filout, bfactor, start, end) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each residue (one key per residue) in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    for res in dPDB["reslist"] :
        if dPDB[res]["resnum"] >= start + 1 and dPDB[res]["resnum"] <= end + 1:
            for atom in dPDB[res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(dPDB[res][atom]["id"], atom, dPDB[res]["resname"],"A", res,dPDB[res][atom]["x"], dPDB[res][atom]["y"],dPDB[res][atom]["z"],dPDB[res][atom]["bfactor"] ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[res][atom]["id"], atom, dPDB[res]["resname"],"A", res,dPDB[res][atom]["x"], dPDB[res][atom]["y"],dPDB[res][atom]["z"] ))

    fout.close()
