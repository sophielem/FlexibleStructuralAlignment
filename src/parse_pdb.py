#!/usr/bin/env python

"""
Parse and write PDB files
"""

import math, string


def parsePDBMultiChains(infile, charge = 1, chargeFromInfile = False, bfactor = True, CG = False) :
    """ purpose: to parse a pdb file (infile)
        input: PDB file
        output: a dico dPDB which contains for each atom of each residue of
        each chain, its corresponding 3D coordinates. Please take a look to
        the code to understand the structure of the dico.
    """
    # Read the pdb file
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # Initialization variables
    dPDB = {}
    dPDB["reslist"] = []
    current_res = 1
    prev_true_id = None

    for line in lines :
        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and ( (line[17:20].strip() == "MET") or  (line[17:20].strip() == "MSE") )) :
            true_id = line[22:27].strip()
            # Renumber the pdb to parse it after with protein peeling index
            if true_id != prev_true_id and prev_true_id != None:
                current_res += 1
            code_res = str(current_res) + line[26:27].strip()
            # first time we encounter the index
            if not code_res in dPDB["reslist"] :
                dPDB["reslist"].append(code_res)
                dPDB[code_res] = {}
                dPDB[code_res]["resname"] = line[17:20].strip()
                dPDB[code_res]["atomlist"] = []
                alternateoccupancy = None
                occupancy = "%s"%(line[16:17])
                if occupancy != " " :
                    alternateoccupancy = occupancy
            # this is not a new residue
            else:
                occupancy = "%s"%(line[16:17])

                # means we are in the first alternate location of that residue
                if occupancy != " " and alternateoccupancy == None :
                    alternateoccupancy = occupancy
            # means we are parsing a CG model so we have to treat the CSE
            # atomtypes which can be redondant in terms of name the same res
            if CG :
                atomtype = "%s_%s"%(line[6:11].strip(), line[12:16].strip())
            else:
                atomtype = line[12:16].strip()
            # means this atom corresponds to the first rotamer found
            # in the PDB for this residue
            if occupancy == alternateoccupancy  or occupancy == " " :
                dPDB[code_res]["atomlist"].append(atomtype)
                dPDB[code_res][atomtype] = {}
                dPDB[code_res][atomtype]["x"] = float(line[30:38])
                dPDB[code_res][atomtype]["y"] = float(line[38:46])
                dPDB[code_res][atomtype]["z"] = float(line[46:54])
                dPDB[code_res][atomtype]["id"] = line[6:11].strip()
                if bfactor == True :
                    dPDB[code_res][atomtype]["bfactor"] = float(line[60:67].strip())

            dPDB[code_res]["resnum"] = current_res
            prev_true_id = true_id

    return dPDB


def write(dPDB, bfactor, res, atom, fout):
    """write the line in the pdb file
    """
    fout.write("ATOM  {:5s} {:^4s}{:3s} {:1s}{:>4s}   {:8.3f}{:8.3f}{:8.3f}\
                1.00{:7.3f} X X\n".format(dPDB[res][atom]["id"],atom,
                                      dPDB[res]["resname"], "A", res,
                                      dPDB[res][atom]["x"], dPDB[res][atom]["y"],
                                      dPDB[res][atom]["z"],
                                      bfactor))


def writePDB(dPDB, filout, bfactor, start, end) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each residue (one key per residue) in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """
    fout = open(filout, "w")

    for res in dPDB["reslist"] :
        # Write only the PU
        if dPDB[res]["resnum"] >= start and dPDB[res]["resnum"] <= end :
            for atom in dPDB[res]["atomlist"] :
                if bfactor :
                    write(dPDB, dPDB[res][atom]["bfactor"], res, atom, fout)
                else:
                    write(dPDB, 1.00, res, atom, fout)
    fout.close()
