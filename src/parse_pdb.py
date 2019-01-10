#!/usr/bin/env python3

"""
Parse and write PDB files
"""


def parse_pdb(infile, chain, reindex, bfactor=True):
    """ purpose: to parse a pdb file (infile)
        args:
            PDB file
        return:
            a dico dict_pdb which contains for each atom of each residue of
        each chain, its corresponding 3D coordinates. Please take a look to
        the code to understand the structure of the dico.
    """
    # Read the pdb file
    file_pdb = open(infile, "r")
    lines = file_pdb.readlines()
    file_pdb.close()

    # Initialization variables
    dict_pdb = {}
    dict_pdb["reslist"] = []
    current_res = 1
    prev_true_id = None

    for line in lines:
        if ((line[0:4] == "ATOM") and line[21] == chain):
            true_id = line[22:27].strip()
            # Renumber the pdb to parse it after with protein peeling index
            if true_id != prev_true_id and prev_true_id is not None:
                current_res += 1
            if reindex:
                code_res = str(current_res) + line[26:27].strip()
            else:
                code_res = str(true_id) + line[26:27].strip()
            # First time we encounter the index
            if code_res not in dict_pdb["reslist"]:
                dict_pdb["reslist"].append(code_res)
                dict_pdb[code_res] = {}
                dict_pdb[code_res]["resname"] = line[17:20].strip()
                dict_pdb[code_res]["atomlist"] = []
                alternateoccupancy = None
                occupancy = "%s"%(line[16:17])
                if occupancy != " ":
                    alternateoccupancy = occupancy
            # This is not a new residue
            else:
                occupancy = "%s"%(line[16:17])

                # Means we are in the first alternate location of that residue
                if occupancy != " " and alternateoccupancy is None:
                    alternateoccupancy = occupancy

            atomtype = line[12:16].strip()
            # Means this atom corresponds to the first rotamer found
            # in the PDB for this residue
            if occupancy == alternateoccupancy or occupancy == " ":
                dict_pdb[code_res]["atomlist"].append(atomtype)
                dict_pdb[code_res][atomtype] = {}
                dict_pdb[code_res][atomtype]["x"] = float(line[30:38])
                dict_pdb[code_res][atomtype]["y"] = float(line[38:46])
                dict_pdb[code_res][atomtype]["z"] = float(line[46:54])
                dict_pdb[code_res][atomtype]["id"] = line[6:11].strip()
                if bfactor is True:
                    dict_pdb[code_res][atomtype]["bfactor"] = float(line[60:67].strip())
            if reindex:
                dict_pdb[code_res]["resnum"] = current_res
            else:
                dict_pdb[code_res]["resnum"] = int(true_id)
            prev_true_id = true_id

    return dict_pdb


def write(dict_pdb, bfactor, res, atom, fout):
    """
    write the line in the pdb file
    args:
        the dictionnary containing the pdb
        the res and atom read
        the output file
    """
    fout.write("ATOM  {:5s} {:^4s} {:3s} {:1s}{:>4s}   {:8.3f}{:8.3f}{:8.3f}\
                1.00{:7.3f} X X\n".format(dict_pdb[res][atom]["id"], atom,
                                          dict_pdb[res]["resname"], "A", res,
                                          dict_pdb[res][atom]["x"],
                                          dict_pdb[res][atom]["y"],
                                          dict_pdb[res][atom]["z"],
                                          bfactor))


def write_pdb(dict_pdb, filout, bfactor, start, end):
    """
    purpose: according to the coordinates in dict_pdb, writes the
    corresponding PDB file.
    If bfactor = True, writes also the information corresponding to the key
    bfactor of each residue (one key per residue) in dict_pdb.
    args:
        a dico with the dict_pdb format
    return:
        PDB file.
    """
    with open(filout, "w") as fout:
        for res in dict_pdb["reslist"]:
            # Write only the PU
            if (dict_pdb[res]["resnum"] >= start and
                    dict_pdb[res]["resnum"] <= end):
                for atom in dict_pdb[res]["atomlist"]:
                    if bfactor:
                        write(dict_pdb, dict_pdb[res][atom]["bfactor"], res,
                              atom, fout)
                    else:

                        write(dict_pdb, 1.00, res, atom, fout)


def write_pdb2(dict_pdb, filout, bfactor):
    """
    purpose: according to the coordinates in dict_pdb, writes the
    corresponding PDB file.
    If bfactor = True, writes also the information corresponding to the key
    bfactor of each residue (one key per residue) in dict_pdb.
    args:
        a dico with the dict_pdb format
    return:
        PDB file.
    """
    with open(filout, "w") as fout:
        for res in dict_pdb["reslist"]:
            for atom in dict_pdb[res]["atomlist"]:
                if bfactor:
                    write(dict_pdb, dict_pdb[res][atom]["bfactor"], res,
                          atom, fout)
                else:
                    write(dict_pdb, 1.00, res, atom, fout)
