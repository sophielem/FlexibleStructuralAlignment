#!/usr/bin/env python3

"""
Functions to call soft and parse the output.
"""

import subprocess as sub
import re
import os
import glob
import src.parse_pdb as parse


def call_executabe(cmd_line):
    """
    Call function in shell and retrieve the ouput. Remove files created by
    TMalign and Protein Peeling.
    args:
        the command line
    return:
        the ouput of the soft
    """
    # Command to communicate with the shell and retrieve the output
    out, err = sub.Popen(cmd_line.split(), stdout=sub.PIPE).communicate()

    # To remove all uncessary files created
    for filename in glob.glob("file*"):
        os.remove(filename)
    # To convert the bytes output of suprocess into a str:
    return out.decode()


def parse_protein_peeling(output):
    """
    Parse the ouput of the protein peeling soft to retrieve PUs.
    args:
        the ouput of the protein peeling soft
    return:
        write pdb of PUs
    """
    # Remove lines uncessary
    output = output[15 : ]
    # Remove the new line
    output = output[ : -1]
    # Split the table in list
    # The index 0 correspond to the number of PU and the others the
    # delimitation of each PU
    output = [itx.split()[4: ] for itx in output]
    dict_pdb = parse.parse_pdb("data/1aoh.pdb")
    dict_pu = {}
    for itx in output:
        # Retrieve the number of PU
        nb_pu = int(itx[0])
        dict_pu[nb_pu] = []
        # For each delimitations, create the pdb file
        for i in range(1, nb_pu*2, 2):
            dict_pu[nb_pu].append(i)
            start = int(itx[i])
            end = int(itx[i + 1])
            parse.write_pdb(dict_pdb,
                            ("results/" + str(nb_pu) + "_" + str(i) + ".pdb"),
                            True, start, end)
    return dict_pu


def parse_tmscore(output, soft):
    """
    Parse the ouput of the TMalign and TMscore soft to retrieve the
    TMscore calulated.
    args:
        the ouput of the soft and the soft to print
    return:
        the TMscore
    """
    matches = re.search("TM-score *= (?P<score>[0-9]*\.[0-9]*)", output)
    print("{} Score normalized by length of the first chain \
           {}".format(soft, matches.group("score")))
    return float(matches.group("score"))


def tm_align(dict_pu, nb_pu, input2):
    """
    Execute TMalign for each PU vs the second protein. Keep the best TMscore
    and remove the region aligned of the second protein and then restart the
    loop for others PU.
    args:
        a dictionary of PU which contains a list of the index of each PU
        the number of PU
        the scond protein
    return:

    """
    # All PU are aligned
    if dict_pu[nb_pu] == []: return None
    # Initialize a minus score
    score = [-100, 0]
    for i in dict_pu[nb_pu]:
        # Read the pdb file of the PU
        pdb_file = "results/" + str(nb_pu) + "_" + str(i) + ".pdb"

        # TMalign
        CMD_LINE = ("bin/./TMalign " + input2 + " " + pdb_file + " -o TM.sup")
        OUTPUT = call_executabe(CMD_LINE)
        # Retrieve the TM score
        SC_TMALIGN = parse_tmscore(OUTPUT, "TMalign")
        # Keep the best score
        if SC_TMALIGN > score[0]:
            score[0] = SC_TMALIGN
            score[1] = i
    dict_pu[nb_pu].remove(score[1])
    # EFFACER LE FICHIER PDB ET RECOMMENCER AVEC LES AUTRES PUs
    tm_align(dict_pu, nb_pu, input2)


def remove_aligned_region(input):
    """

    """
    with open("TM.sup_all_atm", "r") as filin:
        lines = filin.readlines()

    for i, line in enumerate(lines):
        # Find the first line
        if line[0:4] == "ATOM" and line[12:16].strip() == "CA":
            break
    ids = []
    amino = []
    # Keep all residues and index to remove
    while lines[i][0:4].strip() != "TER":
        if lines[i][12:16].strip() == "CA":
            ids.append(lines[i][22:27].strip())
            amino.append(lines[i][17:20])
        i += 1

    
