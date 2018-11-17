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
    index = 1 if soft == "TMalign" else 0
    matches = re.findall("TM-score *= (?P<score>[0-9]*\.[0-9]*)", output)
    if matches == []:
        return None
    else:
        print("{} Score normalized by length of the second chain \
               {}".format(soft, matches[index]))
        return float(matches[index])


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
        CMD_LINE = ("bin/./TMalign " + pdb_file + " " + input2 + " -o TM" + str(i) + ".sup")
        OUTPUT = call_executabe(CMD_LINE)
        # Retrieve the TM score
        SC_TMALIGN = parse_tmscore(OUTPUT, "TMalign")
        if SC_TMALIGN is None:
            break
        else:
            # Keep the best score
            if SC_TMALIGN > score[0]:
                score[0] = SC_TMALIGN
                score[1] = i
    # The second protein is smaller than the first so some regions cant be
    # aligned. They are write in the pdb without alignment
    print(score[0], score[1])
    if SC_TMALIGN is None:
        for i in dict_pu[nb_pu]:
            pdb_file = "results/" + str(nb_pu) + "_" + str(i) + ".pdb"
            with open(pdb_file, "r") as filin:
                text = filin.readlines()
            with open("test_aligned.pdb", "a") as filout:
                for line in text:
                    filout.write(line)
    else:
        dict_pu[nb_pu].remove(score[1])
        remove_aligned_region(input2, str(score[1]))
        tm_align(dict_pu, nb_pu, input2)


def remove_aligned_region(input, idx_pu):
    """

    """
    with open("TM" + idx_pu + ".sup_all_atm", "r") as filin:
        lines = filin.readlines()

    for i, line in enumerate(lines):
        # Find the first line concerning pdb
        if line[0:4] == "ATOM":
            break

    with open("test_aligned.pdb", "a") as filout:
        # write the pdb PU aligned
        while lines[i][0:4].strip() != "TER":
            filout.write(lines[i])
            i += 1
    i += 1
    ids = []
    # Keep all residues and index to remove
    while lines[i][0:4].strip() != "TER":
        ids.append(lines[i][30:38].strip())
        i += 1
    k = 0
    with open(input, "r") as file_protein:
        text = file_protein.readlines()
    with open(input, "w") as file_protein:
        for line in text:
            if line[30:38].strip() not in ids:
                file_protein.write(line)

    # file_protein.truncate()
    # file_protein.close()
