#!/usr/bin/env python3

"""
Functions to call soft and parse the output.
"""

import subprocess as sub
import re
import os
import glob
import copy
import matplotlib.pyplot as plt
import src.parse_pdb as parse

DEBUG = False


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
    out = sub.Popen(cmd_line.split(), stdout=sub.PIPE).communicate()[0]

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
    output = output[15:]
    # Remove the new line
    output = output[:-1]
    # Split the table in list
    # The index 0 correspond to the number of PU and the others the
    # delimitation of each PU
    output = [itx.split()[4:] for itx in output]
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
        if DEBUG:
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
    if dict_pu[nb_pu] == []:
        return None
    # Initialize a minus score
    score = [-100, 0]
    for i in dict_pu[nb_pu]:
        # Read the pdb file of the PU
        pdb_file = "results/" + str(nb_pu) + "_" + str(i) + ".pdb"

        # TMalign
        cmd_line = ("bin/./TMalign " + pdb_file + " " + input2 +
                    " -o TM" + str(i) + ".sup")
        output = call_executabe(cmd_line)
        # Retrieve the TM score
        sc_tmalign = parse_tmscore(output, "TMalign")
        if sc_tmalign is None:
            break
        else:
            # Keep the best score
            if sc_tmalign > score[0]:
                score[0] = sc_tmalign
                score[1] = i
    # The second protein is smaller than the first so some regions cant be
    # aligned. They are write in the pdb without alignment
    if sc_tmalign is None:
        for i in dict_pu[nb_pu]:
            pdb_file = "results/" + str(nb_pu) + "_" + str(i) + ".pdb"
            with open(pdb_file, "r") as filin:
                text = filin.readlines()
            with open("PU_aligned.pdb", "a") as filout:
                for line in text:
                    filout.write(line)
    else:
        dict_pu[nb_pu].remove(score[1])
        remove_aligned_region(input2, str(score[1]))
        tm_align(dict_pu, nb_pu, input2)


def remove_aligned_region(input, idx_pu):
    """
    Remove region of the second protein which has been aligned
    to align others PU with the remaining regions.
    args:
        the file where the aligned region must be removed
        the PU with the best score
    return:

    """
    with open("TM" + idx_pu + ".sup_all_atm", "r") as filin:
        lines = filin.readlines()

    for i, line in enumerate(lines):
        # Find the first line concerning pdb
        if line[0:4] == "ATOM":
            break

    with open("PU_aligned.pdb", "a") as filout:
        # write the pdb PU aligned
        while lines[i][0:4].strip() != "TER":
            filout.write(lines[i])
            i += 1
    i += 1
    idx = []
    # Keep all residues and index to remove
    while lines[i][0:4].strip() != "TER":
        idx.append(lines[i][30:38].strip())
        i += 1
    with open(input, "r") as file_protein:
        text = file_protein.readlines()
    with open(input, "w") as file_protein:
        for line in text:
            # write if the line doesnt contain coordinate for an atom which
            # is in the aligned region to erase
            if line[30:38].strip() not in idx:
                file_protein.write(line)


def main_flex_aln(input1, input2):
    """
    The main loop of the flexible strcutural alignment. For each cutting with
    protein peeling, the optimal TMscore is calculated and keep.
    args:
        the two proteins to align
    return:
        the dictionary containing TMscore for each cutting
    """
    os.system("bin/./dssp " + input1 + " > data/input1.dss")

    # Protein peeling
    cmd_line = ("bin/./peeling11_4.1 -pdb " + input1 +
                " -dssp data/input1.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 \
                -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16")
    output = call_executabe(cmd_line).split("\n")
    dict_pu = parse_protein_peeling(output)

    with open(input2, "r") as filin:
        protein_global = filin.readlines()
    # File the rest of the second protein to be aligned
    input_erasable = "results/input_erasable.pdb"
    # File containing PU which habe been aligned
    pu_aligned = "PU_aligned.pdb"
    dict_tmscore = {}
    for nb_pu in dict_pu:
        # Initialize the file with all atoms for the second protein
        with open(input_erasable, "w") as filout:
            for line in protein_global:
                filout.write(line)
        tm_align(copy.deepcopy(dict_pu), nb_pu, input_erasable)
        # TMscore between the second protein and PU aligned
        cmd_line = ("bin/./TMscore  " + pu_aligned + " " +
                    input2 + " -o TM.sup")
        output = call_executabe(cmd_line)
        dict_tmscore[nb_pu] = parse_tmscore(output, "TMscore")
        # Erase the content of the file containing the PU aligned
        with open(pu_aligned, 'r+') as f_pu_aligned:
            f_pu_aligned.truncate(0)

    # Remove uncessary files created
    os.system("rm TM*")
    os.system("rm results/*.pdb")
    os.system("rm " + pu_aligned)
    # Display the score by the number of PU
    for clef in dict_tmscore:
        print("Number of PUs : {:>2d}  \
               Score : {}".format(clef, dict_tmscore[clef]))
    return dict_tmscore


def display_plot(DICT_TMSCORE_1_2, DICT_TMSCORE_2_1, name_1, name_2):

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(list(DICT_TMSCORE_1_2.keys()), list(DICT_TMSCORE_1_2.values()),
               linestyle = "--", marker = "o")
    ax[0].set_title("{} vs {}".format(name_1, name_2))
    ax[0].set_xlabel("Number of PUs")
    ax[0].set_ylabel("TM score")

    ax[1].plot(list(DICT_TMSCORE_2_1.keys()), list(DICT_TMSCORE_2_1.values()),
               linestyle = "--", marker = "o")
    ax[1].set_title("{} vs {}".format(name_2, name_1))
    ax[1].set_xlabel("Number of PUs")
    ax[1].set_ylabel("TM score")
    plt.subplots_adjust(hspace=1)
    plt.show()
