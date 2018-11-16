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
    for filename in glob.glob("TM*"):
        os.remove(filename)
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
    for itx in output:
        # Retrieve the number of PU
        nb_pu = int(itx[0])
        # For each delimitations, create the pdb file
        for i in range(1, nb_pu*2, 2):
            start = int(itx[i])
            end = int(itx[i + 1])
            parse.write_pdb(dict_pdb,
                            ("results/" + str(nb_pu) + "_" + str(i) + ".pdb"),
                            True, start, end)


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
    return matches.group("score")
