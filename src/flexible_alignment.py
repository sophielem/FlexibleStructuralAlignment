#!/usr/bin/env python3

"""
"""

import subprocess as sub
import shlex as shx
import re
import os
import glob
import src.parse_pdb as parse


def call_executabe(cmd_line):
    """
    """
    # The split function of the shlex module is used to generate a list of args:
    out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()

    # to remove all uncessary files created
    for filename in glob.glob("TM.*"):
        os.remove(filename)
    for filename in glob.glob("file*"):
        os.remove(filename)
    # To convert the bytes output of suprocess into a str:
    return out.decode()


def parse_protein_peeling(output):
    """
    """
    output = output[15 : ]
    output = output[ : -1]
    output = [itx.split()[4: ] for itx in output]
    dPDB = parse.parsePDBMultiChains("data/1aoh.pdb")
    for itx in output:
        nb_pu = int(itx[0])
        for i in range(1, nb_pu*2, 2):
            start = int(itx[i])
            end = int(itx[i + 1])
            parse.writePDB(dPDB, ("results/" + str(nb_pu) + "_" + str(i) + ".pdb"), True, start, end)


def parse_tmscore(output, soft):
    """
    """
    m = re.search("TM-score *= (?P<score>[0-9]*\.[0-9]*)", output)
    print("{} Score normalized by length of the first chain {}".format(soft, m.group("score")))
    return m.group("score")
