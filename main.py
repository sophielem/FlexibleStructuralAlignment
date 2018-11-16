#!/usr/bin/env python3

"""Flexible Structure Alignment

Usage:
  main.py -i <file> -p <file2>

Options:
  -h --help                  help
  --version                  version of the script
  -i --input = file          input pdb
  -p --input2 = file2        input2 pdb
"""

import os
import sys
import docopt
import src.flexible_alignment as flex


if __name__ == '__main__':
    ARG = docopt.docopt(__doc__, version='0.1')
    INPUT1 = ARG["--input"]
    INPUT2 = ARG["--input2"]
    if not(os.path.isfile(INPUT1) and os.path.isfile(INPUT2)):
        print("The file doesn't exist !")
        sys.exit(2)

    os.system("dssp " + INPUT1 + " > data/input1.dss")

    # Protein peeling
    CMD_LINE = ("bin/./peeling11_4.1 -pdb " + INPUT1 +
                " -dssp data/input1.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 \
                -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16")
    OUTPUT = flex.call_executabe(CMD_LINE).split("\n")
    flex.parse_protein_peeling(OUTPUT)

    # TMalign
    CMD_LINE = ("bin/./TMalign " + INPUT1 + " " + INPUT2 + " -o TM.sup")
    OUTPUT = flex.call_executabe(CMD_LINE)
    SC_TMALIGN = flex.parse_tmscore(OUTPUT, "TMalign")

    #TMscore
    CMD_LINE = ("bin/./TMscore " + INPUT1 + " " + INPUT2 + " -o TM.sup")
    OUTPUT = flex.call_executabe(CMD_LINE)
    SC_TMSCORE = flex.parse_tmscore(OUTPUT, "TMscore")
