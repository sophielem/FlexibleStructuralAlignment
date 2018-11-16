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

import docopt
import os
import src.flexible_alignment as flex

if __name__ == '__main__':
    arg = docopt.docopt(__doc__, version='0.1')
    input1 = arg["--input"]
    input2 = arg["--input2"]
    if not(os.path.isfile(input1) and os.path.isfile(input2)):
        print("The file doesn't exist !")
        sys.exit(2)

    os.system("dssp " + input1 + " > data/input1.dss")

    # Protein peeling
    cmd_line = ("bin/./peeling11_4.1 -pdb " + input1 + " -dssp data/input1.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16")
    output = flex.call_executabe(cmd_line).split("\n")
    flex.parse_protein_peeling(output)

    # TMalign
    cmd_line = ("bin/./TMalign " + input1 + " " + input2 + " -o TM.sup")
    output = flex.call_executabe(cmd_line)
    sc_tmalign = flex.parse_tmscore(output, "TMalign")

    #TMscore
    cmd_line = ("bin/./TMscore " + input1 + " " + input2 + " -o TM.sup")
    output = flex.call_executabe(cmd_line)
    sc_tmscore = flex.parse_tmscore(output, "TMscore")
