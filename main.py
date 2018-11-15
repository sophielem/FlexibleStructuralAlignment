#!/usr/bin/env python3

import docopt
import os
import src.flexible_alignment as flex

if __name__ == '__main__':
    # os.system("dssp data/1aoh.pdb > data/1aoh.dss")

    # Protein peeling
    cmd_line = ("bin/./peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16")
    output = flex.call_executabe(cmd_line).split("\n")
    flex.parse_protein_peeling(output)

    # TMalign
    cmd_line = ("bin/./TMalign data/1aoh.pdb data/1aoj.pdb -o TM.sup")
    output = flex.call_executabe(cmd_line)
    sc_tmalign = flex.parse_tmscore(output, "TMalign")

    #TMscore
    cmd_line = ("bin/./TMscore data/1aoh.pdb data/1aoj.pdb -o TM.sup")
    output = flex.call_executabe(cmd_line)
    sc_tmscore = flex.parse_tmscore(output, "TMscore")
