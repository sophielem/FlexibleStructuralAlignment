#!/usr/bin/env python3

"""Flexible Structure Alignment

Usage:
  main.py -i <file> -c <chain1> -f <file2> -d <chain2>

Options:
  -h --help                  help
  --version                  version of the script
  -i --input = file          input pdb
  -c --chain1 = letter1      chain1
  -f --input2 = file2        input2 pdb
  -d --chain2 = letter2      chain2
"""

import os
import sys
import re
import csv
import docopt
import src.flexible_alignment as flex
import src.parse_pdb as parse
import src.re_number as rnb

if __name__ == '__main__':
    ARG = docopt.docopt(__doc__, version='0.1')
    INPUT1 = ARG["--input"]
    INPUT2 = ARG["--input2"]
    CHAIN1 = ARG["--chain1"]
    CHAIN2 = ARG["--chain2"]
    if not(os.path.isfile(INPUT1) and os.path.isfile(INPUT2)):
        print("The file doesn't exist !")
        sys.exit(2)

    NAME_1 = re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT1).group("name")
    NAME_2 = re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT2).group("name")

    LEN_INPUT1 = len(parse.parse_pdb(INPUT1, CHAIN1, True))
    LEN_INPUT2 = len(parse.parse_pdb(INPUT2, CHAIN2, True))

    INPUT1_LONGER = LEN_INPUT1 > LEN_INPUT2
    rnb.re_number_pdb(1, INPUT1, INPUT1)
    rnb.re_number_pdb(1, INPUT2, INPUT2)
    # INPUT1 vs INPUT2
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, NAME_1,
                                                      NAME_2, "*"*16))
    DICT_TMSCORE_1_2 = flex.main_flex_aln(INPUT1, INPUT2, INPUT1_LONGER,
                                          CHAIN1)
    # INPUT2 vs INPUT1
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, NAME_2,
                                                      NAME_1, "*"*16))
    DICT_TMSCORE_2_1 = flex.main_flex_aln(INPUT2, INPUT1, not(INPUT1_LONGER),
                                          CHAIN2)

    flex.display_plot(DICT_TMSCORE_1_2, DICT_TMSCORE_2_1, NAME_1, NAME_2)

    SC_TMALIGN = flex.tmalign_simple(INPUT1_LONGER, INPUT1, INPUT2)

    SC_PARMATT = flex.parmatt(INPUT1, INPUT2, INPUT1_LONGER)

    MAX_SCORE_1_2 = max(DICT_TMSCORE_1_2.keys(),
                        key=(lambda key: DICT_TMSCORE_1_2[key]))
    MAX_SCORE_2_1 = max(DICT_TMSCORE_2_1.keys(),
                        key=(lambda key: DICT_TMSCORE_2_1[key]))
    LIST_SCORE = [["Names", "1_2", "2_1", "TMalign", "parMATT"],
                  [NAME_1 + "_" + NAME_2, DICT_TMSCORE_1_2[MAX_SCORE_1_2],
                   DICT_TMSCORE_2_1[MAX_SCORE_2_1], SC_TMALIGN, SC_PARMATT]]
    with open("results/" + NAME_1 + "_" + NAME_2 + "res.txt", "w") as res_file:
        WRITER = csv.writer(res_file)
        WRITER.writerows(LIST_SCORE)
    res_file.close()
