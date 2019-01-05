#!/usr/bin/env python3

"""Flexible Structure Alignment

Usage:
  main.py -i <file> -f <file2>

Options:
  -h --help                  help
  --version                  version of the script
  -i --input = file          input pdb
  -f --input2 = file2        input2 pdb
"""

import os
import sys
import re
import docopt
import src.flexible_alignment as flex
import src.parse_pdb as parse


if __name__ == '__main__':
    ARG = docopt.docopt(__doc__, version='0.1')
    INPUT1 = ARG["--input"]
    INPUT2 = ARG["--input2"]
    if not(os.path.isfile(INPUT1) and os.path.isfile(INPUT2)):
        print("The file doesn't exist !")
        sys.exit(2)

    name_1=re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT1).group("name")
    name_2=re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT2).group("name")

    LEN_INPUT1 = len(parse.parse_pdb(INPUT1, 'A'))
    LEN_INPUT2 = len(parse.parse_pdb(INPUT2, 'A'))

    TEST = LEN_INPUT1 > LEN_INPUT2
    # INPUT1 vs INPUT2
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_1,
                                                      name_2, "*"*16))
    DICT_TMSCORE_1_2 = flex.main_flex_aln(INPUT1, INPUT2, TEST)
    # INPUT2 vs INPUT1
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_2,
                                                      name_1, "*"*16))
    DICT_TMSCORE_2_1 = flex.main_flex_aln(INPUT2, INPUT1, not(TEST))

    flex.display_plot(DICT_TMSCORE_1_2, DICT_TMSCORE_2_1, name_1, name_2)


    if TEST:
        cmd_line = ("bin/./TMalign " + INPUT1 + " " +  INPUT2 +
                    " -o TM.sup")
    else:
        cmd_line = ("bin/./TMalign " + INPUT2 + " " +  INPUT1 +
                    " -o TM.sup")
    output = flex.call_executabe(cmd_line)
    # Retrieve the TM score
    sc_tmalign = flex.parse_tmscore(output, "TMalign")
    print("\n\t\t  {}\n\t\t  * TMalign *\n\t\t  {}".format("*"*11, "*"*11))
    print("Score : {}".format(sc_tmalign))

    cmd_line = "./bin/parMatt " + INPUT1 + " " + INPUT2 + " -o test"
    flex.call_executabe(cmd_line)
    dict_1 = parse.parse_pdb("test.pdb", 'A')
    dict_2 = parse.parse_pdb("test.pdb", 'B')
    parse.write_pdb(dict_1, "input1.pdb", True, 1, LEN_INPUT1)
    parse.write_pdb(dict_2, "input2.pdb", True, 1, LEN_INPUT2)
    if TEST:
        cmd_line = ("bin/./TMscore input2.pdb input1.pdb -o TM.sup")
    else:
        cmd_line = ("bin/./TMscore input1.pdb input2.pdb -o TM.sup")
    output = flex.call_executabe(cmd_line)
    # Retrieve the TM score
    sc_tmalign = flex.parse_tmscore(output, "TMscore")
    print("\n\t\t  {}\n\t\t  * parMATT *\n\t\t  {}".format("*"*11, "*"*11))
    print("Score : {}".format(sc_tmalign))
