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
import re
import docopt
import src.flexible_alignment as flex


if __name__ == '__main__':
    ARG = docopt.docopt(__doc__, version='0.1')
    INPUT1 = ARG["--input"]
    INPUT2 = ARG["--input2"]
    if not(os.path.isfile(INPUT1) and os.path.isfile(INPUT2)):
        print("The file doesn't exist !")
        sys.exit(2)

    name_1=re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT1).group("name")
    name_2=re.search("[.*/|](?P<name>[0-9]*.*)\.pdb", INPUT2).group("name")

    # INPUT1 vs INPUT2
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_1,
                                                      name_2, "*"*16))
    DICT_TMSCORE_1_2 = flex.main_flex_aln(INPUT1, INPUT2)
    # INPUT2 vs INPUT1
    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_2,
                                                      name_1, "*"*16))
    DICT_TMSCORE_2_1 = flex.main_flex_aln(INPUT2, INPUT1)

    flex.display_plot(DICT_TMSCORE_1_2, DICT_TMSCORE_2_1, name_1, name_2)
