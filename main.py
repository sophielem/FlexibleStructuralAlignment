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
import matplotlib.pyplot as plt
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

    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_1, name_2, "*"*16))
    DICT_TMSCORE_1_2 = flex.main_flex_aln(INPUT1, INPUT2)

    print("\n\t\t{}\n\t\t* {} vs {} *\n\t\t{}".format("*"*16, name_2, name_1, "*"*16))
    DICT_TMSCORE_2_1 = flex.main_flex_aln(INPUT2, INPUT1)

    names = list(DICT_TMSCORE_1_2.keys())
    values = list(DICT_TMSCORE_1_2.values())

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(list(DICT_TMSCORE_1_2.keys()), list(DICT_TMSCORE_1_2.values()), linestyle = "--", marker = "o")
    ax[0].set_title("{} vs {}".format(name_1, name_2))
    ax[0].set_xlabel("Number of PUs")
    ax[0].set_ylabel("TM score")

    ax[1].plot(list(DICT_TMSCORE_2_1.keys()), list(DICT_TMSCORE_2_1.values()), linestyle = "--", marker = "o")
    ax[1].set_title("{} vs {}".format(name_2, name_1))
    ax[1].set_xlabel("Number of PUs")
    ax[1].set_ylabel("TM score")
    plt.subplots_adjust(hspace=1)
    plt.show()
