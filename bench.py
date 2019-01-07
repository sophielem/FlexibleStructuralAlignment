#!/usr/bin/env python3

import os


if __name__ == '__main__':
    with open("data/RIPC_dataset.txt", "r") as bench_file:
        bench_file.readline()
        for line in bench_file:
            line_splitted = line.split("\t")
            prot = "data/" + line_splitted[1][1:-2] + ".pdb"
            chain = line_splitted[1][-2].upper()
            if chain == "_":
                chain = "' '"

            prot2 = "data/" + line_splitted[2][1:-2] + ".pdb"
            chain2 = line_splitted[2][-2].upper()
            if chain2 == "_":
                chain2 = "' '"

            os.system("./main.py -i " + prot + " -c " + chain +
                      " -f " + prot2 + " -d " + chain2)
