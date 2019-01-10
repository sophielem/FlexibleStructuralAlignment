#!/usr/bin/env python3

"""
Benchmark from RIPC_dataset text file
"""

import os
import csv


if __name__ == '__main__':
    ALL_RESLT = [["Names", "1_2", "2_1", "TMalign", "parMATT"]]
    with open("data/RIPC_dataset.txt", "r") as bench_file:
        bench_file.readline()
        # For each protein pair
        for line in bench_file:
            line_splitted = line.split("\t")
            # Name and chain of the first protein
            prot = "data/" + line_splitted[1][1:-2] + ".pdb"
            chain = line_splitted[1][-2].upper()
            if chain == "_":
                chain = "' '"
            # Name and chain of the second protein
            prot2 = "data/" + line_splitted[2][1:-2] + ".pdb"
            chain2 = line_splitted[2][-2].upper()
            if chain2 == "_":
                chain2 = "' '"
            # Principal command
            os.system("./main.py -i " + prot + " -c " + chain +
                      " -f " + prot2 + " -d " + chain2)
            # Read the csv result file and append it in a global list
            with open("results/" + line_splitted[1][1:-2] + "_" +
                      line_splitted[2][1:-2] + "res.txt", "r") as csv_file:
                text = csv.reader(csv_file)
                text = list(text)
            ALL_RESLT.append(text[1])
        # write all the result
        with open("results/bench_res.csv", "w") as res_file:
            WRITER = csv.writer(res_file)
            WRITER.writerows(ALL_RESLT)
        res_file.close()

        os.system("Rscript src/plot_bench.R")
