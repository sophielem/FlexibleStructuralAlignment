#!/usr/bin/env python3

import os
import csv

if __name__ == '__main__':
    all_reslt = [["Names", "1_2", "2_1", "TMalign", "parMATT"]]
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
            with open("results/" + line_splitted[1][1:-2] + line_splitted[2][1:-2] + "res.txt", "r") as csv_file:
                text = csv.reader(csv_file)
                text = list(text)
            all_reslt.append(text[1])
        # write all the result
        with open("results/bench_res.txt", "w") as res_file:
            writer = csv.writer(res_file)
            writer.writerows(all_reslt)
        res_file.close()
