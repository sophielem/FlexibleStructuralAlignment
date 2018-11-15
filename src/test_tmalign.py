#!/usr/bin/env python3

import re
import os
import parse_pdb as parse


os.system("bin/./TMalign data/1aoh.pdb data/1aoj.pdb -o TM.sup > results/tmalign_output.txt")
with open("results/tmalign_output.txt", "r") as filin:
    text = filin.read().splitlines()
text = "".join(text)

m = re.search("TM-score= (?P<score>[0-9]*\.[0-9]*)", text)
print("TMalign Score normalized by length of the first chain {}".format(m.group("score")))


os.system("bin/./TMscore data/1aoh.pdb data/1aoj.pdb -o TM.sup > results/tmscore_output.txt")
with open("results/tmscore_output.txt", "r") as filin:
    text = filin.read().splitlines()
text = "".join(text)

m = re.search("TM-score *= (?P<score>[0-9]*\.[0-9]*)", text)
print("TMscore Score normalized by length of the second chain {}".format(m.group("score")))

os.system("dssp data/1aoh.pdb > data/1aoh.dss")
os.system("bin/./peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16 > results/peeling_output.txt")
with open("results/peeling_output.txt", "r") as filin:
    text = filin.read().splitlines()
text = text[15 : ]
text = [itx.split()[4: ] for itx in text]
dPDB = parse.parsePDBMultiChains("data/1aoh.pdb")
for itx in text:
    nb_pu = int(itx[0])
    for i in range(1, nb_pu*2, 2):
        start = int(itx[i])
        end = int(itx[i + 1])
        parse.writePDB(dPDB, ("data/" + str(nb_pu) + "_" + str(i) + ".pdb"), True, start, end)
break
