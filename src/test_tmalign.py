#!/usr/bin/env python3

import re
import os

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
