#!/usr/bin/env python3

import re
import os
import parse_pdb as parse


import subprocess as sub
import shlex as shx

os.system("dssp data/1aoh.pdb > data/1aoh.dss")

cmd_line = ("bin/./peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss -R2 95 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 -oss 1 -p 0 -cp 0 -npu 16")
# The split function of the shlex module is used to generate a list of args:
out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()

# To convert the bytes output of suprocess into a str:
str_out = out.decode().split("\n")

str_out = str_out[15 : ]
str_out = str_out[ : -1]
str_out = [itx.split()[4: ] for itx in str_out]
dPDB = parse.parsePDBMultiChains("data/1aoh.pdb")
for itx in str_out:
    nb_pu = int(itx[0])
    for i in range(1, nb_pu*2, 2):
        start = int(itx[i])
        end = int(itx[i + 1])
        parse.writePDB(dPDB, ("data/" + str(nb_pu) + "_" + str(i) + ".pdb"), True, start, end)


cmd_line = ("bin/./TMalign data/1aoh.pdb data/1aoj.pdb -o TM.sup")
# The split function of the shlex module is used to generate a list of args:
out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()

# To convert the bytes output of suprocess into a str:
output = out.decode()

m = re.search("TM-score= (?P<score>[0-9]*\.[0-9]*)", output)
print("TMalign Score normalized by length of the first chain {}".format(m.group("score")))


cmd_line = ("bin/./TMscore data/1aoh.pdb data/1aoj.pdb -o TM.sup")
# The split function of the shlex module is used to generate a list of args:
out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()

# To convert the bytes output of suprocess into a str:
output = out.decode()

m = re.search("TM-score *= (?P<score>[0-9]*\.[0-9]*)", output)
print("TMscore Score normalized by length of the second chain {}".format(m.group("score")))
