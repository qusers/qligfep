import re
import shlex
from subprocess import check_output
import os
import stat
from xml.etree.ElementInclude import include
import numpy as np

import functions as f
import settings as s
import IO
from itertools import chain


distance= 40
coord1 = 0.919,0.61,-2.808

with open('/home/daan/tiberius_data/fep/dev/data/ISO1_GIFT1215_1_prot.pdb') as infile, \
     open('/home/daan/tiberius_data/fep/dev/data/test.pdb', 'w') as outfile:
    for line in infile:
        #if line.startswith(include):
        line = IO.pdb_parse_in(line)
        if type(line) == list:
            coord2 = tuple(line[8:11])
            #print(coord2, f.euclidian_overlap(coord1,coord2,distance))

            # then back to string for writing
            outline = IO.pdb_parse_out(line)

        else:
            outline = line

        # block of code that checks whether it is in

#this writes out the pdb file that gets created
        outfile.write(outline  + '\n')


with open('/home/daan/tiberius_data/fep/dev/data/ISO1_GIFT1215_1_prot.pdb') as infile, \
     open('/home/daan/tiberius_data/fep/dev/data/test_resnum.pdb', 'w') as outfile: 
    for line in infile:
        #if line.startswith(include):
        line = IO.pdb_parse_in(line)
        if type(line) == list:
            res_coord = [line[6]] + line[8:11]
            #coord2 = chain(line[6], line[8:11])
            
            #print(coord2, f.euclidian_overlap(coord1,coord2,distance))
            print(res_coord)

            # then back to string for writing
            outline = IO.pdb_parse_out(line)

        else:
            outline = line

        # block of code that checks whether it is in

#this writes out the pdb file that gets created
        outfile.write(outline  + '\n')

