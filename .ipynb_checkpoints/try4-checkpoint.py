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

with open('/home/daan/tiberius_data/fep/dev/data/ISO1_GIFT1215_1_prot.pdb') as infile, \
 open('/home/daan/tiberius_data/fep/dev/data/test.pdb', 'w') as outfile:
    for line in infile:
        #if line.startswith(include):
        line = IO.pdb_parse_in(line)
        outline = IO.pdb_parse_out(line)
        outfile.write(outline  + '\n')
        #outfile.write(line)
        
               


coord1 = 0.919,0.61,-2.808

#coord2 = 0.920,0.75,-2.904
distance= 40


def euclidian_overlap(coord1, coord2, distance):
    """
    Calculates whether two points in space overlap within a certain distance
    Returns:
        Boolean
    """
    if ((coord1[0]-coord2[0])**2 + 
        (coord1[1]-coord2[1])**2 + 
        (coord1[2]-coord2[2])**2) < distance**2:
        return True
        

    else:
        return False




if euclidian_overlap(coord1, coord2, distance) == True:
    print("it works")
else: print("it also works but false")


