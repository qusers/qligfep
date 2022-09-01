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


distance= 5
coord1 = 0.919,0.61,-2.808

with open('/home/daan/tiberius_data/fep/dev/data/ISO1_GIFT1215_1_prot.pdb') as infile, \
     open('/home/daan/tiberius_data/fep/dev/data/try5.pdb', 'w') as outfile:
     for line in infile:
        #if line.startswith(include):
        line = IO.pdb_parse_in(line)
        
        if type(line) == list:
            coord2 = tuple(line[8:11])
            #chain_id = tuple(line[5])
            #print(line)
            to_keep = {}
            
                    #print(line)
            if f.euclidian_overlap(coord1, coord2, distance) == True:
                lst = []
                lst.append(line)
                #print(line)
                if line[5] not in to_keep:
                    to_keep[line[5]]=[]
                to_keep[line[5]].append(line[6])   
                
                #print(to_keep)
                #print(to_keep[line[6]])
            #print(line[6])
            
            #print(line)   
            end_list=[]
            if line[5] not in to_keep:
                    to_keep[line[5]]=[]
            
            if [line[6]]==to_keep[line[5]]:
                #print(line[6])
                #print(to_keep[line[5]])
                end_list.append(line)
                print(end_list)
                    

            
            # then back to string for writing
              #  outline = IO.pdb_parse_out(line)

            #else:
             #   outline = line

        # block of code that checks whether it is in

#this writes out the pdb file that gets created
           # outfile.write(outline  + '\n')