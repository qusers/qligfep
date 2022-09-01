from calendar import c
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
to_keep = {}
end_list=[]

with open('/home/daan/tiberius_data/fep/dev/data/ISO1_GIFT1215_1_prot.pdb') as infile, \
     open('/home/daan/tiberius_data/fep/dev/data/try5.pdb', 'w') as outfile:
     for line in infile:
        #if line.startswith(include):
        line = IO.pdb_parse_in(line)
        
        if type(line) == list:
            coord2 = tuple(line[8:11])
            chain_id = line[5]
            resnum = line[6]
            combo = line[5:7]
        
            
                    #print(line)
            if f.euclidian_overlap(coord1, coord2, distance) == True:
                #lst = []
                #lst.append(line)
                #print(lst)
                #for lst[6] in line[5]:
                 #   end_list.append[line]
                  #  print(end_list)
                if chain_id not in to_keep:
                    to_keep[chain_id]=[]
                to_keep[chain_id].append(resnum)
                #print(to_keep)   
                
                    
                #mylist = list(to_keep)
                #print(mylist)
                #print(line) 
                #for mylist in chain_id:
                 #   end_list.append(line)  
                  #  print(end_list)
               # print(to_keep)
            
            #print(to_keep)
            if chain_id not in to_keep:
                to_keep[chain_id]=[]

            #print(line[6])
          #  print(to_keep[chain_id])
            #print(line)
            print([resnum])

            if resnum in to_keep[chain_id]:
                end_list.append(line)
                print(end_list)
                    

            
            # then back to string for writing
              #  outline = IO.pdb_parse_out(line)

            #else:
             #   outline = line

        # block of code that checks whether it is in

#this writes out the pdb file that gets created
           # outfile.write(outline  + '\n')