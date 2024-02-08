#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:40:43 2020

"""
import argparse
import math
from operator import itemgetter, attrgetter
import fileinput
import sys
import os
import os.path


#parsers for input file, center residue and radius
parser = argparse.ArgumentParser(description = "neutralize sphere")
parser.add_argument("-p",dest= "prot", type = str ,required = True, help = "protein PDB file")
parser.add_argument("-c",dest= "center",type = str, required = True, help = "center residue and chain, RESN:chain")
parser.add_argument("-r",dest= "radius", required = True, type= int, help = "sphere radius")
parser.add_argument("-cc",dest= "change_charge", required = True, type = int, help = "units of charge to change in the sphere")
arg = parser.parse_args()
pdbfile = arg.prot
center_coordinates = arg.center
center_coordinates = center_coordinates.strip('[]').split(':')
sphere_radius = arg.radius
change_charge = arg.change_charge


#input_pdb = open("/Users/k8/Documents/5n2s-gromacs.pdb", "r")
with open(pdbfile, "r") as input_pdb:
    pos_aa = open("pos_aa.txt", "w")
    neg_aa = open("neg_aa.txt", "w")
    for res in input_pdb:
        if res[17:20] == "LYS":
            pos_aa.write(res)
        elif res[17:20] == "ARG":
            pos_aa.write(res)
        elif res[17:20] == "GLU":
            neg_aa.write(res)
        elif res[17:20] == "ASP":
            neg_aa.write(res)
input_pdb.close()
neg_aa.close()
pos_aa.close()

coord1 = []
coord2 = []
with open("pos_aa.txt", "r") as pos_in:
    for line1 in pos_in:
        if line1[13:16] == "NH1" or line1[13:16] == "NZ ":
            coord1.append([line1[17:20]+line1[23:26],[float(line1[32:38]),
                            float(line1[40:46]),
                            float(line1[47:54])
                           ]])        
with open("neg_aa.txt", "r") as neg_in:
    for line2 in neg_in:
        if line2[13:16] == "OE2" or line2[13:16] == "OD2":
            coord2.append([line2[17:20]+line2[23:26],[float(line2[32:38]),
                            float(line2[40:46]),
                            float(line2[47:54])
                           ]])
# Make list of charged residues that participate in a salt bridge
salt_bridge = []
for n in range(0, len(coord2)):
    for j in range(0,len(coord1)):
        dist = math.sqrt(
            ((coord2[n][1][0] - coord1[j][1][0])**2) + ((coord2[n][1][1] - coord1[j][1][1])**2) + ((coord2[n][1][2] - coord1[j][1][2])**2))
        if dist <= 4:
            salt_bridge.append(coord1[j][0])
            salt_bridge.append(coord2[n][0])
neg_in.close()
pos_in.close()

res_salt_bridge = []
with open(pdbfile, "r") as inp:
    for s in inp:
        if s[17:20]+s[23:26] in salt_bridge:
            res_salt_bridge.append(s[17:20]+s[23:26])
inp.close()
res_salt_bridge = list(dict.fromkeys(res_salt_bridge))

# Make list of charged residues within sphere radius
too_close = []
with open("pos_aa.txt", "r") as charged_pos:
    distance_pos = {}
    for charged in charged_pos:
        distance = (math.sqrt(((float(charged[32:38]) - float(center_coordinates[0])) ** 2) + ((float(charged[40:46]) - float(center_coordinates[1])) ** 2) + ((float(charged[48:54]) - float(center_coordinates[2])) ** 2)))
        if charged[13:16] == "NH1" or charged[13:16] == 'NZ ':
            if distance < sphere_radius:
                distance_pos[charged[17:20] + (charged[23:26])] = distance
                if distance < 10:
                    too_close.append(charged[17:20] + charged[23:26])

distance_pos_l = []
for key, value in distance_pos.items():
    temp = [value,key]
    distance_pos_l.append(temp)

with open("neg_aa.txt", "r") as charged_neg:
    distance_neg = {}
    for charged2 in charged_neg:
        if charged2[0:3] != "TER":
            distance = (math.sqrt(((float(charged2[32:38]) - float(center_coordinates[0])) ** 2) + ((float(charged2[40:46]) - float(center_coordinates[1])) ** 2) + ((float(charged2[48:54]) - float(center_coordinates[2])) ** 2)))
            if charged2[13:16] == "OD2" or charged2[13:16] == "OE2":
                if distance < (sphere_radius):
                    distance_neg[charged2[17:20] + (charged2[23:26])] = distance
                    if distance < 10:
                        too_close.append(charged2[17:20] + charged2[23:26])

distance_neg_l = []
for key, value in distance_neg.items():
    temp2 = [value,key]
    distance_neg_l.append(temp2)

# Sort distance list per charged atom in charged residues on distance furthest to closest from the sphere center
sorted_neg = sorted(distance_neg_l, key = itemgetter(0), reverse = True)
sorted_pos = sorted(distance_pos_l, key = itemgetter(0), reverse = True)

# Define top n residues
n = abs(change_charge)
residues_pos = [i[1] for i in sorted_pos]
for ref in residues_pos:
    if ref in too_close:
        residues_pos.remove(ref)
    elif ref in res_salt_bridge:
        residues_pos.remove(ref)
        
residues_neg = [i[1] for i in sorted_neg]
for re in residues_neg:
    if re in too_close:
        residues_neg.remove(re)
    elif re in res_salt_bridge:
        residues_neg.remove(re)
residues_to_decharge_pos = residues_pos[:n]
residues_to_decharge_neg = residues_neg[:n]

# Write new pdb file with decharged residues
def write_pdb(list1):
    with open(pdbfile, "r") as input_pdb:
        with open("output.pdb","w") as output:
            for d in input_pdb:
                if change_charge < 0:
                    if (d[17:20]+d[23:26]) in list1:
                        if d[13:14] != "H":
                            output.write(d[0:19]+"N "+d[21:78]+"\n")
                    elif d[17:20] != "SOL":
                            output.write(d)
                    else:
                        output.write(d)
                elif change_charge > 0:
                    if (d[17:20]+d[23:26]) in list1:
                        if d[13:14] != "H":
                            output.write(d[0:19]+"H "+d[21:78]+"\n")
                    elif d[17:20] != "SOL":
                            output.write(d)
                    else:
                        output.write(d)
    output.close()
    input_pdb.close()

# Clean up files
os.remove("neg_aa.txt")
os.remove("pos_aa.txt")

# Ddefine which residues will be decharged
print("Decharged residues   r")
if change_charge > 0:
    for res in residues_to_decharge_neg:
        for item in sorted_neg:
            if res in item:
                print(item[1],'             ',round(item[0],3))
    list1 = residues_to_decharge_neg
elif change_charge < 0:
    for res in residues_to_decharge_pos:
        for item in sorted_pos:
            if res in item:
                print(item[1],'             ',round(item[0],3))
    list1 = residues_to_decharge_pos

write_pdb(list1)