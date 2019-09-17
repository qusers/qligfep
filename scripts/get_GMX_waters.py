import re

def keep_waters_sphere(pdb, center, r, include=('ATOM','HETATM')): 
    """
    Generate water sphere around provided center coordinates and keep waters that are 
    within the sphere.
    Returns:
        List of waters that should be kept
    """

    xyz_ok= 0
    watname = ['SOL', 'HOH']
    water = []
    water_sphere = []

    #Search for water atoms that are within the sphere radius
    with open(pdb, 'r') as infile:
        for line in infile:
            if line.startswith(include) and line[17:20] in watname:
                water.append(line)
                if ((float(line[30:38])-center[0])**2 + \
                    (float(line[38:46])-center[1])**2 + \
                    (float(line[46:54])-center[2])**2)  \
                    < r**2: 
                    xyz_ok = 1
            if line[12:16] == ' HW2':
                if xyz_ok == 1:
                    water_sphere.append(water)
                water = []
                xyz_ok = 0
                    
    return water_sphere

def replace(string, replacements):
    pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string

def write_waters(waters, outfile, radius):
    # replacements for Q atomtyping
    replacements = {'OW':'O ', 'HW1':'H1 ', 'HW2':'H2 ', 'SOL':'HOH'}
    resid = 0
    atid = 0
    with open(outfile, 'w') as outfile:
        for water in waters:
            resid += 1
            for atom in water:
                atom = replace(atom, replacements)
                atid += 1
                outfile.write('{}{:>5}{}{:>4}{}'.format(atom[0:6], atid, atom[11:22], resid, atom[26:]))
#2YDV-GPCR
#infile  = '2YDV-GPCR_pymemdyn-result.pdb'         # Change filename to gromacs .pdb file -- pymemdyn output file
#outfile = 'water_sphere-2YDV-GPCR.pdb'     # Name of the outputfile for the water sphere used for Q
#center  = [16.477,49.869,43.450] # list of coordinates 2YDV-GPCR

#5g53-GPCR
#infile  = '5g53-GPCR_pymemdyn-result.pdb'         # Change filename to gromacs .pdb file
#outfile = 'water_sphere-5g53-GPCR.pdb'     # Name of the outputfile for the water sphere used for Q
#center  = [20.312,47.776,37.917] # list of coordinates 5g53-GPCR

#5g53-GPCR-GS
infile  = '4EIY-GPCR_gromacs-out.pdb'         # Change filename to gromacs .pdb file
outfile = 'water.pdb'     # Name of the outputfile for the water sphere used for Q
center  = [25.734,30.710,39.767] # list of coordinates 5g53-GPCR-GSalpha

radius  = 25                     # Spheresize for Q

# Get waters and write them out
waters =  keep_waters_sphere(infile, center, radius)
write_waters(waters, outfile, radius)
