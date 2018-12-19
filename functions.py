import IO
import numpy as np

def sigmoid(t, k):
    return ( k * t / (k -t + 1.0))

def COG(pdbfile, include='ATOM,HETATM'):
    """
    Calculates center of geometry of a protein and/or ligand structure.
    Returns:
        center (list): List of float coordinates [x,y,z] that represent the
        center of geometry (precision 3).
    """

    center = [None, None, None]
    include = tuple(include.split(','))

    with open(pdbfile) as pdb:

        # extract coordinates [ [x1,y1,z1], [x2,y2,z2], ... ]
        coordinates = []
        for line in pdb:
            if line.startswith(include):
                coordinates.append([float(line[30:38]),    # x_coord
                                    float(line[38:46]),    # y_coord
                                    float(line[46:54])     # z_coord
                                   ])

        # calculate center of geometry
        center = [sum([coordinates[i][j]/(len(coordinates))
              for i in range(len(coordinates))]) for j in range(3)]
        center = [round(center[i], 3) for i in range(3)]
    return center

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

def overlapping_pairs(pdbfile, reslist, include=('ATOM', 'HETATM')):
    """
    Calculates whether input pdb has overlaying atoms, based on provided residue names
    Returns:
        dictionary of overlaying atoms based on atom number
    """
    coordinates = []
    overlapping_atoms = []
    atomlist = []
    index = 0
    # Parse the input pdbfile
    with open(pdbfile) as infile:
        for line in infile:
            if line.startswith(include):
                line_parse = IO.pdb_parse_in(line)
                if line_parse[4] in reslist:
                    coordinates.append([line_parse[1],
                                        line_parse[8],
                                        line_parse[9],
                                        line_parse[10],
                                        line_parse[13]
                                       ])
    for at1 in coordinates:
        for at2 in coordinates:
            if at1[0] != at2[0]:
                if ((at1[1]-at2[1])**2 + 
                    (at1[2]-at2[2])**2 + 
                    (at1[3]-at2[3])**2) < 0.8:
                    if at1[4] == at2[4] and at1[4].strip() != 'H':
                        overlapping_atoms.append([at1[0], at2[0]])
                        
    total = len(overlapping_atoms)
    for i in range (0, (total/2)):
        atomlist.append(overlapping_atoms[i])
                        
    return atomlist

def kT(T):
    k = 0.0019872041 # kcal/(mol.K)
    kT = k * T
    kT = '{:.3f}'.format(kT)
    return kT

def avg_sem(array):
    FEP_sum = array.sum(axis = 0)
    dG = np.nanmean(FEP_sum)
    cnt = len(FEP_sum)
    sem = np.nanstd(FEP_sum, ddof =1)/np.sqrt(cnt)

    return [dG, sem]
