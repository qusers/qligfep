import argparse

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, ipdb, *args, **kwargs):
        self.ipdb = ipdb
        
    def COG(self, include='ATOM,HETATM'):
        """
        Calculates center of geometry of a protein and/or ligand structure.
        Returns:
            center (list): List of float coordinates [x,y,z] that represent the
            center of geometry (precision 3).
        """

        center = [None, None, None]
        include = tuple(include.split(','))

        with open(self.ipdb) as pdb:

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
            center = [str(round(center[i], 3)) for i in range(3)]
            center = ':'.join(center)
        return center

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='calcCOG',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Calculate Centre of Geometry (COG) for a given .pdb == ')


    parser.add_argument('-i', '--ipdb',
                        dest = "ipdb",
                        required = True,
                        help = "name of pdbfile")

    args = parser.parse_args()
    run = Run(ipdb = args.ipdb)
    print run.COG()
