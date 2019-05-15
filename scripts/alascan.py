import settings as s
from subprocess import check_output
import shlex
import math
import re
import sys
import os
import argparse

import IO

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig, prot, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.lig = lig
        self.prot = prot

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='ala-scan',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate QresFEP inputscript for Ala scan == ')

    
    parser.add_argument('-l', '--ligand',
                        dest = "lig",
                        required = True,
                        help = "name of the ligand")
    
    parser.add_argument('-p', '--protein',
                        dest = "prot",
                        required = True,
                        help = "name of the protein")
    
    parser.add_argument('-log', '--logfile',
                        dest = "log",
                        default = True,
                        help = "set to t")
    
    args = parser.parse_args()
    run = Run(lig = args.lig,
              prot = args.prot,
             )
    
    