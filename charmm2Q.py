#! /usr/bin/env python

__author__ = "Geir V. Isaksen"
__date__ = "2015"

import sys
import os
from string import maketrans
from copy import deepcopy
import numpy as np
import argparse
import settings as s


class Run(object):
    """
    Convert CGenFF parameter files to Q
    """
    def __init__(self, lig, merge, forcefield):
        self.lig = lig.split('.')[0]
        self.charmm_prm = self.lig + '.prm'
        self.charmm_top = self.lig + '.rtf'
        self.FF = forcefield
        self.merge = merge    
    
    def read_Q_prm(self, filename):
        """
        Reads original Q prm file and collect existing atomtypes for parameters. Returns dictionary.
        This is used to pick out non-existing parameters in charmm parameter file.
        """
        qprm = {'[atom_types]': list(),
                '[bonds]': list(),
                '[angles]': list(),
                '[torsions]': list(),
                '[impropers]': list()}

        sections = {'[bonds]': 2,
                        '[angles]': 3,
                        '[torsions]': 4,
                        '[impropers]': 4,
                        '[atom_types]': 1}

        section = 'SECTION_'
        print(filename)
        with open(filename, 'r') as prm:

            for line in prm:
                #Add atomtypes for current section
                if section in sections.keys():
                    if '[' not in line:
                        if not line.startswith('*') and len(line.split()) > (sections[section]):
                            atomtypes = ' '.join(line.split()[0:sections[section]])
                            if '!' not in atomtypes:
                                qprm[section].append(atomtypes)

                #Check if new section starts (bond/angle/torsion/improper)
                if len(line.split()) > 0:
                    if line.split()[0] in sections.keys():
                        section = line.split()[0]

        return qprm


    def read_charmm_par(self, filename, atom_mass, existing_prm, convert_to_opls):
        """
        Reads a charmmm parmater file and returns a dictionary with parameters suitable for Q.
        If a Q parameter file is specified, non-existing parameters will be returned.
        Else all parameters in file will be returned.
        """

        params = {'[bonds]': list(),
                  '[angles]': list(),
                  '[torsions]': list(),
                  '[impropers]': list(),
                  '[atom_types]': list()}

        selector = {'found_bonds': False,
                    'found_angles': False,
                    'found_torsions': False,
                    'found_impropers': False,
                    'found_atoms': False}



        prev_tors = ''

        scale_14 = '1'

        with open(filename, 'r') as prm:
            for line in prm:
                data = line.split('!')[0]
                try:
                    comment = '!%s' % line.split('!')[1].strip('\n')
                except:
                    comment = ' '
                
                if len(line.split()) > 3:
                    line1 = line.split()
                    if line1[0] == 'MASS':
                        atom_mass[line1[2]] = float(line1[3])
                
                if selector['found_bonds']:
                    if not data.startswith('!') and len(data.split()) > 3:
                        types = ' '.join(data.split()[0:2])
                        if types not in existing_prm['[bonds]']:
                            types = types.translate(maketrans('X', '?'))
                            a1, a2 = types.split()[0], types.split()[1]
                            k = 2 * float(line.split()[2])
                            r = float(line.split()[3])
                            params['[bonds]'].append('%6s %6s %9.3f %9.3f %s' % (a1.ljust(6), a2.ljust(6), k, r, comment))

                if selector['found_angles']:
                    if not data.startswith('!') and len(data.split()) > 4:
                        types = ' '.join(data.split()[0:3])
                        if types not in existing_prm['[angles]']:
                            types = types.translate(maketrans('X', '?'))
                            a1, a2, a3 = types.split()[0], types.split()[1], types.split()[2]
                            k = 2 * float(line.split()[3])
                            theta = float(line.split()[4])
                            if len(data.split()) == 7:
                                Kub = float(line.split()[5])
                                Stheta = float(line.split()[6].translate(None, '!'))

                            else:
                                Kub = 0.00
                                Stheta = 0.0000
                            params['[angles]'].append('%6s %6s %6s %9.2f %8.1f %9.2f %9.2f %s' % (a1.ljust(6), a2.ljust(6),
                                                                                           a3.ljust(6), k, theta, Kub,
                                                                                           Stheta, comment))

                if selector['found_torsions']:
                    if not data.startswith('!') and len(data.split()) > 6:
                        types = ' '.join(data.split()[0:4])
                        try:
                            if types not in existing_prm['[torsions]']:
                                if types == prev_tors:
                                    t1, t2, t3, t4, k, minima, phase = params['[torsions]'][-1].split()[0:7]
                                    minima = -float(minima)
                                    prev_comment = ' '.join(params['[torsions]'][-1].split()[8:])
                                    params['[torsions]'][-1] = '%6s %6s %6s %6s %9.3f %9.3f %9.3f     1 %s' \
                                                                     % (t1.ljust(6), t2.ljust(6), t3.ljust(6), t4.ljust(6),
                                                                        float(k), float(minima), float(phase), prev_comment)
                                types = types.translate(maketrans('X', '?'))
                                a1, a2, a3, a4 = types.split()[0], types.split()[1], types.split()[2], types.split()[3]
                                k, minima, phase = map(float, line.split()[4:7])
                                params['[torsions]'].append('%6s %6s %6s %6s %9.3f %9.3f %9.3f     1 %s' %
                                                            (a1.ljust(6), a2.ljust(6), a3.ljust(6), a4.ljust(6), k, minima,
                                                             phase, comment))
                                prev_tors = types
                        except:
                            pass

                #In the charmm OPLS/M (2015) this section is
                if selector['found_impropers']:
                    if not data.startswith('!') and len(data.split()) > 6:
                        try:
                            types = ' '.join(data.split()[0:4])
                            if types not in existing_prm['[impropers]']:
                                types = types.translate(maketrans('X', '?'))
                                a1, a2, a3, a4 = types.split()[0], types.split()[1], types.split()[2], types.split()[3]
                                k = float(line.split()[4])
                                if not convert_to_opls:
                                    k *= 2.
                                imp = float(line.split()[6])
                                params['[impropers]'].append('%6s %6s %6s %6s %9.3f %9.3f %s' %
                                                             (a1.ljust(6), a2.ljust(6), a3.ljust(6), a4.ljust(6), k, imp,
                                                              comment))
                        except:
                            pass

                if selector['found_atoms']:
                    if not data.startswith('!') and len(data.split()) > 3:
                        if ' e14fac ' in line:
                            scale_14 = line.split('e14fac')[-1].split()[0]
                            continue
                        #try:
                        atomtype = line.split()[0]
                        if atomtype not in existing_prm['[atom_types]']:
                            epsilon = float(line.split()[2])
                            rmin = float(line.split()[3])
                            if len(data.split()) > 6:
                                try:
                                    epsilon_1_4 = float(line.split()[5])
                                    rmin_1_4 = float(line.split()[6])
                                except:
                                    epsilon_1_4 = epsilon
                                    rmin_1_4 = rmin
                            else:
                                epsilon_1_4 = epsilon
                                rmin_1_4 = rmin

                            avdw1, avdw2 = rmin, rmin
                            bvdw1 = abs(epsilon)
                            avdw3 = rmin_1_4
                            bvdw2_3 = abs(epsilon_1_4)
                            mass = atom_mass[atomtype]

                            if convert_to_opls:
                                epsilon = abs(epsilon)
                                sigma = (2.**(5./6.) * rmin)
                                avdw1 = float(np.sqrt(4. * epsilon * sigma**12.))
                                avdw2 = avdw1
                                avdw3 = np.sqrt(0.5) * avdw1

                                bvdw1 = np.sqrt(4 * epsilon * sigma**6.)
                                bvdw2_3 = np.sqrt(0.5) * bvdw1

                            params['[atom_types]'].append('%6s %10.4f %10.4f %9.4f %10.4f %9.4f %11.4f 1 %s' %
                                                          (atomtype.ljust(6), avdw1, avdw2, bvdw1,
                                                           avdw3, bvdw2_3, mass, comment))
                        #except:
                        #    pass

                if not line.startswith('!'):
                    if 'BOND' in line:
                        print('Found BONDS')
                        for section in selector.keys():
                            selector[section] = False
                        selector['found_bonds'] = True

                    if 'ANGLE' in line:
                        print('Found ANGLES')
                        for section in selector.keys():
                            selector[section] = False
                        selector['found_angles'] = True

                    if 'DIHEDRAL' in line:
                        print('Found DIHEDRALS')
                        for section in selector.keys():
                            selector[section] = False
                        selector['found_torsions'] = True

                    if 'IMPROPER' in line:
                        print('Found IMPROPER')
                        for section in selector.keys():
                            selector[section] = False
                        selector['found_impropers'] = True

                    if 'CMAP' in line or 'HBOND' in line:
                        for section in selector.keys():
                            selector[section] = False

                    if 'NONBONDED' in line:
                        print('Found NONBONDED')
                        for section in selector.keys():
                            selector[section] = False
                        selector['found_atoms'] = True

        return params, scale_14


    def find_ff_name(self, filename):
        """
        Look in the header of the ff parameter file for ff name and version
        """
        ff_name = 'CHARMM'

        found_name = False
        with open(filename, 'r') as prmfile:
            for line in prmfile:
                if ff_name in line:
                    ff_name = line
                    found_name = True
                    break

        if not found_name:
            ff_name += ' %s' % filename.split('.')[0]
        else:
            ff_name = ff_name.translate(None, '*></\\')

        return ff_name.strip('\n')


    def get_atom_masses(self, top_rtf):
        """
        Reads atomtype atom masses from the topology .rtf file (CHARMM)
        """
        atom_masses = dict()

        with open(top_rtf, 'r') as top:
            for line in top:
                if line.startswith('MASS'):
                    atomtype = line.split()[2]
                    mass = float(line.split()[3])
                    atom_masses[atomtype] = mass
                if line.startswith('RESI'):
                    break

        return atom_masses


    def write_file_header(self, filename, ff_name, scale_14, convert_to_opls):
        """
        Write the header section to the new Q parameter file.
        """
        qprm = open(filename, 'w')
        improper_potential = 'harmonic'
        switch_atoms = 'off'
        vdw_rule = 'arithmetic'
        force_field = 'CHARMM'
        if convert_to_opls:
            improper_potential = 'periodic'
            switch_atoms = 'on'
            vdw_rule = 'geometric'
            force_field = 'AMBER'
        qprm.write('*------------------------------------------------\n'
        #           '*\n'
        #           '*Q-FF parameters: %s parameters\n'
        #           '*\n'
        #           '*------------------------------------------------\n'
        #           '*%s\n'
        #           '*Parameters translated for Q with charmm_Q_par.py (Geir V. Isaksen, 2015)\n'
        #           '*________________________________________________\n'
                   '[options]\n')
        #           'name            Q%s\n'
        #           'vdw_rule    %s   ! vdW combination rule (geometric or arithmetic )\n'
        #           'scale_14         %s         ! electrostatic 1-4 scaling factor\n'
        #           'switch_atoms  %s       ! on = use switch atom; off = use charge group\n'
        #           'improper_definition explicit ! improper representation by 2 or four atoms\n'
        #           'improper_potential      %s\n'
        #           'coulomb_constant 332.0716    ! Constant in electrostatic energy calculation; default = 332.\n'
        #           'force_field %s           ! Force Field Type (GROMOS (default), AMBER or CHARMM)\n\n' %
        #           (ff_name, ff_name, ff_name, vdw_rule, scale_14, switch_atoms, improper_potential, force_field))


    def write_Q_prm(self, filename, new_prm):
        """
        Writes new CHARMM parameters (new_prm) to filename
        """
        #Read original file (that is, if new parameters are to be merged with existing file)
        orig_data = open(filename, 'r').readlines()

        new_data = open(filename, 'w')

        to_add = False

        #Write new parameters to file:
        #If this file is empty, this loop will not produce anything.
        for line in orig_data:
            if to_add:
                if len(new_prm[to_add]) > 0:
                    if not line.startswith('*'):
                        for i in range(len(new_prm[to_add]) - 1, -1, -1):
                            new_data.write('%s\n' % (new_prm[to_add].pop(i)))
                        to_add = False

            new_data.write(line)

            if len(line.split()) > 0:
                if line.split()[0] in new_prm.keys():
                    to_add = line.split()[0]

        order = {1: '[atom_types]',
                 2: '[bonds]',
                 3: '[angles]',
                 4: '[torsions]',
                 5: '[impropers]'}

        headers = {1: '*-iac------Avdw1------Avdw2-----Bvdw1------Avdw3-----Bvdw2&3----mass---SYBYL-name-old-comment',
                   2: '*iaci  iacj     force.c.    dist.\n*------------------------------------------------',
                   3: '*iaci  iacj   iack      forceK   angle0     forceUB    Stheta\n'
                      '*------------------------------------------------------------',
                   4: '*iaci  iacj   iack   iacl      forceK   #minima     phase   #paths\n'
                      '*-----------------------------------------------------------------',
                   5: '*iaci  iacj   iack   iacl      forceK     imp0\n'
                      '*----------------------------------------------'}

        #Check if any sections are not written:
        for nr in sorted(order.keys()):
            section = order[nr]
            if len(new_prm[section]) > 0:
                new_data.write('\n%s\n' % section)
                #new_data.write('%s\n' % headers[nr])
                for line in new_prm[section]:
                    new_data.write('%s\n' % line)


    def read_charmm_top(self, filename):
        """
        Reads CHARMM top (.rtf) file and returns a dictionary with all information for Q lib
        """
        #Q lib friendly stuff will be stored in this dictionary (see further below for content):
        residues = dict()

        #Store patches for residues in different dictionary
        presidues = dict()

        #Boolen that controls if we are collecting lib data or not, in if so, the residue name
        res = False

        #Charge groups list:
        group = list()

        #Add entry to RESI (residues) or PRES (presidues):
        add_too = residues

        #Change charmm residue name of some residue to match Q users expectations:

        translate = {'HSE': 'HIE',
                     'HSP': 'HIP',
                     'HSD': 'HID'}

        #Open top .rtf and start extracting sensible data:
        with open(filename, 'r') as top:
            for line in top:
                if len(line.split()) > 0:
                    #Add new residue if 'RESI' is encountered
                    if len(line.split()) > 2:
                        if line.split()[0].lower() == 'resi' or line.split()[0].lower() == 'pres':
                            if res:
                                if len(group) > 0:
                                    add_too[res]['charge_groups'].append(group)

                            get_lib_entry = True

                            if line.split()[0].lower() == 'resi':
                                add_too = residues

                            elif line.split()[0].lower() == 'pres':
                                add_too = presidues

                        else:
                            get_lib_entry = False

                        if get_lib_entry:
                            res = line.split()[1]
                            if res in translate.keys():
                                res = translate[res]

                            res_charge = line.split()[2]

                            comment = ' '
                            if '!' in line:
                                comment = line.split('!')[-1].strip('\n')
                            add_too[res] = dict()
                            add_too[res]['charge'] = res_charge
                            add_too[res]['comment'] = comment
                            add_too[res]['atoms'] = dict()
                            add_too[res]['bonds'] = list()
                            add_too[res]['head'] = False
                            add_too[res]['tail'] = False
                            add_too[res]['impropers'] = list()
                            add_too[res]['charge_groups'] = list()
                            add_too[res]['delete'] = list()
                            group = list()

                    #Add info from .rtf file to current residue (res)
                    if res:
                        #If PRES, check for atoms to be deleted:
                        if add_too == presidues:
                            if len(line.split()) > 2:
                                if line.split()[0].lower() == 'delete' and line.split()[1].lower() == 'atom':
                                    add_too[res]['delete'].append(line.split()[2])

                        #Check if a new group is defined and add previous group to charge_groups
                        if line.split()[0] == 'GROUP':
                            if len(group) > 0:
                                add_too[res]['charge_groups'].append(group)
                                group = list()

                        #Get atom information (name, type and charge)
                        if line.split()[0] == 'ATOM':
                            atomname, atomtype, charge = line.split()[1:4]
                            charge = charge.translate(None, '!')
                            group.append(atomname)
                            comment = ' '
                            if '!' in line:
                                comment = '!' + line.split('!')[-1]
                            #Add atoms with nr to dictionary in order to write them out in the correct order later.
                            nr = len(add_too[res]['atoms']) + 1
                            add_too[res]['atoms'][nr] = dict()
                            add_too[res]['atoms'][nr]['name'] = atomname
                            add_too[res]['atoms'][nr]['type'] = atomtype
                            add_too[res]['atoms'][nr]['charge'] = charge
                            add_too[res]['atoms'][nr]['comment'] = comment

                        #Get all bonds (pairs of atoms)
                        if line.split()[0] == 'BOND' or line.split()[0] == 'DOUBLE':
                            pair = list()
                            for atom in line.split('!')[0].split()[1:]:
                                if '+' in atom:
                                    print(atom)
                                atom = atom.translate(None, '+-')
                                pair.append(atom)
                                if len(pair) == 2:
                                    add_too[res]['bonds'].append(pair)
                                    pair = list()


                        #Get improper definitions
                        if line.split()[0] == 'IMPR':
                            quad = list()
                            for atom in line.split('!')[0].split()[1:]:
                                quad.append(atom)
                                if len(quad) == 4:
                                    add_too[res]['impropers'].append(quad)
                                    quad = list()

                        #Then, here comes the "tricky" part. Guess the correct head and tail for the Q lib file!
                        #amino acids:   IC -C CA *N HN  <-- N is head
                        #               IC +N CA *C O   <-- C is tail
                        #               If only *, look in the comment to check if it is N- or C-terminus.
                        #DNA:           BILD -O3' O5'  *P   O1P <-- P is head O3' is tail ... not very robust, but will do!

                        #Proteins
                        if line.split()[0] == 'IC':
                            #If not terminal, check if connections head/tail can be extracted:
                            if '-' in line.split()[1]:
                                if not add_too[res]['head']:
                                    for atom in line.split()[2:5]:
                                        if '*' in atom:
                                            add_too[res]['head'] = atom.translate(None, '*')
                            elif '+' in line.split()[1]:
                                if not add_too[res]['tail']:
                                    for atom in line.split()[2:5]:
                                        if '*' in atom:
                                            add_too[res]['tail'] = atom.translate(None, '*')

                        #DNA:
                        if line.split()[0] == 'BILD':
                            if not add_too[res]['tail'] or not residues[res]['head']:
                                if '-' in line.split()[1]:
                                    add_too[res]['tail'] = line.split()[1].translate(None, '-')
                                    for atom in line.split()[2:5]:
                                        if '*' in atom:
                                            add_too[res]['head'] = atom.translate(None, '*')

        #Return residues and patches
        return residues, presidues


    def write_Q_lib(self, filename, residues, topname):
        """
        Write out a readable Q lib file from the residues dictionary collected from the CHARMM top .rtf file
        """
        libfile = open(filename, 'w')

        #Get charged residues with their charge for the header:
        charged_res = '*'
        res_charges = '*'

        for res in sorted(residues.keys()):
            if residues[res]['head'] and residues[res]['tail']:
                if float(residues[res]['charge']) != 0:
                    charged_res += ' %5s,' % res.ljust(4)
                    res_charges += ' %5s ' % residues[res]['charge']
        #Write header information
        libfile.write('*------------------------------------------------------------------\n')
        libfile.write('* Q fragment library file extracted from:\n* %s\n' % topname)
        libfile.write('* Converted with charmm_Q_par.py (Geir V. Isaksen, 2015)')
        libfile.write('\n*\n')
        libfile.write('*%s\n*%s\n' % (charged_res, res_charges))
        libfile.write('*------------------------------------------------------------------\n')

        for res in sorted(residues.keys()):
            #Write residue name, charge and additional comments from top .rtf
            libfile.write('{%s}            !charge: %s || %s\n' % (res, residues[res]['charge'], residues[res]['comment']))

            #Write atoms sections:
            libfile.write('    [atoms]\n')
            for nr in sorted(residues[res]['atoms'].keys()):
                atomname = residues[res]['atoms'][nr]['name']
                atomtype = residues[res]['atoms'][nr]['type']
                atomcharge = residues[res]['atoms'][nr]['charge']
                libfile.write('        %2d %6s %6s %10s\n' % (nr, atomname.ljust(6), atomtype.ljust(6), atomcharge))

            #Write bonds
            if len(residues[res]['bonds']) > 0:
                libfile.write('    [bonds]\n')
                for pair in residues[res]['bonds']:
                    libfile.write('       %6s %6s\n' % (pair[0].ljust(6), pair[1].ljust(6)) )

            #Write connections if they exist:
            if residues[res]['head'] or residues[res]['tail']:
                libfile.write('    [connections]\n')
                if residues[res]['head']:
                    libfile.write('        head %6s\n' % residues[res]['head'].ljust(6))
                if residues[res]['tail']:
                    libfile.write('        tail %6s\n' % residues[res]['tail'].ljust(6))

            #Write impropers
            if len(residues[res]['impropers']) > 0:
                libfile.write('    [impropers]\n')
                for improper in residues[res]['impropers']:
                    libfile.write('         %6s %6s %6s %6s\n' % (improper[0].ljust(6), improper[1].ljust(6),
                                                                  improper[2].ljust(6), improper[3].ljust(6)))

            #Write charge groups
            if len(residues[res]['charge_groups']) > 0:
                libfile.write('    [charge_groups]\n')
                for group in residues[res]['charge_groups']:
                    charge_group = '        '
                    for atom in group:
                        charge_group += '%6s' % atom.ljust(6)
                    libfile.write('%s\n' % charge_group)

            libfile.write('*------------------------------------------------------------------\n')


    def patch_to_residue(self, residues, patches):
        """
        Q lib does not support patches, so we need to generate all residues. Too bad, but ok, here goes
        """

        #Patch translation dictionary (copy residue from original and modify according to patch):
        patch_res = {'ASH': 'ASP',
                     'GLH': 'GLU',
                     'LYN' : 'LYS',
                     'SED': 'SER',
                     'NGLY': 'GLY',
                     'NPRO': 'PRO',
                     'CYX' : 'CYS'}

        #Get rid of these 4-letter residue names. Currently using 4-letter for N- and C-terminals.
        rename = {'ASPP': 'ASH',
                  'GLUP': 'GLH',
                  'LSN': 'LYN',
                  'SERD': 'SED'}

        for presname in rename.keys():
            if presname in patches.keys():
                patches[rename[presname]] = deepcopy(patches[presname])
                del patches[presname]

        #GLY and PRO have special N-terminals, so the do not use standard N-terminal for these two.
        #Prepare charged N terminals for GLY and PRO (these are different from the others)
        if 'GLYP' in patches.keys():
            patches['NGLY'] = deepcopy(patches['GLYP'])
            del patches['GLYP']

        if 'PROP' in patches.keys():
            patches['NPRO'] = deepcopy(patches['PROP'])
            del patches['PROP']


        #Fix DISU so that it is possible to create CYX from the patch:
        if 'DISU' in patches.keys():
            print(patches['DISU'])
            patches['CYX'] = deepcopy(patches['DISU'])
            del patches['DISU']

            atoms = list()
            for nr in sorted(patches['CYX']['atoms'].keys()):
                atomname = patches['CYX']['atoms'][nr]['name']
                if atomname[0].isdigit():
                    atomname = atomname[1:]
                    if atomname not in atoms:
                        atoms.append(atomname)
                        patches['CYX']['atoms'][nr]['name'] = atomname
                    else:
                        del patches['CYX']['atoms'][nr]

            patches['CYX']['bonds'] = list()
            patches['CYX']['charge_groups'] = list()
            #patches['CYX']['charge_groups'].append(atoms)

            atoms = list()
            for atom in patches['CYX']['delete']:
                if atom[0].isdigit():
                    atomname = atom[1:]
                    if not atomname in atoms:
                        atoms.append(atomname)
                        patches['CYX']['delete'][ patches['CYX']['delete'].index(atom)] = atomname
                    else:
                        del patches['CYX']['delete'][ patches['CYX']['delete'].index(atom)]

        for pres in patches.keys():
            if pres in patch_res.keys():
                residues[pres] = deepcopy(residues[patch_res[pres]])
                print('\n--> Creating %s from %s' % (pres, patch_res[pres]))

                #Delete atoms from original residue that is not present in the patch
                if len(patches[pres]['delete']) > 0:

                    for delete_atom in patches[pres]['delete']:
                        print(' - Deleting atom %s from %s to create %s' % (delete_atom, patch_res[pres], pres))
                        #delete from atoms:
                        for nr in residues[pres]['atoms'].keys():
                            if residues[pres]['atoms'][nr]['name'] == delete_atom:
                                del residues[pres]['atoms'][nr]

                        #Delete bonds with atom deleted:
                        delete_indexes = list()

                        for i in range(0, len(residues[pres]['bonds'])):
                            if delete_atom in residues[pres]['bonds'][i]:
                                delete_indexes.append(i)

                        for i in delete_indexes:
                            print(' - Deleting bond %s-%s from %s to create %s' % (residues[pres]['bonds'][i][0],
                                                                               residues[pres]['bonds'][i][1],
                                                                               patch_res[pres], pres))

                            del residues[pres]['bonds'][i]

                        #Delete from charge_groups
                        for i in range(len(residues[pres]['charge_groups'])):
                            charge_group = residues[pres]['charge_groups'][i]

                            if delete_atom in charge_group:
                                print(' - Deleting %s from charge_group: ' % delete_atom + ' '.join(charge_group))
                                del charge_group[charge_group.index(delete_atom)]

                                residues[pres]['charge_groups'][i] = charge_group

                        #Delete impropers containing atom
                        for i in range(len(residues[pres]['impropers'])):
                            if delete_atom in residues[pres]['impropers'][i]:
                                print('Deleting improper with atom %s: ' % delete_atom + \
                                      ' '.join(residues[pres]['impropers'][i]))

                #Modify atomtype, charges and comments:
                patch_atoms = list()
                name_number = dict()
                for nr in sorted(patches[pres]['atoms'].keys()):
                    patch_atoms.append(patches[pres]['atoms'][nr]['name'])
                    name_number[patches[pres]['atoms'][nr]['name']] = nr

                for nr in sorted(residues[pres]['atoms'].keys()):
                    atomname = residues[pres]['atoms'][nr]['name']

                    if atomname in patch_atoms:
                        pnr = name_number[atomname]

                        print(' * Atomtype %s with charge %s changed to atomtype %s with charge %s' % \
                              (residues[pres]['atoms'][nr]['type'], residues[pres]['atoms'][nr]['charge'],
                               patches[pres]['atoms'][pnr]['type'], patches[pres]['atoms'][pnr]['charge']))

                        residues[pres]['atoms'][nr]['type'] = patches[pres]['atoms'][pnr]['type']
                        residues[pres]['atoms'][nr]['charge'] = patches[pres]['atoms'][pnr]['charge']
                        residues[pres]['atoms'][nr]['comment'] = patches[pres]['atoms'][pnr]['comment']

                        del patch_atoms[patch_atoms.index(atomname)]

                #Add new atoms from patch
                if len(patch_atoms) > 0:
                    for atom in patch_atoms:
                        pnr = name_number[atom]
                        nr = max(residues[pres]['atoms']) + 1
                        residues[pres]['atoms'][nr] = dict()
                        print(' + Adding atom %s to %s to create %s' % (patches[pres]['atoms'][pnr]['name'],
                                                                       patch_res[pres], pres))
                        residues[pres]['atoms'][nr]['name'] = patches[pres]['atoms'][pnr]['name']
                        residues[pres]['atoms'][nr]['type'] = patches[pres]['atoms'][pnr]['type']
                        residues[pres]['atoms'][nr]['charge'] = patches[pres]['atoms'][pnr]['charge']
                        residues[pres]['atoms'][nr]['comment'] = patches[pres]['atoms'][pnr]['comment']

                #Add new bonds from patch
                for bond in patches[pres]['bonds']:
                    if bond not in residues[pres]['bonds']:
                        print(' + Adding bond %s-%s to %s to create %s' % (bond[0], bond[1], patch_res[pres], pres))
                        residues[pres]['bonds'].append(bond)

                #add new impropers
                for improper in patches[pres]['impropers']:
                    if improper not in residues[pres]['impropers']:
                        a1, a2, a3, a4 = improper[0:]
                        print(' + Adding improper %s %s %s %s to %s to create %s' % (a1, a2, a3, a4, patch_res[pres], pres))
                        residues[pres]['impropers'].append(improper)

                #add new charge_groups
                for new_group in patches[pres]['charge_groups']:
                    for old_group in residues[pres]['charge_groups']:
                        if len(set(new_group) & set(old_group)) > 0:
                            print(' * Charge group modified from [%s] to [%s]' % ( ' '.join(old_group),' '.join(new_group)))
                            residues[pres]['charge_groups'][residues[pres]['charge_groups'].index(old_group)] = new_group

                #Remove head connection  and improper if N-terminal GLY or PRO patch
                del_imp = list()
                if pres == 'NGLY' or pres == 'NPRO':
                    residues[pres]['head'] = False
                    for i in range(len(residues[pres]['impropers'])):
                        if '-C' in residues[pres]['impropers'][i]:
                            del_imp.append(i)
                for i in del_imp:
                    del residues[pres]['impropers'][i]


                #Total charge of new residue:
                tot_charge = 0.0
                for nr in residues[pres]['atoms'].keys():
                    tot_charge += float(residues[pres]['atoms'][nr]['charge'])
                residues[pres]['charge'] = str(round(tot_charge, 1))

        return residues


    def patch_terminals(self, residues, patches):
        """
        Generate charged and neutral N- and C-terminals for all residues
        """
        orig_residues = deepcopy(residues)

        #Charge N-terminal: NRES
        #Neutral N-terminal: nRES
        #Charged C-terminal: CRES
        #Neutral C-terminal: cRES

        for pres in patches.keys():
            first_letter = False
            res_comment = ' '
            if 'standard n-term' in patches[pres]['comment'].lower():
                first_letter = 'N'
                res_comment = 'N-terminal (standard)'
            if 'neutral n-term' in patches[pres]['comment'].lower():
                first_letter = 'n'
                res_comment = 'N-terminal (neutral)'
            if 'standard c-term' in patches[pres]['comment'].lower():
                first_letter = 'C'
                res_comment = 'C-terminal (standard)'
            if 'neutral c-term' in patches[pres]['comment'].lower():
                first_letter = 'c'
                res_comment = 'C-terminal (neutral)'

            if first_letter:
                for res in orig_residues.keys():

                    #Check that this is a proper head and tail fragment:
                    if orig_residues[res]['head'] and residues[res]['tail']:

                        #Generate new name for terminal residue
                        tres = '%s%s' % (first_letter, res)
                        residues[tres] = deepcopy(residues[res])
                        residues[tres]['comment'] ='%s %s' % (res, res_comment)

                        #Remove head or tail and the spanning improper
                        del_imp = list()
                        if first_letter.lower() == 'n':
                            residues[tres]['head'] = False
                            for i in range(len(residues[tres]['impropers'])):
                                if '-C' in residues[tres]['impropers'][i]:
                                    del_imp.append(i)
                        elif first_letter.lower() == 'c':
                            residues[tres]['tail'] = False
                            for i in range(len(residues[tres]['impropers'])):
                                if '+N' in residues[tres]['impropers'][i]:
                                    del_imp.append(i)

                        for i in del_imp:
                            del residues[tres]['impropers'][i]

                        print('\n--> Creating %s from %s' % (tres, res))

                        #Delete atoms from original residue that is not tresent in the patch
                        if len(patches[pres]['delete']) > 0:

                            for delete_atom in patches[pres]['delete']:
                                print( ' - Deleting atom %s from %s to create %s' % (delete_atom, res, tres))
                                #delete from atoms:
                                for nr in residues[tres]['atoms'].keys():
                                    if residues[tres]['atoms'][nr]['name'] == delete_atom:
                                        del residues[tres]['atoms'][nr]

                                #Delete bonds with atom deleted:
                                delete_indexes = list()

                                for i in range(0, len(residues[tres]['bonds'])):
                                    if delete_atom in residues[tres]['bonds'][i]:
                                        delete_indexes.append(i)

                                for i in delete_indexes:
                                    print( ' - Deleting bond %s-%s from %s to create %s' % (residues[tres]['bonds'][i][0],
                                                                                       residues[tres]['bonds'][i][1],
                                                                                       res, tres))

                                    del residues[tres]['bonds'][i]

                                #Delete from charge_groups
                                for i in range(len(residues[tres]['charge_groups'])):
                                    charge_group = residues[tres]['charge_groups'][i]

                                    if delete_atom in charge_group:
                                        print( ' - Deleting %s from charge_group: ' % delete_atom + ' '.join(charge_group))
                                        del charge_group[charge_group.index(delete_atom)]

                                        residues[tres]['charge_groups'][i] = charge_group

                                #Delete impropers containing atom
                                for i in range(len(residues[tres]['impropers'])):
                                    if delete_atom in residues[tres]['impropers'][i]:
                                        print( 'Deleting improper with atom %s: ' % delete_atom + \
                                              ' '.join(residues[tres]['impropers'][i]))

                        #Modify atomtype, charges and comments:
                        patch_atoms = list()
                        name_number = dict()
                        for nr in sorted(patches[pres]['atoms'].keys()):
                            patch_atoms.append(patches[pres]['atoms'][nr]['name'])
                            name_number[patches[pres]['atoms'][nr]['name']] = nr

                        for nr in sorted(residues[tres]['atoms'].keys()):
                            atomname = residues[tres]['atoms'][nr]['name']

                            if atomname in patch_atoms:
                                pnr = name_number[atomname]

                                print( ' * Atomtype %s with charge %s changed to atomtype %s with charge %s' % \
                                      (residues[tres]['atoms'][nr]['type'], residues[tres]['atoms'][nr]['charge'],
                                       patches[pres]['atoms'][pnr]['type'], patches[pres]['atoms'][pnr]['charge']))

                                residues[tres]['atoms'][nr]['type'] = patches[pres]['atoms'][pnr]['type']
                                residues[tres]['atoms'][nr]['charge'] = patches[pres]['atoms'][pnr]['charge']
                                residues[tres]['atoms'][nr]['comment'] = patches[pres]['atoms'][pnr]['comment']

                                del patch_atoms[patch_atoms.index(atomname)]

                        #Add new atoms from patch
                        if len(patch_atoms) > 0:
                            for atom in patch_atoms:
                                pnr = name_number[atom]
                                nr = max(residues[tres]['atoms']) + 1
                                residues[tres]['atoms'][nr] = dict()
                                print( ' + Adding atom %s to %s to create %s' % (patches[pres]['atoms'][pnr]['name'],
                                                                               res, tres))
                                residues[tres]['atoms'][nr]['name'] = patches[pres]['atoms'][pnr]['name']
                                residues[tres]['atoms'][nr]['type'] = patches[pres]['atoms'][pnr]['type']
                                residues[tres]['atoms'][nr]['charge'] = patches[pres]['atoms'][pnr]['charge']
                                residues[tres]['atoms'][nr]['comment'] = patches[pres]['atoms'][pnr]['comment']

                        #Add new bonds from patch
                        for bond in patches[pres]['bonds']:
                            if bond not in residues[tres]['bonds']:
                                print( ' + Adding bond %s-%s to %s to create %s' % (bond[0], bond[1], res, tres))
                                residues[tres]['bonds'].append(bond)

                        #add new impropers
                        for improper in patches[pres]['impropers']:
                            if improper not in residues[tres]['impropers']:
                                a1, a2, a3, a4 = improper[0:]
                                print( ' + Adding improper %s %s %s %s to %s to create %s' % (a1, a2, a3, a4, res, tres))
                                residues[tres]['impropers'].append(improper)

                        #add new charge_groups
                        for new_group in patches[pres]['charge_groups']:
                            for old_group in residues[tres]['charge_groups']:
                                if len(set(new_group) & set(old_group)) > 0:
                                    print( ' * Charge group modified from [%s] to [%s]' % ( ' '.join(old_group),' '.join(new_group)))
                                    residues[tres]['charge_groups'][residues[tres]['charge_groups'].index(old_group)] = new_group

                        #Total charge of new residue:
                        tot_charge = 0.0
                        for nr in residues[tres]['atoms'].keys():
                            tot_charge += float(residues[tres]['atoms'][nr]['charge'])
                        residues[tres]['charge'] = str(round(tot_charge, 1))

        return residues
    

    def toplevel(self):
        #Try to get force field name and version:
        ff_name = run.find_ff_name(self.charmm_prm)
        convert_to_opls = False
        #Append new parameters to existing Q parameter file?
        qprm = {'[atom_types]': list(),
                '[bonds]': list(),
                '[angles]': list(),
                '[torsions]': list(),
                '[impropers]': list()}

        #Get atom weights from top .rtf file
        atom_mass = run.get_atom_masses(self.charmm_top)

        #Get new parameters from charmm .par file
        new_prm, scale_14 = run.read_charmm_par(self.charmm_prm, atom_mass, qprm, convert_to_opls)

        if self.merge == True:
            print('NOT WORKING YET')
            old_q_prm = s.FF_DIR + self.FF + '.prm'
            
        else:
            old_q_prm = 'Q{}.prm'.format(self.lig)
        
        run.write_file_header(old_q_prm, ff_name, scale_14, convert_to_opls)

        print('')

        print('%6d NEW ATOM TYPES PARAMETERS FOUND' % len(new_prm['[atom_types]']))
        print('%6d NEW BOND PARAMETERS FOUND' % len(new_prm['[bonds]']))
        print('%6d NEW ANGLE PARAMETERS FOUND' % len(new_prm['[angles]']))
        print('%6d NEW TORSION PARAMETERS FOUND' % len(new_prm['[torsions]']))
        print('%6d NEW IMPROPER PARAMETERS FOUND' % len(new_prm['[impropers]']))

        #Write new parameters to Q .prm file
        run.write_Q_prm(old_q_prm, new_prm)

        residues, patches = run.read_charmm_top(self.charmm_top)

        print('\n%6d RESIDUES FOUND' % len(residues.keys()))
        print('%6d PATCHES FOUND' % len(patches.keys()))

        #Use relevant patches and generate suitable residues for Q lib:
        residues = run.patch_to_residue(residues, patches)

        #Generate all N- and C-terminals from patches
        residues = run.patch_terminals(residues, patches)

        #Write Q library file
        q_libfile = 'Q{}.lib'.format(self.lig)
        top_name = run.find_ff_name(self.charmm_top)

        run.write_Q_lib(q_libfile, residues, top_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='charmm2Q',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Convert CHARMM ligand parameters to Q files == ')

    
    parser.add_argument('-l', '--lig',
                        dest = "lig",
                        required = True,
                        help = "name of the ligand")
    
    parser.add_argument('-m', '--merge',
                        dest = "merge",
                        default = False,
                        action = 'store_true',
                        help = "define if ligand prm file needs to be merged to ref")
    
    parser.add_argument('-FF', '--forcefield',
                        dest = "forcefield",
                        required = True,
                        help = "name of CHARMM type forcefield, e.g. CHARMM36")
    args = parser.parse_args()
    
    run = Run(lig = args.lig,
              merge = args.merge,
              forcefield = args.forcefield)
    
    run.toplevel()
