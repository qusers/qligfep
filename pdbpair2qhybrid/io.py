from dataclasses import dataclass
import IO as qio
import re

_BACKBONE_NAMES = {'N', 'CA', 'C', 'O', 'OXT', 'H', 'HA', 'HA2', 'HA3'}

@dataclass(frozen=True)
class AtomRef:
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resi: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    charge: str
    record: str

    def id(self) -> str:
        ins = self.icode.strip()
        return f"{self.chain}:{self.resi}{ins}:{self.name}"

    def residue_id(self) -> str:
        ins = self.icode.strip()
        return f"{self.chain}:{self.resi}{ins}:{self.resname}"

    def is_backbone(self) -> bool:
        return self.name in _BACKBONE_NAMES


def parse_pdb_atom(line: str) -> AtomRef:
    atom = qio.pdb_parse_in(line)
    if not isinstance(atom, list) or len(atom) < 14:
        raise ValueError(f"Unable to parse PDB atom line: {line!r}")

    element = atom[13].strip() if atom[13].strip() else _infer_element(atom[2])
    return AtomRef(
        serial=atom[1],
        name=atom[2].strip(),
        altloc=atom[3],
        resname=atom[4].strip(),
        chain=atom[5],
        resi=atom[6],
        icode=atom[7],
        x=atom[8],
        y=atom[9],
        z=atom[10],
        occupancy=atom[11],
        bfactor=atom[12],
        element=element,
        charge=atom[14].strip(),
        record=atom[0].strip(),
    )


def _infer_element(atom_name: str) -> str:
    atom_name = atom_name.strip()
    if not atom_name:
        return ''
    # PDB atom names are right- or left-justified depending on element
    if len(atom_name) == 1:
        return atom_name[0].upper()
    if atom_name[0].isalpha() and atom_name[1].islower():
        return atom_name[:2].title()
    return atom_name[0].upper()


def read_pdb_atoms(path: str):
    atoms = []
    with open(path) as infile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')):
                atoms.append(parse_pdb_atom(line))
    return atoms
