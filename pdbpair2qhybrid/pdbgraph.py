from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple
import numpy as np
from .io import AtomRef, read_pdb_atoms

_COVALENT_RADII = {
    'H': 0.31,
    'C': 0.76,
    'N': 0.71,
    'O': 0.66,
    'S': 1.05,
    'P': 1.07,
    'F': 0.57,
    'CL': 0.99,
    'BR': 1.14,
    'I': 1.33,
}

@dataclass
class MolecularGraph:
    atoms: Dict[str, AtomRef]
    neighbors: Dict[str, Set[str]] = field(default_factory=dict)
    residue_atoms: Dict[str, List[str]] = field(default_factory=dict)

    @classmethod
    def from_pdb(cls, path: str) -> 'MolecularGraph':
        atoms_list = read_pdb_atoms(path)
        atoms = {atom.id(): atom for atom in atoms_list}
        residue_atoms = {}
        for atom in atoms_list:
            residue_atoms.setdefault(atom.residue_id(), []).append(atom.id())

        graph = cls(atoms=atoms, residue_atoms=residue_atoms)
        graph.neighbors = graph._infer_bonds()
        return graph

    def _infer_bonds(self) -> Dict[str, Set[str]]:
        atom_ids = list(self.atoms.keys())
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms.values()])
        neighbors = {atom_id: set() for atom_id in atom_ids}
        for i, atom_id in enumerate(atom_ids):
            atom_i = self.atoms[atom_id]
            for j in range(i + 1, len(atom_ids)):
                atom_j = self.atoms[atom_ids[j]]
                dist = np.linalg.norm(coords[i] - coords[j])
                if self._is_covalent_pair(atom_i, atom_j, dist):
                    neighbors[atom_id].add(atom_ids[j])
                    neighbors[atom_ids[j]].add(atom_id)
        return neighbors

    def _is_covalent_pair(self, atom_a: AtomRef, atom_b: AtomRef, distance: float) -> bool:
        if distance < 0.1:
            return False
        r_a = _COVALENT_RADII.get(atom_a.element.upper(), 0.75)
        r_b = _COVALENT_RADII.get(atom_b.element.upper(), 0.75)
        cutoff = r_a + r_b + 0.45
        return distance <= cutoff

    def transform(self, rotation: np.ndarray, translation: np.ndarray) -> 'MolecularGraph':
        new_atoms = {}
        for atom_id, atom in self.atoms.items():
            coord = np.dot(rotation, np.array([atom.x, atom.y, atom.z])) + translation
            new_atoms[atom_id] = AtomRef(
                serial=atom.serial,
                name=atom.name,
                altloc=atom.altloc,
                resname=atom.resname,
                chain=atom.chain,
                resi=atom.resi,
                icode=atom.icode,
                x=float(coord[0]),
                y=float(coord[1]),
                z=float(coord[2]),
                occupancy=atom.occupancy,
                bfactor=atom.bfactor,
                element=atom.element,
                charge=atom.charge,
                record=atom.record,
            )
        return MolecularGraph(atoms=new_atoms, neighbors=self.neighbors.copy(), residue_atoms=self.residue_atoms.copy())

    def atom_ids(self):
        return list(self.atoms.keys())

    def coordinates_for(self, atom_ids: List[str]) -> np.ndarray:
        return np.array([[self.atoms[atom_id].x,
                          self.atoms[atom_id].y,
                          self.atoms[atom_id].z] for atom_id in atom_ids])


def read_pdb_as_graph(path: str) -> MolecularGraph:
    return MolecularGraph.from_pdb(path)
