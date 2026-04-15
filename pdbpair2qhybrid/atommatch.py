from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
import numpy as np
from .pdbgraph import MolecularGraph

@dataclass
class AtomMatch:
    a_id: str
    b_id: str
    provenance: str

@dataclass
class AtomMapping:
    pairs: List[AtomMatch] = field(default_factory=list)

    def common_a_ids(self) -> Set[str]:
        return {pair.a_id for pair in self.pairs}

    def common_b_ids(self) -> Set[str]:
        return {pair.b_id for pair in self.pairs}

    def as_dict(self):
        return [vars(pair) for pair in self.pairs]


def match_atoms(graphA: MolecularGraph, graphB: MolecularGraph) -> AtomMapping:
    mapping = AtomMapping()
    unmatched_a = set(graphA.atom_ids())
    unmatched_b = set(graphB.atom_ids())

    # Exact label match by chain:resid:atom name
    for atom_id in list(unmatched_a):
        if atom_id in unmatched_b:
            atom_a = graphA.atoms[atom_id]
            atom_b = graphB.atoms[atom_id]
            if atom_a.element.upper() == atom_b.element.upper():
                mapping.pairs.append(AtomMatch(atom_id, atom_id, 'label'))
                unmatched_a.remove(atom_id)
                unmatched_b.remove(atom_id)

    # Same residue and same element fallback
    for a_id in list(unmatched_a):
        atom_a = graphA.atoms[a_id]
        candidates = [b_id for b_id in unmatched_b
                      if graphB.atoms[b_id].chain == atom_a.chain
                      and graphB.atoms[b_id].resi == atom_a.resi
                      and graphB.atoms[b_id].element.upper() == atom_a.element.upper()]
        if not candidates:
            continue
        coords_a = np.array([atom_a.x, atom_a.y, atom_a.z])
        closest = min(candidates,
                      key=lambda b_id: np.linalg.norm(coords_a - np.array([graphB.atoms[b_id].x,
                                                                            graphB.atoms[b_id].y,
                                                                            graphB.atoms[b_id].z])))
        distance = np.linalg.norm(coords_a - np.array([graphB.atoms[closest].x,
                                                      graphB.atoms[closest].y,
                                                      graphB.atoms[closest].z]))
        if distance <= 1.5:
            mapping.pairs.append(AtomMatch(a_id, closest, 'residue'))
            unmatched_a.remove(a_id)
            unmatched_b.remove(closest)

    # final spatial nearest-neighbor fallback for same element
    for a_id in list(unmatched_a):
        atom_a = graphA.atoms[a_id]
        candidates = [b_id for b_id in unmatched_b
                      if graphB.atoms[b_id].element.upper() == atom_a.element.upper()]
        if not candidates:
            continue
        coords_a = np.array([atom_a.x, atom_a.y, atom_a.z])
        closest = min(candidates,
                      key=lambda b_id: np.linalg.norm(coords_a - np.array([graphB.atoms[b_id].x,
                                                                            graphB.atoms[b_id].y,
                                                                            graphB.atoms[b_id].z])))
        distance = np.linalg.norm(coords_a - np.array([graphB.atoms[closest].x,
                                                      graphB.atoms[closest].y,
                                                      graphB.atoms[closest].z]))
        if distance <= 2.0:
            mapping.pairs.append(AtomMatch(a_id, closest, 'spatial'))
            unmatched_a.remove(a_id)
            unmatched_b.remove(closest)

    return mapping
