from dataclasses import dataclass, field
from typing import List, Set, Dict
from .pdbgraph import MolecularGraph
from .atommatch import AtomMapping

@dataclass
class RegionSplit:
    common_atoms: List[str]
    a_only_atoms: List[str]
    b_only_atoms: List[str]
    boundary_atoms: List[str]
    residues: List[str]
    connected: bool
    warnings: List[str] = field(default_factory=list)


def extract_changed_region(graphA: MolecularGraph, graphB: MolecularGraph, mapping: AtomMapping) -> RegionSplit:
    common_a = mapping.common_a_ids()
    common_b = mapping.common_b_ids()
    a_only = [atom_id for atom_id in graphA.atom_ids() if atom_id not in common_a]
    b_only = [atom_id for atom_id in graphB.atom_ids() if atom_id not in common_b]
    boundary = []

    common_atoms = [atom_id for atom_id in common_a]
    a_set = set(a_only)
    b_set = set(b_only)

    for atom_id in common_atoms:
        a_neighbors = graphA.neighbors.get(atom_id, set())
        b_neighbors = graphB.neighbors.get(atom_id, set())
        if a_neighbors & a_set or b_neighbors & b_set:
            boundary.append(atom_id)

    residues = sorted({graphA.atoms[atom_id].residue_id() for atom_id in a_only if atom_id in graphA.atoms} |
                      {graphB.atoms[atom_id].residue_id() for atom_id in b_only if atom_id in graphB.atoms})

    combined = list(a_only) + list(b_only) + boundary
    connected = _is_connected(combined, graphA, graphB)
    warnings = []
    if not connected and combined:
        warnings.append('Changed region is not fully connected across endpoints')
    if len(residues) > 1 and len(boundary) == 0:
        warnings.append('Changed region spans multiple residues with no boundary atoms')

    return RegionSplit(
        common_atoms=sorted(common_atoms),
        a_only_atoms=sorted(a_only),
        b_only_atoms=sorted(b_only),
        boundary_atoms=sorted(boundary),
        residues=residues,
        connected=connected,
        warnings=warnings,
    )


def _is_connected(atom_ids, graphA: MolecularGraph, graphB: MolecularGraph) -> bool:
    if not atom_ids:
        return True
    atom_set = set(atom_ids)
    neighbors = {atom_id: set() for atom_id in atom_ids}
    for atom_id in atom_ids:
        if atom_id in graphA.neighbors:
            neighbors[atom_id].update(graphA.neighbors[atom_id] & atom_set)
        if atom_id in graphB.neighbors:
            neighbors[atom_id].update(graphB.neighbors[atom_id] & atom_set)

    visited = set()
    stack = [atom_ids[0]]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        stack.extend(neighbors[current] - visited)
    return visited == atom_set


def classify_mode(region_split: RegionSplit, graphA: MolecularGraph, graphB: MolecularGraph) -> str:
    if not region_split.common_atoms:
        return 'adduct'
    if len(region_split.residues) == 1:
        changed_backbone = any(graphA.atoms.get(atom_id, graphB.atoms.get(atom_id)).is_backbone()
                               for atom_id in region_split.a_only_atoms + region_split.b_only_atoms)
        if not changed_backbone:
            return 'residue'
    return 'adduct'
