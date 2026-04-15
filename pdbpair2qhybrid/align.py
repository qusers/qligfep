from typing import Dict, Tuple
import numpy as np
from .pdbgraph import MolecularGraph
from .atommatch import match_atoms, AtomMapping


def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # Compute rotation and translation that aligns Q onto P
    P_center = P.mean(axis=0)
    Q_center = Q.mean(axis=0)
    P_centered = P - P_center
    Q_centered = Q - Q_center
    covariance = np.dot(P_centered.T, Q_centered)
    U, S, Vt = np.linalg.svd(covariance)
    d = np.linalg.det(np.dot(U, Vt))
    D = np.diag([1.0, 1.0, d])
    rotation = np.dot(U, np.dot(D, Vt))
    translation = P_center - np.dot(rotation, Q_center)
    return rotation, translation


def _rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    diff = P - Q
    return float(np.sqrt((diff * diff).sum() / len(P)))


def align_endpoints(graphA: MolecularGraph, graphB: MolecularGraph, mapping: AtomMapping = None):
    if mapping is None:
        mapping = match_atoms(graphA, graphB)

    matched_pairs = [pair for pair in mapping.pairs if pair.provenance in ('label', 'residue')]
    if len(matched_pairs) < 3:
        raise ValueError('Need at least 3 matched atoms for reliable alignment')

    atom_ids_a = [pair.a_id for pair in matched_pairs]
    atom_ids_b = [pair.b_id for pair in matched_pairs]
    coords_a = graphA.coordinates_for(atom_ids_a)
    coords_b = graphB.coordinates_for(atom_ids_b)
    rotation, translation = _kabsch(coords_a, coords_b)
    graphB_aligned = graphB.transform(rotation, translation)
    rmsd = _rmsd(coords_a, graphB_aligned.coordinates_for(atom_ids_b))
    metadata = {
        'matched_atoms': len(matched_pairs),
        'rmsd': rmsd,
    }
    return graphB_aligned, rmsd, metadata
