import json
from pathlib import Path
from typing import List
from .atommatch import AtomMapping
from .diffregion import RegionSplit


def write_atom_map(mapping: AtomMapping, region_split: RegionSplit, anchors: List[str], mode: str, path: str):
    payload = {
        'mode': mode,
        'common': [vars(pair) for pair in mapping.pairs],
        'a_only': sorted(region_split.a_only_atoms),
        'b_only': sorted(region_split.b_only_atoms),
        'boundary': sorted(region_split.boundary_atoms),
        'residues': sorted(region_split.residues),
        'anchors': anchors,
        'warnings': sorted(region_split.warnings),
    }
    Path(path).write_text(json.dumps(payload, indent=2))


def write_build_report(stateA: str, stateB: str, mode: str, rmsd: float, mapping: AtomMapping, region_split: RegionSplit, output_path: str):
    lines = [
        '# pdbpair2qhybrid Build Report',
        '',
        f'* State A: `{stateA}`',
        f'* State B: `{stateB}`',
        f'* Inferred mode: **{mode}**',
        f'* Alignment RMSD: `{rmsd:.4f}` Å',
        '',
        '## Atom counts',
        f'* common atoms: `{len(mapping.pairs)}`',
        f'* A-only atoms: `{len(region_split.a_only_atoms)}`',
        f'* B-only atoms: `{len(region_split.b_only_atoms)}`',
        f'* boundary atoms: `{len(region_split.boundary_atoms)}`',
        '',
        '## Changed region residues',
        f'* `{len(region_split.residues)}` residues: {", ".join(region_split.residues) or "none"}',
        '',
        '## Warnings',
    ]

    if region_split.warnings:
        for warning in region_split.warnings:
            lines.append(f'* {warning}')
    else:
        lines.append('* None')

    lines.extend([
        '',
        '## Match provenance',
    ])

    for pair in mapping.pairs:
        lines.append(f'* `{pair.a_id}` ↔ `{pair.b_id}` ({pair.provenance})')

    Path(output_path).write_text('\n'.join(lines) + '\n')
