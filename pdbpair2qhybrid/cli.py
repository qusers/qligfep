import argparse
import os
from .pdbgraph import read_pdb_as_graph
from .atommatch import match_atoms
from .align import align_endpoints
from .diffregion import extract_changed_region, classify_mode
from .report import write_atom_map, write_build_report


def run(stateA: str, stateB: str, outdir: str, mode: str = 'auto', verbose: bool = False):
    graphA = read_pdb_as_graph(stateA)
    graphB = read_pdb_as_graph(stateB)

    initial_mapping = match_atoms(graphA, graphB)
    graphB_aligned, rmsd, align_meta = align_endpoints(graphA, graphB, mapping=initial_mapping)

    mapping = match_atoms(graphA, graphB_aligned)
    region_split = extract_changed_region(graphA, graphB_aligned, mapping)

    chosen_mode = mode if mode != 'auto' else classify_mode(region_split, graphA, graphB_aligned)
    os.makedirs(outdir, exist_ok=True)

    atom_map_path = os.path.join(outdir, 'atom_map.json')
    report_path = os.path.join(outdir, 'build_report.md')
    write_atom_map(mapping, region_split, [], chosen_mode, atom_map_path)
    write_build_report(stateA, stateB, chosen_mode, rmsd, mapping, region_split, report_path)

    if verbose:
        print(f'Wrote: {atom_map_path}')
        print(f'Wrote: {report_path}')
        print(f'common atoms: {len(mapping.pairs)}')
        print(f'A-only atoms: {len(region_split.a_only_atoms)}')
        print(f'B-only atoms: {len(region_split.b_only_atoms)}')
        print(f'alignment RMSD: {rmsd:.4f}')

    return {
        'stateA': stateA,
        'stateB': stateB,
        'outdir': outdir,
        'mode': chosen_mode,
        'rmsd': rmsd,
        'common': len(mapping.pairs),
        'a_only': len(region_split.a_only_atoms),
        'b_only': len(region_split.b_only_atoms),
    }


def main():
    parser = argparse.ArgumentParser(description='Build a paired-endpoint hybrid mapping for Q-style dual-topology setup.')
    parser.add_argument('--stateA', required=True, help='Endpoint A PDB file')
    parser.add_argument('--stateB', required=True, help='Endpoint B PDB file')
    parser.add_argument('--outdir', required=True, help='Output directory for atom_map.json and build_report.md')
    parser.add_argument('--mode', choices=['auto', 'residue', 'adduct'], default='auto', help='Mode selection')
    parser.add_argument('--verbose', action='store_true', help='Print progress messages')

    args = parser.parse_args()
    run(args.stateA, args.stateB, args.outdir, mode=args.mode, verbose=args.verbose)
