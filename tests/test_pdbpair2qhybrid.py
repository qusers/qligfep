"""
Tests for the new pdbpair2qhybrid paired-endpoint hybrid builder.
"""

import json
from pathlib import Path
import sys
import tempfile
import shutil

sys.path.insert(0, str(Path(__file__).parent.parent))

from pdbpair2qhybrid import (
    read_pdb_as_graph,
    align_endpoints,
    match_atoms,
    extract_changed_region,
    classify_mode,
    write_atom_map,
    write_build_report,
    cli,
)


def _write_pdb(path, lines):
    Path(path).write_text('\n'.join(lines) + '\n')


def test_read_pdb_as_graph_simple_protein(tmp_path):
    pdb_text = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.936  -1.020   1.228  1.00  0.00           C',
        'TER'
    ]
    pdb_file = tmp_path / 'stateA.pdb'
    _write_pdb(pdb_file, pdb_text)

    graph = read_pdb_as_graph(str(pdb_file))
    assert len(graph.atoms) == 5
    assert 'A:1:CB' in graph.atoms
    assert graph.neighbors['A:1:CB']


def test_align_endpoints_with_translation(tmp_path):
    pdb_text = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
    ]
    a_path = tmp_path / 'A.pdb'
    b_path = tmp_path / 'B.pdb'
    _write_pdb(a_path, pdb_text)
    translated = []
    for line in pdb_text:
        if line.startswith('ATOM'):
            x = float(line[30:38]) + 2.0
            y = float(line[38:46]) + 3.0
            z = float(line[46:54]) + 4.0
            translated.append(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}".rstrip())
        else:
            translated.append(line)
    _write_pdb(b_path, translated)

    graphA = read_pdb_as_graph(str(a_path))
    graphB = read_pdb_as_graph(str(b_path))
    aligned, rmsd, _ = align_endpoints(graphA, graphB)

    assert rmsd < 1e-6
    assert abs(aligned.atoms['A:1:CA'].x - graphA.atoms['A:1:CA'].x) < 1e-6


def test_match_atoms_and_extract_changed_region(tmp_path):
    stateA = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.936  -1.020   1.228  1.00  0.00           C',
    ]
    stateB = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
        'ATOM      5  CB  VAL A   1       1.936  -1.020   1.228  1.00  0.00           C',
        'ATOM      6 CG1  VAL A   1       2.912  -1.102   0.411  1.00  0.00           C',
        'ATOM      7 CG2  VAL A   1       1.918  -2.294   1.183  1.00  0.00           C',
    ]
    a_path = tmp_path / 'stateA.pdb'
    b_path = tmp_path / 'stateB.pdb'
    _write_pdb(a_path, stateA)
    _write_pdb(b_path, stateB)

    graphA = read_pdb_as_graph(str(a_path))
    graphB = read_pdb_as_graph(str(b_path))
    mapping = match_atoms(graphA, graphB)
    region = extract_changed_region(graphA, graphB, mapping)

    assert 'A:1:CB' in region.common_atoms
    assert len(region.a_only_atoms) == 0
    assert len(region.b_only_atoms) == 2
    assert region.connected
    assert classify_mode(region, graphA, graphB) == 'residue'


def test_cli_writes_atom_map_and_report(tmp_path):
    stateA = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.936  -1.020   1.228  1.00  0.00           C',
    ]
    stateB = [
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       1.942   1.412   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.138   2.286   0.000  1.00  0.00           O',
        'ATOM      5  CB  VAL A   1       1.936  -1.020   1.228  1.00  0.00           C',
        'ATOM      6 CG1  VAL A   1       2.912  -1.102   0.411  1.00  0.00           C',
        'ATOM      7 CG2  VAL A   1       1.918  -2.294   1.183  1.00  0.00           C',
    ]
    a_path = tmp_path / 'stateA.pdb'
    b_path = tmp_path / 'stateB.pdb'
    outdir = tmp_path / 'build'
    _write_pdb(a_path, stateA)
    _write_pdb(b_path, stateB)

    result = cli.run(str(a_path), str(b_path), str(outdir), mode='auto', verbose=False)
    atom_map_file = outdir / 'atom_map.json'
    report_file = outdir / 'build_report.md'

    assert atom_map_file.exists()
    assert report_file.exists()

    payload = json.loads(atom_map_file.read_text())
    assert payload['mode'] == 'residue'
    assert payload['a_only'] == []
    assert len(payload['b_only']) == 2
    assert 'A:1:CB' in [entry['a_id'] for entry in payload['common']]
