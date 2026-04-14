---
name: io-py
description: "Use when: working with PDB parsing/writing, extending amino acid compatibility, fixing PDB format issues, or implementing new I/O features in IO.py. Reference for PDB v3.3 format specifications."
applyTo: "**/IO.py"
---

# IO.py PDB Parsing & I/O Guide

## Overview

`IO.py` implements a **custom PDB parser** following the **Protein Data Bank v3.3 format specification**. It does not use external libraries (like BioPython) to avoid dependencies. This file handles:
- PDB line parsing (`pdb_parse_in`)
- PDB line writing (`pdb_parse_out`)
- Amino acid conversions (1-letter ↔ 3-letter codes)
- Charged residue state management (HIS, GLU, ASP variants)
- Atom name/type utilities

---

## PDB v3.3 Format Specification

Each ATOM/HETATM line follows fixed-width column positions:

```
Columns    Data Type       Description
-----------------------------------
1-6        Record name     "ATOM  " or "HETATM"
7-11       Integer         Atom serial number
13-16      Atom name       (e.g., "CA", "CB", "N")
17         Character       Alternate location indicator (usually blank)
18-21      Residue name    (e.g., "ALA", "GLY")
22         Character       Chain identifier (usually "A", "B", etc.)
23-26      Integer         Residue sequence number
27         Character       Code for insertion of residue (usually blank)
31-38      Real(8.3)       X orthogonal coordinate (Ångströms)
39-46      Real(8.3)       Y orthogonal coordinate
47-54      Real(8.3)       Z orthogonal coordinate
55-60      Real(6.2)       Occupancy (optional, defaults to 0.0)
61-66      Real(6.2)       Temperature factor (optional, defaults to 0.0)
77-78      LString(2)      Element symbol (optional, e.g., "C", "N")
79-80      LString(2)      Charge on atom (optional, e.g., "+1", "-1")
```

---

## Function Reference

### `pdb_parse_in(line, include=('ATOM','HETATM'))`

**Purpose**: Parse a single PDB line into a structured list

**Returns**: List of 15 elements:
```python
[
    0: "ATOM" or "HETATM",           # Record type
    1: 1001,                          # Atom serial (int)
    2: "CA",                          # Atom name (str)
    3: "",                            # Alternate location (str)
    4: "ALA",                         # Residue name (str)
    5: "A",                           # Chain ID (str)
    6: 1,                             # Residue number (int)
    7: "",                            # Insertion code (str)
    8: 1.234,                         # X coordinate (float)
    9: 5.678,                         # Y coordinate (float)
    10: 9.012,                        # Z coordinate (float)
    11: 1.00,                         # Occupancy (float, default 0.0)
    12: 25.34,                        # B-factor (float, default 0.0)
    13: "C",                          # Element symbol (str, default "  ")
    14: "",                           # Charge (str, default "  ")
]
```

**Example**:
```python
pdb_line = "ATOM   1001  CA  AALA A   1       1.234   5.678   9.012  1.00 25.34           C"
parsed = pdb_parse_in(pdb_line)
print(parsed[2])  # "CA"
print(parsed[6])  # 1 (residue number)
```

**Edge Cases**:
- Occupancy/B-factor/element missing → defaults to 0.0 or empty string
- Atom names can be 3 or 4 characters (e.g., "C", "CA", "HG21")
- Chain ID defaults to "A" if blank

---

### `pdb_parse_out(parsed_list)`

**Purpose**: Convert a parsed list back to a PDB-formatted line

**Input**: 15-element list from `pdb_parse_in()` (or manually constructed)

**Returns**: PDB-formatted string (fixed-width columns)

**Example**:
```python
parsed_atom = [
    "ATOM", 1001, "CA", "", "ALA", "A", 1, "",
    1.234, 5.678, 9.012, 1.00, 25.34, "C", ""
]
pdb_line = pdb_parse_out(parsed_atom)
# Output: "ATOM   1001  CA  AALA A   1       1.234   5.678   9.012  1.00 25.34           C"
```

**Formatting Rules**:
- **3-char atom name** (len ≤ 3): Space-padding pattern `  {name}{space}`
- **4-char atom name** (len == 4): Space-padding pattern ` {name}`
- Coordinates formatted as `8.3f` (8 chars, 3 decimals)
- Occupancy/B-factor as `6.2f` (6 chars, 2 decimals)

---

### `AA(residue_code)`

**Purpose**: Convert between 1-letter, 3-letter, and 4-letter amino acid codes

**Returns**: Converted code as string

**Supported Conversions**:
```python
AA('ALA')    # → 'A'    (3-letter to 1-letter)
AA('A')      # → 'ALA'  (1-letter to 3-letter)
AA('CALA')   # → 'A'    (4-letter to 1-letter, e.g., C-terminal)
```

**Special Cases** (charged residues):
- `HIS` → `HID`, `HIE`, `HIP` (different protonation states)
- `GLU` → `GLH` (protonated carboxyl)
- `ASP` → `ASH` (protonated carboxyl)
- `LYS` → `LYN` (deprotonated amine)

**Example**:
```python
AA('GLH')   # → 'E'  (protonated GLU)
AA('HID')   # → 'H'  (HIS with delta proton)
```

---

### Standard Amino Acids Registry

The module maintains comprehensive atom name lists and charged residue mappings:

```python
atoms = ['N', 'H', 'C', 'O', 'CA', 'HA', ..., 'HH21', 'HH22']

charged_res = {
    'HIS': {'HD1': 'HID', 'HE2': 'HIE'},
    'GLU': {'HE2': 'GLH'},
    'ASP': {'HD2': 'ASH'}
}
```

**Edit when**:
- Adding new non-standard amino acids (e.g., selenocysteine, glyceraldehyde)
- Supporting different protonation states
- Adding new atom types for extended forcefields

---

## Common Patterns

### Reading a PDB File

```python
import IO

with open('protein.pdb', 'r') as f:
    atoms = []
    for line in f:
        parsed = IO.pdb_parse_in(line)
        if isinstance(parsed, list) and len(parsed) == 15:
            atoms.append(parsed)

print(f"Read {len(atoms)} atoms")
```

### Writing a Modified PDB File

```python
import IO

modified_atoms = []
for atom in atoms:
    # Modify coordinates (index 8-10)
    atom[8] += 0.1  # Shift X
    # Modify occupancy (index 11)
    atom[11] = 1.0
    modified_atoms.append(atom)

with open('protein_modified.pdb', 'w') as f:
    for atom in modified_atoms:
        f.write(IO.pdb_parse_out(atom) + '\n')
```

### Filtering Atoms by Residue

```python
# Get all CA atoms
ca_atoms = [atom for atom in atoms if atom[2] == 'CA']

# Get atoms in residue 42
res42_atoms = [atom for atom in atoms if atom[6] == 42]

# Get atoms in chain B
chain_b = [atom for atom in atoms if atom[5] == 'B']
```

### Converting Residue Names

```python
import IO

# Convert all ALA to 3-letter code
for atom in atoms:
    if atom[1] == 'A':  # 1-letter code
        atom[4] = IO.AA('A')  # → 'ALA'

# Handle charged variants
if atom[4] == 'GLU':
    # Check if this should be protonated (GLH)
    if should_protonate:
        atom[4] = IO.AA('GLH')
```

---

## Troubleshooting

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| `pdb_parse_in` returns raw line (not list) | Line doesn't start with "ATOM" or "HETATM" | Check for HEADER, REMARK, or non-standard records |
| Coordinates are `None` or incorrect | Column positions misaligned | Verify PDB file uses v3.3 format; check for unusual whitespace |
| Atom names not recognized | Atom type not in `atoms` list | Add to `atoms` list or check if PDB uses non-standard naming |
| AA conversion fails with KeyError | Unsupported residue code | Check charged_res dict and oneAA/threeAA/fourAA mappings |
| Output PDB doesn't parse in PyMOL/Q | Format spacing wrong | Ensure `pdb_parse_out` receives correctly structured 15-element list |

---

## When to Modify IO.py

**Do**:
- Add support for new residue types (e.g., post-translational modifications)
- Extend charged residue mappings (e.g., different HIS protonation states)
- Add helper functions for common PDB manipulations
- Support alternative file formats (GRO, MOL2) as new functions

**Don't**:
- Change fixed-width column positions (this breaks the PDB format)
- Introduce external dependencies (BioPython, MDTraj) as core imports
- Assume all residues follow standard naming (validate input)

---

## PDB v3.3 Specification Reference

- Official PDB Format Guide: https://www.wwpdb.org/documentation/file-format
- Column Details: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
- ATOM/HETATM Records: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

