---
name: forcefield-development
description: "Use when: adding a new forcefield variant or completely new forcefield system. Multi-step workflow from parameter extraction → converter script → registration → testing. Includes templates and validation."
---

# Forcefield Development Workflow

Complete guide for integrating a new forcefield into QligFEP/QresFEP.

## Overview

This workflow covers **three scenarios**:
1. **New variant** of existing FF (e.g., OPLS2005 → OPLS2020)
2. **New FF from external source** (e.g., CHARMM, AMBER variants)
3. **Custom parametrization** (e.g., specialized proteins, ligands)

---

## Phase 1: Acquire Parameters

### Option A: From Files

If you have a `.prm`, `.ff`, or `.par` file:

```bash
# Copy to FF/ directory
cp /path/to/yourff.prm FF/YOURFF.prm

# Verify format (should look like FF/OPLS2015.prm):
head -20 FF/YOURFF.prm
```

Expected format:
```
atom  C_opls    12.011    0.5
bond  C   H    340.0    1.090
angle CAr  CAr  CAr  70.0  120.0
dihedral  CA  CA  CA  HA  3.625  2  180.0
vdw  C_opls  3.905  0.066
```

### Option B: From Schrodinger (ffld_server)

If using CHARMM or other forcefields via Schrodinger:

```bash
# Generate ligand parameters via Schrodinger
$SCHROD/utilities/ffld_server -iMACMINI -a CHARMM_ALL your_ligand.mol2 -o your_ligand.mae

# Extract forcefield parameters
# (Requires OPLS/CHARMM/AMBER ffld_server module)
```

### Option C: From Literature / Ab Initio

If deriving from publications or QM calculations:

1. **Extract from paper tables** or supplementary materials
2. **Format to Q format** (see template converter below)
3. **Implement in converter script** (Phase 2)

---

## Phase 2: Create Converter Script

Create `<FFNAME>2Q.py` in workspace root if converting from external format.

### Template: Convert External Format to Q

```python
#!/usr/bin/env python
"""
Convert <FFNAME> forcefield parameters to Q format.

Usage:
    python <FFNAME>2Q.py input.ff -o output.prm
    
Reference:
    Input format: <FFNAME> specification (version X.X)
    Output format: Q .prm format (see FF/OPLS2015.prm)
"""

import argparse
import re
from pathlib import Path


class FFConverter:
    """Parse <FFNAME> parameters and write Q .prm format"""
    
    def __init__(self, input_file):
        self.input_file = Path(input_file)
        self.atoms = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.vdws = {}
    
    def parse_input(self):
        """Extract atom types, bonds, angles, dihedrals, VdW from input file"""
        with open(self.input_file, 'r') as f:
            lines = f.readlines()
        
        # Parse by section (customize based on input format)
        section = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if line.startswith('[ATOMS]'):
                section = 'atoms'
                continue
            elif line.startswith('[BONDS]'):
                section = 'bonds'
                continue
            elif line.startswith('[ANGLES]'):
                section = 'angles'
                continue
            elif line.startswith('[DIHEDRALS]'):
                section = 'dihedrals'
                continue
            elif line.startswith('[VDW]'):
                section = 'vdw'
                continue
            
            # Parse based on current section
            if section == 'atoms':
                parts = line.split()
                if len(parts) >= 3:
                    atom_type, mass, charge = parts[0], float(parts[1]), float(parts[2])
                    self.atoms[atom_type] = (mass, charge)
            
            elif section == 'bonds':
                parts = line.split()
                if len(parts) >= 4:
                    t1, t2, fc, r_eq = parts[0], parts[1], float(parts[2]), float(parts[3])
                    self.bonds[(t1, t2)] = (fc, r_eq)
            
            elif section == 'angles':
                parts = line.split()
                if len(parts) >= 5:
                    t1, t2, t3, fc, theta = parts[0], parts[1], parts[2], float(parts[3]), float(parts[4])
                    self.angles[(t1, t2, t3)] = (fc, theta)
            
            elif section == 'dihedrals':
                parts = line.split()
                if len(parts) >= 6:
                    t1, t2, t3, t4, fc, mult = parts[0], parts[1], parts[2], parts[3], float(parts[4]), int(parts[5])
                    phase = float(parts[6]) if len(parts) > 6 else 0.0
                    self.dihedrals[(t1, t2, t3, t4)] = (fc, mult, phase)
            
            elif section == 'vdw':
                parts = line.split()
                if len(parts) >= 3:
                    atom_type, r_min, epsilon = parts[0], float(parts[1]), float(parts[2])
                    self.vdws[atom_type] = (r_min, epsilon)
    
    def write_Q_format(self, output_file):
        """Write parameters in Q .prm format"""
        with open(output_file, 'w') as f:
            # Header
            f.write("! <FFNAME> forcefield parameters for Q\n")
            f.write("! Auto-converted from <FFNAME>2Q.py\n\n")
            
            # Atom types
            f.write("! Atom types\n")
            for atom_type, (mass, charge) in self.atoms.items():
                f.write(f"atom {atom_type:8s} {mass:8.3f}  {charge:8.5f}\n")
            f.write("\n")
            
            # Bonds
            f.write("! Bonds\n")
            for (t1, t2), (fc, r_eq) in self.bonds.items():
                f.write(f"bond {t1:8s} {t2:8s} {fc:8.2f}  {r_eq:8.4f}\n")
            f.write("\n")
            
            # Angles
            f.write("! Angles\n")
            for (t1, t2, t3), (fc, theta) in self.angles.items():
                f.write(f"angle {t1:8s} {t2:8s} {t3:8s} {fc:8.2f}  {theta:8.2f}\n")
            f.write("\n")
            
            # Dihedrals
            f.write("! Dihedrals\n")
            for (t1, t2, t3, t4), (fc, mult, phase) in self.dihedrals.items():
                f.write(f"dihedral {t1:8s} {t2:8s} {t3:8s} {t4:8s} {fc:8.3f}  {mult:2d}  {phase:8.2f}\n")
            f.write("\n")
            
            # VdW
            f.write("! Van der Waals\n")
            for atom_type, (r_min, epsilon) in self.vdws.items():
                f.write(f"vdw {atom_type:8s} {r_min:8.4f}  {epsilon:8.5f}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input <FFNAME> forcefield file')
    parser.add_argument('-o', '--output', default='output.prm', help='Output Q .prm file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    try:
        converter = FFConverter(args.input)
        converter.parse_input()
        
        if args.verbose:
            print(f"✓ Parsed {len(converter.atoms)} atom types")
            print(f"✓ Parsed {len(converter.bonds)} bond parameters")
            print(f"✓ Parsed {len(converter.angles)} angle parameters")
            print(f"✓ Parsed {len(converter.dihedrals)} dihedral parameters")
            print(f"✓ Parsed {len(converter.vdws)} VdW parameters")
        
        converter.write_Q_format(args.output)
        print(f"\n✓ Converted to {args.output}")
        
    except Exception as e:
        print(f"✗ Conversion failed: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
```

### Run converter

```bash
python YOURFF2Q.py input.ff -o FF/YOURFF.prm -v
```

---

## Phase 3: Register in settings.py

### Add to FORCEFIELDS dict

```python
# In settings.py, find FORCEFIELDS dict:
FORCEFIELDS = {
    'OPLS2015': 'OPLS2015.prm',
    'OPLS2005': 'OPLS2005.prm',
    'CHARMM36': 'CHARMM36.prm',
    'AMBER14sb': 'AMBER14sb.prm',
    'YOURFF': 'YOURFF.prm',  # ← ADD HERE
}
```

### Add to CONVERTERS dict (if converter script used)

```python
CONVERTERS = {
    'opls': 'opls2Q.py',
    'charmm': 'charmm2Q.py',
    'yourff': 'YOURFF2Q.py',  # ← ADD HERE
}
```

### Verify load function

```python
def load_forcefield(ff_name):
    """Load forcefield parameters from FF/ directory"""
    if ff_name not in FORCEFIELDS:
        raise KeyError(f"Unknown forcefield: {ff_name}. Available: {list(FORCEFIELDS.keys())}")
    
    prm_file = os.path.join(FF_DIR, FORCEFIELDS[ff_name])
    if not os.path.exists(prm_file):
        raise FileNotFoundError(f"Forcefield file not found: {prm_file}")
    
    # Load and parse prm file (implementation details)
    return prm_file
```

---

## Phase 4: Update CLI Tools

### QligFEP.py, QresFEP.py, QLIE.py

Add your FF to argparse choices:

```python
# In each tool's argument parser:
parser.add_argument(
    '-f', '--forcefield',
    choices=['OPLS2015', 'OPLS2005', 'CHARMM36', 'AMBER14sb', 'YOURFF'],
    default='OPLS2015',
    help='Forcefield variant'
)
```

---

## Phase 5: Validation & Testing

### 1. Test Converter Script (if applicable)

```bash
python YOURFF2Q.py test_input.ff -o FF/YOURFF.prm -v
head -50 FF/YOURFF.prm  # Inspect output
```

### 2. Test FF Loading

```python
import settings as s

# Verify registration
print("Available FF:", list(s.FORCEFIELDS.keys()))

# Test loading
ff_path = s.load_forcefield('YOURFF')
print(f"Loaded from: {ff_path}")
```

### 3. Test CLI Recognition

```bash
# Check that FF appears in help
python QligFEP.py -h | grep -A 5 "forcefield"
python QresFEP.py -h | grep -A 5 "forcefield"
python QLIE.py -h | grep -A 5 "forcefield"
```

### 4. Run Small Test Case

```bash
# Create minimal test structures (or use tutorials)
cd /tmp
python QligFEP.py \
  -l1 tutorials/1.QligFEP_CDK2/1.ligprep/1h1q.pdb \
  -l2 tutorials/1.QligFEP_CDK2/1.ligprep/1h1r.pdb \
  -f YOURFF \
  -S tutorials/1.QligFEP_CDK2/2.protprep/complex_processed.pdb \
  -C localhost \
  -n 3 \
  --no-submit

# Inspect generated files
ls test_setup/
cat test_setup/FEP_submit.sh  # Should show YOURFF parameters
cat test_setup/1_2_FEP/inputfiles/qfep.inp  # Check forcefield refs
```

### 5. Verify Physical Parameters

Compare against reference values:

```bash
# Example: Bond lengths should be reasonable (1.0 - 2.0 Å)
grep "^bond" FF/YOURFF.prm | awk '{print $5}' | sort -n

# Angles should be 90-180°
grep "^angle" FF/YOURFF.prm | awk '{print $6}' | sort -n

# Dihedrals should have periodic multiplicity (1-6)
grep "^dihedral" FF/YOURFF.prm | awk '{print $7}' | sort -n
```

---

## Phase 6: Documentation

### Update README.md

Add entry to forcefield list:

```markdown
## Supported Forcefields

- **OPLS2015**: OPLS2015 parameters
- **OPLS2005**: OPLS2005 parameters  
- **CHARMM36**: CHARMM36 all-atom parameters
- **AMBER14sb**: AMBER 14 SB forcefield
- **YOURFF**: [Reference/Source] - [Key features]
```

### Document in settings.py

Add comment:

```python
FORCEFIELDS = {
    # ...
    'YOURFF': 'YOURFF.prm',  # Your FF v1.0, optimized for [application]
                             # Reference: [Citation]
                             # Added: [Date]
}
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Unknown forcefield" error | Check FORCEFIELDS dict in settings.py |
| Converter script fails | Verify input file format matches expectations; check parsing logic |
| Missing atom types in parameters | Add missing entries to .prm file or converter output |
| Wrong parameters in Q output | Verify .prm file format (see FF/OPLS2015.prm example) |
| CLI doesn't recognize FF | Update argparse `choices=` in QligFEP/QresFEP/QLIE.py |
| FEP simulation diverges | Check parameter ranges (bonds, angles, dihedrals) are physical |

---

## Next Steps

1. **Benchmark against known systems** (e.g., CDK2 from tutorials)
2. **Compare energies** with original forcefield implementation
3. **Share with community** (open PR or discussion)
4. **Document lessons learned** in project wiki

---

## References

- **Q Forcefield Format**: See FF/*.prm examples
- **OPLS Manual**: https://www.schrodinger.com/
- **CHARMM Documentation**: https://www.charmm.org/
- **AMBER Manual**: http://ambermd.org/

