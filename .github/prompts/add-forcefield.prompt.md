---
name: add-forcefield
description: "Use when: creating a new forcefield converter (opls, charmm, amber, openff variants). Guides step-by-step creation of converter script, parameter file registration, and CLI integration."
---

# Add a New Forcefield to QligFEP

Follow this workflow to integrate a new forcefield (new variant or completely new FF):

## Step 1: Create the Converter Script

1. **Examine existing converters** for your source format:
   - OPLS → `opls2Q.py` 
   - CHARMM → `charmm2Q.py`
   - OpenFF → `openff2Q.py`

2. **Create `<FFNAME>2Q.py`** (e.g., `amber14sb2Q.py`) in the workspace root:
   ```python
   """
   Convert <FF_NAME> forcefield parameters to Q format.
   
   Usage: python <FFNAME>2Q.py <input_file> [options]
   """
   import argparse
   import numpy as np
   
   class Converter:
       """Parse <FF_NAME> params and write Q .prm file"""
       
       def __init__(self, input_file):
           self.input_file = input_file
           self.params = {}
       
       def parse(self):
           """Extract atom types, bonds, angles, dihedrals, VdW"""
           # TODO: Implement parsing logic
           pass
       
       def write_Q_format(self, output_file):
           """Write params in Q .prm format (see FF/OPLS2015.prm)"""
           # TODO: Format as:
           # atom type definitions
           # bond parameters
           # angle parameters
           # dihedral parameters
           # VdW parameters (vdw atom_type r_min epsilon)
           pass
   
   if __name__ == "__main__":
       parser = argparse.ArgumentParser()
       parser.add_argument("input", help="Input FF file")
       parser.add_argument("-o", "--output", default="converted.prm")
       args = parser.parse_args()
       
       conv = Converter(args.input)
       conv.parse()
       conv.write_Q_format(args.output)
       print(f"✓ Converted to {args.output}")
   ```

3. **Reference Q parameter format** in `FF/OPLS2015.prm`:
   - Atom types: `atom <name> <mass> <charge>...`
   - Bonds: `bond <type1> <type2> <force_const> <eq_dist>`
   - Angles: `angle <t1> <t2> <t3> <fc> <eq_angle>`
   - Dihedrals: `dihedral <t1> <t2> <t3> <t4> <fc> <mult> <phase>`

## Step 2: Add Parameter File

1. **Copy or generate the `.prm` file**:
   ```bash
   # Example: for AMBER14SB variant
   cp FF/AMBER14sb.prm FF/AMBER14sb-NEW.prm  # or use converter script
   ```

2. **Verify format**:
   - Check against existing `.prm` files in `FF/`
   - Ensure all atom types are defined
   - Validate bond/angle/dihedral coverage

## Step 3: Register in `settings.py`

1. **Add to the `FORCEFIELDS` dict**:
   ```python
   FORCEFIELDS = {
       'OPLS2015': 'OPLS2015.prm',
       'CHARMM36': 'CHARMM36.prm',
       'AMBER14sb': 'AMBER14sb.prm',
       'YOUR_FF': 'YOUR_FF.prm',  # ← Add here
   }
   ```

2. **If it requires special converter**, add to converter registry:
   ```python
   CONVERTERS = {
       'opls': 'opls2Q.py',
       'charmm': 'charmm2Q.py',
       'your_format': 'your_format2Q.py',  # ← Add here
   }
   ```

## Step 4: Add CLI Option

1. **In `QligFEP.py`, `QresFEP.py`, `QLIE.py` argparse setup**:
   ```python
   parser.add_argument(
       '-f', '--forcefield',
       choices=['OPLS2015', 'CHARMM36', 'AMBER14sb', 'YOUR_FF'],  # ← Add here
       default='OPLS2015',
       help='Forcefield variant'
   )
   ```

2. **Verify that `load_forcefield()` in settings.py** loads your .prm file correctly

## Step 5: Test Integration

1. **Test the converter script**:
   ```bash
   python YOUR_FF2Q.py test_input.prm -o FF/YOUR_FF.prm
   ```

2. **Test CLI recognition**:
   ```bash
   python QligFEP.py -h  # Should show your FF in choices
   python QligFEP.py -l1 lig1.pdb -l2 lig2.pdb -f YOUR_FF -S prot.pdb -C localhost
   ```

3. **Verify generated input files** use correct parameters from your .prm file

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Unknown forcefield" | Check `FORCEFIELDS` dict in settings.py |
| Wrong parameters in output | Verify .prm file format matches OPLS example |
| Converter script fails | Check input format matches expected structure |
| Missing atom types | Add missing entries to .prm file or converter output |

## Quick Reference

- **Parameter file location**: `FF/`
- **Converter script location**: Root directory
- **Settings registration**: `settings.py` (FORCEFIELDS & CONVERTERS dicts)
- **CLI validation**: Updated in argparse `choices=` for each tool
- **Documentation**: Add entry to README.md showing new FF availability
