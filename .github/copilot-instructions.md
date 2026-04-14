# QligFEP & QresFEP Workspace Guidelines

## Overview

QligFEP/QresFEP is a Python package that automates **FEP (Free Energy Perturbation) calculation setup** for the **Q software** package. It provides three main CLI tools for generating ligand mutations, protein mutations, and LIE simulations.

**Key Papers**:
- QligFEP: Jespers et al. https://doi.org/10.1186/s13321-019-0348-5
- QresFEP: Jespers et al. https://doi.org/10.1021/acs.jctc.9b00538 and Koenekoop et al. https://doi.org/10.1038/s42004-025-01771-0

---

## Code Organization

### Main Entry Points (CLI)
- **QligFEP.py** (`Run` class): Dual-topology ligand FEP from 2 ligand structures → generates Q input files + HPC submission scripts
- **QresFEP.py** (`Run` class): Protein mutations (single/dual topology) → similar output pipeline
- **QLIE.py** (`Run` class): LIE (Linear Interaction Energy) MD simulations

### Core Utilities
- **IO.py**: PDB parsing/writing (follows PDB v3.3 format), amino acid conversions, charged residue utilities
- **functions.py**: Lambda spacing (linear, sigmoid, sigmoidal, exponential), geometric overlays, center-of-geometry
- **settings.py**: Runtime config template (Q paths, HPC clusters, forcefields)—populated by `qligfep_init.sh`

### Supporting Tools
- **Forcefield Converters**: `opls2Q.py`, `charmm2Q.py`, `openff2Q.py` (converts external FF → Q parameter format)
- **Analysis Tools**: `analyze_FEP.py`, `analyze_LIE.py` (post-processing, energy extraction)
- **Protein Prep**: `protprep.py` (prepares systems for spherical boundary simulations)
- **Utility Scripts** in `scripts/`: Alascan, counter-ions, center-of-geometry, renumbering, etc.

### Directory Structure
```
FF/              # Forcefield parameter files (AMBER, CHARMM, OPLS variants)
INPUTS/          # Template Q input files (MD, equilibration, FEP) + submission scripts
tutorials/       # Two end-to-end examples: CDK2 ligands (QligFEP) + T4L mutants (QresFEP)
src/             # (Currently minimal—package module)
```

---

## Build & Setup

### Initial Setup
```bash
bash qligfep_init.sh  # Interactive setup (prompts for Q path, Schrodinger path, default HPC cluster)
```
This script:
1. Sets `QLIGFEP` environment variable
2. Adds repo to `$PATH` and `$PYTHONPATH`
3. Populates `settings.py` with Q install location and HPC cluster configs

### Running the Tools
```bash
python QligFEP.py -l1 ligand1.pdb -l2 ligand2.pdb -f OPLS2015 -S protein.pdb -C cluster_name [options]
python QresFEP.py -m A99V -mc A -S protein.pdb -C cluster_name [options]
python QLIE.py -l ligand.pdb -f OPLS2015 -S water.pdb -C cluster_name [options]
```

**Each tool uses argparse.** Run with `-h` for full option documentation.

### Dependencies
- **Python 3.7+** (officially Python 3 only; Python 2 branch archived)
- **External**: Q software (FEP engine), Schrodinger suite (PyMOL, Protein Prep Wizard, ffld_server), CGenFF
- **Python packages**: numpy, mdtraj (optional), matplotlib (optional), rdkit (optional, for LIE SMARTS)

**Note**: No CI/test infrastructure exists—add testing infrastructure for new features.

---

## Conventions

### Naming & Style
- **Files**: kebab-case (`analyze_FEP.py`)
- **Functions/Variables**: snake_case
- **Constants**: CAPS
- **Classes**: PascalCase

### Module Pattern
Each tool (QligFEP, QresFEP, QLIE) follows a common structure:
1. `Run` class with modular methods for each pipeline step
2. Argparse CLI with option validation
3. Instantiates `Run` class → executes pipeline → outputs Q scripts

### Lambda Spacing
Built-in functions: `linear()`, `sigmoid()`, `sigmoidal()`, `exponential()`  
Used to distribute FEP intermediate states—each has trade-offs for QM/MM sampling.

### I/O Patterns
- **PDB Parsing**: Custom parser (does not rely on external libraries)—see `IO.py` for format spec
- **Batch Operations**: Heavy use of `glob()` for processing multiple files
- **Configuration**: Templated inputs in `INPUTS/` → filled with parameters → saved per-system

---

## Development Workflow

### Adding a New Forcefield
1. Create converter script: `<ffname>2Q.py` (follow `opls2Q.py`, `charmm2Q.py` patterns)
2. Add `.prm` file to `FF/` directory
3. Register in `settings.py` and CLI option validation

### Adding a New Analysis Tool
1. Create `analyze_<name>.py` with energy parsing logic
2. Follow `analyze_FEP.py` structure (import Q trajectory, parse energies, plot results)

### Adding HPC Cluster Support
1. Edit `settings.py`: Add cluster config dict with keys: `Q_DIR`, `MODULES`, `NTASKS`, `TIME`, `PARTITION`, command paths
2. Verify submission script generation (see `INPUTS/FEP_submit*.sh` templates)

### Extending the Tutorial
- `tutorials/1.QligFEP_CDK2/`: Ligand mutation example (OPLS FF, protein prep, FEP setup, analysis)
- `tutorials/2.QresFEP_T4L/`: Residue mutation example (different topologies, thermal stability)

---

## Common Pitfalls

1. **Missing Q installation**: Tools fail silently if `settings.py` not initialized. Run `qligfep_init.sh` before first use.
2. **Schrodinger path issues**: Some converters (`ffld_server`, `Protein Prep`) require correct Schrodinger location in `settings.py`.
3. **PDB format assumptions**: Custom PDB parser expects v3.3 format—molecular dynamics output varies; use PyMOL or `protprep.py` to standardize.
4. **Lambda state count**: Ensure number of lambda states matches Q simulation capacity (typically 51-101 states).
5. **Spherical boundary systems**: `protprep.py` is required when using `qprep_protprep.inp`; Cartesian systems don't need it.

---

## Key Files for Agents

| File | Purpose | Edit When |
|------|---------|-----------|
| `settings.py` | HPC/Q configuration, FF paths | Adding clusters or forcefields |
| `functions.py` | Lambda spacing algorithms | Adding new sampling strategies |
| `IO.py` | PDB I/O, standard residues | Extending residue compatibility |
| `QligFEP.py`, `QresFEP.py`, `QLIE.py` | Main workflows | Adding new options or pipeline stages |
| `FF/*.prm` | Forcefield parameters | Updating FF versions |
| `INPUTS/*.inp` | Q template input files | Changing MD/FEP defaults |

---

## External Requirements

- **Q Software**: https://github.com/esguerra/Q6 (FEP energy calculations)
- **Schrodinger Suite**: PyMOL, Protein Preparation Wizard, ffld_server (structure prep, parameter generation)
- **CGenFF**: For CHARMM forcefield parameterization

These are configured in `settings.py` during initialization.

---

## Further Reading

- **Full README**: [README.md](../../README.md)
- **QligFEP Tutorial**: [tutorials/1.QligFEP_CDK2/README.md](../../tutorials/1.QligFEP_CDK2/README.md)
- **QresFEP Tutorial**: [tutorials/2.QresFEP_T4L/README.md](../../tutorials/2.QresFEP_T4L/README.md)
- **Contact**: w.jespers@rug.nl
